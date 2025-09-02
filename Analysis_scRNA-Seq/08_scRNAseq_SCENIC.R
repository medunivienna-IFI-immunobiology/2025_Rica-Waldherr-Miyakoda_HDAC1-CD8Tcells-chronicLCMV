# title: analysis using SCENIC (Single Cell rEgulatory Network Inference and Clustering)
# author: Monika Waldherr
# input: scRNA-seq expression matrix with geneIDs as gene-symbols as rownames and cells as columns

# using the old R version from https://github.com/aertslab/SCENIC

#### install required packages ####
BiocManager::install(c("AUCell"))
# special for arrow which is needed for RcisTarget
Sys.setenv("NOT_CRAN" = "true")
install.packages("arrow")

BiocManager::install("RcisTarget")
BiocManager::install(c("GRNBoost2")) # faster alternative (slow: GENIE3)
BiocManager::install(c("GENIE3")) # installed bc SCENIC required it

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools"))
remotes::install_github("bokeh/rbokeh")
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# To support paralell execution:
BiocManager::install(c("doMC", "doRNG"))
# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

## install SCENIC
devtools::install_github("aertslab/SCENIC")
packageVersion("SCENIC")

## get species-specific databases
# By default, SCENIC uses the databases that score the motifs 
# in the promoter of the genes (up to 500bp upstream the TSS), 
# and in the 20kb around the TSS (+/-10kbp).
# for mouse (use mm10):
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/old/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/old/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
# mc9nr: Motif collection version 9: 24k motifs
dir.create("cisTarget_databases")
setwd("cisTarget_databases")
for(featherURL in dbFiles){
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
}
setwd("..")

#### create working directory ####
dir.create("SCENIC_cd8t")
setwd("SCENIC_cd8t")
dir.create("int")

#### load data from our Seurat object ####
library(SeuratObject)
library(Seurat)
cd8t_cluster <- readRDS("../cd8t_cluster_ownIdentsRes0.6.rds")
exprMat_df <- data.frame(cd8t_cluster@assays$SCT@counts)
exprMat <- as.matrix.data.frame(exprMat_df, rownames.force = T)
write.table(exprMat, file = "int/exprMat.tab", row.names = T)

# to use Seurat clusters as cell annotation (e.g. for visualization):
cellInfo <- data.frame(seuratCluster=Idents(cd8t_cluster))
cellInfo$nGene <- colSums(exprMat>0)
cellInfo$genotype <- cd8t_cluster@meta.data$genotype
metadata <- cd8t_cluster@meta.data
metadata$TexGroup <- metadata$Res0.6Names
metadata$TexGroup <- gsub("TexProg", "early-TexProg", metadata$TexGroup)
metadata$TexGroup <- gsub("TexProl", "early-non-TexProg", metadata$TexGroup)
metadata$TexGroup <- gsub("TexExh", "early-non-TexProg", metadata$TexGroup)
metadata$TexGroup <- gsub("TexEarly", "early-non-TexProg", metadata$TexGroup)
metadata$TexGroup <- gsub("TexInt", "early-non-TexProg", metadata$TexGroup)
metadata$TexGroup <- gsub("TexCX3CR1", "early-non-TexProg", metadata$TexGroup)
metadata$TexGroup <- gsub("TexCyt", "early-non-TexProg", metadata$TexGroup)
cellInfo$TexGroup <- metadata$TexGroup
head(cellInfo)
cbind(table(cellInfo$seuratCluster))
cbind(table(cellInfo$genotype))
saveRDS(cellInfo, file="int/cellInfo.Rds")
cellInfo <- readRDS("int/cellInfo.Rds")
# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(seuratCluster = c("Naive" = "seashell3",
                                  "TexProg" = "#0084B3",
                                  "TexProl" = "#00B1B5",
                                  "TexExh" = "palevioletred1",
                                  "TexEarly" = "#A71B4B",
                                  "TexInt" = "#F7B347",
                                  "TexCX3CR1" = "khaki3",
                                  "TexCyt" = "darkseagreen"))
colVars$seuratCluster <- colVars$seuratCluster[intersect(names(colVars$seuratCluster), cellInfo$seuratCluster)]
saveRDS(colVars, file="int/colVars.Rds")
colVars <- readRDS("int/colVars.Rds")
#plot.new(); legend(0,1, fill = colVars$seuratCluster, legend = names(colVars$seuratCluster))

#### initialize SCENIC settings ####
library(SCENIC)
library(data.table)
library(AUCell)

org <- "mgi" # specify organism
dbDir <- "../cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "SCENIC analysis CD8+ Tcells" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)
motifAnnotations_mgi <- motifAnnotations # bug fix, object is created but with wrong name
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)

write.table(motifAnnotations_mgi, file = "int/motifAnno.tab")
motifAnnotations_mgi <- read.table("int/motifAnno.tab", header = T)


# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# Databases:
# scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
# scenicOptions@settings$db_mcVersion <- "v8"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
scenicOptions <- readRDS("int/scenicOptions.Rds")

#### Co-expression network ####
# first step of SCENIC is to infer potential TF targets based on the expression data
# for this GENIE3 is used if datasets are small (runtime of several days for 3-5K cells)
# for larger datasets GRNboost was implemented and should run much faster
# both allow to identify also non-linear relationships

## gene filter/selection
genesKept <- geneFiltering(exprMat, scenicOptions = scenicOptions,
                           minCountsPerGene = 2*.01*ncol(exprMat), # genes have at least a count of 2 in at least 1% of all cells
                           minSamples = ncol(exprMat)*.01) # genes are detected in at least 1% of all cells
interestingGenes <- c("Gzma", "Eomes", "Prf1", "Tox", "Gzmk", "Klre1", "Tnfrsf9", "Havcr2",
                      "Hdac1", "Cx3cr1",
                      "Hdac2", "Runx1", "Runx2", "Runx3", "Zeb2", "Prdm1", "Tbx21", "Irf4", "Id2",
                      "Ccl5", "Cd52", "Ctla2a", "Malat1", "Ms4a4b", "Hcst", "S100a4", "Fau", "Ifi214", "Ifi206",
                      "Ptma", "Rps17", "Hsp90ab1", "Anxa2")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)] # Ifi214 and Ifi206 are not in the database
# filter expression matrix
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered) # 6,002 x 16,393

## correlation
# in order to distinguish potential activation from repression, targets are split into positive and negative correlated targets
# i.e. Spearman correlation between the TF and the potential target)

# log the expression
exprMat_filtered <- log2(exprMat_filtered+1)
write.table(exprMat_filtered, file = "int/exprMat_filtered.tab", row.names = T)

# exprMat_filtered <- as.matrix.data.frame(read.table("int/exprMat_filtered.tab", header = T), rownames.force = T)

# calculate the correlation
runCorrelation(exprMat_filtered, scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

## Run GENIE3 (on the server since it uses a lot of RAM, did run ~30 hours)
# use script 08a to run on the server, it saves scenicOptions as Rds file afterwards

#### Build and score the GRN ####
# use script 08b to run on a server with enough RAM, it saves scenicOptions as Rds file afterwards

#### Visualize regulators ####
cellInfo <- readRDS("int/cellInfo.Rds")
rownames(cellInfo) <- gsub("-", ".", rownames(cellInfo))
colVars <- readRDS("int/colVars.Rds")
scenicOptions <- readRDS("int/scenicOptions_afterRunSCENIC_3.Rds")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

#### generate an AUC dotplot
library(reshape2)
library(dplyr)
library(ggplot2)

cellInfo_woNaive <- subset(cellInfo, TexGroup != "Naive")
cellInfo_woNaive$TexGroup <- paste(cellInfo_woNaive$TexGroup, cellInfo_woNaive$genotype, sep = " - ")
regulonAUC_woNaive <- regulonAUC[,rownames(cellInfo_woNaive)]
aucMat <- getAUC(regulonAUC_woNaive)

df <- melt(aucMat)
colnames(df) <- c("Regulon", "Cell", "AUC")
df$Cluster <- cellInfo_woNaive[df$Cell, "TexGroup"]

df_summary <- df %>%
  group_by(Cluster, Regulon) %>%
  summarize(MeanAUC = mean(AUC), .groups = "drop")

# log2FC dotplot KO vs WT AUC score for all regulons in TexProg and nonTexProg
#### TexProg ####
# get means of AUC per regulon for both groups seperately
df_TexProg_WT <- subset(df_summary, Cluster == "early-TexProg - WT")
df_TexProg_KO <- subset(df_summary, Cluster == "early-TexProg - HDAC1cKO")
TexProg_WT_mn1 <- round(df_TexProg_WT$MeanAUC, 3)
TexProg_KO_mn2 <- round(df_TexProg_KO$MeanAUC, 3)

df_pval <- df %>%
  filter(Cluster %in% c("early-TexProg - WT", "early-TexProg - HDAC1cKO")) %>%
  group_by(Regulon) %>%
  summarize(
    p_value = wilcox.test(AUC ~ Cluster)$p.value,
    .groups = "drop"
  )
e <- min(setdiff(c(TexProg_WT_mn1, TexProg_KO_mn2), 0))/10
log2_avgFC <- round(log2(TexProg_WT_mn1 + e) - log2(TexProg_KO_mn2 + e), 2)

myregulons <- data.frame(regulon = df_TexProg_WT$Regulon, log2_avgFC = log2_avgFC)
myregulons$adj_pval = p.adjust(df_pval$p_value, method = 'fdr')
myregulons_sorted <- myregulons[order(myregulons[[2]]), ]
myregulons_sorted$regulon <- factor(myregulons_sorted$regulon, 
                                    levels = myregulons_sorted$regulon[order(myregulons_sorted$log2_avgFC)])
myregulons_sorted_selected <- myregulons_sorted[abs(myregulons_sorted$log2_avgFC) > 0.25 & myregulons$adj_pval < 0.05, ]

# plot regulon activity
pdf(file = "DotPlot_RegulonActivity_TexProg_WT_vs_KO.pdf", width = 6, height = 10)
ggplot(myregulons_sorted_selected, aes(x = log2_avgFC, y = regulon, color = log2_avgFC, size = -log10(adj_pval))) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "lightgrey", high = "red", midpoint = 0, limits = c(-1.2, 1.2)) +
  theme_bw() +
  labs(
    title = "Regulon Activity in early-TexProg",
    subtitle = "WT vs HDAC1cKO",
    x = "average log2FC",
    y = "Regulon",
    color = "average log2FC"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # rotate x labels if many clusters
  ) +
  xlim(c(-1.2, 1.2))
dev.off()

#### nonTexProg ####
# get means of AUC per regulon for both groups seperately
df_nonTexProg_WT <- subset(df_summary, Cluster == "early-non-TexProg - WT")
df_nonTexProg_KO <- subset(df_summary, Cluster == "early-non-TexProg - HDAC1cKO")
nonTexProg_WT_mn1 <- round(df_nonTexProg_WT$MeanAUC, 3)
nonTexProg_KO_mn2 <- round(df_nonTexProg_KO$MeanAUC, 3)

df_pval <- df %>%
  filter(Cluster %in% c("early-non-TexProg - WT", "early-non-TexProg - HDAC1cKO")) %>%
  group_by(Regulon) %>%
  summarize(
    p_value = wilcox.test(AUC ~ Cluster)$p.value,
    .groups = "drop"
  )
e <- min(setdiff(c(nonTexProg_WT_mn1, nonTexProg_KO_mn2), 0))/10
log2_avgFC <- round(log2(nonTexProg_WT_mn1 + e) - log2(nonTexProg_KO_mn2 + e), 2)

myregulons <- data.frame(regulon = df_nonTexProg_WT$Regulon, log2_avgFC = log2_avgFC)
myregulons$adj_pval = p.adjust(df_pval$p_value, method = 'fdr')
myregulons$adj_pval[myregulons$adj_pval == 0] <- 10^-300
myregulons_sorted <- myregulons[order(myregulons[[2]]), ]
myregulons_sorted$regulon <- factor(myregulons_sorted$regulon, 
                                    levels = myregulons_sorted$regulon[order(myregulons_sorted$log2_avgFC)])
myregulons_sorted_selected <- myregulons_sorted[abs(myregulons_sorted$log2_avgFC) > 0.25 & myregulons$adj_pval < 0.05, ]

# plot regulon activity
pdf(file = "DotPlot_RegulonActivity_nonTexProg_WT_vs_KO.pdf", width = 6, height = 10)
ggplot(myregulons_sorted_selected, aes(x = log2_avgFC, y = regulon, color = log2_avgFC, size = -log10(adj_pval))) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "lightgrey", high = "red", midpoint = 0, limits = c(-1.3, 1.3)) +
  theme_bw() +
  labs(
    title = "Regulon Activity in early-non-TexProg",
    subtitle = "WT vs HDAC1cKO",
    x = "average log2FC",
    y = "Regulon",
    color = "average log2FC"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # rotate x labels if many clusters
  ) +
  xlim(c(-1.3, 1.3))
dev.off()




