# title: cluster analysis
# author: Monika Waldherr
# input: Seurat object saved after cluster characterization with script 05_scRNAseq_Clustering_characterization.R
#         marker gene file?

#### load libraries ####
library(SeuratObject)
library(Seurat)
library(EnhancedVolcano)

#### load Seurat object ####
cd8t_cluster <- readRDS("cd8t_cluster_ownIdentsRes0.6.rds")

#### differential gene expression analysis for KO vs WT in early-non-TexProg ####
# subset only early-non-Texprog
cd8t_cluster_nonTexprog <- subset(cd8t_cluster, Res0.6Names != "Naive" & Res0.6Names != "TexProg")

# find all DEGs without cutoffs for KO vs WT
Idents(cd8t_cluster_nonTexprog) <- "genotype"
degs <- FindMarkers(cd8t_cluster_nonTexprog, ident.1 = "HDAC1cKO",
                    ident.2 = "WT",
                    logfc.threshold = 0, min.pct = 0.01, min.diff.pct = 0,
                    only.pos = F)

# add gene names as column
degs = data.frame(gene=row.names(degs), degs)

# save as csv
write_csv2(degs, "../nonTexprogDEGs.csv")

# define genes of interest for plotting
mylabs <- c("Gzma", "Eomes", "Prf1", "Tox", "Gzmk", "Klre1", "Tnfrsf9", "Havrc2",
            "Hdac1", "Cx3cr1")

# plot DEGs as volcano plot and save as pdf
pdf("../VolcanoPlot_KOvsWT_earlynonTexprog_Res0.6.pdf", height = 12, width = 12)
myplot <- EnhancedVolcano(degs, lab = degs$gene,
                          x = "avg_log2FC", 
                          y = "p_val_adj",
                          selectLab = mylabs,
                          labSize = 6,
                          pointSize = 4,
                          FCcutoff = 0.25, 
                          pCutoff = 0.05,
                          title = "DEGs HDAC1cKO vs. WT in early-non-Texprog cells",
                          subtitle = "",
                          drawConnectors = T,
                          widthConnectors = 0.75,
                          typeConnectors = "open",
                          endsConnectors = "last",
                          arrowheads = F,
                          max.overlaps = Inf)
myplot <- myplot + ggplot2::coord_cartesian(xlim = c(-2,2))
print(myplot)
dev.off()

#### run GO term (BP) enrichment analysis on all clusters ####

library("dplyr")
library("clusterProfiler")
library("org.Mm.eg.db")
library("AnnotationHub")
library("readr")

#load marker genes for all clusters resolution 0.6
allMarkers <- read.csv2("../Allmarkers_default_Res.0.6.csv")

allMarkersSorted <- allMarkers %>% group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = T) # make sure genes are grouped by cluster and sorted by descending log2FC

df <- allMarkersSorted[,7:6] # take only gene and cluster columns
dfsample <- split(df$gene, df$cluster) # split data frame into one gene list per cluster
length(dfsample) # should be 8 for resolution 0.6

dfsample <- dfsample[c(1,7,8,5,4,6,2,3)] # reorder gene lists to order used before

# change gene symbols to entrez ids, some might not map and therefore be lost
for (i in 1:length(dfsample)) {
  dfsample[[i]] = bitr(dfsample[[i]], fromType = "SYMBOL", 
                       toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
}

# generate list of genes per cluster
genelist <- list()
for (i in 1:length(dfsample)) {
  genelist[i] <- dfsample[[i]][2]
  names(genelist)[i] <- names(dfsample)[i]
}

# compare clusters regarding their enriched GO terms (Biological Processes)
GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", 
                                OrgDb = "org.Mm.eg.db", ont = "BP",
                                pool = T)
write_csv2(GOclusterplot@compareClusterResult,
           "../GOBPenriched_DotplotData_ALLclustersRes0.6_AllSignMarkers.csv")

# take top 10 terms per cluster
top10 <- GOclusterplot@compareClusterResult %>% group_by(Cluster) %>%
  top_n(n = 10, wt = -p.adjust)
top10IDpadj <- top10[,c(2,7)]

# save top 10 terms per cluster with corresponding adjusted p-value
write_delim(top10IDpadj, "../GOBPenriched_DotPlotData_ALLclustersRes0.6_AllSignMarkers_Top10forRevigo.txt",
            delim = "\t", col_names = F)

# read in reduced GO terms after usage of webbased Revigo tool
revigoterms <- read_delim("../Revigo_BP_Table_ALLclustersRes0.6_AllSignMarkers_Top10.tsv")
revigoterms <- subset(revigoterms, Representative == "null")
revigonames <- revigoterms$Name

# plot enriched GO terms for all clusters which belong to the revigo-reduced categories
pdf("../GOBPenriched_DotPlot_ALLclustersRes0.6.pdf", 
    width = 16, height = 10)
enrichplot::dotplot(GOclusterplot, font.size = 8, showCategory = revigonames)
dev.off()
