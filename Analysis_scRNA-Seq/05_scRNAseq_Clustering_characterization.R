# title: cluster characterization
# author: Monika Waldherr
# input: Seurat object saved after clustering with script 04_scRNAseq_Clustering.R
#         CSV file containing published cluster signatures from Daniel et al. (2021)

#### load libraries ####
library(SeuratObject)
library(Seurat)
library(reshape)
library(RColorBrewer)

#### load Seurat object ####
cd8t_cluster <- readRDS("cd8t_cluster.rds")

# Assign identity of clusters and assay, normalize, scale and 
# visually verify that anticipated clustering was set
Idents(object = cd8t_cluster) <- "integrated_snn_res.0.6"
DefaultAssay(cd8t_cluster) <- "RNA"
cd8t_cluster <- NormalizeData(cd8t_cluster)
cd8t_cluster <- ScaleData(cd8t_cluster)
DimPlot(cd8t_cluster, label = T)

#### calculation of enrichment scores for published gene sets from Daniel et al. ####
# load published signatures
tcellSigs <- read.csv("../Daniel2021_TcellSignatures.csv", sep = ";")

# get top 100 genes for each cluster sorted by average log2FC
tcellSigs100 <- tcellSigs %>% group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC)

# get the names of the different clusters as factors in a vector
clusterNames <- levels(as.factor(tcellSigs100$cluster))

# remove clusters which are not of interest
clusterNames <- clusterNames[-c(7)]

# reorder cluster names
clusterNames <- clusterNames[c(10,7,3,
                               4,6,5,
                               8,1,2,9)]

# initialize list for all gene names
geneLists <- list()

# save gene names for every cluster sorted by log2FC in descending order
for (i in clusterNames) {
  geneLists[[i]] <- subset(tcellSigs100, cluster == i) %>% 
    arrange(desc(avg_log2FC)) %>%
    pull(gene)
}

# add module score for each tcell subset (cluster) to Seurat object
for (i in clusterNames) {
  cd8t_cluster <- AddModuleScore(cd8t_cluster,
                                 features = list(geneLists[[i]]),
                                 name = as.character(i))
}
colnames(cd8t_cluster@meta.data)[c(37:46)] <- clusterNames

# print feature plot showing average expression of published cluster signatures in data
pdf("FeaturePlot_ClusterComparedtoDanieletal_Res.0.6.pdf", height = 8, width = 10)
FeaturePlot(cd8t_cluster,
            features = clusterNames, min.cutoff = "q1", max.cutoff = "q99",
            cols = rev(brewer.pal(n = 7, name = "RdBu")),
            order = T, pt.size = 0.1)
dev.off()

# print dotplot showing average expression of published cluster signatures in data;
# order of identities was changed later for publication according to found similarities with published gene sets
# and identified marker genes, but here are still ordered and named as in original clustering 
pdf("DotPlot_ClusterComparedtoDanieletal_Res.0.6.pdf", height = 8, width = 10)
DotPlot(cd8t_cluster, features = clusterNames, cols = "RdBu", col.min = -1.7, col.max = 1.7) + 
  RotatedAxis() + xlab("Signatures")
dev.off()

#### find marker genes for all clusters ####
allmarkers <- FindAllMarkers(object = cd8t_cluster,
                             only.pos = T,
                             min.pct = 0.1,
                             logfc.threshold = 0.25)

allmarkers <- subset(allmarkers, p_val_adj < 0.05)

write_csv2(allmarkers, "../Allmarkers_default_Res.0.6.csv")

#### rename and reorder clusters of resolution 0.6 ####
newclusterids <- c("TexInt", "TexEarly",
                   "TexCyt", "TexCX3CR1",
                   "TexProg", "Naive",
                   "TexExh", "TexProl") # new IDs in order of currently set identifiers
names(newclusterids) <- levels(cd8t_cluster) # name vector of new IDs using names of current identifiers
cd8t_cluster <- RenameIdents(cd8t_cluster, newclusterids) # rename identifiers to new IDs

reorderedlevels <- c("Naive", "TexProg",
                     "TexProl", "TexExh",
                     "TexEarly",
                     "TexInt", "TexCX3CR1",
                     "TexCyt") # vector with reordered cluster names
levels(cd8t_cluster) <- reorderedlevels # reorder levels in Seurat object

# define colors for future plotting
mypal <- c("seashell3", "#0084B3", "#00B1B5", "palevioletred1", "#A71B4B",
           "#F7B347", "khaki3", "darkseagreen")

# plot overall UMAP
pdf("../UMAP_ALL_Res0.6.pdf", width = 10, height = 8)
DimPlot(cd8t_cluster, label = F, pt.size = 0.8,
        cols = mypal)
dev.off()

# add genotype column to metadata for more convenient separation
metadata <- cd8t_cluster@meta.data
metadata <- metadata %>% mutate(genotype = ifelse(grepl("WT", 
                                                        HTO_classification),
                                                  "1WT", "2HDAC1cKO"))
metadata$genotype <- as.factor(metadata$genotype)
levels(metadata$genotype) <- c("WT", "HDAC1cKO")

# add column with new cluster names of resolution 0.6 clustering to metadata
metadata$Res0.6Names <- Idents(cd8t_cluster)

# save updated metadata back into Seurat object
cd8t_cluster@meta.data <- metadata

# plot UMAP split by genotype
pdf("../UMAP_ALLsplit_Res0.6.pdf", width = 15, height = 8)
DimPlot(cd8t_cluster, label = F, pt.size = 0.8,
        cols = mypal, split.by = "genotype")
dev.off()

# save Seurat object with updated Idents and metadata
saveRDS(cd8t_cluster, "../cd8t_cluster_ownIdentsRes0.6.rds")

#### plot cell number and frequencies of cluster as barplot ####
pdf("../CellNumberBarplot_Res0.6.pdf", width = 12, height = 8)
plt1 <- table(Idents(cd8t_cluster), cd8t_cluster$genotype)
plt1 <- as.data.frame(plt1)
ggplot(plt1, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.5) +
  xlab("") +
  ylab("Cell Number") +
  scale_fill_manual(values = mypal[c(2,8)]) +
  theme(legend.title = element_blank())
dev.off()

pdf("../FrequencyBarplot_Res0.6.pdf", width = 12, height = 8)
plt1 <- table(Idents(cd8t_cluster), cd8t_cluster$genotype)
plt1 <- as.data.frame(plt1)
ggplot(plt1, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("") +
  ylab("Frequency") +
  scale_fill_manual(values = mypal[1:8]) +
  theme(legend.title = element_blank())
dev.off()

#### find unique marker genes and plot heatmap ####
mymarkers <- read_csv2("../Allmarkers_default_Res.0.6.csv")
mymarkers$cluster <- factor(mymarkers$cluster, levels = reorderedlevels)
mymarkers <- mymarkers[order(mymarkers$cluster),]
mymarkers %>%
  group_by(gene) %>%
  dplyr::filter(n() == 1) %>%
  ungroup() -> mymarkers_unique
write_csv2(mymarkers_unique, "../Alluniquemarkers_Res.0.6.csv")

cd8t_cluster_avg <- AverageExpression(cd8t_cluster, return.seurat = T)
pdf("../Heatmap_Alluniquemarkers_Res.0.6.pdf", width = 6, height = 15)
DoHeatmap(cd8t_cluster_avg,
          features = mymarkers_unique$gene, 
          group.bar = T, 
          group.colors = mypal,
          label = T,
          size = 3,
          draw.lines = F,
          raster = F) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  guides(color = "none") +
  theme(axis.text.y=element_blank())
dev.off()

#### plot expression of genes of interest as feature plots ####
genesofinterest <- c("Tcf7", "Slamf6", "Havcr2", "Cx3cr1", "Tox", "Gzma")
pdf("../FeaturePlots_GenesofInterest.pdf")
FeaturePlot(cd8t_cluster, features = genesofinterest, max.cutoff = 6.5, min.cutoff = 1,
            pt.size = 0.1, order = T)
dev.off()









