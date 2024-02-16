# title: regression of unwanted variation from cell cycle and integration of replicates
# author: Monika Waldherr
# input: Seurat object saved after gene-based QC with script 02_scRNAseq_QC_genes.R

#### load libraries ####
library(SeuratObject)
library(Seurat)
library(tidyverse)
library(celldex)
library(SingleR)
library(RCurl)
library(cowplot)

#### load Seurat object ####
filtered_cd8t <- readRDS("filtered_cd8t.rds")

#### regress unwanted variation from cell cycle ####

# start by normalizing
filtered_cd8t <- NormalizeData(filtered_cd8t)

# assign cell cycle scores and visualize phase separation
s.genes <- str_to_title(cc.genes.updated.2019$s.genes) # load s phase genes
g2m.genes <- str_to_title(cc.genes.updated.2019$g2m.genes) # load g2m phase genes

filtered_cd8t <- CellCycleScoring(filtered_cd8t,
                                  g2m.features = g2m.genes,
                                  s.features = s.genes)

# find variable features and scale data again, just to be sure
filtered_cd8t <- FindVariableFeatures(filtered_cd8t, 
                                      selection.method = "vst",
                                      nfeatures = 2000, 
                                      verbose = FALSE)
filtered_cd8t <- ScaleData(filtered_cd8t)

# run PCA and visualize impact of cell cycle phases and replicates
# before regression or integration respectively
filtered_cd8t <- RunPCA(filtered_cd8t)
filtered_cd8t <- RunUMAP(filtered_cd8t, dims = 1:40,
                         reduction = "pca")

pdf("../CellCyclePCAbeforeReg.pdf", width = 12, height = 8)
DimPlot(filtered_cd8t,
        reduction = "pca",
        group.by= "Phase",
        split.by = "sample")
dev.off()

pdf("../ReplicatesUMAPbeforeInt.pdf", width = 12, height = 8)
DimPlot(filtered_cd8t, reduction = "umap", group.by = "sample")
dev.off()

# SCTransform to regress out cell cycle variation in each replicate
cd8t_phase <- filtered_cd8t
split_cd8t_phase <- SplitObject(cd8t_phase, split.by = "sample")
split_cd8t_phase <- split_cd8t_phase[c("rep1", "rep2")]

options(future.globals.maxSize = 4000 * 1024^2)

for (i in 1:length(split_cd8t_phase)) {
  split_cd8t_phase[[i]] <- SCTransform(split_cd8t_phase[[i]],
                                       vars.to.regress = c("S.Score",
                                                           "G2M.Score"))
  split_cd8t_phase[[i]] <- RunPCA(split_cd8t_phase[[i]],
                                  features = c(s.genes, g2m.genes),
                                  verbose = F)
}

# run PCA on both replicates and visualize impact of cell cycle genes after regression
split_cd8t_rep1 <- RunPCA(split_cd8t_phase$rep1)
split_cd8t_rep2 <- RunPCA(split_cd8t_phase$rep2)

DimPlot(split_cd8t_rep1,
        reduction = "pca",
        group.by= "Phase",
        split.by = "sample") -> p.rep1
DimPlot(split_cd8t_rep2,
        reduction = "pca",
        group.by= "Phase",
        split.by = "sample") -> p.rep2

pdf("../CellCyclePCAafterReg.pdf", width = 12, height = 8)
p.rep1 | p.rep2
dev.off()

#### integration of replicates ####

# select 3000 genes to be used for integration
integ_features <- SelectIntegrationFeatures(object.list = split_cd8t_phase,
                                            nfeatures = 3000)
# define selected integration genes as anchor genes
split_cd8t_phase <- PrepSCTIntegration(object.list = split_cd8t_phase,
                                       anchor.features = integ_features)
# find pairwise integration anchors (takes some time)
integ_anchors <- FindIntegrationAnchors(object.list = split_cd8t_phase,
                                        normalization.method = "SCT",
                                        anchor.features = integ_features)
# output: found 19430 anchors, retained 12281 after filtering

# integrate replicates along defined anchors
cd8t_integrated <- IntegrateData(anchorset = integ_anchors,
                                 normalization.method = "SCT")

cd8t_integrated <- RunPCA(cd8t_integrated)
PCAPlot(cd8t_integrated, split.by = "sample", group.by = "sample")
cd8t_integrated <- RunUMAP(cd8t_integrated, dims = 1:40,
                           reduction = "pca")
pdf("../ReplicatesUMAPafterInt.pdf", width = 8, height = 6)
DimPlot(cd8t_integrated, group.by = "sample")
dev.off()

#### save Seurat object were data is filtered, cell cycle variation regressed out, and replicates are integrated ####
saveRDS(cd8t_integrated, "cd8t_integrated.rds")

#### from here go on to clustering script ####