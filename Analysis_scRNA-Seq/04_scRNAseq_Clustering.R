# title: clustering
# author: Monika Waldherr
# input: Seurat object saved after cell cycle variation regression and replicate integration
# with script 03_scRNAseq_QC_regression_integration.R

#### load libraries ####
library(SeuratObject)
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

#### load Seurat object ####
cd8t_cluster <- readRDS("cd8t_integrated.rds")

#### define number of principal components (PCs) to use ####
# determine percent of variation associated with each PC
pct <- cd8t_cluster[["pca"]]@stdev / sum(cd8t_cluster[["pca"]]@stdev) * 100

# calculate cumulative percents for each PC
cumu <- cumsum(pct)

# determine which PC exhibits cumulative percent greater than 90% and 
# % variation associated with the PC is less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# determine the difference between variation of PC and subsequent PC and
# define the last point where change of % of variation is more than 0.1%
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), 
            decreasing = T)[1] + 1

# minimum of the two calculations will be the number of PCs to use
pcs <- min(co1, co2)

# plot the elbow plot with the calculated cut-off as a vertical line
pdf("../EllbowPCstoUse.pdf", width = 8, height = 8)
ElbowPlot(object = cd8t_cluster, 
          ndims = 40) + geom_vline(xintercept = 14)
dev.off()

#### run clustering ####
# determine the K-nearest neighbor graph
cd8t_cluster <- FindNeighbors(object = cd8t_cluster, 
                              dims = 1:pcs)

# determine the clusters for various resolutions (the higher the resolution,
# the more clusters will be calculated)                            
cd8t_cluster <- FindClusters(object = cd8t_cluster,
                             resolution = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                                            1.1, 1.2, 1.3, 1.4))

# explore the different resolutions to find the most useful one

# Assign identity of clusters
Idents(object = cd8t_cluster) <- "integrated_snn_res.0.6" # resolution of 0.6 was used in publication
DimPlot(cd8t_cluster)

#### save Seurat object with clustering information for several resolutions ####
saveRDS(cd8t_cluster, "cd8t_cluster.rds")

#### from here go on to the cluster characterization script ####

