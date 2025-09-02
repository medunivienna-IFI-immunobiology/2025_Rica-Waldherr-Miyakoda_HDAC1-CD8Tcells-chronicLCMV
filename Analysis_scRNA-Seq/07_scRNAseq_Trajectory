# title: trajectory analysis with Monocle3
# author: Monika Waldherr
# input: Seurat object

# install dependencies
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
# install monocle 3
devtools::install_github('cole-trapnell-lab/monocle3')

remotes::install_github('satijalab/seurat-wrappers@02754e1517e51a5ee093e8c2f247b49f85ee2d40')

# load libraries
library("SeuratObject")
library("Seurat")
library("SeuratWrappers")
library("monocle3")
library("ggplot2")

# load Seurat object
cd8t_cluster <- readRDS("./cd8t_cluster_ownIdentsRes0.6.rds")

# run monocle on Seurat object
monocle_object <- as.cell_data_set(cd8t_cluster)
rm(cd8t_cluster)
monocle_object <- cluster_cells(cds = monocle_object, reduction_method = "UMAP")
monocle_object <- learn_graph(monocle_object, use_partition = TRUE)
monocle_object <- order_cells(monocle_object, reduction_method = "UMAP") # defines start nodes for pseudotime
saveRDS(monocle_object, "./A_10_scRNAseq_Revision_Trajectory/Monocle3object.Rds")

library(wesanderson)

pdf("./Plot_Trajectory_Monocle3_Pseudotime.pdf", width = 10, height = 8)
plot_cells(monocle_object,
           color_cells_by = "pseudotime",
           cell_size = 1.5,
           trajectory_graph_color = "grey4",
           show_trajectory_graph = TRUE,
           alpha = 0.5,
           trajectory_graph_segment_size = 1.0,
           label_branch_points = T,
           label_roots = F,
           label_leaves = F) +
  scale_color_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous"))
dev.off()
