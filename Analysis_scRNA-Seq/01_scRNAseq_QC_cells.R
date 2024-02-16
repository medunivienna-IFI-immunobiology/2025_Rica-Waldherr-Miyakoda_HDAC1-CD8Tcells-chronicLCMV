# title: Cell-based quality control of scRNA-Seq data
# author: Monika Waldherr
# input: path to directory containing Cell Ranger output (raw count matrix with corresponding feature and barcode lists)

#### load libraries ####
library(SeuratObject)
library(Seurat)
library(tidyverse)
library(celldex)
library(SingleR)

#### read in data and explore ####
mywd <- getwd() # set working directory to the directory containing the input data before assigning mywd

# create a Seurat object for each replicate
for (file in c("rep1_raw_feature_bc_matrix", 
               "rep2_raw_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0(mywd, "/", file))
  
  seurat_obj <- CreateSeuratObject(counts = seurat_data[["Gene Expression"]], # read in gene expression/count data
                                   min.features = 50, # set to 50 to get enough negative cells for decent demultiplexing
                                   project = file)
  seurat_obj[["HTO"]] <- CreateAssayObject(seurat_data[["Antibody Capture"]][, # read in hashtag/genotype information
                                                                             colnames(x = seurat_obj)], 
                                           min.cells = 1) # set to 1 to remove hashtags/antibodies that didn't work
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, 
                                                     pattern = "^mt-") # add percentage of mitochondrial genes, 
                                                                        # in mouse they start with "mt-"
  seurat_obj[["mitoRatio"]] <- seurat_obj@meta.data$percent.mt / 100 # add the mitochondrial content also as ratio
  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) /
    log10(seurat_obj$nCount_RNA) # add novelty score as additional QC metric to evaluate complexity
  assign(file, seurat_obj)
}

# plot general information of data
rep1_raw_feature_bc_matrix # output: 32289 features across 11813 samples within 2 assays

rep2_raw_feature_bc_matrix # output: 32289 features across 18363 samples within 2 assays

# merge the two Seurat objects into one without integrating them
merged_seurat <- merge(x = rep1_raw_feature_bc_matrix,
                       y = rep2_raw_feature_bc_matrix,
                       add.cell.id = c("rep1", "rep2")) # add replicate information

#### demultiplex samples and assign HTO identities ####
cd8t_norm <- NormalizeData(merged_seurat) # normalize RNA data with log normalization
cd8t_norm <- FindVariableFeatures(cd8t_norm,
                                  selection.method = "vst") # find and then scale variable features
cd8t_norm <- ScaleData(cd8t_norm,
                       features = VariableFeatures(cd8t_norm))
cd8t_norm <- NormalizeData(cd8t_norm, assay = "HTO",
                           normalization.method = "CLR") # normalize HTO data with centered-log ratio transformation

cd8t_demux <- HTODemux(cd8t_norm, assay = "HTO", positive.quantile = 0.99) # demultiplex with a threshold of 0.99
# output:
# Cutoff for HTO-WT-naive : 5 reads
# Cutoff for HTO-HDAC1cKO-naive : 5 reads
# Cutoff for HTO-WT-gp33pos : 45 reads
# Cutoff for HTO-HDAC1cKO-gp33pos : 33 reads

# plot general information of Seurat object
cd8t_demux
# output:
# 32289 features across 30176 samples within 2 assays 
# Active assay: RNA (32285 features, 2000 variable features)

# plot the number of Doublets, Negatives and Singlets
table(cd8t_demux$HTO_classification.global) # 2398 Doublets, 7108 Negatives, 20670 Singlets

# plot demux results
pdf("../DemuxRidgePlot.pdf", width = 11, height = 6)
Idents(cd8t_demux) <- "HTO_maxID"
RidgePlot(cd8t_demux, assay = "HTO", 
          features = rownames(cd8t_demux[["HTO"]])[1:4], ncol = 2)
dev.off()

pdf("../DemuxClassificationGlobalPlot.pdf")
Idents(cd8t_demux) <- "HTO_classification.global"
VlnPlot(cd8t_demux, features = "nCount_RNA", pt.size = 0.1, log = T)
dev.off()

pdf("../DemuxHeatmapPlot.pdf", width = 6, height = 3)
HTOHeatmap(cd8t_demux, assay = "HTO", ncells = 2000, raster = F)
dev.off()

#### remove doublets and negatives ####
Idents(cd8t_demux) <- "HTO_classification.global" # set idents to global HTO classification
cd8t_demux_s <- subset(cd8t_demux, idents = "Singlet") # only keep cells classified as Singlet (= 20670)
summary(cd8t_demux_s@meta.data) # plot summary of metadata

# redo finding variable features and scaling for cells of interest
cd8t_demux_s <- FindVariableFeatures(cd8t_demux_s,
                                     selection.method = "vst")
cd8t_demux_s <- ScaleData(cd8t_demux_s,
                          features = VariableFeatures(cd8t_demux_s))

# plot general information of Seurat object only containing singlets
cd8t_demux_s
# output:
# 32289 features across 20670 samples within 2 assays 
# Active assay: RNA (32285 features, 2000 variable features)

#### add cell type annotation using SingleR ####
ref <- celldex::MouseRNAseqData() # get reference
pred <- SingleR(test = cd8t_demux_s@assays$RNA@data, ref = ref, 
                labels = ref$label.fine) # predict cell type (takes some time)
pred.labels <- data.frame(pred[,"labels", drop = F]) # only get labels from predicition
colnames(pred.labels) <- "SingleR_labels" # rename column
cd8t_demux_s <- AddMetaData(cd8t_demux_s, pred.labels) # add cell type labels as metadata to Seurat object

write.csv(table(cd8t_demux_s@meta.data$SingleR_labels), 
          file = "../TableCellTypesPredictedRaw.csv") # save number of cells falling into different cell types

# only keep Tcells and NKcells (as they have a very similar gene profile)
cd8t_demux_s <- subset(cd8t_demux_s, 
                       subset = SingleR_labels == "T cells" | 
                         SingleR_labels == "NK cells")

# plot general information of Seurat object
cd8t_demux_s # 32289 features across 20171 samples within 2 assays

#### add sample identifier and rename metadata columns ####
metadata <- cd8t_demux_s@meta.data
metadata$cells <- rownames(metadata) # add rownames as cells column

metadata$sample <- NA # add column for sample identifier
metadata$sample[which(str_detect(metadata$cells, "^rep1_"))] <- "rep1"
metadata$sample[which(str_detect(metadata$cells, "^rep2_"))] <- "rep2"

metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA) # rename columns to seq_folder, nUMI and nGene

cd8t_demux_s@meta.data <- metadata # add changed metadata back to Seurat object

#### visualize QC metrics of raw singlet data ####
# as violin plot
pdf("../RawViolinMetadataQC.pdf", width = 16, height = 8)
VlnPlot(cd8t_demux_s, 
        features = c("nGene", "nUMI", "percent.mt", "log10GenesPerUMI"),
        ncol = 4)
dev.off()

metadata <- cd8t_demux_s@meta.data # for easier usage fetch metadata from object

# visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Number of Cells per Replicate") -> p.nCells

pdf("../RawNCellsReps.pdf")
p.nCells
dev.off()

# visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500) -> p.nUMI

# visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300) -> p.nGene

# visualize the overall complexity of the gene expression by visualizing the 
# genes detected per UMI (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) -> p.complexity

# visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 5) -> p.mito

pdf("../RawQCmetrics.pdf", width = 10, height = 8)
require(gridExtra)
grid.arrange(p.nUMI, p.nGene, p.complexity, p.mito, ncol = 2)
dev.off()

# visualize the correlation between genes detected and number of UMIs and 
# determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 300) +
  facet_wrap(~sample) -> p.UMIGenesMito

pdf("../RawQC_UMIGenesMitoCorr.pdf", width = 10, height = 8)
p.UMIGenesMito
dev.off()

#### save Seurat object with still raw data (but filtered for singlets and T+NK cells) ####
saveRDS(cd8t_demux_s, file = "cd8t_demux_s.rds") 

#### from here go on to gene-based QC script ####











