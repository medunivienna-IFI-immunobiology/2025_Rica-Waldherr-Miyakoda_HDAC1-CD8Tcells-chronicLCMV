# title: Gene-based quality control of scRNA-Seq data
# author: Monika Waldherr
# input: Seurat object saved after cell-based QC with script 01_scRNAseq_QC_cells.R

#### load libraries ####
library(SeuratObject)
library(Seurat)
library(tidyverse)
library(celldex)
library(SingleR)
library(RCurl)
library(cowplot)

#### load Seurat object ####
cd8t_demux_s <- readRDS("cd8t_demux_s.rds")

#### filter number of detected UMI, genes, complexity and % mitochondrial genes & visualize ####
filtered_cd8t <- subset(x = cd8t_demux_s, 
                        subset= (nUMI >= 500) & 
                          (nGene >= 300) & 
                          (nGene <= 3000) &
                          (log10GenesPerUMI > 0.80) & 
                          (percent.mt <= 5))

# plot general information of Seurat object
filtered_cd8t # 32289 features across 16393 samples within 2 assay

# visualize QC metrics as violin plots
pdf("../FiltViolinMetadataQC.pdf", width = 16, height = 8)
VlnPlot(filtered_cd8t, 
        features = c("nGene", "nUMI", "percent.mt", "log10GenesPerUMI"),
        ncol = 4)
dev.off()

# get filtered metadata
metadata_clean <- filtered_cd8t@meta.data

# visualize the number of cell counts per sample
metadata_clean %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Number of Cells per Replicate") -> p.nCells

pdf("../FiltNCellsReps.pdf")
p.nCells
dev.off()

# visualize the number UMIs/transcripts per cell
metadata_clean %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500) -> pf.nUMI

# visualize the distribution of genes detected per cell via histogram
metadata_clean %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300) -> pf.nGene

# visualize the overall complexity of the gene expression by visualizing the 
# genes detected per UMI (novelty score)
metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) -> pf.complexity

# visualize the distribution of mitochondrial gene expression detected per cell
metadata_clean %>% 
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 5) -> pf.mito

pdf("../FiltQCmetrics.pdf", width = 10, height = 8)
require(gridExtra)
grid.arrange(pf.nUMI, pf.nGene, pf.complexity, pf.mito, ncol = 2)
dev.off()

# visualize the correlation between genes detected and number of UMIs and 
# determine whether strong presence of cells with low numbers of genes/UMIs
metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample) -> pf.UMIGenesMito

pdf("../FiltQC_UMIGenesMitoCorr.pdf", width = 10, height = 8)
pf.UMIGenesMito
dev.off()

#### save filtered Seurat object ####
saveRDS(filtered_cd8t, file = "filtered_cd8t.rds")




