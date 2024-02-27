# Publication: "HDAC1 controls the generation and maintenance of effector-like CD8+ T cells during chronic viral infection"
This README gives an overview of all R scripts used for analysis and visualization of scRNA-Seq and ATAC-Seq data in the study (DOI: XXX). They are available from the respective subfolders in this repository. The raw and processed data used as input will be available in the GEO database. Additionally, relevant metadata generated within the study can be found in the corresponding subfolder.

## scRNA-Seq Analysis
Raw sequencing output was processed to count-feature-barcode matrices as described in the methods section of the publication. Here only the subsequently used R scripts are described and made available. They are split in logically connected units to make versatile usage as easy as possible.

Large parts of the scripts are based on the tutorials provided by the Satija lab and the tutorial from HBC training:
* https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
* https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html
#### 01_scRNAseq_QC_cells.R
R script to run the cell-based quality control of scRNA-Seq data using the Seurat package. As input the path to the directory containing Cell Ranger output (raw count matrix with corresponding feature and barcode lists) has to be specified. The script will read in the data, convert it to a Seurat object, demuliplex and assign the HTO identities, remove doublets and negatives, remove unwanted cell types, add sample identifiers, rename some metadata columns, and visualize the QC metrics of the raw singlet data. Finally, the updated Seurat object will be saved.
#### 02_scRNAseq_QC_genes.R
R script to run the gene-based quality control of scRNA-Seq data using the Seurat package. As input the Seurat object resulting from script 01_scRNAseq_QC_cells.R is required. The data will be filtered by the number of detected UMIs, genes, complexity and % of mitochondrial genes. The same QC metrics as for the raw data will be visualized for the filtered data. The filtered Seurat object will be saved.
#### 03_scRNAseq_QC_regression_integration.R
R script to regress out unwanted variation from cell cycle phases and to integrate the replicates. PCA and UMAP before and after regression will be visualized. As input the Seurat object resulting from script 02_scRNAseq_QC_genes.R is required. The final Seurat object will again be saved.
#### 04_scRNAseq_Clustering.R
R script to determine the number of principal components to use and to perform the clustering with different resolutions. Requires the Seurat object generated in the script 03_scRNAseq_QC_regression_integration.R as input. All clusterings will be saved within the resulting Seurat object.

## ATAC-Seq Analysis
Raw sequencing output was processed to a consensus sequence count matrix as described in the methods section of the publication. Here only the subsequently used R scripts are described and made available.
