# Publication: "HDAC1 controls the generation and maintenance of effector-like CD8+ T cells during chronic viral infection"
This README gives an overview of all R scripts used for analysis and visualization of scRNA-Seq and ATAC-Seq data in the study (DOI: XXX). They are available from the respective subfolders in this repository. The raw and processed data used as input will be available in the GEO database.

## scRNA-Seq Analysis
Large parts of the scripts are based on the tutorials provided by the Satija lab (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) and the tutorial from HBC training (https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html).
#### 01_scRNAseq_QC_cells.R
R script to run the cell-based quality control of scRNA-Seq data using the Seurat package. As input the path to the directory containing Cell Ranger output (raw count matrix with corresponding feature and barcode lists) has to be specified. The script will read in the data, convert it to a Seurat object, demuliplex and assign the HTO identities, remove doublets and negatives, remove unwanted cell types, add sample identifiers, rename some metadata columns, and visualize the QC metrics of the raw singlet data. Finally, the updated Seurat object will be saved.
## ATAC-Seq Analysis
