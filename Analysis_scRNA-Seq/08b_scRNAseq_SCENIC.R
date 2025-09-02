# title: SCENIC (Single Cell rEgulatory Network Inference and Clustering) - three steps
# author: Monika Waldherr
# input: expression matrix and options RDS from SCENIC preparation after GENIE. annotation table

# load libraries
library(SCENIC)

# load scenicOptions_afterGenie.Rds
scenicOptions <- readRDS("int/scenicOptions_afterGenie.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

# reload the expression matrix and get a log transformed one too
exprMat <- as.matrix.data.frame(read.table("int/exprMat.tab", header = T), rownames.force = T)
exprMat_log <- log2(exprMat + 1) # Optional: log expression (for TF expression plot, it does not affect any other calculation)

# reload annotations
motifAnnotations_mgi <- read.table("int/motifAnno.tab", header = T)

# run SCENIC in three steps
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) # needs motifAnnotations_mgi, reqires a lot of RAM -> run on server
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
