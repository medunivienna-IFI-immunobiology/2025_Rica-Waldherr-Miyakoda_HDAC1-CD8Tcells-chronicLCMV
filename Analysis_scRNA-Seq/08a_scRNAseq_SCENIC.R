# title: SCENIC (Single Cell rEgulatory Network Inference and Clustering) - GENIE
# author: Monika Waldherr
# input: filtered expression matrix and options RDS from SCENIC preparation

# load libraries
library(SCENIC)

# load data
motifAnnotations_mgi <- read.table("int/motifAnno.tab", header = T)
exprMat_filtered <- as.matrix.data.frame(read.table("int/exprMat_filtered.tab", header = T), rownames.force = T)
scenicOptions <- readRDS("int/scenicOptions.Rds")

## Run GENIE3
set.seed(400)
runGenie3(exprMat_filtered, scenicOptions, nParts = 30)

saveRDS(scenicOptions, "int/scenicOptions_afterGenie.Rds")
