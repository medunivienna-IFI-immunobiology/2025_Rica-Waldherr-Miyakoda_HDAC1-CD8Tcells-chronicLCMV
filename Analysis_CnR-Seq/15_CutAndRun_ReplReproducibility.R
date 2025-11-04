# title: revision - Cut&Run analysis - Check Replicate Reproducibility
# author: Monika Waldherr
# input: BED files where fragment counts per bin are indicated

# code adapted from https://yezhengstat.github.io/CUTTag_tutorial/#33_Remove_duplicates_[optionalrequired]
library(dplyr)
library(corrplot)

# get seperate correlation matrices for celltypes and genotypes
setwd("./")
cellgenotypes <- c("TprogWT", "TprogKO", "nonTprogWT", "nonTprogKO")

for (type in cellgenotypes){
  sampleList <- list.files(paste0("./fragmentCounts/", type), pattern = "*fragmentsCount.bin500.bed")
  #histList <- c("H3K27ac", "H3K27me3", "HDAC1", "RUNX3")
  histInfoList <- gsub("CnR-20241114-CD8--", "", sampleList)
  histInfoList <- gsub("_RmBlacklist.mapped_sorted.fragmentsCount.bin500.bed", "", histInfoList)
  
  reprod = c()
  fragCount = NULL
  for(i in 1:length(sampleList)){
    histInfo = histInfoList[i]
    if(is.null(fragCount)){
      
      fragCount = read.table(paste0("./fragmentCounts/", type, "/", sampleList[i]), header = FALSE) 
      colnames(fragCount) = c("chrom", "bin", histInfo)
      
    }else{
      
      fragCountTmp = read.table(paste0("./fragmentCounts/", type, "/", sampleList[i]), header = FALSE)
      colnames(fragCountTmp) = c("chrom", "bin", histInfo)
      fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
      
    }
  }
  
  M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs")
  
  pdf(paste0("./ReplicateReproducibility_", type, ".pdf"), width = 10, height = 10)
  corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="alphabet", addrect = 5, rect.col = "black",
           rect.lwd = 3, cl.pos = "b", tl.col = "black", tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2,
           number.cex = 1, col = colorRampPalette(c("midnightblue","white","darkred"))(100))
  dev.off()
}

