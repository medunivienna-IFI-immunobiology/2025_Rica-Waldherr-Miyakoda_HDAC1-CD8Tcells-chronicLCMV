# title: subset consensus region bed files of ATAC-Seq data to only keep the ones of interest for TF motif analysis
# author: Monika Waldherr
# input: conensus_regions.bed file containing all called consensus peaks, CSV files containing names of consensus regions of interest

#### load libraries ####
library(ggplot2)

# read in bed file containing all consensus regions
allcons <- read.delim("consensus_regions.bed", header = F)
colnames(allcons) <- c("chr", "start", "end", "consensus")

# read in consensus regions of interest Tprog
DARsTprog <- read.delim("dfATAC_Tprog_significantDARs.csv", sep = ";")

# split DARs in WT (neg log2(FC)) and KO (pos log2(FC))
DARsTprogWT <- subset(DARsTprog, DARsTprog$log2FoldChange__group_KO_Tprog_vs_WT_Tprog < 0)
DARsTprogKO <- subset(DARsTprog, DARsTprog$log2FoldChange__group_KO_Tprog_vs_WT_Tprog > 0)

# only keep consensus regions in allcons which are found in DARs
WTcons <- allcons[allcons$consensus %in% DARsTprogWT$consensus,]
KOcons <- allcons[allcons$consensus %in% DARsTprogKO$consensus,]

# save as textfile with .bed ending for homer processing
write.table(WTcons, "TprogWTopenCons.bed", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(KOcons, "TprogKOopenCons.bed", col.names = F, row.names = F, sep = "\t", quote = F)

# read in consensus regions of interest Tterm
DARsTterm <- read.delim("dfATAC_Tterm_significantDARs.csv", sep = ";")

# split DARs in WT (neg log2(FC)) and KO (pos log2(FC))
DARsTtermWT <- subset(DARsTterm, DARsTterm$log2FoldChange__group_KO_Tterm_vs_WT_Tterm < 0)
DARsTtermKO <- subset(DARsTterm, DARsTterm$log2FoldChange__group_KO_Tterm_vs_WT_Tterm > 0)

# only keep consensus regions in allcons which are found in DARs
WTcons <- allcons[allcons$consensus %in% DARsTtermWT$consensus,]
KOcons <- allcons[allcons$consensus %in% DARsTtermKO$consensus,]

# save as textfile with .bed ending for homer processing
write.table(WTcons, "TtermWTopenCons.bed", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(KOcons, "TtermKOopenCons.bed", col.names = F, row.names = F, sep = "\t", quote = F)
