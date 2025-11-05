# title: Cut&Run analysis - differential peak analysis using DiffBind
# author: Monika Waldherr
# input: RUNX3 samplesheet with locations of reads and called peaks

# install and load DiffBind
# BiocManager::install("DiffBind")
library("DiffBind")

setwd("./CutAndRun")
mysamples <- read.csv("CnR_DiffBind_samplesheet_RUNX3.csv")

target <- dba(sampleSheet = mysamples)
target # plot metadata associated with created dba object

#target.count <- dba.count(target)

# set summits = T to use peaks as they were called
target.count <- dba.count(target, summits = TRUE)

dba.plotHeatmap(target.count, margin = 20)
mysamples[,DBA_CONDITION]
dba.plotPCA(target.count, attributes=DBA_CONDITION, label=DBA_ID)

target.count <- dba.contrast(target.count, categories = DBA_CONDITION)
target.count

target.analysed <- dba.analyze(target.count, method = DBA_ALL_METHODS)
dba.show(target.analysed, bContrasts = T)

dba.plotVenn(target.analysed, contrast = 1, method=DBA_ALL_METHODS)
dba.plotMA(target.analysed, method=DBA_DESEQ2, bNormalized = FALSE, yrange = c(-5, 5))
dba.plotMA(target.analysed, method=DBA_DESEQ2, bNormalized = TRUE, yrange = c(-5, 5))
dba.plotMA(target.analysed, method=DBA_EDGER, bNormalized = TRUE, yrange = c(-5, 5)) # EDGER normalization actually does nothing

# use DESEQ2 results (default) for further steps

dba.plotVolcano(target.analysed, method=DBA_DESEQ2)
#dba.plotVolcano(target.analysed, method = DBA_EDGER)
dba.plotPCA(target.analysed, contrast = 1, method=DBA_DESEQ2)

report <- dba.report(target.analysed, DataType = DBA_DATA_FRAME, method=DBA_DESEQ2)
report$PeakID <- rownames(report)
report.df <- as.data.frame(report[,c(10,1:9)])
#write.table(report.df, "RUNX3_DiffBind_DESeq2_report.csv", sep="\t", quote=F, row.names=F)
#write.table(report.df[,c(1:4)], "RUNX3_DiffBind_DESeq2_report_forHOMER.tab", sep="\t", quote=F, row.names=F)

write.table(report.df, "RUNX3_DiffBind_DESeq2_report_summitsT.tab", sep="\t", quote=F, row.names=F)
report.df.homer <- report.df[,c(1:4)]
report.df.homer$Strand <- "+"
write.table(report.df.homer, "RUNX3_DiffBind_DESeq2_report_summitsT_forHOMER.tab", sep="\t", quote=F, row.names=F)
