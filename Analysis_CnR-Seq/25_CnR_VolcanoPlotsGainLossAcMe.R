# title: Cut&Run analysis - visualization of gain and loss of H3K27ac and H3K27me3 upon loss of HDAC1
# author: Monika Waldherr
# input: report files from DiffBind analysis with -th set to 1 (= no threshold)

#### Volcano plot showing H3K27ac gain and loss ####
report.df <- read.delim("H3K27ac_DiffBind_DESeq2_report_summitsT.tab", sep="\t")
CnR_Tprog_H3K27ac <- read.delim("H3K27ac_DiffBind_DESeq2_report_summitsT_forHOMER_ann.tab")
colnames(CnR_Tprog_H3K27ac)[1] <- "PeakID"
peakIDFold <- report.df[,c(1,8)]
CnR_Tprog_H3K27ac <- merge(CnR_Tprog_H3K27ac, peakIDFold, by = "PeakID")
length(subset(CnR_Tprog_H3K27ac, Fold < -0.5)$Gene.Name) # 241 regions
length(subset(CnR_Tprog_H3K27ac, Fold > 0.5)$Gene.Name) # 678 regions
length(unique(subset(CnR_Tprog_H3K27ac, Fold < -0.5)$Gene.Name)) # 189 genes
length(unique(subset(CnR_Tprog_H3K27ac, Fold > 0.5)$Gene.Name)) # 571 genes

moreAc <- subset(CnR_Tprog_H3K27ac, Fold > 0.5)
lessAc <- subset(CnR_Tprog_H3K27ac, Fold < -0.5)

length(intersect(unique(moreAc$Gene.Name), unique(hdac_peaks_ann$Gene.Name))) # 510 genes
length(intersect(unique(lessAc$Gene.Name), unique(hdac_peaks_ann$Gene.Name))) # 174 genes

custom_colors <- ifelse(
  report.df$FDR < 0.05 & report.df$Fold > 0.5, 'red',
  ifelse(report.df$FDR < 0.05 & report.df$Fold < -0.5, 'blue', 'grey')
)
names(custom_colors) <- rownames(report.df)

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
pdf("./VolcanoPlot_CnR-KOvsWT_H3K27ac.pdf", height = 12, width = 12)
myplot <- EnhancedVolcano(report.df, lab = NA,
                          x = "Fold", 
                          y = "FDR",
                          pointSize = 4,
                          FCcutoff = 0.5, 
                          pCutoff = 0.05,
                          title = "H3K27ac changes in TexProg",
                          subtitle = "678 regions (1965 total) gained acetylation, 241 regions (714 total) lost acetylation upon loss of HDAC1",
                          colCustom = custom_colors,
                          legendPosition = "none")
myplot <- myplot + ggplot2::coord_cartesian(xlim = c(-2,2))
print(myplot)
dev.off()

#### Volcano plot showing H3K27me3 gain and loss ####
report.df <- read.delim("H3K27me3_DiffBind_DESeq2_report_summitsT.tab", sep="\t")

report.df$Strand <- rep("+", nrow(report.df))

write.table(report.df[,c(1:4,11)], "H3K27me3_DiffBind_DESeq2_report_summitsT_forHOMER.tab", sep="\t", quote=F, row.names=F)

CnR_Tprog_H3K27me3 <- read.delim("H3K27me3_DiffBind_DESeq2_report_summitsT_forHOMER_ann.tab")
colnames(CnR_Tprog_H3K27me3)[1] <- "PeakID"
peakIDFold <- report.df[,c(1,8)]
CnR_Tprog_H3K27me3 <- merge(CnR_Tprog_H3K27me3, peakIDFold, by = "PeakID")
length(subset(CnR_Tprog_H3K27me3, Fold < -0.5)$Gene.Name) # 954 regions
length(subset(CnR_Tprog_H3K27me3, Fold > 0.5)$Gene.Name) # 537 regions
length(unique(subset(CnR_Tprog_H3K27me3, Fold < -0.5)$Gene.Name)) # 765 genes
length(unique(subset(CnR_Tprog_H3K27me3, Fold > 0.5)$Gene.Name)) # 432 genes

moreMe <- subset(CnR_Tprog_H3K27me3, Fold > 0.5)
lessMe <- subset(CnR_Tprog_H3K27me3, Fold < -0.5)

length(intersect(unique(moreMe$Gene.Name), unique(hdac_peaks_ann$Gene.Name))) # 334 genes
length(intersect(unique(lessMe$Gene.Name), unique(hdac_peaks_ann$Gene.Name))) # 675 genes

length(intersect(intersect(unique(moreAc$Gene.Name), unique(hdac_peaks_ann$Gene.Name)), unique(lessMe$Gene.Name))) # 112 genes show more Ac and less Me + HDAC1 binding

custom_colors <- ifelse(
  report.df$FDR < 0.05 & report.df$Fold > 0.5, 'red',
  ifelse(report.df$FDR < 0.05 & report.df$Fold < -0.5, 'blue', 'grey')
)
names(custom_colors) <- rownames(report.df)

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
pdf("./VolcanoPlot_CnR-KOvsWT_H3K27me3.pdf", height = 12, width = 12)
myplot <- EnhancedVolcano(report.df, lab = NA,
                          x = "Fold", 
                          y = "FDR",
                          pointSize = 4,
                          FCcutoff = 0.5, 
                          pCutoff = 0.05,
                          title = "H3K27me3 changes in TexProg",
                          subtitle = "537 regions (2598 total) gained methylation, 954 regions (3133 total) lost methylation upon loss of HDAC1",
                          colCustom = custom_colors,
                          legendPosition = "none")
myplot <- myplot + ggplot2::coord_cartesian(xlim = c(-2,2))
print(myplot)
dev.off()
