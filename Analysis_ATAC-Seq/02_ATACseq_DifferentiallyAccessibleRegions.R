# title: differential accessibility analysis of ATAC-Seq data
# author: Monika Waldherr
# input: dds object from script 01_ATACseq_QC.R

#### load libraries ####
library(ggplot2)
library(pander)
library(DESeq2)
library(ashr)

#### load dds object ###
dds <- readRDS("../ATACseq_dds.rds")

#### pairwise differential accessibility testing ####
plot_list <- list() # initiate plot_list

# pairwise comparison
pairwise_comp = list(c("KO_Tprog","WT_Tprog"),
                     c("KO_Tterm","WT_Tterm"))

# go over each specified pairwise comparison to be conducted
for (i in 1:length(pairwise_comp)){
  # extract current pairwise comparison
  comp_i <- pairwise_comp[[i]]
  
  # pairwise comparison by specification of contrast
  coef_i <- paste0("group_", comp_i[1], "_vs_", comp_i[2])
  dds_res <- DESeq(dds)
  
  # calculate FC shrinkage
  dds_res_shrink_i <- lfcShrink(dds_res, contrast = c("group", comp_i[1],  comp_i[2]), type = "ashr")
  
  # add DE testing results to df (shrinked FC)
  names_res <- names(dds_res_shrink_i)
  newnames_res <- paste0(names_res, "__", coef_i)
  names(dds_res_shrink_i) <- newnames_res
  df <- cbind(df, dds_res_shrink_i)
}

write.csv2(df, file = "../dfATACcountswithsignificances.csv",
           quote = F, row.names = F)

# extract significant DARs for TexProg and non-TexProg
First <- df[, c(1:8)]
TexProg <- df[grep("Tprog",names(resdf))]
nonTexProg <- df[grep("Tterm",names(resdf))]

dfTexProg <- cbind(First, TexProg)
dfnonTexProg <- cbind(First, nonTexProg)

# only keep significant results, order by log2 fold change and save
dfTexProg <- subset(dfTexProg, padj__group_KO_Tprog_vs_WT_Tprog < 0.05)
dfTexProg <- dfTexProg[order(dfTexProg$log2FoldChange__group_KO_Tprog_vs_WT_Tprog),]
dfnonTexProg <- subset(dfnonTexProg, padj__group_KO_Tterm_vs_WT_Tterm < 0.05)
dfnonTexProg <- dfnonTexProg[order(dfnonTexProg$log2FoldChange__group_KO_Tterm_vs_WT_Tterm),]

write.csv2(dfTexProg, file = "../dfATAC_Texprog_significantDARs.csv",
           quote = F, row.names = F)
write.csv2(dfnonTexProg, file = "../dfATAC_nonTexprog_significantDARs.csv",
           quote = F, row.names = F)

#### plot region annotations for all significant DAR groups ####
dfTexProg$homerAnno <- sapply(str_split(dfTexProg$homer_Annotation, " "),"[[",1) # add short annotation to data frame
dfnonTexProg$homerAnno <- sapply(str_split(dfnonTexProg$homer_Annotation, " "),"[[",1)

## Texprog KO vs WT ##
# get only necessary infos from df
openWT_Texprog <- subset(dfTexProg, dfTexProg$log2FoldChange__group_KO_Tprog_vs_WT_Tprog < 0)
WTpienums <- as.data.frame(table(openWT_Texprog$homerAnno))
WTpienums$Freq <- WTpienums$Freq*(-1)
WTpienums$group <- rep("open in WT", nrow(WTpienums))

openKO_Texprog <- subset(dfTexProg, dfTexProg$log2FoldChange__group_KO_Tprog_vs_WT_Tprog > 0)
KOpienums <- as.data.frame(table(openKO_Texprog$homerAnno))
KOpienums$group <- rep("open in KO", nrow(KOpienums))

Prognums <- rbind(WTpienums, KOpienums)
Prognums$Var1 <- factor(Prognums$Var1, levels =  c("intron", "Intergenic", "promoter-TSS",
                                                   "exon", "5'", "TTS", "3'", "non-coding"))
Prognums <- rbind(Prognums, c("5'", 0, "open in KO"), c("5'", 0, "open in WT"))
Prognums$Freq <- as.numeric(Prognums$Freq)

# plot colored barplots
pdf("../Texprog_signDARsKOvsWT_annotation.pdf", width = 10, height = 6)
ggplot(Prognums, aes(Var1, Freq, fill = group)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#1F78B4", "#A6CEE3")) +
  xlab("") +
  ylab("# of DARs") +
  ggtitle("Texprog KO vs WT") + 
  ylim(-(max(abs(Prognums$Freq))+5), (max(abs(Prognums$Freq))+5)) +
  theme_bw()
dev.off()

## Tterm KO vs WT ##
# get only necessary infos from df
openWT_nonTexprog <- subset(dfnonTexProg, dfnonTexProg$log2FoldChange__group_KO_Tterm_vs_WT_Tterm < 0)
WTpienums <- as.data.frame(table(openWT_nonTexprog$homerAnno))
WTpienums$Freq <- WTpienums$Freq*(-1)
WTpienums$group <- rep("open in WT", nrow(WTpienums))

openKO_nonTexprog <- subset(dfnonTexProg, dfnonTexProg$log2FoldChange__group_KO_Tterm_vs_WT_Tterm > 0)
KOpienums <- as.data.frame(table(openKO_nonTexprog$homerAnno))
KOpienums$group <- rep("open in KO", nrow(KOpienums))

Termnums <- rbind(WTpienums, KOpienums)
Termnums$Var1 <- factor(Termnums$Var1, levels =  c("intron", "Intergenic", "promoter-TSS",
                                                   "exon", "5'", "TTS", "3'", "non-coding"))

# plot colored barplots
pdf("../nonTexprog_signDARsKOvsWT_annotation.pdf", width = 10, height = 6)
ggplot(Termnums, aes(Var1, Freq, fill = group)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("mediumseagreen", "#B2DF8A")) +
  xlab("") +
  ylab("# of DARs") +
  ggtitle("nonTexprog KO vs WT") + 
  ylim(-(max(abs(Termnums$Freq))+5), (max(abs(Termnums$Freq))+5)) +
  theme_bw()
dev.off()





