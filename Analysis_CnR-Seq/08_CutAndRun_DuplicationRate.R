# title: revision - Cut&Run analysis - Check Duplication Rate and Library Size
# author: Monika Waldherr
# input: summary TXT files from Picard MarkDuplicates

# code adapted from https://yezhengstat.github.io/CUTTag_tutorial/#33_Remove_duplicates_[optionalrequired]
library(dplyr)
library(ggplot2)
library(viridis)

setwd("./")
sampleList <- list.files("./picard_summary", pattern = "*mDup.txt") # dupMark.txt or rmDup.txt
histList <- c("IgG", "H3K27ac", "H3K27me3", "HDAC1", "RUNX3")
histInfoList <- gsub("CnR-20241114-CD8--", "", sampleList)
histInfoList <- gsub("_picard_dupMark.txt", "", histInfoList)

## Summarize the duplication information from the picard summary outputs.
dupResult = c()
for(i in 1:length(sampleList)){
  dupRes = read.table(paste0("./picard_summary/", sampleList[i]), header = TRUE, fill = TRUE)
  histInfo = histInfoList[i]
  histInfo_split = strsplit(histInfo, "-")[[1]]
  dupResult = data.frame(SampleName = histInfo, Histone = rev(histInfo_split)[2],
                         Replicate = rev(histInfo_split)[1], Genotype = rev(histInfo_split)[3], CellType = histInfo_split[1],
                         MappedFragNum_hg38 = dupRes$READ_PAIRS_EXAMINED[1] %>% as.character %>% as.numeric,
                         DuplicationRate = dupRes$PERCENT_DUPLICATION[1] %>% as.character %>% as.numeric * 100,
                         EstimatedLibrarySize = dupRes$ESTIMATED_LIBRARY_SIZE[1] %>% as.character %>% as.numeric) %>% 
    mutate(UniqueFragNum = MappedFragNum_hg38 * (1-DuplicationRate/100)) %>% 
    rbind(dupResult, .)
}
dupResult$Histone = factor(dupResult$Histone, levels = histList)
dupResult$Genotype = factor(dupResult$Genotype, levels = c("WT", "KO"))
dupResult$CellType = factor(dupResult$CellType, levels = c("Ly108", "Tim3"))

## generate boxplot figure for the  duplication rate
#install.packages("ggpubr")
library(ggpubr)

fig4A = dupResult %>% ggplot(aes(x = Histone, y = DuplicationRate, fill = Histone)) +
  geom_boxplot(outlier.shape = NA) +
  guides(fill="none") +
  geom_jitter(aes(color = CellType, shape = Genotype), position = position_jitter(0.15), size = 5) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "B", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "H") +
  theme_bw(base_size = 18) +
  ylab("Duplication Rate (%)") +
  xlab("") +
  scale_x_discrete(guide = guide_axis(angle = 45))

fig4B = dupResult %>% ggplot(aes(x = Histone, y = EstimatedLibrarySize, fill = Histone)) +
  geom_boxplot(outlier.shape = NA) +
  guides(fill="none") +
  geom_jitter(aes(color = CellType, shape = Genotype), position = position_jitter(0.15), size = 5) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "B", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "H") +
  theme_bw(base_size = 18) +
  ylab("Estimated Library Size") +
  xlab("") +
  scale_x_discrete(guide = guide_axis(angle = 45))

fig4C = dupResult %>% ggplot(aes(x = Histone, y = UniqueFragNum, fill = Histone)) +
  geom_boxplot(outlier.shape = NA) +
  guides(fill="none") +
  geom_jitter(aes(color = CellType, shape = Genotype), position = position_jitter(0.15), size = 5) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "B", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "H") +
  theme_bw(base_size = 18) +
  ylab("# of Unique Fragments") +
  xlab("") +
  scale_x_discrete(guide = guide_axis(angle = 45))

pdf("./PicardDuplication.pdf", width = 20, height = 10)
ggpubr::ggarrange(fig4A, fig4B, fig4C, ncol = 3, common.legend = TRUE, legend = "top")
dev.off()
