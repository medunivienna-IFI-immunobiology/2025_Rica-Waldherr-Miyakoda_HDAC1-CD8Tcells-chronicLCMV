# title: revision - Cut&Run analysis - Check Fragment Lengths
# author: Monika Waldherr
# input: TXT files with first column size, second column count (9th column of bam file is size)

# code adapted from https://yezhengstat.github.io/CUTTag_tutorial/#33_Remove_duplicates_[optionalrequired]
library(dplyr)
library(ggplot2)
library(viridis)

setwd("./")
sampleList <- list.files("./fragmentLen", pattern = "*fragmentLen.txt")
histList <- c("IgG", "H3K27ac", "H3K27me3", "HDAC1", "RUNX3")
histInfoList <- gsub("CnR-20241114-CD8--", "", sampleList)
histInfoList <- gsub("_fragmentLen.txt", "", histInfoList)

fragLen = c()
for(i in 1:length(sampleList)){
  histInfo = histInfoList[i]
  histInfo_split = strsplit(histInfo, "-")[[1]]
  fragLen = read.table(paste0("./fragmentLen/", sampleList[i]), header = FALSE) %>%
    mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)),
           Histone = rev(histInfo_split)[2], Replicate = rev(histInfo_split)[1],
           Genotype = rev(histInfo_split)[3], CellType = histInfo_split[1],
           SampleName = histInfo) %>%
    rbind(fragLen, .) 
}
fragLen$Histone = factor(fragLen$Histone, levels = histList)
fragLen$Genotype = factor(fragLen$Genotype, levels = c("WT", "KO"))
fragLen$CellType = factor(fragLen$CellType, levels = c("Ly108", "Tim3"))
fragLen$HistoneGenotype = factor(paste0(fragLen$Histone, " ", fragLen$Genotype),
                                 levels = c("IgG WT", "IgG KO",
                                            "H3K27ac WT", "H3K27ac KO", 
                                            "H3K27me3 WT", "H3K27me3 KO",
                                            "HDAC1 WT", "HDAC1 KO", 
                                            "RUNX3 WT", "RUNX3 KO"))

## Generate the fragment size density plot (violin plot)
fig5A = fragLen %>% ggplot(aes(x = HistoneGenotype, y = fragLen, weight = Weight, fill = Histone)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("")

# HDAC1
fragLen_HDAC1_Tprog <- subset(fragLen, Histone == "HDAC1" & CellType == "Ly108")
fig5B = fragLen_HDAC1_Tprog %>% ggplot(aes(x = fragLen, y = fragCount, color = Replicate, group = SampleName, linetype = Genotype)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("Tprog - HDAC1")

fragLen_HDAC1_nonTprog <- subset(fragLen, Histone == "HDAC1" & CellType == "Tim3")
fig5C = fragLen_HDAC1_nonTprog %>% ggplot(aes(x = fragLen, y = fragCount, color = Replicate, group = SampleName, linetype = Genotype)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("nonTprog - HDAC1")

# RUNX3
fragLen_RUNX3_Tprog <- subset(fragLen, Histone == "RUNX3" & CellType == "Ly108")
fig5D = fragLen_RUNX3_Tprog %>% ggplot(aes(x = fragLen, y = fragCount, color = Replicate, group = SampleName, linetype = Genotype)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("Tprog - RUNX3")

fragLen_RUNX3_nonTprog <- subset(fragLen, Histone == "RUNX3" & CellType == "Tim3")
fig5E = fragLen_RUNX3_nonTprog %>% ggplot(aes(x = fragLen, y = fragCount, color = Replicate, group = SampleName, linetype = Genotype)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("nonTprog - RUNX3")

# IgG
fragLen_IgG_Tprog <- subset(fragLen, Histone == "IgG" & CellType == "Ly108")
fig5F = fragLen_IgG_Tprog %>% ggplot(aes(x = fragLen, y = fragCount, color = Replicate, group = SampleName, linetype = Genotype)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("Tprog - IgG")

fragLen_IgG_nonTprog <- subset(fragLen, Histone == "IgG" & CellType == "Tim3")
fig5G = fragLen_IgG_nonTprog %>% ggplot(aes(x = fragLen, y = fragCount, color = Replicate, group = SampleName, linetype = Genotype)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("nonTprog - IgG")

# H3K27ac
fragLen_H3K27ac_Tprog <- subset(fragLen, Histone == "H3K27ac" & CellType == "Ly108")
fig5H = fragLen_H3K27ac_Tprog %>% ggplot(aes(x = fragLen, y = fragCount, color = Replicate, group = SampleName, linetype = Genotype)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("Tprog - H3K27ac")

fragLen_H3K27ac_nonTprog <- subset(fragLen, Histone == "H3K27ac" & CellType == "Tim3")
fig5I = fragLen_H3K27ac_nonTprog %>% ggplot(aes(x = fragLen, y = fragCount, color = Replicate, group = SampleName, linetype = Genotype)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("nonTprog - H3K27ac")

# H3K27me3
fragLen_H3K27me3_Tprog <- subset(fragLen, Histone == "H3K27me3" & CellType == "Ly108")
fig5J = fragLen_H3K27me3_Tprog %>% ggplot(aes(x = fragLen, y = fragCount, color = Replicate, group = SampleName, linetype = Genotype)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("Tprog - H3K27me3")

fragLen_H3K27me3_nonTprog <- subset(fragLen, Histone == "H3K27me3" & CellType == "Tim3")
fig5K = fragLen_H3K27me3_nonTprog %>% ggplot(aes(x = fragLen, y = fragCount, color = Replicate, group = SampleName, linetype = Genotype)) +
  geom_line(linewidth = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("nonTprog - H3K27me3")

library(ggpubr)
pdf("./FragmentLengthViolinPlots.pdf", width = 15, height = 10)
fig5A
dev.off()

pdf("./FragmentLengthHDAC1.pdf", width = 20, height = 10)
ggpubr::ggarrange(fig5B, fig5C, ncol = 2, common.legend = TRUE, legend = "bottom")
dev.off()

pdf("./FragmentLengthRUNX3.pdf", width = 20, height = 10)
ggpubr::ggarrange(fig5D, fig5E, ncol = 2, common.legend = TRUE, legend = "bottom")
dev.off()

pdf("./FragmentLengthIgG.pdf", width = 20, height = 10)
ggpubr::ggarrange(fig5F, fig5G, ncol = 2, common.legend = TRUE, legend = "bottom")
dev.off()

pdf("./FragmentLengthH3K27ac.pdf", width = 20, height = 10)
ggpubr::ggarrange(fig5H, fig5I, ncol = 2, common.legend = TRUE, legend = "bottom")
dev.off()

pdf("./FragmentLengthH3K27me3.pdf", width = 20, height = 10)
ggpubr::ggarrange(fig5J, fig5K, ncol = 2, common.legend = TRUE, legend = "bottom")
dev.off()







