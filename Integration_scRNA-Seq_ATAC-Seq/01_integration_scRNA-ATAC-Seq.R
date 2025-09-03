library(Seurat)
library(RColorBrewer)
library(ggpubr)
library(svglite)
library(pheatmap)
library(ggrepel)

# read in clustered Seurat object

cd8t_cluster <- readRDS("cd8t_cluster_ownIdentsRes0.6.rds")

# Assign identity of clusters and assay, normalize and scale
Idents(object = cd8t_cluster) <- "Res0.6Names"
DefaultAssay(cd8t_cluster) <- "RNA"

mypal <- c("seashell3", "#0084B3", "#00B1B5", "palevioletred1",
           "#D44D35", "#F7B347", "khaki3", "darkseagreen")

mypal <- c("seashell3", "#0084B3", "#00B1B5", "palevioletred1",
           "#A71B4B", "#D44D35",
           "#F7B347", "khaki3",
           "darkseagreen", "#BAEEAE")

# read in gene list of interest

progdf <- read.csv2("../dfATAC_Texprog_significantDARs.csv")

progdfWT <- subset(progdf, progdf$log2FoldChange__group_KO_Tprog_vs_WT_Tprog < 0)
progdfWT[progdfWT == "Septin11"] <- "Sept11"
progdfKO <- subset(progdf, progdf$log2FoldChange__group_KO_Tprog_vs_WT_Tprog > 0)

termdf <- read.csv2("../dfATAC_nonTexprog_significantDARs.csv")

termdfWT <- subset(termdf, termdf$log2FoldChange__group_KO_Tterm_vs_WT_Tterm < 0)
termdfWT[termdfWT == "Marchf2"] <- "March2"
termdfWT[termdfWT == "Mars1"] <- "Mars"
termdfKO <- subset(termdf, termdf$log2FoldChange__group_KO_Tterm_vs_WT_Tterm > 0)
termdfKO[termdfKO == "Septin11"] <- "Sept11"
termdfKO[termdfKO == "H4f16"] <- "His4h4"
termdfKO[termdfKO == "Mtarc1"] <- "Marc1"
termdfKO[termdfKO == "Klrh1"] <- "Gm156"
termdfKO[termdfKO == "H1f8"] <- "H1foo"

# add module score

cd8t_cluster <- AddModuleScore(cd8t_cluster,
                               features = list(progdfWT$gencode_gene_name),
                               name = "WT",
                               ctrl = 100, #1000
                               nbin = 24) #25
cd8t_cluster <- AddModuleScore(cd8t_cluster,
                               features = list(progdfKO$gencode_gene_name),
                               name = "KO",
                               ctrl = 100,
                               nbin = 24)

cd8t_cluster <- AddModuleScore(cd8t_cluster,
                               features = list(termdfWT$gencode_gene_name),
                               name = "WT",
                               ctrl = 100, #1000
                               nbin = 24) #25
cd8t_cluster <- AddModuleScore(cd8t_cluster,
                               features = list(termdfKO$gencode_gene_name),
                               name = "KO",
                               ctrl = 100,
                               nbin = 24)

# get metadata (where module scores are stored)
# change column indices according to needed data, and exchange Tprog and Tterm where necessary
metadata <- cd8t_cluster@meta.data
colnames(metadata)[40:41] <- c("WT", "KO")
cd8t_cluster@meta.data <- metadata

FeaturePlot(cd8t_cluster,
            features = "WT", min.cutoff = "q5", max.cutoff = "q99",
            cols = rev(brewer.pal(n = 7, name = "RdBu")), pt.size = 0.4,
            order = T) +
  ggtitle("in WT nonTprog") -> pwt

FeaturePlot(cd8t_cluster,
            features = "KO", min.cutoff = "q5", max.cutoff = "q99",
            cols = rev(brewer.pal(n = 7, name = "RdBu")), pt.size = 0.4,
            order = T) +
  ggtitle("in KO nonTprog") -> pko

myplots <- list(pwt, pko)
myarrangedplot <- ggarrange(plotlist =  myplots, ncol = 2, 
                            common.legend = T, legend = "right")

pdf("../TtermDARsOnCluster.pdf", width = 10, height = 5)
annotate_figure(myarrangedplot, 
                top = text_grob("Genes with more accessible regions", 
                                face = "bold", size = 20))
dev.off()

# get module scores and cluster annotations for each cell in separate df

modulescore_df <- subset(metadata, select = c("Res0.6Names", "genotype", "WT", "KO"))


ggplot(modulescore_df, aes(x = Res0.6Names, y = WT)) +
  geom_hline(yintercept = 0, color = "grey30", linetype = "dashed") +
  geom_violin(aes(fill = Res0.6Names), draw_quantiles = 0.5, trim = T,
              scale = "count") +
  xlab("") +
  ylab("module score") +
  scale_fill_manual(values = mypal) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Genes with more accessible regions in WT nonTprog") +
  ylim(-0.05, 0.15) -> pwtmod

ggplot(modulescore_df, aes(x = Res0.6Names, y = KO)) +
  geom_hline(yintercept = 0, color = "grey30", linetype = "dashed") +
  geom_violin(aes(fill = Res0.6Names), draw_quantiles = 0.5, trim = T,
              scale = "count") +
  xlab("") +
  ylab("module score") +
  scale_fill_manual(values = mypal) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Genes with more accessible regions in H1-cKO nonTprog") +
  ylim(-0.05, 0.15) -> pkomod

mymodplots <- list(pwtmod, pkomod)
mymodarrangedplot <- ggarrange(plotlist =  mymodplots, nrow = 2)

pdf("../TtermDARsModuleScores_scRNACluster.pdf", width = 12, height = 6)
mymodarrangedplot
dev.off()

# get statistics on violin plots
library(multcomp)
library(emmeans)

WTmoduleAnova <- lm(WT ~ Res0.6Names, data = modulescore_df)
anova(WTmoduleAnova)
summary(WTmoduleAnova)

m1 <- emmeans(WTmoduleAnova, ~ Res0.6Names)
m1
pairs(m1)
capture.output(list(summary(WTmoduleAnova), m1, pairs(m1)), 
               file = "../WTTtermDARs_ModuleScores_AnovaTukeyResults.txt")


KOmoduleAnova <- lm(KO ~ Res0.6Names, data = modulescore_df)
anova(KOmoduleAnova)
summary(KOmoduleAnova)

m2 <- emmeans(KOmoduleAnova, ~ Res0.6Names)
m2
pairs(m2)
capture.output(list(summary(KOmoduleAnova), m2, pairs(m2)), 
               file = "../KOTtermDARs_ModuleScores_AnovaTukeyResults.txt")
