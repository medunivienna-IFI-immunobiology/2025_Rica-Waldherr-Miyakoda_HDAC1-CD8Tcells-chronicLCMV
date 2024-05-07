# title: gene set enrichment analysis of ATAC-Seq data
# author: Monika Waldherr
# input: CSV files with open DARs from script 02_ATACseq_DifferentiallyAccessibleRegions.R

#### load libraries ####
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationHub)
library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(tibble)
library(genekitr)
library(org.Mm.eg.db)

#### read in data from CSV files and subset ####
dfTexProg <- read.csv2("../dfATAC_Texprog_significantDARs.csv")
dfnonTexProg <- read.csv2("../dfATAC_nonTexprog_significantDARs.csv")

openWT_Texprog <- subset(dfTexProg,
                         dfTexProg$log2FoldChange__group_KO_Tprog_vs_WT_Tprog < 0) # DARs open in WT Texprog
openKO_Texprog <- subset(dfTexProg,
                         dfTexProg$log2FoldChange__group_KO_Tprog_vs_WT_Tprog > 0) # DARs open in KO Texprog
openWT_nonTexprog <- subset(dfnonTexProg,
                            dfnonTexProg$log2FoldChange__group_KO_Tterm_vs_WT_Tterm < 0) # DARs open in WT nonTexprog
openKO_nonTexprog <- subset(dfnonTexProg,
                            dfnonTexProg$log2FoldChange__group_KO_Tterm_vs_WT_Tterm > 0) # DARs open in KO nonTexprog

# prepare data for gene set enrichment analysis function
df1 <- subset(openWT_Texprog, select = c("gencode_gene_name", "log2FoldChange__group_KO_Tprog_vs_WT_Tprog"))
colnames(df1) <- c("gene", "log2FC") # select only gene name and log2FC column and rename
df1$group <- rep("openWT_Texprog", nrow(df1)) # add column defining group

df2 <- subset(openKO_Texprog, select = c("gencode_gene_name", "log2FoldChange__group_KO_Tprog_vs_WT_Tprog"))
colnames(df2) <- c("gene", "log2FC")
df2$group <- rep("openKO_Texprog", nrow(df2))

df3 <- subset(openWT_nonTexprog, select = c("gencode_gene_name", "log2FoldChange__group_KO_Tterm_vs_WT_Tterm"))
colnames(df3) <- c("gene", "log2FC")
df3$group <- rep("openWT_nonTexprog", nrow(df3))

df4 <- subset(openKO_nonTexprog, select = c("gencode_gene_name", "log2FoldChange__group_KO_Tterm_vs_WT_Tterm"))
colnames(df4) <- c("gene", "log2FC")
df4$group <- rep("openKO_nonTexprog", nrow(df4))

df <- rbind(df1, df2, df3, df4) # combine all four data frames

dfsample <- split(df$gene, df$group) # split into gene lists per group

# add ENTREZ IDs
dfsample$openWT_Texprog = bitr(dfsample$openWT_Texprog, fromType = "SYMBOL",
                             toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
dfsample$openKO_Texprog = bitr(dfsample$openKO_Texprog, fromType = "SYMBOL",
                             toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
dfsample$openWT_nonTexprog = bitr(dfsample$openWT_nonTexprog, fromType = "SYMBOL",
                             toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
dfsample$openKO_nonTexprog = bitr(dfsample$openKO_nonTexprog, fromType = "SYMBOL",
                             toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

# get list of named character vectors containing ENTREZ IDs
genelist <- list("openWT_Texprog" = dfsample$openWT_Texprog$ENTREZID,
                 "openKO_Texprog" = dfsample$openKO_Texprog$ENTREZID,
                 "openWT_nonTexprog" = dfsample$openWT_nonTexprog$ENTREZID,
                 "openKO_nonTexprog" = dfsample$openKO_nonTexprog$ENTREZID)

# run enrichment analysis for GO terms (Biological Processes) on all groups
GOclusterplotBP <- compareCluster(geneCluster = genelist, fun = "enrichGO", 
                                  OrgDb = "org.Mm.eg.db", ont = "BP",
                                  pool = T)

# get table of all GO terms, group by cluster and sort by adjusted p-value, then save
allGOtermsBP <- GOclusterplotBP@compareClusterResult %>% group_by(Cluster) %>%
  top_n(n = nrow(GOclusterplotBP@compareClusterResult), wt = -p.adjust)

write.csv2(allGOtermsBP, "../GOBPenriched_allTerms.csv")

# get top 10 GO terms per group and corresponding adjusted p-values to use with Reviogo (web-based)
allGOtermsBPpadj <- GOclusterplotBP@compareClusterResult %>% group_by(Cluster) %>%
  top_n(n = 10, wt = -p.adjust)
allGOtermsBPpadj <- allGOtermsBPpadj[,c(2,7)]
write_delim(allGOtermsBPpadj, "../GOBPenriched_allTerms_Top10forRevigo.txt",
            delim = "\t", col_names = F)
## run Revigo on http://revigo.irb.hr/ ##

# read in results from Revigo (web-based)
revigoterms <- read_delim("../GOBPenriched_all_Revigo.tsv")
revigoterms <- subset(revigoterms, Representative == "null")
revigonames <- revigoterms$Name

# plot dotplot of enriched GO terms in all groups (reduced with Revigo)
pdf("../DotPlot_GOBPenriched_signDARs.pdf")
dotplot(GOclusterplotBP, font.size = 10, showCategory = revigonames) + 
  ggtitle("GO Biological Processes") +
  xlab("")
dev.off()

#### gsea analysis of DARs compared to effector like gene sets ####

# load immune signature db sets for mouse
m_df <- msigdbr(species = "Mus musculus", category = "C7", subcategory = "IMMUNESIGDB")
m_df_eff_tcell <- filter(m_df, grepl('CD8_TCELL', gs_name)) %>%
  filter(grepl('EFF', gs_name)) %>%
  filter(grepl('NAIVE|EXH', gs_name))

# prepare gene sets to have a list of gene sets
fgsea_sets <- m_df_eff_tcell %>% split(x = .$gene_symbol, f = .$gs_name)

# get gene lists of interest and rank them
allDARs <- read.csv2("../dfATACcountswithsignificances.csv")
allDARs$padj__group_KO_Tprog_vs_WT_Tprog <- as.numeric(allDARs$padj__group_KO_Tprog_vs_WT_Tprog)
allDARs$log2FoldChange__group_KO_Tprog_vs_WT_Tprog <- as.numeric(allDARs$log2FoldChange__group_KO_Tprog_vs_WT_Tprog)
myDARsTprog <- subset(allDARs, select = c("homer_Gene.Name",
                                          "log2FoldChange__group_KO_Tprog_vs_WT_Tprog",
                                          "padj__group_KO_Tprog_vs_WT_Tprog",
                                          "homer_Annotation"))
colnames(myDARsTprog) <- c("gene", "log2FC_DARs", "padj_DARs", "annotation")

# make sure to have unique gene list for ranking
myDARsTprog_uniq <- myDARsTprog %>%
  group_by(gene) %>%
  top_n(n = 1, wt = -padj_DARs)

myDARsTprog_uniq <- myDARsTprog_uniq %>% filter(!gene %in% myDARsTprog_uniq$gene[duplicated(myDARsTprog_uniq$gene)])

myDARsTprog_uniq$neglog_padj_DARs <- -log10(myDARsTprog_uniq$padj_DARs)
myDARsTprog_uniq$rank <- myDARsTprog_uniq$log2FC_DARs*myDARsTprog_uniq$neglog_padj_DARs

mydf <- myDARsTprog_uniq[,c(1,2)]
mydf <- mydf[order(mydf$log2FC_DARs, decreasing = F),]
ranks <- deframe(mydf)

# run enrichment analysis
fgseaRes <- fgsea(fgsea_sets, stats = ranks, scoreType = "std")

# clean up results to only have significant terms
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(padj < 0.05)

# plot normalized enrichment scores for all significant terms
pdf("../GSEA_EnrichmentPlot_ImmuneSigDb_CD8TcellEffectorExhausted_DARs.pdf", width = 12, height = 8)
ggplot(fgseaResTidy, 
       aes(pathway, NES)) +
  geom_col(aes(fill= NES)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="ImmuneSigDb CD8 Tcell Pathways (restricted to Effector & (Naive or Exhausted)") + 
  theme_minimal()
dev.off()

# plot typical GSEA curve for "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN - Exhausted Signature"
pdf("../GSEA_EnrichmentPlot_ExhaustedSignature_DARs.pdf", width = 10, height = 6)
plotEnrichment(fgsea_sets$GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN, ranks) +
  ggtitle("GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN - Exhausted Signature")
dev.off()

# plot typical GSEA curve for "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP - Effector Signature"
pdf("../GSEA_EnrichmentPlot_EffectorSignature_DARs.pdf", width = 10, height = 6)
plotEnrichment(fgsea_sets$GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP, ranks) +
  ggtitle("GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP - Effector Signature")
dev.off()

# save table with all significant results
for (i in 1:nrow(fgseaResTidy)){
  df_fgseaResTidsy$leadingEdge[i] <- paste(unlist(fgseaResTidy$leadingEdge[i]), collapse="/")
}
write.csv2(df_fgseaResTidsy, "../GSEA_EnrichmentPlot_Table.csv")





















