# title: quality control of ATAC-Seq data
# author: Monika Waldherr
# input: path to directory containing three files:
#       .) consensus region count matrix for all samples
#       .) consensus region annotation file
#       .) sample annotation files

#### load libraries ####
library(ggplot2)
library(pander)
library(DESeq2)

#### load count matrix, consensus info and sample info ####
countsdf <- read.table("all_counts.csv", header = T, sep = ",")
colnames(countsdf)[1] <- "consensus"
consdf <- read.table("consensus_regions_annotation.csv", header = T, sep = ",",
                     fill = T, quote = "\"")
impconsdf <- consdf[,c(1:4,11:13,21)] # use unique id (peak/consensus id), chromosome, start, end, gene name and consensus characterization/annotation (here from gencode and homer)
colnames(impconsdf)[1] <- "consensus"
si <- read.table("all_annotation.csv", header = T, sep = ",")
si <- si[,c(1,3,4)]
colnames(si) = c("sample", "genotype", "celltype")
si$group <- paste0(si$genotype, "_", si$celltype)
rownames(si) = si$sample
si$sample = factor(si$sample)
si$genotype = factor(si$genotype)
si$celltype = factor(si$celltype)
si$group = factor(si$group)

df <- merge(countsdf, impconsdf) # combine all df to one containing all needed infos
df <- df[, c(1,14:ncol(df),2:13)] # resort columns

#### remove missing/low count peaks ####
# only keep rows where at least in one sample more than 50 reads were found #
df <- df[apply(df[,9:ncol(df)], 1, max) > 50,] # reduces from 121,367 to 58,045 consensus seqs

# only keep rows where no sample has less than 10 counts
df <- df[apply(df[,9:ncol(df)], 1, min) >= 10,] # 56,821 consensus seqs left

#### data exploration ####
rowsummary = data.frame(rowmeans = apply(df[, 9:20], 1, mean),
                        rowsds = apply(df[, 9:20], 1, sd))

# plot peak means vs peak sds to see how well they correlate
pdf("../PeakMeansSDs_all.pdf")
ggplot(data=rowsummary, aes(x=rowmeans, y=rowsds)) + geom_point() + 
  xlab("Peak means") + ylab("Peak SDs")
dev.off()

# remove H2-K1 (promoter-TSS) outlier
df <- df[apply(df[,9:ncol(df)], 1, mean) < 3000,]

# remove CD4 gene since this cannot be distinguished from artifact of KO generation
df <- df[df$homer_Gene.Name != "Cd4",]

# plot peak means vs peak sds to see how well they correlate after removal of two data points
rowsummary_woout = data.frame(rowmeans = apply(df[, 9:20], 1, mean),
                        rowsds = apply(df[, 9:20], 1, sd))

pdf("../PeakMeansSDs_woout.pdf")
ggplot(data=rowsummary_woout, aes(x=rowmeans, y=rowsds)) + geom_point() + 
  xlab("Peak means") + ylab("Peak SDs")
dev.off()

#### data normalization ####
counts = df[,9:ncol(df)]
dds = DESeqDataSetFromMatrix(countData = counts[,order(colnames(counts))], 
                             colData = si, design = ~ group)
dds = DESeq(dds) # calculate default deseq model

# get count table before normalization and plot
cm = data.frame(counts(dds, normalized=F))
rownames(cm) = df$consensus

pdf("../BoxplotCountsBeforeNorm.pdf")
par(mar=c(8,4,4,2))
boxplot(log2(cm+1),
        las=2,
        main="Boxplots (before normalization)",
        xaxt="n",
        yaxt="n",
        ylab="log2 counts",
        lwd=1.5, cex.main=1)
axis(side=1, at= 1:length(colnames(cm)),las=2, labels=colnames(cm), cex.axis = 0.75)
axis(side=2)
dev.off()

# get count table after normalization and plot
cmnorm = data.frame(counts(dds, normalized=TRUE))
rownames(cmnorm) = df$consensus

pdf("../BoxplotCountsAfterNorm.pdf")
par(mar=c(8,4,4,2))
boxplot(log2(cmnorm+1),
        las=2,
        main="Boxplots (after normalization)",
        xaxt="n",
        yaxt="n",
        ylab="log2 counts",
        lwd=1.5, cex.main=1)
axis(side=1, at= 1:length(colnames(cmnorm)),las=2, labels=colnames(cmnorm), cex.axis = 0.75)
axis(side=2)
dev.off()

#### data visualization ####

# add normalized counts to df and save
names(cmnorm) <- paste0("normalized_", names(cmnorm))
df <- cbind(df, cmnorm)
write.csv2(df, file = "../dfATACcountswithinfosandnormalization.csv",
           quote = F, row.names = F)

#### pca ####
m <- log2(cmnorm + 1)
names(m) <- names(cm)
pca = prcomp(t(m))

pcaData = as.data.frame(pca$x)
pcaData$sample=rownames(pcaData)
pcaData=merge(pcaData, si)
percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))
p=ggplot(data=pcaData, aes(x = PC1, y = PC2, color=group)) + geom_point(size=3)
p=p+xlab(paste0("PC1: ", percentVar[1], "% variance"))
p=p+ylab(paste0("PC2: ", percentVar[2], "% variance"))

pdf("../PCAgroup.pdf")
print(p)
dev.off()

#### save DESeq2 object and move on to 02_ATACseq_DifferentiallyAccessibleRegions ####
saveRDS(dds, "../ATACseq_dds.rds")







