# title: Cut&Run analysis - visualization of Hdac1 and Runx3 binding overlapping with more open DARs
# author: Monika Waldherr
# input: numbers from intersections

#### make pie charts of how many open regions show HDAC1 and/or RUNX3 binding ####
library(RColorBrewer)
mygreens <- brewer.pal(9, "Greens")
myYlGn <- brewer.pal(9, "YlGn")
mypal <- c(myYlGn[2], myYlGn[3], mygreens[4], "grey90")

# WT open, HDAC1 WT, RUNX3 more WT
prop <- c(27, 17, 21, (102-27-17-21))
mylabs <- c("27 HDAC1 bound", "17 HDAC1 and more RUNX3 bound in WT", "21 more RUNX3 bound in WT",
            "37 open regions w/o HDAC1 or more RUNX3 binding in WT")

pdf("PieChart_TexProg_openWT_HDAC1RUNX3-WT.pdf", width = 15, height = 8)
pie(prop, labels = mylabs, col = mypal, clockwise = T, main = "HDAC1 and more RUNX3 binding in WT to 102 DARs in TexProg WT")
dev.off()

# WT open, HDAC1 KO, RUNX3 more KO
prop <- c(1, (102-0-1))
mylabs <- c("1 more RUNX3 bound in KO",
            "101 regions w/o HDAC1 or more RUNX3 binding in KO")

pdf("PieChart_TexProg_openWT_HDAC1RUNX3-KO.pdf", width = 15, height = 8)
pie(prop, labels = mylabs, col = mypal, clockwise = T, main = "no HDAC1 and more RUNX3 binding in KO to 102 DARs in TexProg WT")
dev.off()

# KO open, HDAC1 KO, RUNX3 KO
prop <- c(118, (217-118))
mylabs <- c("118 more RUNX3 bound in KO",
            "99 open regions w/o more RUNX3 bound in KO")

pdf("PieChart_TexProg_openKO_HDAC1RUNX3-KO.pdf", width = 15, height = 8)
pie(prop, labels = mylabs, col = mypal, clockwise = T,
    main = "no HDAC1 and more RUNX3 binding in KO to 217 DARs in TexProg KO")
dev.off()

# KO open, HDAC1 WT, RUNX3 WT
prop <- c(21, (217-21))
mylabs <- c("21 HDAC1 bound in WT",
            "196 open regions w/o HDAC1 bound in WT or more RUNX3 bound in WT")

pdf("PieChart_TexProg_openKO_HDAC1RUNX3-WT.pdf", width = 15, height = 8)
pie(prop, labels = mylabs, col = mypal, clockwise = T,
    main = "HDAC1 and more RUNX3 binding in WT to 217 DARs in TexProg KO")
dev.off()

# final table for pie charts
# for DARs open in WT TexProg
intersect_HDAC1_ATAC <- read.delim("../intersections/intersect_TexProg-ATACWTopen_TexProg-HDAC1.bed",
                                   header = F)
colnames(intersect_HDAC1_ATAC) <- c("chr", "start", "end", "consensus")
intersect_RUNX3_ATAC <- read.delim("../intersections/intersect_TexProg-ATACWTopen_TexProg-RUNX3.bed",
                                   header = F)
colnames(intersect_RUNX3_ATAC) <- c("chr", "start", "end", "consensus")

mytableWT <- subset(ATAC_Tprog_WT,
                    select = c("consensus", "gencode_chr", "gencode_start", "gencode_end", "homer_Annotation", "homer_Gene.Name"))
mytableWT$HDAC1binding <- "NA"
mytableWT[mytableWT$consensus %in% intersect_HDAC1_ATAC$consensus,]$HDAC1binding <- "TRUE"

mytableWT$RUNX3moreWT <- "NA"
mytableWT[mytableWT$consensus %in% intersect_RUNX3_ATAC$consensus,]$RUNX3moreWT <- "TRUE"

write.table(mytableWT, "TexProgWT_openDARs_HDAC1binding_RUNX3binding_comb.tab", sep = "\t", quote = F, row.names = F)

# for DARs open in KO TexProg
intersect_HDAC1_ATAC <- read.delim("../intersections/intersect_TexProg-ATACKOopen_TexProg-HDAC1.bed",
                                   header = F)
colnames(intersect_HDAC1_ATAC) <- c("chr", "start", "end", "consensus")
intersect_RUNX3_ATAC <- read.delim("../intersections/intersect_TexProg-ATACKOopen_TexProg-RUNX3.bed",
                                   header = F)
colnames(intersect_RUNX3_ATAC) <- c("chr", "start", "end", "consensus")

mytableKO <- subset(ATAC_Tprog_KO,
                    select = c("consensus", "gencode_chr", "gencode_start", "gencode_end", "homer_Annotation", "homer_Gene.Name"))
mytableKO$HDAC1binding <- "NA"
mytableKO[mytableKO$consensus %in% intersect_HDAC1_ATAC$consensus,]$HDAC1binding <- "TRUE"

mytableKO$RUNX3moreKO <- "NA"
mytableKO[mytableKO$consensus %in% intersect_RUNX3_ATAC$consensus,]$RUNX3moreKO <- "TRUE"

write.table(mytableKO, "TexProgKO_openDARs_HDAC1binding_RUNX3binding_comb.tab", sep = "\t", quote = F, row.names = F)
