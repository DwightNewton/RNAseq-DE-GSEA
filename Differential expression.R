library(DESeq2)
library(Rsamtools);
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
library(pheatmap)
library(RColorBrewer)
library(variancePartition)
library(vsn)
options(stringsAsFactors = FALSE)
options(scipen = 999)

#PV
load("seFullData_proteincoding.rData")

#Subset by cell-type
se_PV_exons_groupCT_age_sex_group <- seFullData_pc[,seFullData_pc$Cell.Type == "PVALB"]
se_PV_exons_groupCT_age_sex_group$Subject.Group <- relevel(as.factor(se_PV_exons_groupCT_age_sex_group$Subject.Group), "Control")

#Define design for DEseq2, and filter genes based on cell-type specific expression
dds0_PV_exons_groupCT_age_sex_group <- DESeqDataSet(se_PV_exons_groupCT_age_sex_group, design = ~ Age + Sex + Subject.Group)
isexprPV <- rowSums(counts(dds0_PV_exons_groupCT_age_sex_group)) > 30 & rowSums(counts(dds0_PV_exons_groupCT_age_sex_group) == 0) <= 60
sum(isexprPV)
dds0_PV_exons_groupCT_age_sex_group <- dds0_PV_exons_groupCT_age_sex_group[isexprPV,]

#call parallel cores - perform DE
register(MulticoreParam(workers=10))
dds_PV_exons_groupCT_age_sex_group <- DESeq(dds0_PV_exons_groupCT_age_sex_group, parallel=TRUE)
save(dds_PV_exons_groupCT_age_sex_group, file="dds_PV_exons_groupCT_age_sex_group.rData")

#Pull results for each contrast - no IF since FDRtool will be used for FDR correction later on
res_PV_exons_groupCT_age_sex_group_MDD <- results(dds_PV_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "MDD", "Control"),  independentFiltering=FALSE, parallel = TRUE)
res_PV_exons_groupCT_age_sex_group_Bipolar <- results(dds_PV_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "Bipolar", "Control"),  independentFiltering=FALSE, parallel = TRUE)
res_PV_exons_groupCT_age_sex_group_SCHIZ <- results(dds_PV_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "SCHIZ", "Control"),  independentFiltering=FALSE, parallel = TRUE)

save(res_PV_exons_groupCT_age_sex_group_MDD, res_PV_exons_groupCT_age_sex_group_Bipolar, res_PV_exons_groupCT_age_sex_group_SCHIZ, file="res_PV_exons_groupCT_age_sex_group.rData")

#Post-DE analysis using FDRtool
#p-value distribution and MA plots - Uncorrected
library(ggplot2)
options(stringsAsFactors = FALSE)

##MDD
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_PV_MDD.pdf")
ggplot(as.data.frame(res_PV_exons_groupCT_age_sex_group_MDD), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()

##Bipolar
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_PV_Bipolar.pdf")
ggplot(as.data.frame(res_PV_exons_groupCT_age_sex_group_Bipolar), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()

###SCHIZ
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_PV_SCHIZ.pdf")
ggplot(as.data.frame(res_PV_exons_groupCT_age_sex_group_SCHIZ), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()


###MA plots - uncorrected
###MDD
pdf("MAplot_FDR_MDD_PV_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_PV_exons_groupCT_age_sex_group_MDD, alpha=0.05)
dev.off()

###Bipolar
pdf("MAplot_FDR_Bipolar_PV_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_PV_exons_groupCT_age_sex_group_Bipolar, alpha=0.05)
dev.off()

###SCHIZ
pdf("MAplot_FDR_SCHIZ_PV_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_PV_exons_groupCT_age_sex_group_SCHIZ, alpha=0.05)
dev.off()


#FDRtool - CORRECTED pvalue histograms and MA plots
library(fdrtool)
#Correct overestimation of null  variance using fdrtool

#MDD
#PV
fdrT_exons_PV_MDD <- res_PV_exons_groupCT_age_sex_group_MDD
fdrT_exons_PV_MDD <- fdrT_exons_PV_MDD[ !is.na(fdrT_exons_PV_MDD$padj), ]
fdrT_exons_PV_MDD <- fdrT_exons_PV_MDD[ !is.na(fdrT_exons_PV_MDD$pvalue), ]
fdrT_exons_PV_MDD <- fdrT_exons_PV_MDD[, -which(names(fdrT_exons_PV_MDD) == "padj")]
FDR.fdrT_exons_PV_MDD <- fdrtool(fdrT_exons_PV_MDD$stat, statistic= "normal", plot = T)
fdrT_exons_PV_MDD[,"padj"]  <- p.adjust(FDR.fdrT_exons_PV_MDD$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_PV_MDD.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_PV_MDD), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_MDD_PV_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_PV_MDD, alpha=0.05)
dev.off()

#Bipolar
#PV
fdrT_exons_PV_Bipolar <- res_PV_exons_groupCT_age_sex_group_Bipolar
fdrT_exons_PV_Bipolar <- fdrT_exons_PV_Bipolar[ !is.na(fdrT_exons_PV_Bipolar$padj), ]
fdrT_exons_PV_Bipolar <- fdrT_exons_PV_Bipolar[ !is.na(fdrT_exons_PV_Bipolar$pvalue), ]
fdrT_exons_PV_Bipolar <- fdrT_exons_PV_Bipolar[, -which(names(fdrT_exons_PV_Bipolar) == "padj")]
FDR.fdrT_exons_PV_Bipolar <- fdrtool(fdrT_exons_PV_Bipolar$stat, statistic= "normal", plot = T)
fdrT_exons_PV_Bipolar[,"padj"]  <- p.adjust(FDR.fdrT_exons_PV_Bipolar$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_PV_Bipolar.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_PV_Bipolar), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_Bipolar_PV_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_PV_Bipolar, alpha=0.05)
dev.off()

#SCHIZ
#PV
fdrT_exons_PV_SCHIZ <- res_PV_exons_groupCT_age_sex_group_SCHIZ
fdrT_exons_PV_SCHIZ <- fdrT_exons_PV_SCHIZ[ !is.na(fdrT_exons_PV_SCHIZ$padj), ]
fdrT_exons_PV_SCHIZ <- fdrT_exons_PV_SCHIZ[ !is.na(fdrT_exons_PV_SCHIZ$pvalue), ]
fdrT_exons_PV_SCHIZ <- fdrT_exons_PV_SCHIZ[, -which(names(fdrT_exons_PV_SCHIZ) == "padj")]
FDR.fdrT_exons_PV_SCHIZ <- fdrtool(fdrT_exons_PV_SCHIZ$stat, statistic= "normal", plot = T)
fdrT_exons_PV_SCHIZ[,"padj"]  <- p.adjust(FDR.fdrT_exons_PV_SCHIZ$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_PV_SCHIZ.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_PV_SCHIZ), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_SCHIZ_PV_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_PV_SCHIZ, alpha=0.05)
dev.off()


##################Venn diagrams - output data for later merger
metadata <- as.data.frame(rowData(dds0_PV_exons_groupCT_age_sex_group))
metadata <- metadata[!duplicated(metadata$Gene.name),]
#FDR
FDR_exons_groupCT_sex_age_MDD_PV  <- as.data.frame(fdrT_exons_PV_MDD[,c("log2FoldChange", "padj")])
FDR_exons_groupCT_sex_age_Bipolar_PV  <- as.data.frame(fdrT_exons_PV_Bipolar[,c("log2FoldChange", "padj")])
FDR_exons_groupCT_sex_age_SCHIZ_PV  <- as.data.frame(fdrT_exons_PV_SCHIZ[,c("log2FoldChange", "padj")])

#Add metadata
FDR_exons_groupCT_sex_age_MDD_PV  <- merge(FDR_exons_groupCT_sex_age_MDD_PV, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")
FDR_exons_groupCT_sex_age_Bipolar_PV  <- merge(FDR_exons_groupCT_sex_age_Bipolar_PV, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")
FDR_exons_groupCT_sex_age_SCHIZ_PV  <- merge(FDR_exons_groupCT_sex_age_SCHIZ_PV, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")

#Get lists for venns
vecFDR_exons_groupCT_sex_age_MDD_PV  <- as.character(FDR_exons_groupCT_sex_age_MDD_PV$Gene.name[FDR_exons_groupCT_sex_age_MDD_PV$padj < 0.05])
vecFDR_exons_groupCT_sex_age_Bipolar_PV  <- as.character(FDR_exons_groupCT_sex_age_Bipolar_PV$Gene.name[FDR_exons_groupCT_sex_age_Bipolar_PV$padj < 0.05])
vecFDR_exons_groupCT_sex_age_SCHIZ_PV  <- as.character(FDR_exons_groupCT_sex_age_SCHIZ_PV$Gene.name[FDR_exons_groupCT_sex_age_SCHIZ_PV$padj < 0.05])

save(vecFDR_exons_groupCT_sex_age_MDD_PV, vecFDR_exons_groupCT_sex_age_Bipolar_PV, vecFDR_exons_groupCT_sex_age_SCHIZ_PV, file="CT_VennData_FDRtooled_PV_exons_groupCT_sex_age.rData")



#PYR23

load("seFullData_proteincoding.rData")
se_PYR23_exons_groupCT_age_sex_group <- seFullData_pc[,seFullData_pc$Cell.Type == "Pyr_L2n3"]
se_PYR23_exons_groupCT_age_sex_group$Subject.Group <- relevel(as.factor(se_PYR23_exons_groupCT_age_sex_group$Subject.Group), "Control")
dds0_PYR23_exons_groupCT_age_sex_group <- DESeqDataSet(se_PYR23_exons_groupCT_age_sex_group, design = ~ Age + Sex + Subject.Group)
isexprPYR23 <- rowSums(counts(dds0_PYR23_exons_groupCT_age_sex_group)) > 30 & rowSums(counts(dds0_PYR23_exons_groupCT_age_sex_group) == 0) <= 60
sum(isexprPYR23)
dds0_PYR23_exons_groupCT_age_sex_group <- dds0_PYR23_exons_groupCT_age_sex_group[isexprPYR23,]
#call parallel cores
register(MulticoreParam(workers=10))
dds_PYR23_exons_groupCT_age_sex_group <- DESeq(dds0_PYR23_exons_groupCT_age_sex_group, parallel=TRUE)
save(dds_PYR23_exons_groupCT_age_sex_group, file="dds_PYR23_exons_groupCT_age_sex_group.rData")

res_PYR23_exons_groupCT_age_sex_group_MDD <- results(dds_PYR23_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "MDD", "Control"),  independentFiltering=FALSE, parallel = TRUE)
res_PYR23_exons_groupCT_age_sex_group_Bipolar <- results(dds_PYR23_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "Bipolar", "Control"),  independentFiltering=FALSE, parallel = TRUE)
res_PYR23_exons_groupCT_age_sex_group_SCHIZ <- results(dds_PYR23_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "SCHIZ", "Control"),  independentFiltering=FALSE, parallel = TRUE)

save(res_PYR23_exons_groupCT_age_sex_group_MDD, res_PYR23_exons_groupCT_age_sex_group_Bipolar, res_PYR23_exons_groupCT_age_sex_group_SCHIZ, file="res_PYR23_exons_groupCT_age_sex_group.rData")

#Post-DE analysis
#p-value distribution and MA plots - Uncorrected
library(ggplot2)
options(stringsAsFactors = FALSE)

##MDD
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_PYR23_MDD.pdf")
ggplot(as.data.frame(res_PYR23_exons_groupCT_age_sex_group_MDD), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()

##Bipolar
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_PYR23_Bipolar.pdf")
ggplot(as.data.frame(res_PYR23_exons_groupCT_age_sex_group_Bipolar), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()

###SCHIZ
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_PYR23_SCHIZ.pdf")
ggplot(as.data.frame(res_PYR23_exons_groupCT_age_sex_group_SCHIZ), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()


###MA plots - uncorrected
###MDD
pdf("MAplot_FDR_MDD_PYR23_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_PYR23_exons_groupCT_age_sex_group_MDD, alpha=0.05)
dev.off()

###Bipolar
pdf("MAplot_FDR_Bipolar_PYR23_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_PYR23_exons_groupCT_age_sex_group_Bipolar, alpha=0.05)
dev.off()

###SCHIZ
pdf("MAplot_FDR_SCHIZ_PYR23_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_PYR23_exons_groupCT_age_sex_group_SCHIZ, alpha=0.05)
dev.off()


#FDRtool - CORRECTED pvalue histograms and MA plots
library(fdrtool)
#Correct overestimation of null  variance using fdrtool

#MDD
#PYR23
fdrT_exons_PYR23_MDD <- res_PYR23_exons_groupCT_age_sex_group_MDD
fdrT_exons_PYR23_MDD <- fdrT_exons_PYR23_MDD[ !is.na(fdrT_exons_PYR23_MDD$padj), ]
fdrT_exons_PYR23_MDD <- fdrT_exons_PYR23_MDD[ !is.na(fdrT_exons_PYR23_MDD$pvalue), ]
fdrT_exons_PYR23_MDD <- fdrT_exons_PYR23_MDD[, -which(names(fdrT_exons_PYR23_MDD) == "padj")]
FDR.fdrT_exons_PYR23_MDD <- fdrtool(fdrT_exons_PYR23_MDD$stat, statistic= "normal", plot = T)
fdrT_exons_PYR23_MDD[,"padj"]  <- p.adjust(FDR.fdrT_exons_PYR23_MDD$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_PYR23_MDD.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_PYR23_MDD), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_MDD_PYR23_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_PYR23_MDD, alpha=0.05)
dev.off()

#Bipolar
#PYR23
fdrT_exons_PYR23_Bipolar <- res_PYR23_exons_groupCT_age_sex_group_Bipolar
fdrT_exons_PYR23_Bipolar <- fdrT_exons_PYR23_Bipolar[ !is.na(fdrT_exons_PYR23_Bipolar$padj), ]
fdrT_exons_PYR23_Bipolar <- fdrT_exons_PYR23_Bipolar[ !is.na(fdrT_exons_PYR23_Bipolar$pvalue), ]
fdrT_exons_PYR23_Bipolar <- fdrT_exons_PYR23_Bipolar[, -which(names(fdrT_exons_PYR23_Bipolar) == "padj")]
FDR.fdrT_exons_PYR23_Bipolar <- fdrtool(fdrT_exons_PYR23_Bipolar$stat, statistic= "normal", plot = T)
fdrT_exons_PYR23_Bipolar[,"padj"]  <- p.adjust(FDR.fdrT_exons_PYR23_Bipolar$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_PYR23_Bipolar.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_PYR23_Bipolar), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_Bipolar_PYR23_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_PYR23_Bipolar, alpha=0.05)
dev.off()

#SCHIZ
#PYR23
fdrT_exons_PYR23_SCHIZ <- res_PYR23_exons_groupCT_age_sex_group_SCHIZ
fdrT_exons_PYR23_SCHIZ <- fdrT_exons_PYR23_SCHIZ[ !is.na(fdrT_exons_PYR23_SCHIZ$padj), ]
fdrT_exons_PYR23_SCHIZ <- fdrT_exons_PYR23_SCHIZ[ !is.na(fdrT_exons_PYR23_SCHIZ$pvalue), ]
fdrT_exons_PYR23_SCHIZ <- fdrT_exons_PYR23_SCHIZ[, -which(names(fdrT_exons_PYR23_SCHIZ) == "padj")]
FDR.fdrT_exons_PYR23_SCHIZ <- fdrtool(fdrT_exons_PYR23_SCHIZ$stat, statistic= "normal", plot = T)
fdrT_exons_PYR23_SCHIZ[,"padj"]  <- p.adjust(FDR.fdrT_exons_PYR23_SCHIZ$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_PYR23_SCHIZ.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_PYR23_SCHIZ), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_SCHIZ_PYR23_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_PYR23_SCHIZ, alpha=0.05)
dev.off()


##################Venn diagrams
metadata <- as.data.frame(rowData(dds0_PYR23_exons_groupCT_age_sex_group))
metadata <- metadata[!duplicated(metadata$Gene.name),]
#FDR
FDR_exons_groupCT_sex_age_MDD_PYR23  <- as.data.frame(fdrT_exons_PYR23_MDD[,c("log2FoldChange", "padj")])
FDR_exons_groupCT_sex_age_Bipolar_PYR23  <- as.data.frame(fdrT_exons_PYR23_Bipolar[,c("log2FoldChange", "padj")])
FDR_exons_groupCT_sex_age_SCHIZ_PYR23  <- as.data.frame(fdrT_exons_PYR23_SCHIZ[,c("log2FoldChange", "padj")])

#Add metadata
FDR_exons_groupCT_sex_age_MDD_PYR23  <- merge(FDR_exons_groupCT_sex_age_MDD_PYR23, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")
FDR_exons_groupCT_sex_age_Bipolar_PYR23  <- merge(FDR_exons_groupCT_sex_age_Bipolar_PYR23, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")
FDR_exons_groupCT_sex_age_SCHIZ_PYR23  <- merge(FDR_exons_groupCT_sex_age_SCHIZ_PYR23, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")

#Get lists for venns
vecFDR_exons_groupCT_sex_age_MDD_PYR23  <- as.character(FDR_exons_groupCT_sex_age_MDD_PYR23$Gene.name[FDR_exons_groupCT_sex_age_MDD_PYR23$padj < 0.05])
vecFDR_exons_groupCT_sex_age_Bipolar_PYR23  <- as.character(FDR_exons_groupCT_sex_age_Bipolar_PYR23$Gene.name[FDR_exons_groupCT_sex_age_Bipolar_PYR23$padj < 0.05])
vecFDR_exons_groupCT_sex_age_SCHIZ_PYR23  <- as.character(FDR_exons_groupCT_sex_age_SCHIZ_PYR23$Gene.name[FDR_exons_groupCT_sex_age_SCHIZ_PYR23$padj < 0.05])

save(vecFDR_exons_groupCT_sex_age_MDD_PYR23, vecFDR_exons_groupCT_sex_age_Bipolar_PYR23, vecFDR_exons_groupCT_sex_age_SCHIZ_PYR23, file="CT_VennData_FDRtooled_PYR23_exons_groupCT_sex_age.rData")



#PYR56

load("seFullData_proteincoding.rData")
se_PYR56_exons_groupCT_age_sex_group <- seFullData_pc[,seFullData_pc$Cell.Type == "Pyr_L5n6"]
se_PYR56_exons_groupCT_age_sex_group$Subject.Group <- relevel(as.factor(se_PYR56_exons_groupCT_age_sex_group$Subject.Group), "Control")
dds0_PYR56_exons_groupCT_age_sex_group <- DESeqDataSet(se_PYR56_exons_groupCT_age_sex_group, design = ~ Age + Sex + Subject.Group)
isexprPYR56 <- rowSums(counts(dds0_PYR56_exons_groupCT_age_sex_group)) > 30 & rowSums(counts(dds0_PYR56_exons_groupCT_age_sex_group) == 0) <= 60
sum(isexprPYR56)
dds0_PYR56_exons_groupCT_age_sex_group <- dds0_PYR56_exons_groupCT_age_sex_group[isexprPYR56,]
#call parallel cores
register(MulticoreParam(workers=10))
dds_PYR56_exons_groupCT_age_sex_group <- DESeq(dds0_PYR56_exons_groupCT_age_sex_group, parallel=TRUE)
save(dds_PYR56_exons_groupCT_age_sex_group, file="dds_PYR56_exons_groupCT_age_sex_group.rData")

res_PYR56_exons_groupCT_age_sex_group_MDD <- results(dds_PYR56_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "MDD", "Control"),  independentFiltering=FALSE, parallel = TRUE)
res_PYR56_exons_groupCT_age_sex_group_Bipolar <- results(dds_PYR56_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "Bipolar", "Control"),  independentFiltering=FALSE, parallel = TRUE)
res_PYR56_exons_groupCT_age_sex_group_SCHIZ <- results(dds_PYR56_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "SCHIZ", "Control"),  independentFiltering=FALSE, parallel = TRUE)

save(res_PYR56_exons_groupCT_age_sex_group_MDD, res_PYR56_exons_groupCT_age_sex_group_Bipolar, res_PYR56_exons_groupCT_age_sex_group_SCHIZ, file="res_PYR56_exons_groupCT_age_sex_group.rData")

#Post-DE analysis
#p-value distribution and MA plots - Uncorrected
library(ggplot2)
options(stringsAsFactors = FALSE)

##MDD
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_PYR56_MDD.pdf")
ggplot(as.data.frame(res_PYR56_exons_groupCT_age_sex_group_MDD), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()

##Bipolar
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_PYR56_Bipolar.pdf")
ggplot(as.data.frame(res_PYR56_exons_groupCT_age_sex_group_Bipolar), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()

###SCHIZ
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_PYR56_SCHIZ.pdf")
ggplot(as.data.frame(res_PYR56_exons_groupCT_age_sex_group_SCHIZ), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()


###MA plots - uncorrected
###MDD
pdf("MAplot_FDR_MDD_PYR56_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_PYR56_exons_groupCT_age_sex_group_MDD, alpha=0.05)
dev.off()

###Bipolar
pdf("MAplot_FDR_Bipolar_PYR56_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_PYR56_exons_groupCT_age_sex_group_Bipolar, alpha=0.05)
dev.off()

###SCHIZ
pdf("MAplot_FDR_SCHIZ_PYR56_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_PYR56_exons_groupCT_age_sex_group_SCHIZ, alpha=0.05)
dev.off()


#FDRtool - CORRECTED pvalue histograms and MA plots
library(fdrtool)
#Correct overestimation of null  variance using fdrtool

#MDD
#PYR56
fdrT_exons_PYR56_MDD <- res_PYR56_exons_groupCT_age_sex_group_MDD
fdrT_exons_PYR56_MDD <- fdrT_exons_PYR56_MDD[ !is.na(fdrT_exons_PYR56_MDD$padj), ]
fdrT_exons_PYR56_MDD <- fdrT_exons_PYR56_MDD[ !is.na(fdrT_exons_PYR56_MDD$pvalue), ]
fdrT_exons_PYR56_MDD <- fdrT_exons_PYR56_MDD[, -which(names(fdrT_exons_PYR56_MDD) == "padj")]
FDR.fdrT_exons_PYR56_MDD <- fdrtool(fdrT_exons_PYR56_MDD$stat, statistic= "normal", plot = T)
fdrT_exons_PYR56_MDD[,"padj"]  <- p.adjust(FDR.fdrT_exons_PYR56_MDD$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_PYR56_MDD.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_PYR56_MDD), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_MDD_PYR56_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_PYR56_MDD, alpha=0.05)
dev.off()

#Bipolar
#PYR56
fdrT_exons_PYR56_Bipolar <- res_PYR56_exons_groupCT_age_sex_group_Bipolar
fdrT_exons_PYR56_Bipolar <- fdrT_exons_PYR56_Bipolar[ !is.na(fdrT_exons_PYR56_Bipolar$padj), ]
fdrT_exons_PYR56_Bipolar <- fdrT_exons_PYR56_Bipolar[ !is.na(fdrT_exons_PYR56_Bipolar$pvalue), ]
fdrT_exons_PYR56_Bipolar <- fdrT_exons_PYR56_Bipolar[, -which(names(fdrT_exons_PYR56_Bipolar) == "padj")]
FDR.fdrT_exons_PYR56_Bipolar <- fdrtool(fdrT_exons_PYR56_Bipolar$stat, statistic= "normal", plot = T)
fdrT_exons_PYR56_Bipolar[,"padj"]  <- p.adjust(FDR.fdrT_exons_PYR56_Bipolar$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_PYR56_Bipolar.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_PYR56_Bipolar), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_Bipolar_PYR56_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_PYR56_Bipolar, alpha=0.05)
dev.off()

#SCHIZ
#PYR56
fdrT_exons_PYR56_SCHIZ <- res_PYR56_exons_groupCT_age_sex_group_SCHIZ
fdrT_exons_PYR56_SCHIZ <- fdrT_exons_PYR56_SCHIZ[ !is.na(fdrT_exons_PYR56_SCHIZ$padj), ]
fdrT_exons_PYR56_SCHIZ <- fdrT_exons_PYR56_SCHIZ[ !is.na(fdrT_exons_PYR56_SCHIZ$pvalue), ]
fdrT_exons_PYR56_SCHIZ <- fdrT_exons_PYR56_SCHIZ[, -which(names(fdrT_exons_PYR56_SCHIZ) == "padj")]
FDR.fdrT_exons_PYR56_SCHIZ <- fdrtool(fdrT_exons_PYR56_SCHIZ$stat, statistic= "normal", plot = T)
fdrT_exons_PYR56_SCHIZ[,"padj"]  <- p.adjust(FDR.fdrT_exons_PYR56_SCHIZ$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_PYR56_SCHIZ.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_PYR56_SCHIZ), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_SCHIZ_PYR56_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_PYR56_SCHIZ, alpha=0.05)
dev.off()


##################Venn diagrams
metadata <- as.data.frame(rowData(dds0_PYR56_exons_groupCT_age_sex_group))
metadata <- metadata[!duplicated(metadata$Gene.name),]
#FDR
FDR_exons_groupCT_sex_age_MDD_PYR56  <- as.data.frame(fdrT_exons_PYR56_MDD[,c("log2FoldChange", "padj")])
FDR_exons_groupCT_sex_age_Bipolar_PYR56  <- as.data.frame(fdrT_exons_PYR56_Bipolar[,c("log2FoldChange", "padj")])
FDR_exons_groupCT_sex_age_SCHIZ_PYR56  <- as.data.frame(fdrT_exons_PYR56_SCHIZ[,c("log2FoldChange", "padj")])

#Add metadata
FDR_exons_groupCT_sex_age_MDD_PYR56  <- merge(FDR_exons_groupCT_sex_age_MDD_PYR56, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")
FDR_exons_groupCT_sex_age_Bipolar_PYR56  <- merge(FDR_exons_groupCT_sex_age_Bipolar_PYR56, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")
FDR_exons_groupCT_sex_age_SCHIZ_PYR56  <- merge(FDR_exons_groupCT_sex_age_SCHIZ_PYR56, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")

#Get lists for venns
vecFDR_exons_groupCT_sex_age_MDD_PYR56  <- as.character(FDR_exons_groupCT_sex_age_MDD_PYR56$Gene.name[FDR_exons_groupCT_sex_age_MDD_PYR56$padj < 0.05])
vecFDR_exons_groupCT_sex_age_Bipolar_PYR56  <- as.character(FDR_exons_groupCT_sex_age_Bipolar_PYR56$Gene.name[FDR_exons_groupCT_sex_age_Bipolar_PYR56$padj < 0.05])
vecFDR_exons_groupCT_sex_age_SCHIZ_PYR56  <- as.character(FDR_exons_groupCT_sex_age_SCHIZ_PYR56$Gene.name[FDR_exons_groupCT_sex_age_SCHIZ_PYR56$padj < 0.05])

save(vecFDR_exons_groupCT_sex_age_MDD_PYR56, vecFDR_exons_groupCT_sex_age_Bipolar_PYR56, vecFDR_exons_groupCT_sex_age_SCHIZ_PYR56, file="CT_VennData_FDRtooled_PYR56_exons_groupCT_sex_age.rData")


#SST

load("seFullData_proteincoding.rData")
se_SST_exons_groupCT_age_sex_group <- seFullData_pc[,seFullData_pc$Cell.Type == "SST"]
se_SST_exons_groupCT_age_sex_group$Subject.Group <- relevel(as.factor(se_SST_exons_groupCT_age_sex_group$Subject.Group), "Control")
dds0_SST_exons_groupCT_age_sex_group <- DESeqDataSet(se_SST_exons_groupCT_age_sex_group, design = ~ Age + Sex + Subject.Group)
isexprSST <- rowSums(counts(dds0_SST_exons_groupCT_age_sex_group)) > 30 & rowSums(counts(dds0_SST_exons_groupCT_age_sex_group) == 0) <= 60
sum(isexprSST)
dds0_SST_exons_groupCT_age_sex_group <- dds0_SST_exons_groupCT_age_sex_group[isexprSST,]
#call parallel cores
register(MulticoreParam(workers=10))
dds_SST_exons_groupCT_age_sex_group <- DESeq(dds0_SST_exons_groupCT_age_sex_group, parallel=TRUE)
save(dds_SST_exons_groupCT_age_sex_group, file="dds_SST_exons_groupCT_age_sex_group.rData")

res_SST_exons_groupCT_age_sex_group_MDD <- results(dds_SST_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "MDD", "Control"),  independentFiltering=FALSE, parallel = TRUE)
res_SST_exons_groupCT_age_sex_group_Bipolar <- results(dds_SST_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "Bipolar", "Control"),  independentFiltering=FALSE, parallel = TRUE)
res_SST_exons_groupCT_age_sex_group_SCHIZ <- results(dds_SST_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "SCHIZ", "Control"),  independentFiltering=FALSE, parallel = TRUE)

save(res_SST_exons_groupCT_age_sex_group_MDD, res_SST_exons_groupCT_age_sex_group_Bipolar, res_SST_exons_groupCT_age_sex_group_SCHIZ, file="res_SST_exons_groupCT_age_sex_group.rData")

#Post-DE analysis
#p-value distribution and MA plots - Uncorrected
library(ggplot2)
options(stringsAsFactors = FALSE)

##MDD
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_SST_MDD.pdf")
ggplot(as.data.frame(res_SST_exons_groupCT_age_sex_group_MDD), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()

##Bipolar
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_SST_Bipolar.pdf")
ggplot(as.data.frame(res_SST_exons_groupCT_age_sex_group_Bipolar), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()

###SCHIZ
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_SST_SCHIZ.pdf")
ggplot(as.data.frame(res_SST_exons_groupCT_age_sex_group_SCHIZ), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()


###MA plots - uncorrected
###MDD
pdf("MAplot_FDR_MDD_SST_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_SST_exons_groupCT_age_sex_group_MDD, alpha=0.05)
dev.off()

###Bipolar
pdf("MAplot_FDR_Bipolar_SST_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_SST_exons_groupCT_age_sex_group_Bipolar, alpha=0.05)
dev.off()

###SCHIZ
pdf("MAplot_FDR_SCHIZ_SST_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_SST_exons_groupCT_age_sex_group_SCHIZ, alpha=0.05)
dev.off()


#FDRtool - CORRECTED pvalue histograms and MA plots
library(fdrtool)
#Correct overestimation of null  variance using fdrtool

#MDD
#SST
fdrT_exons_SST_MDD <- res_SST_exons_groupCT_age_sex_group_MDD
fdrT_exons_SST_MDD <- fdrT_exons_SST_MDD[ !is.na(fdrT_exons_SST_MDD$padj), ]
fdrT_exons_SST_MDD <- fdrT_exons_SST_MDD[ !is.na(fdrT_exons_SST_MDD$pvalue), ]
fdrT_exons_SST_MDD <- fdrT_exons_SST_MDD[, -which(names(fdrT_exons_SST_MDD) == "padj")]
FDR.fdrT_exons_SST_MDD <- fdrtool(fdrT_exons_SST_MDD$stat, statistic= "normal", plot = T)
fdrT_exons_SST_MDD[,"padj"]  <- p.adjust(FDR.fdrT_exons_SST_MDD$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_SST_MDD.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_SST_MDD), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_MDD_SST_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_SST_MDD, alpha=0.05)
dev.off()

#Bipolar
#SST
fdrT_exons_SST_Bipolar <- res_SST_exons_groupCT_age_sex_group_Bipolar
fdrT_exons_SST_Bipolar <- fdrT_exons_SST_Bipolar[ !is.na(fdrT_exons_SST_Bipolar$padj), ]
fdrT_exons_SST_Bipolar <- fdrT_exons_SST_Bipolar[ !is.na(fdrT_exons_SST_Bipolar$pvalue), ]
fdrT_exons_SST_Bipolar <- fdrT_exons_SST_Bipolar[, -which(names(fdrT_exons_SST_Bipolar) == "padj")]
FDR.fdrT_exons_SST_Bipolar <- fdrtool(fdrT_exons_SST_Bipolar$stat, statistic= "normal", plot = T)
fdrT_exons_SST_Bipolar[,"padj"]  <- p.adjust(FDR.fdrT_exons_SST_Bipolar$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_SST_Bipolar.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_SST_Bipolar), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_Bipolar_SST_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_SST_Bipolar, alpha=0.05)
dev.off()

#SCHIZ
#SST
fdrT_exons_SST_SCHIZ <- res_SST_exons_groupCT_age_sex_group_SCHIZ
fdrT_exons_SST_SCHIZ <- fdrT_exons_SST_SCHIZ[ !is.na(fdrT_exons_SST_SCHIZ$padj), ]
fdrT_exons_SST_SCHIZ <- fdrT_exons_SST_SCHIZ[ !is.na(fdrT_exons_SST_SCHIZ$pvalue), ]
fdrT_exons_SST_SCHIZ <- fdrT_exons_SST_SCHIZ[, -which(names(fdrT_exons_SST_SCHIZ) == "padj")]
FDR.fdrT_exons_SST_SCHIZ <- fdrtool(fdrT_exons_SST_SCHIZ$stat, statistic= "normal", plot = T)
fdrT_exons_SST_SCHIZ[,"padj"]  <- p.adjust(FDR.fdrT_exons_SST_SCHIZ$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_SST_SCHIZ.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_SST_SCHIZ), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_SCHIZ_SST_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_SST_SCHIZ, alpha=0.05)
dev.off()


##################Venn diagrams
metadata <- as.data.frame(rowData(dds0_SST_exons_groupCT_age_sex_group))
metadata <- metadata[!duplicated(metadata$Gene.name),]
#FDR
FDR_exons_groupCT_sex_age_MDD_SST  <- as.data.frame(fdrT_exons_SST_MDD[,c("log2FoldChange", "padj")])
FDR_exons_groupCT_sex_age_Bipolar_SST  <- as.data.frame(fdrT_exons_SST_Bipolar[,c("log2FoldChange", "padj")])
FDR_exons_groupCT_sex_age_SCHIZ_SST  <- as.data.frame(fdrT_exons_SST_SCHIZ[,c("log2FoldChange", "padj")])

#Add metadata
FDR_exons_groupCT_sex_age_MDD_SST  <- merge(FDR_exons_groupCT_sex_age_MDD_SST, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")
FDR_exons_groupCT_sex_age_Bipolar_SST  <- merge(FDR_exons_groupCT_sex_age_Bipolar_SST, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")
FDR_exons_groupCT_sex_age_SCHIZ_SST  <- merge(FDR_exons_groupCT_sex_age_SCHIZ_SST, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")

#Get lists for venns
vecFDR_exons_groupCT_sex_age_MDD_SST  <- as.character(FDR_exons_groupCT_sex_age_MDD_SST$Gene.name[FDR_exons_groupCT_sex_age_MDD_SST$padj < 0.05])
vecFDR_exons_groupCT_sex_age_Bipolar_SST  <- as.character(FDR_exons_groupCT_sex_age_Bipolar_SST$Gene.name[FDR_exons_groupCT_sex_age_Bipolar_SST$padj < 0.05])
vecFDR_exons_groupCT_sex_age_SCHIZ_SST  <- as.character(FDR_exons_groupCT_sex_age_SCHIZ_SST$Gene.name[FDR_exons_groupCT_sex_age_SCHIZ_SST$padj < 0.05])

save(vecFDR_exons_groupCT_sex_age_MDD_SST, vecFDR_exons_groupCT_sex_age_Bipolar_SST, vecFDR_exons_groupCT_sex_age_SCHIZ_SST, file="CT_VennData_FDRtooled_SST_exons_groupCT_sex_age.rData")


#VIP

load("seFullData_proteincoding.rData")
se_VIP_exons_groupCT_age_sex_group <- seFullData_pc[,seFullData_pc$Cell.Type == "VIP"]
se_VIP_exons_groupCT_age_sex_group$Subject.Group <- relevel(as.factor(se_VIP_exons_groupCT_age_sex_group$Subject.Group), "Control")
dds0_VIP_exons_groupCT_age_sex_group <- DESeqDataSet(se_VIP_exons_groupCT_age_sex_group, design = ~ Age + Sex + Subject.Group)
isexprVIP <- rowSums(counts(dds0_VIP_exons_groupCT_age_sex_group)) > 30 & rowSums(counts(dds0_VIP_exons_groupCT_age_sex_group) == 0) <= 60
sum(isexprVIP)
dds0_VIP_exons_groupCT_age_sex_group <- dds0_VIP_exons_groupCT_age_sex_group[isexprVIP,]
#call parallel cores
register(MulticoreParam(workers=10))
dds_VIP_exons_groupCT_age_sex_group <- DESeq(dds0_VIP_exons_groupCT_age_sex_group, parallel=TRUE)
save(dds_VIP_exons_groupCT_age_sex_group, file="dds_VIP_exons_groupCT_age_sex_group.rData")

res_VIP_exons_groupCT_age_sex_group_MDD <- results(dds_VIP_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "MDD", "Control"),  independentFiltering=FALSE, parallel = TRUE)
res_VIP_exons_groupCT_age_sex_group_Bipolar <- results(dds_VIP_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "Bipolar", "Control"),  independentFiltering=FALSE, parallel = TRUE)
res_VIP_exons_groupCT_age_sex_group_SCHIZ <- results(dds_VIP_exons_groupCT_age_sex_group, contrast=c("Subject.Group", "SCHIZ", "Control"),  independentFiltering=FALSE, parallel = TRUE)

save(res_VIP_exons_groupCT_age_sex_group_MDD, res_VIP_exons_groupCT_age_sex_group_Bipolar, res_VIP_exons_groupCT_age_sex_group_SCHIZ, file="res_VIP_exons_groupCT_age_sex_group.rData")

#Post-DE analysis
#p-value distribution and MA plots - Uncorrected
library(ggplot2)
options(stringsAsFactors = FALSE)

##MDD
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_VIP_MDD.pdf")
ggplot(as.data.frame(res_VIP_exons_groupCT_age_sex_group_MDD), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()

##Bipolar
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_VIP_Bipolar.pdf")
ggplot(as.data.frame(res_VIP_exons_groupCT_age_sex_group_Bipolar), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()

###SCHIZ
pdf("pvalue_hist_uncorrected_exons_groupCT_sex_age_VIP_SCHIZ.pdf")
ggplot(as.data.frame(res_VIP_exons_groupCT_age_sex_group_SCHIZ), aes(x=pvalue)) +
  geom_histogram(binwidth=0.05)
dev.off()


###MA plots - uncorrected
###MDD
pdf("MAplot_FDR_MDD_VIP_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_VIP_exons_groupCT_age_sex_group_MDD, alpha=0.05)
dev.off()

###Bipolar
pdf("MAplot_FDR_Bipolar_VIP_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_VIP_exons_groupCT_age_sex_group_Bipolar, alpha=0.05)
dev.off()

###SCHIZ
pdf("MAplot_FDR_SCHIZ_VIP_uncorrected_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(res_VIP_exons_groupCT_age_sex_group_SCHIZ, alpha=0.05)
dev.off()


#FDRtool - CORRECTED pvalue histograms and MA plots
library(fdrtool)
#Correct overestimation of null  variance using fdrtool

#MDD
#VIP
fdrT_exons_VIP_MDD <- res_VIP_exons_groupCT_age_sex_group_MDD
fdrT_exons_VIP_MDD <- fdrT_exons_VIP_MDD[ !is.na(fdrT_exons_VIP_MDD$padj), ]
fdrT_exons_VIP_MDD <- fdrT_exons_VIP_MDD[ !is.na(fdrT_exons_VIP_MDD$pvalue), ]
fdrT_exons_VIP_MDD <- fdrT_exons_VIP_MDD[, -which(names(fdrT_exons_VIP_MDD) == "padj")]
FDR.fdrT_exons_VIP_MDD <- fdrtool(fdrT_exons_VIP_MDD$stat, statistic= "normal", plot = T)
fdrT_exons_VIP_MDD[,"padj"]  <- p.adjust(FDR.fdrT_exons_VIP_MDD$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_VIP_MDD.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_VIP_MDD), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_MDD_VIP_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_VIP_MDD, alpha=0.05)
dev.off()

#Bipolar
#VIP
fdrT_exons_VIP_Bipolar <- res_VIP_exons_groupCT_age_sex_group_Bipolar
fdrT_exons_VIP_Bipolar <- fdrT_exons_VIP_Bipolar[ !is.na(fdrT_exons_VIP_Bipolar$padj), ]
fdrT_exons_VIP_Bipolar <- fdrT_exons_VIP_Bipolar[ !is.na(fdrT_exons_VIP_Bipolar$pvalue), ]
fdrT_exons_VIP_Bipolar <- fdrT_exons_VIP_Bipolar[, -which(names(fdrT_exons_VIP_Bipolar) == "padj")]
FDR.fdrT_exons_VIP_Bipolar <- fdrtool(fdrT_exons_VIP_Bipolar$stat, statistic= "normal", plot = T)
fdrT_exons_VIP_Bipolar[,"padj"]  <- p.adjust(FDR.fdrT_exons_VIP_Bipolar$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_VIP_Bipolar.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_VIP_Bipolar), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_Bipolar_VIP_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_VIP_Bipolar, alpha=0.05)
dev.off()

#SCHIZ
#VIP
fdrT_exons_VIP_SCHIZ <- res_VIP_exons_groupCT_age_sex_group_SCHIZ
fdrT_exons_VIP_SCHIZ <- fdrT_exons_VIP_SCHIZ[ !is.na(fdrT_exons_VIP_SCHIZ$padj), ]
fdrT_exons_VIP_SCHIZ <- fdrT_exons_VIP_SCHIZ[ !is.na(fdrT_exons_VIP_SCHIZ$pvalue), ]
fdrT_exons_VIP_SCHIZ <- fdrT_exons_VIP_SCHIZ[, -which(names(fdrT_exons_VIP_SCHIZ) == "padj")]
FDR.fdrT_exons_VIP_SCHIZ <- fdrtool(fdrT_exons_VIP_SCHIZ$stat, statistic= "normal", plot = T)
fdrT_exons_VIP_SCHIZ[,"padj"]  <- p.adjust(FDR.fdrT_exons_VIP_SCHIZ$pval, method = "BH")

#pvalue histogram

pdf("pvalue_hist_FDRtooled_exons_groupCT_sex_age_VIP_SCHIZ.pdf")
ggplot(as.data.frame(FDR.fdrT_exons_VIP_SCHIZ), aes(x=pval)) +
  geom_histogram(binwidth=0.05)
dev.off()

#MA plot
pdf("MAplot_FDR_SCHIZ_VIP_FDRtooled_unshrunk_exons_groupCT_sex_age_group_.pdf")
DESeq2::plotMA(fdrT_exons_VIP_SCHIZ, alpha=0.05)
dev.off()


##################Venn diagrams
metadata <- as.data.frame(rowData(dds0_VIP_exons_groupCT_age_sex_group))
metadata <- metadata[!duplicated(metadata$Gene.name),]
#FDR
FDR_exons_groupCT_sex_age_MDD_VIP  <- as.data.frame(fdrT_exons_VIP_MDD[,c("log2FoldChange", "padj")])
FDR_exons_groupCT_sex_age_Bipolar_VIP  <- as.data.frame(fdrT_exons_VIP_Bipolar[,c("log2FoldChange", "padj")])
FDR_exons_groupCT_sex_age_SCHIZ_VIP  <- as.data.frame(fdrT_exons_VIP_SCHIZ[,c("log2FoldChange", "padj")])

#Add metadata
FDR_exons_groupCT_sex_age_MDD_VIP  <- merge(FDR_exons_groupCT_sex_age_MDD_VIP, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")
FDR_exons_groupCT_sex_age_Bipolar_VIP  <- merge(FDR_exons_groupCT_sex_age_Bipolar_VIP, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")
FDR_exons_groupCT_sex_age_SCHIZ_VIP  <- merge(FDR_exons_groupCT_sex_age_SCHIZ_VIP, metadata[,c(1,6)], by.x="row.names", by.y="Gene.stable.ID")

#Get lists for venns
vecFDR_exons_groupCT_sex_age_MDD_VIP  <- as.character(FDR_exons_groupCT_sex_age_MDD_VIP$Gene.name[FDR_exons_groupCT_sex_age_MDD_VIP$padj < 0.05])
vecFDR_exons_groupCT_sex_age_Bipolar_VIP  <- as.character(FDR_exons_groupCT_sex_age_Bipolar_VIP$Gene.name[FDR_exons_groupCT_sex_age_Bipolar_VIP$padj < 0.05])
vecFDR_exons_groupCT_sex_age_SCHIZ_VIP  <- as.character(FDR_exons_groupCT_sex_age_SCHIZ_VIP$Gene.name[FDR_exons_groupCT_sex_age_SCHIZ_VIP$padj < 0.05])

save(vecFDR_exons_groupCT_sex_age_MDD_VIP, vecFDR_exons_groupCT_sex_age_Bipolar_VIP, vecFDR_exons_groupCT_sex_age_SCHIZ_VIP, file="CT_VennData_FDRtooled_VIP_exons_groupCT_sex_age.rData")




