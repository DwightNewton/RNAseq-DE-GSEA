library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(fdrtool)
options(stringsAsFactors = FALSE)
options(scipen = 999)

#Generate heatmaps for all genes, 5percent FDR, and 25percent FDR
filelist <- list.files(pattern=".rData")
for(i in filelist){
  load(i)
}

#FDRtooling
########MDD
fdrT_CT_PV_MDD <- res_PV_exons_groupCT_age_sex_group_MDD
fdrT_CT_PV_MDD <- fdrT_CT_PV_MDD[ !is.na(fdrT_CT_PV_MDD$padj), ]
fdrT_CT_PV_MDD <- fdrT_CT_PV_MDD[ !is.na(fdrT_CT_PV_MDD$pvalue), ]
fdrT_CT_PV_MDD <- fdrT_CT_PV_MDD[, -which(names(fdrT_CT_PV_MDD) == "padj")]
FDR.fdrT_CT_PV_MDD <- fdrtool(fdrT_CT_PV_MDD$stat, statistic= "normal", plot = T)
fdrT_CT_PV_MDD[,"padj"]  <- p.adjust(FDR.fdrT_CT_PV_MDD$pval, method = "BH")


fdrT_CT_PYR23_MDD <- res_PYR23_exons_groupCT_age_sex_group_MDD
fdrT_CT_PYR23_MDD <- fdrT_CT_PYR23_MDD[ !is.na(fdrT_CT_PYR23_MDD$padj), ]
fdrT_CT_PYR23_MDD <- fdrT_CT_PYR23_MDD[ !is.na(fdrT_CT_PYR23_MDD$pvalue), ]
fdrT_CT_PYR23_MDD <- fdrT_CT_PYR23_MDD[, -which(names(fdrT_CT_PYR23_MDD) == "padj")]
FDR.fdrT_CT_PYR23_MDD <- fdrtool(fdrT_CT_PYR23_MDD$stat, statistic= "normal", plot = T)
fdrT_CT_PYR23_MDD[,"padj"]  <- p.adjust(FDR.fdrT_CT_PYR23_MDD$pval, method = "BH")

fdrT_CT_PYR56_MDD <- res_PYR56_exons_groupCT_age_sex_group_MDD
fdrT_CT_PYR56_MDD <- fdrT_CT_PYR56_MDD[ !is.na(fdrT_CT_PYR56_MDD$padj), ]
fdrT_CT_PYR56_MDD <- fdrT_CT_PYR56_MDD[ !is.na(fdrT_CT_PYR56_MDD$pvalue), ]
fdrT_CT_PYR56_MDD <- fdrT_CT_PYR56_MDD[, -which(names(fdrT_CT_PYR56_MDD) == "padj")]
FDR.fdrT_CT_PYR56_MDD <- fdrtool(fdrT_CT_PYR56_MDD$stat, statistic= "normal", plot = T)
fdrT_CT_PYR56_MDD[,"padj"]  <- p.adjust(FDR.fdrT_CT_PYR56_MDD$pval, method = "BH")

fdrT_CT_SST_MDD <- res_SST_exons_groupCT_age_sex_group_MDD
fdrT_CT_SST_MDD <- fdrT_CT_SST_MDD[ !is.na(fdrT_CT_SST_MDD$padj), ]
fdrT_CT_SST_MDD <- fdrT_CT_SST_MDD[ !is.na(fdrT_CT_SST_MDD$pvalue), ]
fdrT_CT_SST_MDD <- fdrT_CT_SST_MDD[, -which(names(fdrT_CT_SST_MDD) == "padj")]
FDR.fdrT_CT_SST_MDD <- fdrtool(fdrT_CT_SST_MDD$stat, statistic= "normal", plot = T)
fdrT_CT_SST_MDD[,"padj"]  <- p.adjust(FDR.fdrT_CT_SST_MDD$pval, method = "BH")

fdrT_CT_VIP_MDD <- res_VIP_exons_groupCT_age_sex_group_MDD
fdrT_CT_VIP_MDD <- fdrT_CT_VIP_MDD[ !is.na(fdrT_CT_VIP_MDD$padj), ]
fdrT_CT_VIP_MDD <- fdrT_CT_VIP_MDD[ !is.na(fdrT_CT_VIP_MDD$pvalue), ]
fdrT_CT_VIP_MDD <- fdrT_CT_VIP_MDD[, -which(names(fdrT_CT_VIP_MDD) == "padj")]
FDR.fdrT_CT_VIP_MDD <- fdrtool(fdrT_CT_VIP_MDD$stat, statistic= "normal", plot = T)
fdrT_CT_VIP_MDD[,"padj"]  <- p.adjust(FDR.fdrT_CT_VIP_MDD$pval, method = "BH")



########Bipolar
fdrT_CT_PV_Bipolar <- res_PV_exons_groupCT_age_sex_group_Bipolar
fdrT_CT_PV_Bipolar <- fdrT_CT_PV_Bipolar[ !is.na(fdrT_CT_PV_Bipolar$padj), ]
fdrT_CT_PV_Bipolar <- fdrT_CT_PV_Bipolar[ !is.na(fdrT_CT_PV_Bipolar$pvalue), ]
fdrT_CT_PV_Bipolar <- fdrT_CT_PV_Bipolar[, -which(names(fdrT_CT_PV_Bipolar) == "padj")]
FDR.fdrT_CT_PV_Bipolar <- fdrtool(fdrT_CT_PV_Bipolar$stat, statistic= "normal", plot = T)
fdrT_CT_PV_Bipolar[,"padj"]  <- p.adjust(FDR.fdrT_CT_PV_Bipolar$pval, method = "BH")


fdrT_CT_PYR23_Bipolar <- res_PYR23_exons_groupCT_age_sex_group_Bipolar
fdrT_CT_PYR23_Bipolar <- fdrT_CT_PYR23_Bipolar[ !is.na(fdrT_CT_PYR23_Bipolar$padj), ]
fdrT_CT_PYR23_Bipolar <- fdrT_CT_PYR23_Bipolar[ !is.na(fdrT_CT_PYR23_Bipolar$pvalue), ]
fdrT_CT_PYR23_Bipolar <- fdrT_CT_PYR23_Bipolar[, -which(names(fdrT_CT_PYR23_Bipolar) == "padj")]
FDR.fdrT_CT_PYR23_Bipolar <- fdrtool(fdrT_CT_PYR23_Bipolar$stat, statistic= "normal", plot = T)
fdrT_CT_PYR23_Bipolar[,"padj"]  <- p.adjust(FDR.fdrT_CT_PYR23_Bipolar$pval, method = "BH")

fdrT_CT_PYR56_Bipolar <- res_PYR56_exons_groupCT_age_sex_group_Bipolar
fdrT_CT_PYR56_Bipolar <- fdrT_CT_PYR56_Bipolar[ !is.na(fdrT_CT_PYR56_Bipolar$padj), ]
fdrT_CT_PYR56_Bipolar <- fdrT_CT_PYR56_Bipolar[ !is.na(fdrT_CT_PYR56_Bipolar$pvalue), ]
fdrT_CT_PYR56_Bipolar <- fdrT_CT_PYR56_Bipolar[, -which(names(fdrT_CT_PYR56_Bipolar) == "padj")]
FDR.fdrT_CT_PYR56_Bipolar <- fdrtool(fdrT_CT_PYR56_Bipolar$stat, statistic= "normal", plot = T)
fdrT_CT_PYR56_Bipolar[,"padj"]  <- p.adjust(FDR.fdrT_CT_PYR56_Bipolar$pval, method = "BH")

fdrT_CT_SST_Bipolar <- res_SST_exons_groupCT_age_sex_group_Bipolar
fdrT_CT_SST_Bipolar <- fdrT_CT_SST_Bipolar[ !is.na(fdrT_CT_SST_Bipolar$padj), ]
fdrT_CT_SST_Bipolar <- fdrT_CT_SST_Bipolar[ !is.na(fdrT_CT_SST_Bipolar$pvalue), ]
fdrT_CT_SST_Bipolar <- fdrT_CT_SST_Bipolar[, -which(names(fdrT_CT_SST_Bipolar) == "padj")]
FDR.fdrT_CT_SST_Bipolar <- fdrtool(fdrT_CT_SST_Bipolar$stat, statistic= "normal", plot = T)
fdrT_CT_SST_Bipolar[,"padj"]  <- p.adjust(FDR.fdrT_CT_SST_Bipolar$pval, method = "BH")

fdrT_CT_VIP_Bipolar <- res_VIP_exons_groupCT_age_sex_group_Bipolar
fdrT_CT_VIP_Bipolar <- fdrT_CT_VIP_Bipolar[ !is.na(fdrT_CT_VIP_Bipolar$padj), ]
fdrT_CT_VIP_Bipolar <- fdrT_CT_VIP_Bipolar[ !is.na(fdrT_CT_VIP_Bipolar$pvalue), ]
fdrT_CT_VIP_Bipolar <- fdrT_CT_VIP_Bipolar[, -which(names(fdrT_CT_VIP_Bipolar) == "padj")]
FDR.fdrT_CT_VIP_Bipolar <- fdrtool(fdrT_CT_VIP_Bipolar$stat, statistic= "normal", plot = T)
fdrT_CT_VIP_Bipolar[,"padj"]  <- p.adjust(FDR.fdrT_CT_VIP_Bipolar$pval, method = "BH")


########SCHIZ
fdrT_CT_PV_SCHIZ <- res_PV_exons_groupCT_age_sex_group_SCHIZ
fdrT_CT_PV_SCHIZ <- fdrT_CT_PV_SCHIZ[ !is.na(fdrT_CT_PV_SCHIZ$padj), ]
fdrT_CT_PV_SCHIZ <- fdrT_CT_PV_SCHIZ[ !is.na(fdrT_CT_PV_SCHIZ$pvalue), ]
fdrT_CT_PV_SCHIZ <- fdrT_CT_PV_SCHIZ[, -which(names(fdrT_CT_PV_SCHIZ) == "padj")]
FDR.fdrT_CT_PV_SCHIZ <- fdrtool(fdrT_CT_PV_SCHIZ$stat, statistic= "normal", plot = T)
fdrT_CT_PV_SCHIZ[,"padj"]  <- p.adjust(FDR.fdrT_CT_PV_SCHIZ$pval, method = "BH")


fdrT_CT_PYR23_SCHIZ <- res_PYR23_exons_groupCT_age_sex_group_SCHIZ
fdrT_CT_PYR23_SCHIZ <- fdrT_CT_PYR23_SCHIZ[ !is.na(fdrT_CT_PYR23_SCHIZ$padj), ]
fdrT_CT_PYR23_SCHIZ <- fdrT_CT_PYR23_SCHIZ[ !is.na(fdrT_CT_PYR23_SCHIZ$pvalue), ]
fdrT_CT_PYR23_SCHIZ <- fdrT_CT_PYR23_SCHIZ[, -which(names(fdrT_CT_PYR23_SCHIZ) == "padj")]
FDR.fdrT_CT_PYR23_SCHIZ <- fdrtool(fdrT_CT_PYR23_SCHIZ$stat, statistic= "normal", plot = T)
fdrT_CT_PYR23_SCHIZ[,"padj"]  <- p.adjust(FDR.fdrT_CT_PYR23_SCHIZ$pval, method = "BH")

fdrT_CT_PYR56_SCHIZ <- res_PYR56_exons_groupCT_age_sex_group_SCHIZ
fdrT_CT_PYR56_SCHIZ <- fdrT_CT_PYR56_SCHIZ[ !is.na(fdrT_CT_PYR56_SCHIZ$padj), ]
fdrT_CT_PYR56_SCHIZ <- fdrT_CT_PYR56_SCHIZ[ !is.na(fdrT_CT_PYR56_SCHIZ$pvalue), ]
fdrT_CT_PYR56_SCHIZ <- fdrT_CT_PYR56_SCHIZ[, -which(names(fdrT_CT_PYR56_SCHIZ) == "padj")]
FDR.fdrT_CT_PYR56_SCHIZ <- fdrtool(fdrT_CT_PYR56_SCHIZ$stat, statistic= "normal", plot = T)
fdrT_CT_PYR56_SCHIZ[,"padj"]  <- p.adjust(FDR.fdrT_CT_PYR56_SCHIZ$pval, method = "BH")

fdrT_CT_SST_SCHIZ <- res_SST_exons_groupCT_age_sex_group_SCHIZ
fdrT_CT_SST_SCHIZ <- fdrT_CT_SST_SCHIZ[ !is.na(fdrT_CT_SST_SCHIZ$padj), ]
fdrT_CT_SST_SCHIZ <- fdrT_CT_SST_SCHIZ[ !is.na(fdrT_CT_SST_SCHIZ$pvalue), ]
fdrT_CT_SST_SCHIZ <- fdrT_CT_SST_SCHIZ[, -which(names(fdrT_CT_SST_SCHIZ) == "padj")]
FDR.fdrT_CT_SST_SCHIZ <- fdrtool(fdrT_CT_SST_SCHIZ$stat, statistic= "normal", plot = T)
fdrT_CT_SST_SCHIZ[,"padj"]  <- p.adjust(FDR.fdrT_CT_SST_SCHIZ$pval, method = "BH")

fdrT_CT_VIP_SCHIZ <- res_VIP_exons_groupCT_age_sex_group_SCHIZ
fdrT_CT_VIP_SCHIZ <- fdrT_CT_VIP_SCHIZ[ !is.na(fdrT_CT_VIP_SCHIZ$padj), ]
fdrT_CT_VIP_SCHIZ <- fdrT_CT_VIP_SCHIZ[ !is.na(fdrT_CT_VIP_SCHIZ$pvalue), ]
fdrT_CT_VIP_SCHIZ <- fdrT_CT_VIP_SCHIZ[, -which(names(fdrT_CT_VIP_SCHIZ) == "padj")]
FDR.fdrT_CT_VIP_SCHIZ <- fdrtool(fdrT_CT_VIP_SCHIZ$stat, statistic= "normal", plot = T)
fdrT_CT_VIP_SCHIZ[,"padj"]  <- p.adjust(FDR.fdrT_CT_VIP_SCHIZ$pval, method = "BH")



###MDD
PV_MDD <- as.data.frame(res_PV_exons_groupCT_age_sex_group_MDD)
PYR23_MDD <- as.data.frame(res_PYR23_exons_groupCT_age_sex_group_MDD)
PYR56_MDD <- as.data.frame(res_PYR56_exons_groupCT_age_sex_group_MDD)
SST_MDD <- as.data.frame(res_SST_exons_groupCT_age_sex_group_MDD)
VIP_MDD <- as.data.frame(res_VIP_exons_groupCT_age_sex_group_MDD)

PV_MDD$Gene.stable.ID <- rownames(PV_MDD)
PYR23_MDD$Gene.stable.ID <- rownames(PYR23_MDD)
PYR56_MDD$Gene.stable.ID <- rownames(PYR56_MDD)
SST_MDD$Gene.stable.ID <- rownames(SST_MDD)
VIP_MDD$Gene.stable.ID <- rownames(VIP_MDD)

PV_MDD <- PV_MDD[,c("log2FoldChange", "Gene.stable.ID")]
PYR23_MDD <- PYR23_MDD[,c("log2FoldChange", "Gene.stable.ID")]
PYR56_MDD <- PYR56_MDD[,c("log2FoldChange", "Gene.stable.ID")]
SST_MDD <- SST_MDD[,c("log2FoldChange", "Gene.stable.ID")]
VIP_MDD <- VIP_MDD[,c("log2FoldChange", "Gene.stable.ID")]

names(PV_MDD) <- c("PV_LFC", "Gene.stable.ID")
names(PYR23_MDD) <- c("PYR23_LFC", "Gene.stable.ID")
names(PYR56_MDD) <- c("PYR56_LFC", "Gene.stable.ID")
names(SST_MDD) <- c("SST_LFC", "Gene.stable.ID")
names(VIP_MDD) <- c("VIP_LFC", "Gene.stable.ID")

merged_MDD <- merge(PV_MDD,merge(PYR23_MDD, merge(PYR56_MDD,merge(SST_MDD, VIP_MDD, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)



rownames(merged_MDD) <- merged_MDD$Gene.stable.ID
merged_MDD <- merged_MDD[,-1]
merged_MDD[is.na(merged_MDD)] <- 0
merged_MDD <- merged_MDD[,c(2,3,1,4,5)]
merged_MDD <- as.matrix(merged_MDD)


pdf("MDD gene heatmap all genes.pdf")
Heatmap(merged_MDD,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()

#Genes detected
library(Vennerable)
library(gplots)
library(venn)


genes_Detected <- list(PV_MDD$Gene.stable.ID, PYR23_MDD$Gene.stable.ID, PYR56_MDD$Gene.stable.ID, SST_MDD$Gene.stable.ID, VIP_MDD$Gene.stable.ID)

par(cex=2)
# pdf("Genes detected.pdf")
venn(genes_Detected, snames=c("PV", "PYR23", "PYR56", "SST", "VIP"), zcolor = c("#7CAE00", "#F8766D", "#b01005", "#00BFC4", "#C77CFF"), par=TRUE, box = FALSE)
# dev.off()


####Bipolar
PV_Bipolar <- as.data.frame(res_PV_exons_groupCT_age_sex_group_Bipolar)
PYR23_Bipolar <- as.data.frame(res_PYR23_exons_groupCT_age_sex_group_Bipolar)
PYR56_Bipolar <- as.data.frame(res_PYR56_exons_groupCT_age_sex_group_Bipolar)
SST_Bipolar <- as.data.frame(res_SST_exons_groupCT_age_sex_group_Bipolar)
VIP_Bipolar <- as.data.frame(res_VIP_exons_groupCT_age_sex_group_Bipolar)

PV_Bipolar$Gene.stable.ID <- rownames(PV_Bipolar)
PYR23_Bipolar$Gene.stable.ID <- rownames(PYR23_Bipolar)
PYR56_Bipolar$Gene.stable.ID <- rownames(PYR56_Bipolar)
SST_Bipolar$Gene.stable.ID <- rownames(SST_Bipolar)
VIP_Bipolar$Gene.stable.ID <- rownames(VIP_Bipolar)

PV_Bipolar <- PV_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
PYR23_Bipolar <- PYR23_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
PYR56_Bipolar <- PYR56_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
SST_Bipolar <- SST_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
VIP_Bipolar <- VIP_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]

names(PV_Bipolar) <- c("PV_LFC", "Gene.stable.ID")
names(PYR23_Bipolar) <- c("PYR23_LFC", "Gene.stable.ID")
names(PYR56_Bipolar) <- c("PYR56_LFC", "Gene.stable.ID")
names(SST_Bipolar) <- c("SST_LFC", "Gene.stable.ID")
names(VIP_Bipolar) <- c("VIP_LFC", "Gene.stable.ID")

merged_Bipolar <- merge(PV_Bipolar,merge(PYR23_Bipolar, merge(PYR56_Bipolar,merge(SST_Bipolar, VIP_Bipolar, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)



rownames(merged_Bipolar) <- merged_Bipolar$Gene.stable.ID
merged_Bipolar <- merged_Bipolar[,-1]
merged_Bipolar[is.na(merged_Bipolar)] <- 0
merged_Bipolar <- merged_Bipolar[,c(2,3,1,4,5)]
merged_Bipolar <- as.matrix(merged_Bipolar)


pdf("Bipolar gene heatmap all genes.pdf")
Heatmap(merged_Bipolar,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


#SCHIZ

PV_SCHIZ <- as.data.frame(res_PV_exons_groupCT_age_sex_group_SCHIZ)
PYR23_SCHIZ <- as.data.frame(res_PYR23_exons_groupCT_age_sex_group_SCHIZ)
PYR56_SCHIZ <- as.data.frame(res_PYR56_exons_groupCT_age_sex_group_SCHIZ)
SST_SCHIZ <- as.data.frame(res_SST_exons_groupCT_age_sex_group_SCHIZ)
VIP_SCHIZ <- as.data.frame(res_VIP_exons_groupCT_age_sex_group_SCHIZ)

PV_SCHIZ$Gene.stable.ID <- rownames(PV_SCHIZ)
PYR23_SCHIZ$Gene.stable.ID <- rownames(PYR23_SCHIZ)
PYR56_SCHIZ$Gene.stable.ID <- rownames(PYR56_SCHIZ)
SST_SCHIZ$Gene.stable.ID <- rownames(SST_SCHIZ)
VIP_SCHIZ$Gene.stable.ID <- rownames(VIP_SCHIZ)

PV_SCHIZ <- PV_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
PYR23_SCHIZ <- PYR23_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
PYR56_SCHIZ <- PYR56_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
SST_SCHIZ <- SST_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
VIP_SCHIZ <- VIP_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]

names(PV_SCHIZ) <- c("PV_LFC", "Gene.stable.ID")
names(PYR23_SCHIZ) <- c("PYR23_LFC", "Gene.stable.ID")
names(PYR56_SCHIZ) <- c("PYR56_LFC", "Gene.stable.ID")
names(SST_SCHIZ) <- c("SST_LFC", "Gene.stable.ID")
names(VIP_SCHIZ) <- c("VIP_LFC", "Gene.stable.ID")

merged_SCHIZ <- merge(PV_SCHIZ,merge(PYR23_SCHIZ, merge(PYR56_SCHIZ,merge(SST_SCHIZ, VIP_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
rownames(merged_SCHIZ) <- merged_SCHIZ$Gene.stable.ID
merged_SCHIZ <- merged_SCHIZ[,-1]
merged_SCHIZ[is.na(merged_SCHIZ)] <- 0
merged_SCHIZ <- merged_SCHIZ[,c(2,3,1,4,5)]
merged_SCHIZ <- as.matrix(merged_SCHIZ)


pdf("SCHIZ gene heatmap all genes.pdf")
Heatmap(merged_SCHIZ,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


#All disorders
merged_MDD <- merge(PV_MDD,merge(PYR23_MDD, merge(PYR56_MDD,merge(SST_MDD, VIP_MDD, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
rownames(merged_MDD) <- merged_MDD$Gene.stable.ID
merged_MDD[is.na(merged_MDD)] <- 0
merged_MDD <- merged_MDD[,c(1,3,4,2,5,6)]

merged_Bipolar <- merge(PV_Bipolar,merge(PYR23_Bipolar, merge(PYR56_Bipolar,merge(SST_Bipolar, VIP_Bipolar, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
rownames(merged_Bipolar) <- merged_Bipolar$Gene.stable.ID
merged_Bipolar[is.na(merged_Bipolar)] <- 0
merged_Bipolar <- merged_Bipolar[,c(1,3,4,2,5,6)]

merged_SCHIZ <- merge(PV_SCHIZ,merge(PYR23_SCHIZ, merge(PYR56_SCHIZ,merge(SST_SCHIZ, VIP_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
rownames(merged_SCHIZ) <- merged_SCHIZ$Gene.stable.ID
merged_SCHIZ[is.na(merged_SCHIZ)] <- 0
merged_SCHIZ <- merged_SCHIZ[,c(1,3,4,2,5,6)]


colnames(merged_MDD)[2:6] <- paste0("MDD_", colnames(merged_MDD)[2:6])
colnames(merged_Bipolar)[2:6] <- paste0("Bipolar_", colnames(merged_Bipolar)[2:6])
colnames(merged_SCHIZ)[2:6] <- paste0("SCHIZ_", colnames(merged_SCHIZ)[2:6])

merged_hmap <- merge(merged_MDD, merge(merged_Bipolar, merged_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
rownames(merged_hmap) <- merged_hmap$Gene.stable.ID
merged_hmap <- merged_hmap[,-1]
merged_hmap <- as.matrix(merged_hmap)

pdf("All_disorders gene heatmap all genes.pdf")
Heatmap(merged_hmap, column_split = c(rep("MDD", 5), rep("BD", 5), rep("SCZ", 5)),
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"), border = TRUE)
dev.off()


#25percent FDR
###MDD
PV_MDD <- as.data.frame(fdrT_CT_PV_MDD)
PYR23_MDD <- as.data.frame(fdrT_CT_PYR23_MDD)
PYR56_MDD <- as.data.frame(fdrT_CT_PYR56_MDD)
SST_MDD <- as.data.frame(fdrT_CT_SST_MDD)
VIP_MDD <- as.data.frame(fdrT_CT_VIP_MDD)

PV_MDD$Gene.stable.ID <- rownames(PV_MDD)
PYR23_MDD$Gene.stable.ID <- rownames(PYR23_MDD)
PYR56_MDD$Gene.stable.ID <- rownames(PYR56_MDD)
SST_MDD$Gene.stable.ID <- rownames(SST_MDD)
VIP_MDD$Gene.stable.ID <- rownames(VIP_MDD)

PV_MDD$log2FoldChange[PV_MDD$padj >= 0.25] <- 0
PYR23_MDD$log2FoldChange[PYR23_MDD$padj >= 0.25] <- 0
PYR56_MDD$log2FoldChange[PYR56_MDD$padj >= 0.25] <- 0
SST_MDD$log2FoldChange[SST_MDD$padj >= 0.25] <- 0
VIP_MDD$log2FoldChange[VIP_MDD$padj >= 0.25] <- 0

PV_MDD <- PV_MDD[,c("log2FoldChange", "Gene.stable.ID")]
PYR23_MDD <- PYR23_MDD[,c("log2FoldChange", "Gene.stable.ID")]
PYR56_MDD <- PYR56_MDD[,c("log2FoldChange", "Gene.stable.ID")]
SST_MDD <- SST_MDD[,c("log2FoldChange", "Gene.stable.ID")]
VIP_MDD <- VIP_MDD[,c("log2FoldChange", "Gene.stable.ID")]

names(PV_MDD) <- c("PV_LFC","Gene.stable.ID")
names(PYR23_MDD) <- c("PYR23_LFC","Gene.stable.ID")
names(PYR56_MDD) <- c("PYR56_LFC", "Gene.stable.ID")
names(SST_MDD) <- c("SST_LFC", "Gene.stable.ID")
names(VIP_MDD) <- c("VIP_LFC", "Gene.stable.ID")

merged_MDD <- merge(PV_MDD,merge(PYR23_MDD, merge(PYR56_MDD,merge(SST_MDD, VIP_MDD, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)



rownames(merged_MDD) <- merged_MDD$Gene.stable.ID
merged_MDD <- merged_MDD[,-1]
merged_MDD[is.na(merged_MDD)] <- 0
merged_MDD <- merged_MDD[!(rowSums(merged_MDD)==0),]
merged_MDD <- merged_MDD[,c(2,3,1,4,5)]
merged_MDD <- as.matrix(merged_MDD)


pdf("MDD gene heatmap 25percent FDR.pdf")
Heatmap(merged_MDD,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


###Bipolar
PV_Bipolar <- as.data.frame(fdrT_CT_PV_Bipolar)
PYR23_Bipolar <- as.data.frame(fdrT_CT_PYR23_Bipolar)
PYR56_Bipolar <- as.data.frame(fdrT_CT_PYR56_Bipolar)
SST_Bipolar <- as.data.frame(fdrT_CT_SST_Bipolar)
VIP_Bipolar <- as.data.frame(fdrT_CT_VIP_Bipolar)

PV_Bipolar$Gene.stable.ID <- rownames(PV_Bipolar)
PYR23_Bipolar$Gene.stable.ID <- rownames(PYR23_Bipolar)
PYR56_Bipolar$Gene.stable.ID <- rownames(PYR56_Bipolar)
SST_Bipolar$Gene.stable.ID <- rownames(SST_Bipolar)
VIP_Bipolar$Gene.stable.ID <- rownames(VIP_Bipolar)

PV_Bipolar$log2FoldChange[PV_Bipolar$padj >= 0.25] <- 0
PYR23_Bipolar$log2FoldChange[PYR23_Bipolar$padj >= 0.25] <- 0
PYR56_Bipolar$log2FoldChange[PYR56_Bipolar$padj >= 0.25] <- 0
SST_Bipolar$log2FoldChange[SST_Bipolar$padj >= 0.25] <- 0
VIP_Bipolar$log2FoldChange[VIP_Bipolar$padj >= 0.25] <- 0

PV_Bipolar <- PV_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
PYR23_Bipolar <- PYR23_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
PYR56_Bipolar <- PYR56_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
SST_Bipolar <- SST_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
VIP_Bipolar <- VIP_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]

names(PV_Bipolar) <- c("PV_LFC","Gene.stable.ID")
names(PYR23_Bipolar) <- c("PYR23_LFC","Gene.stable.ID")
names(PYR56_Bipolar) <- c("PYR56_LFC", "Gene.stable.ID")
names(SST_Bipolar) <- c("SST_LFC", "Gene.stable.ID")
names(VIP_Bipolar) <- c("VIP_LFC", "Gene.stable.ID")

merged_Bipolar <- merge(PV_Bipolar,merge(PYR23_Bipolar, merge(PYR56_Bipolar,merge(SST_Bipolar, VIP_Bipolar, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)



rownames(merged_Bipolar) <- merged_Bipolar$Gene.stable.ID
merged_Bipolar <- merged_Bipolar[,-1]
merged_Bipolar[is.na(merged_Bipolar)] <- 0
merged_Bipolar <- merged_Bipolar[!(rowSums(merged_Bipolar)==0),]
merged_Bipolar <- merged_Bipolar[,c(2,3,1,4,5)]
merged_Bipolar <- as.matrix(merged_Bipolar)


pdf("Bipolar gene heatmap 25percent FDR.pdf")
Heatmap(merged_Bipolar,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


###SCHIZ
PV_SCHIZ <- as.data.frame(fdrT_CT_PV_SCHIZ)
PYR23_SCHIZ <- as.data.frame(fdrT_CT_PYR23_SCHIZ)
PYR56_SCHIZ <- as.data.frame(fdrT_CT_PYR56_SCHIZ)
SST_SCHIZ <- as.data.frame(fdrT_CT_SST_SCHIZ)
VIP_SCHIZ <- as.data.frame(fdrT_CT_VIP_SCHIZ)

PV_SCHIZ$Gene.stable.ID <- rownames(PV_SCHIZ)
PYR23_SCHIZ$Gene.stable.ID <- rownames(PYR23_SCHIZ)
PYR56_SCHIZ$Gene.stable.ID <- rownames(PYR56_SCHIZ)
SST_SCHIZ$Gene.stable.ID <- rownames(SST_SCHIZ)
VIP_SCHIZ$Gene.stable.ID <- rownames(VIP_SCHIZ)

PV_SCHIZ$log2FoldChange[PV_SCHIZ$padj >= 0.25] <- 0
PYR23_SCHIZ$log2FoldChange[PYR23_SCHIZ$padj >= 0.25] <- 0
PYR56_SCHIZ$log2FoldChange[PYR56_SCHIZ$padj >= 0.25] <- 0
SST_SCHIZ$log2FoldChange[SST_SCHIZ$padj >= 0.25] <- 0
VIP_SCHIZ$log2FoldChange[VIP_SCHIZ$padj >= 0.25] <- 0

PV_SCHIZ <- PV_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
PYR23_SCHIZ <- PYR23_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
PYR56_SCHIZ <- PYR56_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
SST_SCHIZ <- SST_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
VIP_SCHIZ <- VIP_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]

names(PV_SCHIZ) <- c("PV_LFC","Gene.stable.ID")
names(PYR23_SCHIZ) <- c("PYR23_LFC","Gene.stable.ID")
names(PYR56_SCHIZ) <- c("PYR56_LFC", "Gene.stable.ID")
names(SST_SCHIZ) <- c("SST_LFC", "Gene.stable.ID")
names(VIP_SCHIZ) <- c("VIP_LFC", "Gene.stable.ID")

merged_SCHIZ <- merge(PV_SCHIZ,merge(PYR23_SCHIZ, merge(PYR56_SCHIZ,merge(SST_SCHIZ, VIP_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)



rownames(merged_SCHIZ) <- merged_SCHIZ$Gene.stable.ID
merged_SCHIZ <- merged_SCHIZ[,-1]
merged_SCHIZ[is.na(merged_SCHIZ)] <- 0
merged_SCHIZ <- merged_SCHIZ[!(rowSums(merged_SCHIZ)==0),]
merged_SCHIZ <- merged_SCHIZ[,c(2,3,1,4,5)]
merged_SCHIZ <- as.matrix(merged_SCHIZ)


pdf("SCHIZ gene heatmap 25percent FDR.pdf")
Heatmap(merged_SCHIZ,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()



#All disorders
merged_MDD <- merge(PV_MDD,merge(PYR23_MDD, merge(PYR56_MDD,merge(SST_MDD, VIP_MDD, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
rownames(merged_MDD) <- merged_MDD$Gene.stable.ID
merged_MDD[is.na(merged_MDD)] <- 0
merged_MDD <- merged_MDD[!(rowSums(merged_MDD[,c(2:6)])==0),]
merged_MDD <- merged_MDD[,c(1,3,4,2,5,6)]

merged_Bipolar <- merge(PV_Bipolar,merge(PYR23_Bipolar, merge(PYR56_Bipolar,merge(SST_Bipolar, VIP_Bipolar, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
rownames(merged_Bipolar) <- merged_Bipolar$Gene.stable.ID
merged_Bipolar[is.na(merged_Bipolar)] <- 0
merged_Bipolar <- merged_Bipolar[!(rowSums(merged_Bipolar[,c(2:6)])==0),]
merged_Bipolar <- merged_Bipolar[,c(1,3,4,2,5,6)]

merged_SCHIZ <- merge(PV_SCHIZ,merge(PYR23_SCHIZ, merge(PYR56_SCHIZ,merge(SST_SCHIZ, VIP_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
rownames(merged_SCHIZ) <- merged_SCHIZ$Gene.stable.ID
merged_SCHIZ[is.na(merged_SCHIZ)] <- 0
merged_SCHIZ <- merged_SCHIZ[!(rowSums(merged_SCHIZ[,c(2:6)])==0),]
merged_SCHIZ <- merged_SCHIZ[,c(1,3,4,2,5,6)]


colnames(merged_MDD)[2:6] <- paste0("MDD_", colnames(merged_MDD)[2:6])
colnames(merged_Bipolar)[2:6] <- paste0("Bipolar_", colnames(merged_Bipolar)[2:6])
colnames(merged_SCHIZ)[2:6] <- paste0("SCHIZ_", colnames(merged_SCHIZ)[2:6])

merged_hmap <- merge(merged_MDD, merge(merged_Bipolar, merged_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
rownames(merged_hmap) <- merged_hmap$Gene.stable.ID
merged_hmap <- merged_hmap[,-1]
merged_hmap[is.na(merged_hmap)] <- 0
merged_hmap <- merged_hmap[!(rowSums(merged_hmap)==0),]
merged_hmap <- as.matrix(merged_hmap)

pdf("All_disorders gene heatmap 25percentFDR.pdf")
Heatmap(merged_hmap, column_split = c(rep("MDD", 5), rep("BD", 5), rep("SCZ", 5)),
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"), border = TRUE)
dev.off()






#5percent FDR
###MDD
PV_MDD <- as.data.frame(fdrT_CT_PV_MDD)
PYR23_MDD <- as.data.frame(fdrT_CT_PYR23_MDD)
PYR56_MDD <- as.data.frame(fdrT_CT_PYR56_MDD)
SST_MDD <- as.data.frame(fdrT_CT_SST_MDD)
VIP_MDD <- as.data.frame(fdrT_CT_VIP_MDD)

PV_MDD$Gene.stable.ID <- rownames(PV_MDD)
PYR23_MDD$Gene.stable.ID <- rownames(PYR23_MDD)
PYR56_MDD$Gene.stable.ID <- rownames(PYR56_MDD)
SST_MDD$Gene.stable.ID <- rownames(SST_MDD)
VIP_MDD$Gene.stable.ID <- rownames(VIP_MDD)

PV_MDD$log2FoldChange[PV_MDD$padj >= 0.05] <- 0
PYR23_MDD$log2FoldChange[PYR23_MDD$padj >= 0.05] <- 0
PYR56_MDD$log2FoldChange[PYR56_MDD$padj >= 0.05] <- 0
SST_MDD$log2FoldChange[SST_MDD$padj >= 0.05] <- 0
VIP_MDD$log2FoldChange[VIP_MDD$padj >= 0.05] <- 0

PV_MDD <- PV_MDD[,c("log2FoldChange", "Gene.stable.ID")]
PYR23_MDD <- PYR23_MDD[,c("log2FoldChange", "Gene.stable.ID")]
PYR56_MDD <- PYR56_MDD[,c("log2FoldChange", "Gene.stable.ID")]
SST_MDD <- SST_MDD[,c("log2FoldChange", "Gene.stable.ID")]
VIP_MDD <- VIP_MDD[,c("log2FoldChange", "Gene.stable.ID")]

names(PV_MDD) <- c("PV_LFC","Gene.stable.ID")
names(PYR23_MDD) <- c("PYR23_LFC","Gene.stable.ID")
names(PYR56_MDD) <- c("PYR56_LFC", "Gene.stable.ID")
names(SST_MDD) <- c("SST_LFC", "Gene.stable.ID")
names(VIP_MDD) <- c("VIP_LFC", "Gene.stable.ID")

merged_MDD <- merge(PV_MDD,merge(PYR23_MDD, merge(PYR56_MDD,merge(SST_MDD, VIP_MDD, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)



rownames(merged_MDD) <- merged_MDD$Gene.stable.ID
merged_MDD <- merged_MDD[,-1]
merged_MDD[is.na(merged_MDD)] <- 0
merged_MDD <- merged_MDD[!(rowSums(merged_MDD)==0),]
merged_MDD <- merged_MDD[,c(2,3,1,4,5)]
merged_MDD <- as.matrix(merged_MDD)


pdf("MDD gene heatmap 5percent FDR.pdf")
Heatmap(merged_MDD,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


###Bipolar
PV_Bipolar <- as.data.frame(fdrT_CT_PV_Bipolar)
PYR23_Bipolar <- as.data.frame(fdrT_CT_PYR23_Bipolar)
PYR56_Bipolar <- as.data.frame(fdrT_CT_PYR56_Bipolar)
SST_Bipolar <- as.data.frame(fdrT_CT_SST_Bipolar)
VIP_Bipolar <- as.data.frame(fdrT_CT_VIP_Bipolar)

PV_Bipolar$Gene.stable.ID <- rownames(PV_Bipolar)
PYR23_Bipolar$Gene.stable.ID <- rownames(PYR23_Bipolar)
PYR56_Bipolar$Gene.stable.ID <- rownames(PYR56_Bipolar)
SST_Bipolar$Gene.stable.ID <- rownames(SST_Bipolar)
VIP_Bipolar$Gene.stable.ID <- rownames(VIP_Bipolar)

PV_Bipolar$log2FoldChange[PV_Bipolar$padj >= 0.05] <- 0
PYR23_Bipolar$log2FoldChange[PYR23_Bipolar$padj >= 0.05] <- 0
PYR56_Bipolar$log2FoldChange[PYR56_Bipolar$padj >= 0.05] <- 0
SST_Bipolar$log2FoldChange[SST_Bipolar$padj >= 0.05] <- 0
VIP_Bipolar$log2FoldChange[VIP_Bipolar$padj >= 0.05] <- 0

PV_Bipolar <- PV_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
PYR23_Bipolar <- PYR23_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
PYR56_Bipolar <- PYR56_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
SST_Bipolar <- SST_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
VIP_Bipolar <- VIP_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]

names(PV_Bipolar) <- c("PV_LFC","Gene.stable.ID")
names(PYR23_Bipolar) <- c("PYR23_LFC","Gene.stable.ID")
names(PYR56_Bipolar) <- c("PYR56_LFC", "Gene.stable.ID")
names(SST_Bipolar) <- c("SST_LFC", "Gene.stable.ID")
names(VIP_Bipolar) <- c("VIP_LFC", "Gene.stable.ID")

merged_Bipolar <- merge(PV_Bipolar,merge(PYR23_Bipolar, merge(PYR56_Bipolar,merge(SST_Bipolar, VIP_Bipolar, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)



rownames(merged_Bipolar) <- merged_Bipolar$Gene.stable.ID
merged_Bipolar <- merged_Bipolar[,-1]
merged_Bipolar[is.na(merged_Bipolar)] <- 0
merged_Bipolar <- merged_Bipolar[!(rowSums(merged_Bipolar)==0),]
merged_Bipolar <- merged_Bipolar[,c(2,3,1,4,5)]
merged_Bipolar <- as.matrix(merged_Bipolar)


pdf("Bipolar gene heatmap 5percent FDR.pdf")
Heatmap(merged_Bipolar,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


###SCHIZ
PV_SCHIZ <- as.data.frame(fdrT_CT_PV_SCHIZ)
PYR23_SCHIZ <- as.data.frame(fdrT_CT_PYR23_SCHIZ)
PYR56_SCHIZ <- as.data.frame(fdrT_CT_PYR56_SCHIZ)
SST_SCHIZ <- as.data.frame(fdrT_CT_SST_SCHIZ)
VIP_SCHIZ <- as.data.frame(fdrT_CT_VIP_SCHIZ)

PV_SCHIZ$Gene.stable.ID <- rownames(PV_SCHIZ)
PYR23_SCHIZ$Gene.stable.ID <- rownames(PYR23_SCHIZ)
PYR56_SCHIZ$Gene.stable.ID <- rownames(PYR56_SCHIZ)
SST_SCHIZ$Gene.stable.ID <- rownames(SST_SCHIZ)
VIP_SCHIZ$Gene.stable.ID <- rownames(VIP_SCHIZ)

PV_SCHIZ$log2FoldChange[PV_SCHIZ$padj >= 0.05] <- 0
PYR23_SCHIZ$log2FoldChange[PYR23_SCHIZ$padj >= 0.05] <- 0
PYR56_SCHIZ$log2FoldChange[PYR56_SCHIZ$padj >= 0.05] <- 0
SST_SCHIZ$log2FoldChange[SST_SCHIZ$padj >= 0.05] <- 0
VIP_SCHIZ$log2FoldChange[VIP_SCHIZ$padj >= 0.05] <- 0

PV_SCHIZ <- PV_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
PYR23_SCHIZ <- PYR23_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
PYR56_SCHIZ <- PYR56_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
SST_SCHIZ <- SST_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
VIP_SCHIZ <- VIP_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]

names(PV_SCHIZ) <- c("PV_LFC","Gene.stable.ID")
names(PYR23_SCHIZ) <- c("PYR23_LFC","Gene.stable.ID")
names(PYR56_SCHIZ) <- c("PYR56_LFC", "Gene.stable.ID")
names(SST_SCHIZ) <- c("SST_LFC", "Gene.stable.ID")
names(VIP_SCHIZ) <- c("VIP_LFC", "Gene.stable.ID")

merged_SCHIZ <- merge(PV_SCHIZ,merge(PYR23_SCHIZ, merge(PYR56_SCHIZ,merge(SST_SCHIZ, VIP_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)



rownames(merged_SCHIZ) <- merged_SCHIZ$Gene.stable.ID
merged_SCHIZ <- merged_SCHIZ[,-1]
merged_SCHIZ[is.na(merged_SCHIZ)] <- 0
merged_SCHIZ <- merged_SCHIZ[!(rowSums(merged_SCHIZ)==0),]
merged_SCHIZ <- merged_SCHIZ[,c(2,3,1,4,5)]
merged_SCHIZ <- as.matrix(merged_SCHIZ)


pdf("SCHIZ gene heatmap 5percent FDR.pdf")
Heatmap(merged_SCHIZ,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()

#All disorders
merged_MDD <- merge(PV_MDD,merge(PYR23_MDD, merge(PYR56_MDD,merge(SST_MDD, VIP_MDD, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
rownames(merged_MDD) <- merged_MDD$Gene.stable.ID
merged_MDD[is.na(merged_MDD)] <- 0
merged_MDD <- merged_MDD[!(rowSums(merged_MDD[,c(2:6)])==0),]
merged_MDD <- merged_MDD[,c(1,3,4,2,5,6)]

merged_Bipolar <- merge(PV_Bipolar,merge(PYR23_Bipolar, merge(PYR56_Bipolar,merge(SST_Bipolar, VIP_Bipolar, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
rownames(merged_Bipolar) <- merged_Bipolar$Gene.stable.ID
merged_Bipolar[is.na(merged_Bipolar)] <- 0
merged_Bipolar <- merged_Bipolar[!(rowSums(merged_Bipolar[,c(2:6)])==0),]
merged_Bipolar <- merged_Bipolar[,c(1,3,4,2,5,6)]

merged_SCHIZ <- merge(PV_SCHIZ,merge(PYR23_SCHIZ, merge(PYR56_SCHIZ,merge(SST_SCHIZ, VIP_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
rownames(merged_SCHIZ) <- merged_SCHIZ$Gene.stable.ID
merged_SCHIZ[is.na(merged_SCHIZ)] <- 0
merged_SCHIZ <- merged_SCHIZ[!(rowSums(merged_SCHIZ[,c(2:6)])==0),]
merged_SCHIZ <- merged_SCHIZ[,c(1,3,4,2,5,6)]


colnames(merged_MDD)[2:6] <- paste0("MDD_", colnames(merged_MDD)[2:6])
colnames(merged_Bipolar)[2:6] <- paste0("Bipolar_", colnames(merged_Bipolar)[2:6])
colnames(merged_SCHIZ)[2:6] <- paste0("SCHIZ_", colnames(merged_SCHIZ)[2:6])

merged_hmap <- merge(merged_MDD, merge(merged_Bipolar, merged_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
rownames(merged_hmap) <- merged_hmap$Gene.stable.ID
merged_hmap <- merged_hmap[,-1]
merged_hmap[is.na(merged_hmap)] <- 0
merged_hmap <- merged_hmap[!(rowSums(merged_hmap)==0),]
merged_hmap <- as.matrix(merged_hmap)

pdf("All_disorders gene heatmap 5percentFDR.pdf")
Heatmap(merged_hmap, column_split = c(rep("MDD", 5), rep("BD", 5), rep("SCZ", 5)),
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"), border = TRUE)
dev.off()





###Within cell-types, across disorders
#All genes
#PV
PV_MDD <- as.data.frame(res_PV_exons_groupCT_age_sex_group_MDD)
PV_MDD$Gene.stable.ID <- rownames(PV_MDD)
PV_MDD <- PV_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(PV_MDD) <- c("MDD_LFC", "Gene.stable.ID")

PV_Bipolar <- as.data.frame(res_PV_exons_groupCT_age_sex_group_Bipolar)
PV_Bipolar$Gene.stable.ID <- rownames(PV_Bipolar)
PV_Bipolar <- PV_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(PV_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

PV_SCHIZ <- as.data.frame(res_PV_exons_groupCT_age_sex_group_SCHIZ)
PV_SCHIZ$Gene.stable.ID <- rownames(PV_SCHIZ)
PV_SCHIZ <- PV_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(PV_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_PV <- merge(PV_MDD, merge(PV_Bipolar, PV_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_PV) <- merged_PV$Gene.stable.ID
merged_PV <- merged_PV[,-1]
merged_PV[is.na(merged_PV)] <- 0
merged_PV <- as.matrix(merged_PV)

pdf("PV gene heatmap all genes.pdf")
Heatmap(merged_PV,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


#PYR23
PYR23_MDD <- as.data.frame(res_PYR23_exons_groupCT_age_sex_group_MDD)
PYR23_MDD$Gene.stable.ID <- rownames(PYR23_MDD)
PYR23_MDD <- PYR23_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR23_MDD) <- c("MDD_LFC", "Gene.stable.ID")

PYR23_Bipolar <- as.data.frame(res_PYR23_exons_groupCT_age_sex_group_Bipolar)
PYR23_Bipolar$Gene.stable.ID <- rownames(PYR23_Bipolar)
PYR23_Bipolar <- PYR23_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR23_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

PYR23_SCHIZ <- as.data.frame(res_PYR23_exons_groupCT_age_sex_group_SCHIZ)
PYR23_SCHIZ$Gene.stable.ID <- rownames(PYR23_SCHIZ)
PYR23_SCHIZ <- PYR23_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR23_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_PYR23 <- merge(PYR23_MDD, merge(PYR23_Bipolar, PYR23_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_PYR23) <- merged_PYR23$Gene.stable.ID
merged_PYR23 <- merged_PYR23[,-1]
merged_PYR23[is.na(merged_PYR23)] <- 0
merged_PYR23 <- as.matrix(merged_PYR23)

pdf("PYR23 gene heatmap all genes.pdf")
Heatmap(merged_PYR23,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


#PYR56
PYR56_MDD <- as.data.frame(res_PYR56_exons_groupCT_age_sex_group_MDD)
PYR56_MDD$Gene.stable.ID <- rownames(PYR56_MDD)
PYR56_MDD <- PYR56_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR56_MDD) <- c("MDD_LFC", "Gene.stable.ID")

PYR56_Bipolar <- as.data.frame(res_PYR56_exons_groupCT_age_sex_group_Bipolar)
PYR56_Bipolar$Gene.stable.ID <- rownames(PYR56_Bipolar)
PYR56_Bipolar <- PYR56_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR56_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

PYR56_SCHIZ <- as.data.frame(res_PYR56_exons_groupCT_age_sex_group_SCHIZ)
PYR56_SCHIZ$Gene.stable.ID <- rownames(PYR56_SCHIZ)
PYR56_SCHIZ <- PYR56_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR56_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_PYR56 <- merge(PYR56_MDD, merge(PYR56_Bipolar, PYR56_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_PYR56) <- merged_PYR56$Gene.stable.ID
merged_PYR56 <- merged_PYR56[,-1]
merged_PYR56[is.na(merged_PYR56)] <- 0
merged_PYR56 <- as.matrix(merged_PYR56)

pdf("PYR56 gene heatmap all genes.pdf")
Heatmap(merged_PYR56,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


#SST
SST_MDD <- as.data.frame(res_SST_exons_groupCT_age_sex_group_MDD)
SST_MDD$Gene.stable.ID <- rownames(SST_MDD)
SST_MDD <- SST_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(SST_MDD) <- c("MDD_LFC", "Gene.stable.ID")

SST_Bipolar <- as.data.frame(res_SST_exons_groupCT_age_sex_group_Bipolar)
SST_Bipolar$Gene.stable.ID <- rownames(SST_Bipolar)
SST_Bipolar <- SST_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(SST_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

SST_SCHIZ <- as.data.frame(res_SST_exons_groupCT_age_sex_group_SCHIZ)
SST_SCHIZ$Gene.stable.ID <- rownames(SST_SCHIZ)
SST_SCHIZ <- SST_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(SST_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_SST <- merge(SST_MDD, merge(SST_Bipolar, SST_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_SST) <- merged_SST$Gene.stable.ID
merged_SST <- merged_SST[,-1]
merged_SST[is.na(merged_SST)] <- 0
merged_SST <- as.matrix(merged_SST)

pdf("SST gene heatmap all genes.pdf")
Heatmap(merged_SST,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


#VIP
VIP_MDD <- as.data.frame(res_VIP_exons_groupCT_age_sex_group_MDD)
VIP_MDD$Gene.stable.ID <- rownames(VIP_MDD)
VIP_MDD <- VIP_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(VIP_MDD) <- c("MDD_LFC", "Gene.stable.ID")

VIP_Bipolar <- as.data.frame(res_VIP_exons_groupCT_age_sex_group_Bipolar)
VIP_Bipolar$Gene.stable.ID <- rownames(VIP_Bipolar)
VIP_Bipolar <- VIP_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(VIP_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

VIP_SCHIZ <- as.data.frame(res_VIP_exons_groupCT_age_sex_group_SCHIZ)
VIP_SCHIZ$Gene.stable.ID <- rownames(VIP_SCHIZ)
VIP_SCHIZ <- VIP_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(VIP_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_VIP <- merge(VIP_MDD, merge(VIP_Bipolar, VIP_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_VIP) <- merged_VIP$Gene.stable.ID
merged_VIP <- merged_VIP[,-1]
merged_VIP[is.na(merged_VIP)] <- 0
merged_VIP <- as.matrix(merged_VIP)

pdf("VIP gene heatmap all genes.pdf")
Heatmap(merged_VIP,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()






#All cell.types
merged_PV <- as.data.frame(merged_PV)
merged_PV$Gene.stable.ID <- rownames(merged_PV)
colnames(merged_PV)[1:3] <- paste0("PV_", colnames(merged_PV)[1:3])

merged_PYR23 <- as.data.frame(merged_PYR23)
merged_PYR23$Gene.stable.ID <- rownames(merged_PYR23)
colnames(merged_PYR23)[1:3] <- paste0("PYR23_", colnames(merged_PYR23)[1:3])

merged_PYR56 <- as.data.frame(merged_PYR56)
merged_PYR56$Gene.stable.ID <- rownames(merged_PYR56)
colnames(merged_PYR56)[1:3] <- paste0("PYR56_", colnames(merged_PYR56)[1:3])

merged_SST <- as.data.frame(merged_SST)
merged_SST$Gene.stable.ID <- rownames(merged_SST)
colnames(merged_SST)[1:3] <- paste0("SST_", colnames(merged_SST)[1:3])

merged_VIP <- as.data.frame(merged_VIP)
merged_VIP$Gene.stable.ID <- rownames(merged_VIP)
colnames(merged_VIP)[1:3] <- paste0("VIP_", colnames(merged_VIP)[1:3])


merged_hmap <- merge(merged_PYR23,merge(merged_PYR56, merge(merged_PV, merge(merged_SST, merged_VIP, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
merged_hmap <- merged_hmap[,-1]
merged_hmap[is.na(merged_hmap)] <- 0
merged_hmap <- merged_hmap[!(rowSums(merged_hmap)==0),]
merged_hmap <- as.matrix(merged_hmap)

pdf("All_genes heatmap all genes.pdf")
Heatmap(merged_hmap, column_split = c(rep("PYR L2/3", 3), rep("PYR L5/6", 3), rep("PV", 3), rep("SST", 3), rep("VIP", 3)),
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"), border = TRUE)
dev.off()













#########25percentFDR
#All genes
#PV
PV_MDD <- as.data.frame(fdrT_CT_PV_MDD)
PV_MDD$Gene.stable.ID <- rownames(PV_MDD)
PV_MDD$log2FoldChange[PV_MDD$padj >= 0.25] <- 0
PV_MDD <- PV_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(PV_MDD) <- c("MDD_LFC", "Gene.stable.ID")

PV_Bipolar <- as.data.frame(fdrT_CT_PV_Bipolar)
PV_Bipolar$Gene.stable.ID <- rownames(PV_Bipolar)
PV_Bipolar$log2FoldChange[PV_Bipolar$padj >= 0.25] <- 0
PV_Bipolar <- PV_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(PV_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

PV_SCHIZ <- as.data.frame(fdrT_CT_PV_SCHIZ)
PV_SCHIZ$Gene.stable.ID <- rownames(PV_SCHIZ)
PV_SCHIZ$log2FoldChange[PV_SCHIZ$padj >= 0.25] <- 0
PV_SCHIZ <- PV_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(PV_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_PV <- merge(PV_MDD, merge(PV_Bipolar, PV_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_PV) <- merged_PV$Gene.stable.ID
merged_PV <- merged_PV[,-1]
merged_PV[is.na(merged_PV)] <- 0
merged_PV <- merged_PV[!(rowSums(merged_PV)==0),]
merged_PV <- as.matrix(merged_PV)

pdf("PV gene heatmap 25percent FDR.pdf")
Heatmap(merged_PV,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


#PYR23
PYR23_MDD <- as.data.frame(fdrT_CT_PYR23_MDD)
PYR23_MDD$Gene.stable.ID <- rownames(PYR23_MDD)
PYR23_MDD$log2FoldChange[PYR23_MDD$padj >= 0.25] <- 0
PYR23_MDD <- PYR23_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR23_MDD) <- c("MDD_LFC", "Gene.stable.ID")

PYR23_Bipolar <- as.data.frame(fdrT_CT_PYR23_Bipolar)
PYR23_Bipolar$Gene.stable.ID <- rownames(PYR23_Bipolar)
PYR23_Bipolar$log2FoldChange[PYR23_Bipolar$padj >= 0.25] <- 0
PYR23_Bipolar <- PYR23_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR23_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

PYR23_SCHIZ <- as.data.frame(fdrT_CT_PYR23_SCHIZ)
PYR23_SCHIZ$Gene.stable.ID <- rownames(PYR23_SCHIZ)
PYR23_SCHIZ$log2FoldChange[PYR23_SCHIZ$padj >= 0.25] <- 0
PYR23_SCHIZ <- PYR23_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR23_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_PYR23 <- merge(PYR23_MDD, merge(PYR23_Bipolar, PYR23_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_PYR23) <- merged_PYR23$Gene.stable.ID
merged_PYR23 <- merged_PYR23[,-1]
merged_PYR23[is.na(merged_PYR23)] <- 0
merged_PYR23 <- merged_PYR23[!(rowSums(merged_PYR23)==0),]
merged_PYR23 <- as.matrix(merged_PYR23)

pdf("PYR23 gene heatmap 25percent FDR.pdf")
Heatmap(merged_PYR23,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()

#PYR56
PYR56_MDD <- as.data.frame(fdrT_CT_PYR56_MDD)
PYR56_MDD$Gene.stable.ID <- rownames(PYR56_MDD)
PYR56_MDD$log2FoldChange[PYR56_MDD$padj >= 0.25] <- 0
PYR56_MDD <- PYR56_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR56_MDD) <- c("MDD_LFC", "Gene.stable.ID")

PYR56_Bipolar <- as.data.frame(fdrT_CT_PYR56_Bipolar)
PYR56_Bipolar$Gene.stable.ID <- rownames(PYR56_Bipolar)
PYR56_Bipolar$log2FoldChange[PYR56_Bipolar$padj >= 0.25] <- 0
PYR56_Bipolar <- PYR56_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR56_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

PYR56_SCHIZ <- as.data.frame(fdrT_CT_PYR56_SCHIZ)
PYR56_SCHIZ$Gene.stable.ID <- rownames(PYR56_SCHIZ)
PYR56_SCHIZ$log2FoldChange[PYR56_SCHIZ$padj >= 0.25] <- 0
PYR56_SCHIZ <- PYR56_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR56_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_PYR56 <- merge(PYR56_MDD, merge(PYR56_Bipolar, PYR56_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_PYR56) <- merged_PYR56$Gene.stable.ID
merged_PYR56 <- merged_PYR56[,-1]
merged_PYR56[is.na(merged_PYR56)] <- 0
merged_PYR56 <- merged_PYR56[!(rowSums(merged_PYR56)==0),]
merged_PYR56 <- as.matrix(merged_PYR56)

pdf("PYR56 gene heatmap 25percent FDR.pdf")
Heatmap(merged_PYR56,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


#SST
SST_MDD <- as.data.frame(fdrT_CT_SST_MDD)
SST_MDD$Gene.stable.ID <- rownames(SST_MDD)
SST_MDD$log2FoldChange[SST_MDD$padj >= 0.25] <- 0
SST_MDD <- SST_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(SST_MDD) <- c("MDD_LFC", "Gene.stable.ID")

SST_Bipolar <- as.data.frame(fdrT_CT_SST_Bipolar)
SST_Bipolar$Gene.stable.ID <- rownames(SST_Bipolar)
SST_Bipolar$log2FoldChange[SST_Bipolar$padj >= 0.25] <- 0
SST_Bipolar <- SST_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(SST_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

SST_SCHIZ <- as.data.frame(fdrT_CT_SST_SCHIZ)
SST_SCHIZ$Gene.stable.ID <- rownames(SST_SCHIZ)
SST_SCHIZ$log2FoldChange[SST_SCHIZ$padj >= 0.25] <- 0
SST_SCHIZ <- SST_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(SST_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_SST <- merge(SST_MDD, merge(SST_Bipolar, SST_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_SST) <- merged_SST$Gene.stable.ID
merged_SST <- merged_SST[,-1]
merged_SST[is.na(merged_SST)] <- 0
merged_SST <- merged_SST[!(rowSums(merged_SST)==0),]
merged_SST <- as.matrix(merged_SST)

pdf("SST gene heatmap 25percent FDR.pdf")
Heatmap(merged_SST,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


#VIP
VIP_MDD <- as.data.frame(fdrT_CT_VIP_MDD)
VIP_MDD$Gene.stable.ID <- rownames(VIP_MDD)
VIP_MDD$log2FoldChange[VIP_MDD$padj >= 0.25] <- 0
VIP_MDD <- VIP_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(VIP_MDD) <- c("MDD_LFC", "Gene.stable.ID")

VIP_Bipolar <- as.data.frame(fdrT_CT_VIP_Bipolar)
VIP_Bipolar$Gene.stable.ID <- rownames(VIP_Bipolar)
VIP_Bipolar$log2FoldChange[VIP_Bipolar$padj >= 0.25] <- 0
VIP_Bipolar <- VIP_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(VIP_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

VIP_SCHIZ <- as.data.frame(fdrT_CT_VIP_SCHIZ)
VIP_SCHIZ$Gene.stable.ID <- rownames(VIP_SCHIZ)
VIP_SCHIZ$log2FoldChange[VIP_SCHIZ$padj >= 0.25] <- 0
VIP_SCHIZ <- VIP_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(VIP_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_VIP <- merge(VIP_MDD, merge(VIP_Bipolar, VIP_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_VIP) <- merged_VIP$Gene.stable.ID
merged_VIP <- merged_VIP[,-1]
merged_VIP[is.na(merged_VIP)] <- 0
merged_VIP <- merged_VIP[!(rowSums(merged_VIP)==0),]
merged_VIP <- as.matrix(merged_VIP)

pdf("VIP gene heatmap 25percent FDR.pdf")
Heatmap(merged_VIP,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()



#All cell.types
merged_PV <- as.data.frame(merged_PV)
merged_PV$Gene.stable.ID <- rownames(merged_PV)
colnames(merged_PV)[1:3] <- paste0("PV_", colnames(merged_PV)[1:3])

merged_PYR23 <- as.data.frame(merged_PYR23)
merged_PYR23$Gene.stable.ID <- rownames(merged_PYR23)
colnames(merged_PYR23)[1:3] <- paste0("PYR23_", colnames(merged_PYR23)[1:3])

merged_PYR56 <- as.data.frame(merged_PYR56)
merged_PYR56$Gene.stable.ID <- rownames(merged_PYR56)
colnames(merged_PYR56)[1:3] <- paste0("PYR56_", colnames(merged_PYR56)[1:3])

merged_SST <- as.data.frame(merged_SST)
merged_SST$Gene.stable.ID <- rownames(merged_SST)
colnames(merged_SST)[1:3] <- paste0("SST_", colnames(merged_SST)[1:3])

merged_VIP <- as.data.frame(merged_VIP)
merged_VIP$Gene.stable.ID <- rownames(merged_VIP)
colnames(merged_VIP)[1:3] <- paste0("VIP_", colnames(merged_VIP)[1:3])


merged_hmap <- merge(merged_PYR23,merge(merged_PYR56, merge(merged_PV, merge(merged_SST, merged_VIP, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
merged_hmap <- merged_hmap[,-1]
merged_hmap[is.na(merged_hmap)] <- 0
merged_hmap <- merged_hmap[!(rowSums(merged_hmap)==0),]
merged_hmap <- as.matrix(merged_hmap)

pdf("All_genes heatmap 25percent FDR.pdf")
Heatmap(merged_hmap, column_split = c(rep("PYR L2/3", 3), rep("PYR L5/6", 3), rep("PV", 3), rep("SST", 3), rep("VIP", 3)),
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"), border = TRUE)
dev.off()






#########5percentFDR
#All genes
#PV
PV_MDD <- as.data.frame(fdrT_CT_PV_MDD)
PV_MDD$Gene.stable.ID <- rownames(PV_MDD)
PV_MDD$log2FoldChange[PV_MDD$padj >= 0.05] <- 0
PV_MDD <- PV_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(PV_MDD) <- c("MDD_LFC", "Gene.stable.ID")

PV_Bipolar <- as.data.frame(fdrT_CT_PV_Bipolar)
PV_Bipolar$Gene.stable.ID <- rownames(PV_Bipolar)
PV_Bipolar$log2FoldChange[PV_Bipolar$padj >= 0.05] <- 0
PV_Bipolar <- PV_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(PV_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

PV_SCHIZ <- as.data.frame(fdrT_CT_PV_SCHIZ)
PV_SCHIZ$Gene.stable.ID <- rownames(PV_SCHIZ)
PV_SCHIZ$log2FoldChange[PV_SCHIZ$padj >= 0.05] <- 0
PV_SCHIZ <- PV_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(PV_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_PV <- merge(PV_MDD, merge(PV_Bipolar, PV_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_PV) <- merged_PV$Gene.stable.ID
merged_PV <- merged_PV[,-1]
merged_PV[is.na(merged_PV)] <- 0
merged_PV <- merged_PV[!(rowSums(merged_PV)==0),]
merged_PV <- as.matrix(merged_PV)

pdf("PV gene heatmap 5percent FDR.pdf")
Heatmap(merged_PV,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


#PYR23
PYR23_MDD <- as.data.frame(fdrT_CT_PYR23_MDD)
PYR23_MDD$Gene.stable.ID <- rownames(PYR23_MDD)
PYR23_MDD$log2FoldChange[PYR23_MDD$padj >= 0.05] <- 0
PYR23_MDD <- PYR23_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR23_MDD) <- c("MDD_LFC", "Gene.stable.ID")

PYR23_Bipolar <- as.data.frame(fdrT_CT_PYR23_Bipolar)
PYR23_Bipolar$Gene.stable.ID <- rownames(PYR23_Bipolar)
PYR23_Bipolar$log2FoldChange[PYR23_Bipolar$padj >= 0.05] <- 0
PYR23_Bipolar <- PYR23_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR23_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

PYR23_SCHIZ <- as.data.frame(fdrT_CT_PYR23_SCHIZ)
PYR23_SCHIZ$Gene.stable.ID <- rownames(PYR23_SCHIZ)
PYR23_SCHIZ$log2FoldChange[PYR23_SCHIZ$padj >= 0.05] <- 0
PYR23_SCHIZ <- PYR23_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR23_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_PYR23 <- merge(PYR23_MDD, merge(PYR23_Bipolar, PYR23_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_PYR23) <- merged_PYR23$Gene.stable.ID
merged_PYR23 <- merged_PYR23[,-1]
merged_PYR23[is.na(merged_PYR23)] <- 0
merged_PYR23 <- merged_PYR23[!(rowSums(merged_PYR23)==0),]
merged_PYR23 <- as.matrix(merged_PYR23)

pdf("PYR23 gene heatmap 5percent FDR.pdf")
Heatmap(merged_PYR23,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()

#PYR56
PYR56_MDD <- as.data.frame(fdrT_CT_PYR56_MDD)
PYR56_MDD$Gene.stable.ID <- rownames(PYR56_MDD)
PYR56_MDD$log2FoldChange[PYR56_MDD$padj >= 0.05] <- 0
PYR56_MDD <- PYR56_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR56_MDD) <- c("MDD_LFC", "Gene.stable.ID")

PYR56_Bipolar <- as.data.frame(fdrT_CT_PYR56_Bipolar)
PYR56_Bipolar$Gene.stable.ID <- rownames(PYR56_Bipolar)
PYR56_Bipolar$log2FoldChange[PYR56_Bipolar$padj >= 0.05] <- 0
PYR56_Bipolar <- PYR56_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR56_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

PYR56_SCHIZ <- as.data.frame(fdrT_CT_PYR56_SCHIZ)
PYR56_SCHIZ$Gene.stable.ID <- rownames(PYR56_SCHIZ)
PYR56_SCHIZ$log2FoldChange[PYR56_SCHIZ$padj >= 0.05] <- 0
PYR56_SCHIZ <- PYR56_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(PYR56_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_PYR56 <- merge(PYR56_MDD, merge(PYR56_Bipolar, PYR56_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_PYR56) <- merged_PYR56$Gene.stable.ID
merged_PYR56 <- merged_PYR56[,-1]
merged_PYR56[is.na(merged_PYR56)] <- 0
merged_PYR56 <- merged_PYR56[!(rowSums(merged_PYR56)==0),]
merged_PYR56 <- as.matrix(merged_PYR56)

pdf("PYR56 gene heatmap 5percent FDR.pdf")
Heatmap(merged_PYR56,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


#SST
SST_MDD <- as.data.frame(fdrT_CT_SST_MDD)
SST_MDD$Gene.stable.ID <- rownames(SST_MDD)
SST_MDD$log2FoldChange[SST_MDD$padj >= 0.05] <- 0
SST_MDD <- SST_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(SST_MDD) <- c("MDD_LFC", "Gene.stable.ID")

SST_Bipolar <- as.data.frame(fdrT_CT_SST_Bipolar)
SST_Bipolar$Gene.stable.ID <- rownames(SST_Bipolar)
SST_Bipolar$log2FoldChange[SST_Bipolar$padj >= 0.05] <- 0
SST_Bipolar <- SST_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(SST_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

SST_SCHIZ <- as.data.frame(fdrT_CT_SST_SCHIZ)
SST_SCHIZ$Gene.stable.ID <- rownames(SST_SCHIZ)
SST_SCHIZ$log2FoldChange[SST_SCHIZ$padj >= 0.05] <- 0
SST_SCHIZ <- SST_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(SST_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_SST <- merge(SST_MDD, merge(SST_Bipolar, SST_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_SST) <- merged_SST$Gene.stable.ID
merged_SST <- merged_SST[,-1]
merged_SST[is.na(merged_SST)] <- 0
merged_SST <- merged_SST[!(rowSums(merged_SST)==0),]
merged_SST <- as.matrix(merged_SST)

pdf("SST gene heatmap 5percent FDR.pdf")
Heatmap(merged_SST,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()


#VIP
VIP_MDD <- as.data.frame(fdrT_CT_VIP_MDD)
VIP_MDD$Gene.stable.ID <- rownames(VIP_MDD)
VIP_MDD$log2FoldChange[VIP_MDD$padj >= 0.05] <- 0
VIP_MDD <- VIP_MDD[,c("log2FoldChange", "Gene.stable.ID")]
names(VIP_MDD) <- c("MDD_LFC", "Gene.stable.ID")

VIP_Bipolar <- as.data.frame(fdrT_CT_VIP_Bipolar)
VIP_Bipolar$Gene.stable.ID <- rownames(VIP_Bipolar)
VIP_Bipolar$log2FoldChange[VIP_Bipolar$padj >= 0.05] <- 0
VIP_Bipolar <- VIP_Bipolar[,c("log2FoldChange", "Gene.stable.ID")]
names(VIP_Bipolar) <- c("BD_LFC", "Gene.stable.ID")

VIP_SCHIZ <- as.data.frame(fdrT_CT_VIP_SCHIZ)
VIP_SCHIZ$Gene.stable.ID <- rownames(VIP_SCHIZ)
VIP_SCHIZ$log2FoldChange[VIP_SCHIZ$padj >= 0.05] <- 0
VIP_SCHIZ <- VIP_SCHIZ[,c("log2FoldChange", "Gene.stable.ID")]
names(VIP_SCHIZ) <- c("SCZ_LFC", "Gene.stable.ID")


merged_VIP <- merge(VIP_MDD, merge(VIP_Bipolar, VIP_SCHIZ, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)

rownames(merged_VIP) <- merged_VIP$Gene.stable.ID
merged_VIP <- merged_VIP[,-1]
merged_VIP[is.na(merged_VIP)] <- 0
merged_VIP <- merged_VIP[!(rowSums(merged_VIP)==0),]
merged_VIP <- as.matrix(merged_VIP)

pdf("VIP gene heatmap 5percent FDR.pdf")
Heatmap(merged_VIP,
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"))
dev.off()

#All cell.types
merged_PV <- as.data.frame(merged_PV)
merged_PV$Gene.stable.ID <- rownames(merged_PV)
colnames(merged_PV)[1:3] <- paste0("PV_", colnames(merged_PV)[1:3])

merged_PYR23 <- as.data.frame(merged_PYR23)
merged_PYR23$Gene.stable.ID <- rownames(merged_PYR23)
colnames(merged_PYR23)[1:3] <- paste0("PYR23_", colnames(merged_PYR23)[1:3])

merged_PYR56 <- as.data.frame(merged_PYR56)
merged_PYR56$Gene.stable.ID <- rownames(merged_PYR56)
colnames(merged_PYR56)[1:3] <- paste0("PYR56_", colnames(merged_PYR56)[1:3])

merged_SST <- as.data.frame(merged_SST)
merged_SST$Gene.stable.ID <- rownames(merged_SST)
colnames(merged_SST)[1:3] <- paste0("SST_", colnames(merged_SST)[1:3])

merged_VIP <- as.data.frame(merged_VIP)
merged_VIP$Gene.stable.ID <- rownames(merged_VIP)
colnames(merged_VIP)[1:3] <- paste0("VIP_", colnames(merged_VIP)[1:3])


merged_hmap <- merge(merged_PYR23,merge(merged_PYR56, merge(merged_PV, merge(merged_SST, merged_VIP, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
merged_hmap <- merged_hmap[,-1]
merged_hmap[is.na(merged_hmap)] <- 0
merged_hmap <- merged_hmap[!(rowSums(merged_hmap)==0),]
merged_hmap <- as.matrix(merged_hmap)

pdf("All_genes heatmap 5percent FDR.pdf")
Heatmap(merged_hmap, column_split = c(rep("PYR L2/3", 3), rep("PYR L5/6", 3), rep("PV", 3), rep("SST", 3), rep("VIP", 3)),
        col=colorRamp2(c(-4, 0, 4),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"), border = TRUE)
dev.off()

