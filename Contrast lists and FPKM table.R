library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(fdrtool)
options(stringsAsFactors = FALSE)
options(scipen = 999)

load("res_PV_exons_groupCT_age_sex_group.rData")
load("res_PYR23_exons_groupCT_age_sex_group.rData")
load("res_PYR56_exons_groupCT_age_sex_group.rData")
load("res_SST_exons_groupCT_age_sex_group.rData")
load("res_VIP_exons_groupCT_age_sex_group.rData")

#FDRtool-ing - run once and save
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

save(fdrT_CT_PV_MDD, fdrT_CT_PYR23_MDD, fdrT_CT_PYR56_MDD, fdrT_CT_SST_MDD, fdrT_CT_VIP_MDD, fdrT_CT_PV_Bipolar, fdrT_CT_PYR23_Bipolar, fdrT_CT_PYR56_Bipolar, fdrT_CT_SST_Bipolar, fdrT_CT_VIP_Bipolar, fdrT_CT_PV_SCHIZ, fdrT_CT_PYR23_SCHIZ, fdrT_CT_PYR56_SCHIZ, fdrT_CT_SST_SCHIZ, fdrT_CT_VIP_SCHIZ,file="FDRtooled_results.rData")


##Merge outputs for each cell-type by contrast, then merge them all
#Need se for rowData
load("seFullData_proteincoding.rData")
load("FDRtooled_results.rData")
meta <- as.data.frame(rowData(seFullData_pc))


PV_MDD <- as.data.frame(fdrT_CT_PV_MDD)
PV_MDD$Gene.stable.ID <- rownames(PV_MDD)
PV_MDD$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% PV_MDD$Gene.stable.ID]
PV_MDD <- PV_MDD[!duplicated(PV_MDD$Gene.name), c(7,8,2,5,6)]
colnames(PV_MDD)[3:5] <- paste0("PV_MDD_", colnames(PV_MDD)[3:5])

PYR23_MDD <- as.data.frame(fdrT_CT_PYR23_MDD)
PYR23_MDD$Gene.stable.ID <- rownames(PYR23_MDD)
PYR23_MDD$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% PYR23_MDD$Gene.stable.ID]
PYR23_MDD <- PYR23_MDD[!duplicated(PYR23_MDD$Gene.name), c(7,8,2,5,6)]
colnames(PYR23_MDD)[3:5] <- paste0("PYR23_MDD_", colnames(PYR23_MDD)[3:5])

PYR56_MDD <- as.data.frame(fdrT_CT_PYR56_MDD)
PYR56_MDD$Gene.stable.ID <- rownames(PYR56_MDD)
PYR56_MDD$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% PYR56_MDD$Gene.stable.ID]
PYR56_MDD <- PYR56_MDD[!duplicated(PYR56_MDD$Gene.name), c(7,8,2,5,6)]
colnames(PYR56_MDD)[3:5] <- paste0("PYR56_MDD_", colnames(PYR56_MDD)[3:5])

SST_MDD <- as.data.frame(fdrT_CT_SST_MDD)
SST_MDD$Gene.stable.ID <- rownames(SST_MDD)
SST_MDD$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% SST_MDD$Gene.stable.ID]
SST_MDD <- SST_MDD[!duplicated(SST_MDD$Gene.name), c(7,8,2,5,6)]
colnames(SST_MDD)[3:5] <- paste0("SST_MDD_", colnames(SST_MDD)[3:5])

VIP_MDD <- as.data.frame(fdrT_CT_VIP_MDD)
VIP_MDD$Gene.stable.ID <- rownames(VIP_MDD)
VIP_MDD$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% VIP_MDD$Gene.stable.ID]
VIP_MDD <- VIP_MDD[!duplicated(VIP_MDD$Gene.name), c(7,8,2,5,6)]
colnames(VIP_MDD)[3:5] <- paste0("VIP_MDD_", colnames(VIP_MDD)[3:5])

PV_Bipolar <- as.data.frame(fdrT_CT_PV_Bipolar)
PV_Bipolar$Gene.stable.ID <- rownames(PV_Bipolar)
PV_Bipolar$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% PV_Bipolar$Gene.stable.ID]
PV_Bipolar <- PV_Bipolar[!duplicated(PV_Bipolar$Gene.name), c(7,8,2,5,6)]
colnames(PV_Bipolar)[3:5] <- paste0("PV_Bipolar_", colnames(PV_Bipolar)[3:5])

PYR23_Bipolar <- as.data.frame(fdrT_CT_PYR23_Bipolar)
PYR23_Bipolar$Gene.stable.ID <- rownames(PYR23_Bipolar)
PYR23_Bipolar$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% PYR23_Bipolar$Gene.stable.ID]
PYR23_Bipolar <- PYR23_Bipolar[!duplicated(PYR23_Bipolar$Gene.name), c(7,8,2,5,6)]
colnames(PYR23_Bipolar)[3:5] <- paste0("PYR23_Bipolar_", colnames(PYR23_Bipolar)[3:5])

PYR56_Bipolar <- as.data.frame(fdrT_CT_PYR56_Bipolar)
PYR56_Bipolar$Gene.stable.ID <- rownames(PYR56_Bipolar)
PYR56_Bipolar$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% PYR56_Bipolar$Gene.stable.ID]
PYR56_Bipolar <- PYR56_Bipolar[!duplicated(PYR56_Bipolar$Gene.name), c(7,8,2,5,6)]
colnames(PYR56_Bipolar)[3:5] <- paste0("PYR56_Bipolar_", colnames(PYR56_Bipolar)[3:5])

SST_Bipolar <- as.data.frame(fdrT_CT_SST_Bipolar)
SST_Bipolar$Gene.stable.ID <- rownames(SST_Bipolar)
SST_Bipolar$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% SST_Bipolar$Gene.stable.ID]
SST_Bipolar <- SST_Bipolar[!duplicated(SST_Bipolar$Gene.name), c(7,8,2,5,6)]
colnames(SST_Bipolar)[3:5] <- paste0("SST_Bipolar_", colnames(SST_Bipolar)[3:5])

VIP_Bipolar <- as.data.frame(fdrT_CT_VIP_Bipolar)
VIP_Bipolar$Gene.stable.ID <- rownames(VIP_Bipolar)
VIP_Bipolar$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% VIP_Bipolar$Gene.stable.ID]
VIP_Bipolar <- VIP_Bipolar[!duplicated(VIP_Bipolar$Gene.name), c(7,8,2,5,6)]
colnames(VIP_Bipolar)[3:5] <- paste0("VIP_Bipolar_", colnames(VIP_Bipolar)[3:5])

PV_SCHIZ <- as.data.frame(fdrT_CT_PV_SCHIZ)
PV_SCHIZ$Gene.stable.ID <- rownames(PV_SCHIZ)
PV_SCHIZ$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% PV_SCHIZ$Gene.stable.ID]
PV_SCHIZ <- PV_SCHIZ[!duplicated(PV_SCHIZ$Gene.name), c(7,8,2,5,6)]
colnames(PV_SCHIZ)[3:5] <- paste0("PV_SCHIZ_", colnames(PV_SCHIZ)[3:5])

PYR23_SCHIZ <- as.data.frame(fdrT_CT_PYR23_SCHIZ)
PYR23_SCHIZ$Gene.stable.ID <- rownames(PYR23_SCHIZ)
PYR23_SCHIZ$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% PYR23_SCHIZ$Gene.stable.ID]
PYR23_SCHIZ <- PYR23_SCHIZ[!duplicated(PYR23_SCHIZ$Gene.name), c(7,8,2,5,6)]
colnames(PYR23_SCHIZ)[3:5] <- paste0("PYR23_SCHIZ_", colnames(PYR23_SCHIZ)[3:5])

PYR56_SCHIZ <- as.data.frame(fdrT_CT_PYR56_SCHIZ)
PYR56_SCHIZ$Gene.stable.ID <- rownames(PYR56_SCHIZ)
PYR56_SCHIZ$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% PYR56_SCHIZ$Gene.stable.ID]
PYR56_SCHIZ <- PYR56_SCHIZ[!duplicated(PYR56_SCHIZ$Gene.name), c(7,8,2,5,6)]
colnames(PYR56_SCHIZ)[3:5] <- paste0("PYR56_SCHIZ_", colnames(PYR56_SCHIZ)[3:5])

SST_SCHIZ <- as.data.frame(fdrT_CT_SST_SCHIZ)
SST_SCHIZ$Gene.stable.ID <- rownames(SST_SCHIZ)
SST_SCHIZ$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% SST_SCHIZ$Gene.stable.ID]
SST_SCHIZ <- SST_SCHIZ[!duplicated(SST_SCHIZ$Gene.name), c(7,8,2,5,6)]
colnames(SST_SCHIZ)[3:5] <- paste0("SST_SCHIZ_", colnames(SST_SCHIZ)[3:5])

VIP_SCHIZ <- as.data.frame(fdrT_CT_VIP_SCHIZ)
VIP_SCHIZ$Gene.stable.ID <- rownames(VIP_SCHIZ)
VIP_SCHIZ$Gene.name <- meta$Gene.name[meta$Gene.stable.ID %in% VIP_SCHIZ$Gene.stable.ID]
VIP_SCHIZ <- VIP_SCHIZ[!duplicated(VIP_SCHIZ$Gene.name), c(7,8,2,5,6)]
colnames(VIP_SCHIZ)[3:5] <- paste0("VIP_SCHIZ_", colnames(VIP_SCHIZ)[3:5])



outputTable_MDD <- merge(PYR23_MDD, merge(PYR56_MDD, merge(PV_MDD, merge(SST_MDD, VIP_MDD, by=c("Gene.stable.ID", "Gene.name"), all=TRUE), by=c("Gene.stable.ID", "Gene.name"), all=TRUE), by=c("Gene.stable.ID", "Gene.name"), all=TRUE), by=c("Gene.stable.ID", "Gene.name"), all=TRUE)

outputTable_BD <- merge(PYR23_Bipolar, merge(PYR56_Bipolar, merge(PV_Bipolar, merge(SST_Bipolar, VIP_Bipolar, by=c("Gene.stable.ID", "Gene.name"), all=TRUE), by=c("Gene.stable.ID", "Gene.name"), all=TRUE), by=c("Gene.stable.ID", "Gene.name"), all=TRUE), by=c("Gene.stable.ID", "Gene.name"), all=TRUE)

outputTable_SCZ <- merge(PYR23_SCHIZ, merge(PYR56_SCHIZ, merge(PV_SCHIZ, merge(SST_SCHIZ, VIP_SCHIZ, by=c("Gene.stable.ID", "Gene.name"), all=TRUE), by=c("Gene.stable.ID", "Gene.name"), all=TRUE), by=c("Gene.stable.ID", "Gene.name"), all=TRUE), by=c("Gene.stable.ID", "Gene.name"), all=TRUE)


outputTable_LFC <- merge(outputTable_MDD, merge(outputTable_BD, outputTable_SCZ, by=c("Gene.stable.ID", "Gene.name"), all=TRUE),by=c("Gene.stable.ID", "Gene.name"), all=TRUE)
  
#Get FPKM for all the genes
ddslist <- list.files(pattern = "dds")
for(i in ddslist){load(i)}


rowData(dds_PV_exons_groupCT_age_sex_group)
RD <- rowData(dds_PV_exons_groupCT_age_sex_group)
RD$basepairs <- RD$RD$Transcript.length..including.UTRs.and.CDS.
rowData(dds_PV_exons_groupCT_age_sex_group) <- RD
PV_FPKM <- as.data.frame(fpkm(dds_PV_exons_groupCT_age_sex_group))
PV_FPKM$Gene.stable.ID <- rownames(PV_FPKM)

rowData(dds_PYR23_exons_groupCT_age_sex_group)
RD <- rowData(dds_PYR23_exons_groupCT_age_sex_group)
RD$basepairs <- RD$RD$Transcript.length..including.UTRs.and.CDS.
rowData(dds_PYR23_exons_groupCT_age_sex_group) <- RD
PYR23_FPKM <- as.data.frame(fpkm(dds_PYR23_exons_groupCT_age_sex_group))
PYR23_FPKM$Gene.stable.ID <- rownames(PYR23_FPKM)

rowData(dds_PYR56_exons_groupCT_age_sex_group)
RD <- rowData(dds_PYR56_exons_groupCT_age_sex_group)
RD$basepairs <- RD$RD$Transcript.length..including.UTRs.and.CDS.
rowData(dds_PYR56_exons_groupCT_age_sex_group) <- RD
PYR56_FPKM <- as.data.frame(fpkm(dds_PYR56_exons_groupCT_age_sex_group))
PYR56_FPKM$Gene.stable.ID <- rownames(PYR56_FPKM)

rowData(dds_SST_exons_groupCT_age_sex_group)
RD <- rowData(dds_SST_exons_groupCT_age_sex_group)
RD$basepairs <- RD$RD$Transcript.length..including.UTRs.and.CDS.
rowData(dds_SST_exons_groupCT_age_sex_group) <- RD
SST_FPKM <- as.data.frame(fpkm(dds_SST_exons_groupCT_age_sex_group))
SST_FPKM$Gene.stable.ID <- rownames(SST_FPKM)

rowData(dds_VIP_exons_groupCT_age_sex_group)
RD <- rowData(dds_VIP_exons_groupCT_age_sex_group)
RD$basepairs <- RD$RD$Transcript.length..including.UTRs.and.CDS.
rowData(dds_VIP_exons_groupCT_age_sex_group) <- RD
VIP_FPKM <- as.data.frame(fpkm(dds_VIP_exons_groupCT_age_sex_group))
VIP_FPKM$Gene.stable.ID <- rownames(VIP_FPKM)


outputTable_counts <- merge(PYR23_FPKM, merge(PYR56_FPKM, merge(PV_FPKM, merge(SST_FPKM, VIP_FPKM, by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE), by="Gene.stable.ID", all=TRUE)
outputTable_counts <- outputTable_counts[outputTable_counts$Gene.stable.ID %in% outputTable_LFC$Gene.stable.ID,]
outputTable_counts$Gene.name <- outputTable_LFC$Gene.name[outputTable_LFC$Gene.stable.ID %in% outputTable_counts$Gene.stable.ID]

totalColData <- rbind(as.data.frame(colData(dds_PV_exons_groupCT_age_sex_group)), as.data.frame(colData(dds_PYR23_exons_groupCT_age_sex_group)), as.data.frame(colData(dds_PYR56_exons_groupCT_age_sex_group)), as.data.frame(colData(dds_SST_exons_groupCT_age_sex_group)), as.data.frame(colData(dds_VIP_exons_groupCT_age_sex_group)))

totalColData$matchvalue <- rownames(totalColData)
totalColData <- totalColData[match(colnames(outputTable_counts)[2:380], totalColData$matchvalue),]


write.csv(outputTable_LFC, "PITT_tetrad_cohort_SCT_contrasts.csv")
write.csv(outputTable_counts, "PITT_tetrad_cohort_SCT_FPKM.csv")
write.csv(totalColData, "FPKMmetadata.csv")
