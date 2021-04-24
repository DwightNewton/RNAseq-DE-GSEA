library(DESeq2)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
library(pheatmap)
library(RColorBrewer)
library(variancePartition)
library(vsn)
options(stringsAsFactors = FALSE)
options(scipen = 999)

#Add "biotype" to rowData, restrict to only genes (protein coding) - then variancepartition analysis
#  biotype refers to genomic element type, i.e. protein-coding gene v.s. miRNA, pseudo-gene, etc

biotype_meta <- read.csv("metadata_with_biotype.txt")
load("seFullData with covariate coding.rData")

#Get CT-specific SEs for DE, so off cell-type effects dont influence baseMean and variance estimates in DESeq2
se_PV <- seFullData[,seFullData$Cell.Type == "PVALB"]
se_PV$Subject.Group <- relevel(as.factor(se_PV$Subject.Group), "Control")

se_PYR23 <- seFullData[,seFullData$Cell.Type == "Pyr_L2n3"]
se_PYR23$Subject.Group <- relevel(as.factor(se_PYR23$Subject.Group), "Control")

se_PYR56 <- seFullData[,seFullData$Cell.Type == "Pyr_L5n6"]
se_PYR56$Subject.Group <- relevel(as.factor(se_PYR56$Subject.Group), "Control")

se_SST <- seFullData[,seFullData$Cell.Type == "SST"]
se_SST$Subject.Group <- relevel(as.factor(se_SST$Subject.Group), "Control")

se_VIP <- seFullData[,seFullData$Cell.Type == "VIP"]
se_VIP$Subject.Group <- relevel(as.factor(se_VIP$Subject.Group), "Control")


prot_coding <- unique(biotype_meta$Gene.stable.ID[biotype_meta$Gene.type == "protein_coding"])

#19961 genes
#restrict genomic elements to protein-coding only
se_PV <- se_PV[rownames(se_PV) %in% prot_coding,]
se_PYR23 <- se_PYR23[rownames(se_PYR23) %in% prot_coding,]
se_PYR56 <- se_PYR56[rownames(se_PYR56) %in% prot_coding,]
se_SST <- se_SST[rownames(se_SST) %in% prot_coding,]
se_VIP <- se_VIP[rownames(se_VIP) %in% prot_coding,]

dds0_PV_pc <- DESeqDataSet(se_PV, design = ~ Subject.Group)
dds0_PYR23_pc <- DESeqDataSet(se_PYR23, design = ~ Subject.Group)
dds0_PYR56_pc <- DESeqDataSet(se_PYR56, design = ~ Subject.Group)
dds0_SST_pc <- DESeqDataSet(se_SST, design = ~ Subject.Group)
dds0_VIP_pc <- DESeqDataSet(se_VIP, design = ~ Subject.Group)

PV_vsd_pc <- vst(dds0_PV_pc,blind=FALSE)
PV_vstplot_pc <- meanSdPlot(assay(PV_vsd_pc))
pdf("PV_vst_pc_unfiltered_variance.pdf")
PV_vstplot_pc
dev.off()

PYR23_vsd_pc <- vst(dds0_PYR23_pc,blind=FALSE)
PYR23_vstplot_pc <- meanSdPlot(assay(PYR23_vsd_pc))
pdf("PYR23_vst_pc_unfiltered_variance.pdf")
PYR23_vstplot_pc
dev.off()

PYR56_vsd_pc <- vst(dds0_PYR56_pc,blind=FALSE)
PYR56_vstplot_pc <- meanSdPlot(assay(PYR56_vsd_pc))
pdf("PYR56_vst_pc_unfiltered_variance.pdf")
PYR56_vstplot_pc
dev.off()

SST_vsd_pc <- vst(dds0_SST_pc,blind=FALSE)
SST_vstplot_pc <- meanSdPlot(assay(SST_vsd_pc))
pdf("SST_vst_pc_unfiltered_variance.pdf")
SST_vstplot_pc
dev.off()

VIP_vsd_pc <- vst(dds0_VIP_pc,blind=FALSE)
VIP_vstplot_pc <- meanSdPlot(assay(VIP_vsd_pc))
pdf("VIP_vst_pc_unfiltered_variance.pdf")
VIP_vstplot_pc
dev.off()

#Save original files for any future use
save(se_PV, se_PYR23, se_PYR56, se_SST, se_VIP, file="CT_specific_SEs_proteinCoding.rData")

###Filtered as before (10 counts, at least 20% non 0-counts)
###########PV
#Use design = ~1 for variance parition purposes
dds0_PV_variancePart <- DESeqDataSet(se_PV, design = ~ 1)
dds0_PV_variancePart <- estimateSizeFactors(dds0_PV_variancePart)
PV_variancePart_Counts <- counts(dds0_PV_variancePart)
isexprPV <- rowSums(PV_variancePart_Counts) > 30 & rowSums(PV_variancePart_Counts == 0) <= 60

quantVST_PV <- vst((dds0_PV_variancePart)[isexprPV,])
vstCounts_PV <- assay(quantVST_PV)

#Visualize if transformation controlled variability in SD over
PV_vstplot <- meanSdPlot(assay(quantVST_PV))
pdf("Genes_only_PV (filtered)_vst_noscale_variance.pdf")
PV_vstplot$gg
dev.off()

#######If normalization looks good - do variance partition, but just 1 go
PVinfo <- as.data.frame(colData(dds0_PV_variancePart))
#Variance partition analyses
formOverall <- ~ Age + pH + RIN + PMI + (1|Sex) + (1|Suicide) + (1|Subject.Group)

PV_varPart_Overall <- fitExtractVarPartModel(vstCounts_PV, formOverall, PVinfo)
PV_vp_Overall <- sortCols(PV_varPart_Overall)

#Plot of gene-wise variance explained per variable
pdf("Variance Parition (PV) all covariates_group_Genesonly(filtered).pdf")
plotVarPart(PV_vp_Overall)
dev.off()

#Quick output of mean/median variance explained per variable
Overall_medians_PV <- apply(PV_varPart_Overall, 2, median)
Overall_means_PV <- apply(PV_varPart_Overall, 2, mean)
Overall_output_PV <- rbind(Overall_medians_PV, Overall_means_PV)
write.csv(Overall_output_PV, "PV_Average variance explained_group.csv")

PV_geneMeta <- as.data.frame(rowData(dds0_PV_variancePart))

#Top 200 genes explaiend by each variable
Overall2_pH_top200 <- head(PV_varPart_Overall[order(PV_varPart_Overall$pH, decreasing=TRUE),], n=200)
Overall2_Age_top200 <- head(PV_varPart_Overall[order(PV_varPart_Overall$Age, decreasing=TRUE),], n=200)
Overall2_PMI_top200 <- head(PV_varPart_Overall[order(PV_varPart_Overall$PMI, decreasing=TRUE),], n=200)
Overall2_RIN_top200 <- head(PV_varPart_Overall[order(PV_varPart_Overall$RIN, decreasing=TRUE),], n=200)
Overall2_Sex_top200 <- head(PV_varPart_Overall[order(PV_varPart_Overall$Sex, decreasing=TRUE),], n=200)
Overall2_Suicide_top200 <- head(PV_varPart_Overall[order(PV_varPart_Overall$Suicide, decreasing=TRUE),], n=200)
Overall2_SubjectGroup_top200 <- head(PV_varPart_Overall[order(PV_varPart_Overall$Subject.Group, decreasing=TRUE),], n=200)

Overall2_pH_top200$pHGenes <- rownames(Overall2_pH_top200)
Overall2_Age_top200$AgeGenes <- rownames(Overall2_Age_top200)
Overall2_PMI_top200$PMIGenes <- rownames(Overall2_PMI_top200)
Overall2_RIN_top200$RINGenes <- rownames(Overall2_RIN_top200)
Overall2_Sex_top200$SexGenes <- rownames(Overall2_Sex_top200)
Overall2_Suicide_top200$SuicideGenes <- rownames(Overall2_Suicide_top200)
Overall2_SubjectGroup_top200$Subject.GroupGenes <- rownames(Overall2_SubjectGroup_top200)

Overall2_pH_top200 <- Overall2_pH_top200[,c("pH", "pHGenes")]
Overall2_Age_top200 <- Overall2_Age_top200[,c("Age", "AgeGenes")]
Overall2_PMI_top200 <- Overall2_PMI_top200[,c("PMI", "PMIGenes")]
Overall2_RIN_top200 <- Overall2_RIN_top200[,c("RIN", "RINGenes")]
Overall2_Sex_top200 <- Overall2_Sex_top200[,c("Sex", "SexGenes")]
Overall2_Suicide_top200 <- Overall2_Suicide_top200[,c("Suicide", "SuicideGenes")]
Overall2_SubjectGroup_top200 <- Overall2_SubjectGroup_top200[,c("Subject.Group", "Subject.GroupGenes")]

#Simple cbind of vectors to get single output file
Overall2_top200 <- cbind(Overall2_pH_top200, Overall2_Age_top200, Overall2_PMI_top200, Overall2_RIN_top200, Overall2_Sex_top200, Overall2_Suicide_top200, Overall2_SubjectGroup_top200)
write.csv(Overall2_top200, "PV_top200 variance contributors (subjectGroup).csv")


###########PYR23
dds0_PYR23_variancePart <- DESeqDataSet(se_PYR23, design = ~ 1)
dds0_PYR23_variancePart <- estimateSizeFactors(dds0_PYR23_variancePart)
PYR23_variancePart_Counts <- counts(dds0_PYR23_variancePart)
isexprPYR23 <- rowSums(PYR23_variancePart_Counts) > 30 & rowSums(PYR23_variancePart_Counts == 0) <= 60

quantVST_PYR23 <- vst((dds0_PYR23_variancePart)[isexprPYR23,])
vstCounts_PYR23 <- assay(quantVST_PYR23)

PYR23_vstplot <- meanSdPlot(assay(quantVST_PYR23))
pdf("Genes_only_PYR23 (filtered)_vst_noscale_variance.pdf")
PYR23_vstplot$gg
dev.off()

#######If normalization looks good - do variance partition, but just 1 go
PYR23info <- as.data.frame(colData(dds0_PYR23_variancePart))
#Variance partition analyses
formOverall <- ~ Age + pH + RIN + PMI + (1|Sex) + (1|Suicide) + (1|Subject.Group)

PYR23_varPart_Overall <- fitExtractVarPartModel(vstCounts_PYR23, formOverall, PYR23info)
PYR23_vp_Overall <- sortCols(PYR23_varPart_Overall)

pdf("Variance Parition (PYR23) all covariates_group_Genesonly(filtered).pdf")
plotVarPart(PYR23_vp_Overall)
dev.off()

Overall_medians_PYR23 <- apply(PYR23_varPart_Overall, 2, median)
Overall_means_PYR23 <- apply(PYR23_varPart_Overall, 2, mean)
Overall_output_PYR23 <- rbind(Overall_medians_PYR23, Overall_means_PYR23)
write.csv(Overall_output_PYR23, "PYR23_Average variance explained_group.csv")

PYR23_geneMeta <- as.data.frame(rowData(dds0_PYR23_variancePart))

Overall2_pH_top200 <- head(PYR23_varPart_Overall[order(PYR23_varPart_Overall$pH, decreasing=TRUE),], n=200)
Overall2_Age_top200 <- head(PYR23_varPart_Overall[order(PYR23_varPart_Overall$Age, decreasing=TRUE),], n=200)
Overall2_PMI_top200 <- head(PYR23_varPart_Overall[order(PYR23_varPart_Overall$PMI, decreasing=TRUE),], n=200)
Overall2_RIN_top200 <- head(PYR23_varPart_Overall[order(PYR23_varPart_Overall$RIN, decreasing=TRUE),], n=200)
Overall2_Sex_top200 <- head(PYR23_varPart_Overall[order(PYR23_varPart_Overall$Sex, decreasing=TRUE),], n=200)
Overall2_Suicide_top200 <- head(PYR23_varPart_Overall[order(PYR23_varPart_Overall$Suicide, decreasing=TRUE),], n=200)
Overall2_SubjectGroup_top200 <- head(PYR23_varPart_Overall[order(PYR23_varPart_Overall$Subject.Group, decreasing=TRUE),], n=200)

Overall2_pH_top200$pHGenes <- rownames(Overall2_pH_top200)
Overall2_Age_top200$AgeGenes <- rownames(Overall2_Age_top200)
Overall2_PMI_top200$PMIGenes <- rownames(Overall2_PMI_top200)
Overall2_RIN_top200$RINGenes <- rownames(Overall2_RIN_top200)
Overall2_Sex_top200$SexGenes <- rownames(Overall2_Sex_top200)
Overall2_Suicide_top200$SuicideGenes <- rownames(Overall2_Suicide_top200)
Overall2_SubjectGroup_top200$Subject.GroupGenes <- rownames(Overall2_SubjectGroup_top200)

Overall2_pH_top200 <- Overall2_pH_top200[,c("pH", "pHGenes")]
Overall2_Age_top200 <- Overall2_Age_top200[,c("Age", "AgeGenes")]
Overall2_PMI_top200 <- Overall2_PMI_top200[,c("PMI", "PMIGenes")]
Overall2_RIN_top200 <- Overall2_RIN_top200[,c("RIN", "RINGenes")]
Overall2_Sex_top200 <- Overall2_Sex_top200[,c("Sex", "SexGenes")]
Overall2_Suicide_top200 <- Overall2_Suicide_top200[,c("Suicide", "SuicideGenes")]
Overall2_SubjectGroup_top200 <- Overall2_SubjectGroup_top200[,c("Subject.Group", "Subject.GroupGenes")]

Overall2_top200 <- cbind(Overall2_pH_top200, Overall2_Age_top200, Overall2_PMI_top200, Overall2_RIN_top200, Overall2_Sex_top200, Overall2_Suicide_top200, Overall2_SubjectGroup_top200)
write.csv(Overall2_top200, "PYR23_top200 variance contributors (subjectGroup).csv")



###########PYR56
dds0_PYR56_variancePart <- DESeqDataSet(se_PYR56, design = ~ 1)
dds0_PYR56_variancePart <- estimateSizeFactors(dds0_PYR56_variancePart)
PYR56_variancePart_Counts <- counts(dds0_PYR56_variancePart)
isexprPYR56 <- rowSums(PYR56_variancePart_Counts) > 30 & rowSums(PYR56_variancePart_Counts == 0) <= 60

quantVST_PYR56 <- vst((dds0_PYR56_variancePart)[isexprPYR56,])
vstCounts_PYR56 <- assay(quantVST_PYR56)

PYR56_vstplot <- meanSdPlot(assay(quantVST_PYR56))
pdf("Genes_only_PYR56 (filtered)_vst_noscale_variance.pdf")
PYR56_vstplot$gg
dev.off()

#######If normalization looks good - do variance partition, but just 1 go
PYR56info <- as.data.frame(colData(dds0_PYR56_variancePart))
#Variance partition analyses
formOverall <- ~ Age + pH + RIN + PMI + (1|Sex) + (1|Suicide) + (1|Subject.Group)

PYR56_varPart_Overall <- fitExtractVarPartModel(vstCounts_PYR56, formOverall, PYR56info)
PYR56_vp_Overall <- sortCols(PYR56_varPart_Overall)

pdf("Variance Parition (PYR56) all covariates_group_Genesonly(filtered).pdf")
plotVarPart(PYR56_vp_Overall)
dev.off()

Overall_medians_PYR56 <- apply(PYR56_varPart_Overall, 2, median)
Overall_means_PYR56 <- apply(PYR56_varPart_Overall, 2, mean)
Overall_output_PYR56 <- rbind(Overall_medians_PYR56, Overall_means_PYR56)
write.csv(Overall_output_PYR56, "PYR56_Average variance explained_group.csv")

PYR56_geneMeta <- as.data.frame(rowData(dds0_PYR56_variancePart))

Overall2_pH_top200 <- head(PYR56_varPart_Overall[order(PYR56_varPart_Overall$pH, decreasing=TRUE),], n=200)
Overall2_Age_top200 <- head(PYR56_varPart_Overall[order(PYR56_varPart_Overall$Age, decreasing=TRUE),], n=200)
Overall2_PMI_top200 <- head(PYR56_varPart_Overall[order(PYR56_varPart_Overall$PMI, decreasing=TRUE),], n=200)
Overall2_RIN_top200 <- head(PYR56_varPart_Overall[order(PYR56_varPart_Overall$RIN, decreasing=TRUE),], n=200)
Overall2_Sex_top200 <- head(PYR56_varPart_Overall[order(PYR56_varPart_Overall$Sex, decreasing=TRUE),], n=200)
Overall2_Suicide_top200 <- head(PYR56_varPart_Overall[order(PYR56_varPart_Overall$Suicide, decreasing=TRUE),], n=200)
Overall2_SubjectGroup_top200 <- head(PYR56_varPart_Overall[order(PYR56_varPart_Overall$Subject.Group, decreasing=TRUE),], n=200)

Overall2_pH_top200$pHGenes <- rownames(Overall2_pH_top200)
Overall2_Age_top200$AgeGenes <- rownames(Overall2_Age_top200)
Overall2_PMI_top200$PMIGenes <- rownames(Overall2_PMI_top200)
Overall2_RIN_top200$RINGenes <- rownames(Overall2_RIN_top200)
Overall2_Sex_top200$SexGenes <- rownames(Overall2_Sex_top200)
Overall2_Suicide_top200$SuicideGenes <- rownames(Overall2_Suicide_top200)
Overall2_SubjectGroup_top200$Subject.GroupGenes <- rownames(Overall2_SubjectGroup_top200)

Overall2_pH_top200 <- Overall2_pH_top200[,c("pH", "pHGenes")]
Overall2_Age_top200 <- Overall2_Age_top200[,c("Age", "AgeGenes")]
Overall2_PMI_top200 <- Overall2_PMI_top200[,c("PMI", "PMIGenes")]
Overall2_RIN_top200 <- Overall2_RIN_top200[,c("RIN", "RINGenes")]
Overall2_Sex_top200 <- Overall2_Sex_top200[,c("Sex", "SexGenes")]
Overall2_Suicide_top200 <- Overall2_Suicide_top200[,c("Suicide", "SuicideGenes")]
Overall2_SubjectGroup_top200 <- Overall2_SubjectGroup_top200[,c("Subject.Group", "Subject.GroupGenes")]

Overall2_top200 <- cbind(Overall2_pH_top200, Overall2_Age_top200, Overall2_PMI_top200, Overall2_RIN_top200, Overall2_Sex_top200, Overall2_Suicide_top200, Overall2_SubjectGroup_top200)
write.csv(Overall2_top200, "PYR56_top200 variance contributors (subjectGroup).csv")




###########SST
dds0_SST_variancePart <- DESeqDataSet(se_SST, design = ~ 1)
dds0_SST_variancePart <- estimateSizeFactors(dds0_SST_variancePart)
SST_variancePart_Counts <- counts(dds0_SST_variancePart)
isexprSST <- rowSums(SST_variancePart_Counts) > 30 & rowSums(SST_variancePart_Counts == 0) <= 60

quantVST_SST <- vst((dds0_SST_variancePart)[isexprSST,])
vstCounts_SST <- assay(quantVST_SST)

SST_vstplot <- meanSdPlot(assay(quantVST_SST))
pdf("Genes_only_SST (filtered)_vst_noscale_variance.pdf")
SST_vstplot$gg
dev.off()

#######If normalization looks good - do variance partition, but just 1 go
SSTinfo <- as.data.frame(colData(dds0_SST_variancePart))
#Variance partition analyses
formOverall <- ~ Age + pH + RIN + PMI + (1|Sex) + (1|Suicide) + (1|Subject.Group)

SST_varPart_Overall <- fitExtractVarPartModel(vstCounts_SST, formOverall, SSTinfo)
SST_vp_Overall <- sortCols(SST_varPart_Overall)

pdf("Variance Parition (SST) all covariates_group_Genesonly(filtered).pdf")
plotVarPart(SST_vp_Overall)
dev.off()

Overall_medians_SST <- apply(SST_varPart_Overall, 2, median)
Overall_means_SST <- apply(SST_varPart_Overall, 2, mean)
Overall_output_SST <- rbind(Overall_medians_SST, Overall_means_SST)
write.csv(Overall_output_SST, "SST_Average variance explained_group.csv")

SST_geneMeta <- as.data.frame(rowData(dds0_SST_variancePart))

Overall2_pH_top200 <- head(SST_varPart_Overall[order(SST_varPart_Overall$pH, decreasing=TRUE),], n=200)
Overall2_Age_top200 <- head(SST_varPart_Overall[order(SST_varPart_Overall$Age, decreasing=TRUE),], n=200)
Overall2_PMI_top200 <- head(SST_varPart_Overall[order(SST_varPart_Overall$PMI, decreasing=TRUE),], n=200)
Overall2_RIN_top200 <- head(SST_varPart_Overall[order(SST_varPart_Overall$RIN, decreasing=TRUE),], n=200)
Overall2_Sex_top200 <- head(SST_varPart_Overall[order(SST_varPart_Overall$Sex, decreasing=TRUE),], n=200)
Overall2_Suicide_top200 <- head(SST_varPart_Overall[order(SST_varPart_Overall$Suicide, decreasing=TRUE),], n=200)
Overall2_SubjectGroup_top200 <- head(SST_varPart_Overall[order(SST_varPart_Overall$Subject.Group, decreasing=TRUE),], n=200)

Overall2_pH_top200$pHGenes <- rownames(Overall2_pH_top200)
Overall2_Age_top200$AgeGenes <- rownames(Overall2_Age_top200)
Overall2_PMI_top200$PMIGenes <- rownames(Overall2_PMI_top200)
Overall2_RIN_top200$RINGenes <- rownames(Overall2_RIN_top200)
Overall2_Sex_top200$SexGenes <- rownames(Overall2_Sex_top200)
Overall2_Suicide_top200$SuicideGenes <- rownames(Overall2_Suicide_top200)
Overall2_SubjectGroup_top200$Subject.GroupGenes <- rownames(Overall2_SubjectGroup_top200)

Overall2_pH_top200 <- Overall2_pH_top200[,c("pH", "pHGenes")]
Overall2_Age_top200 <- Overall2_Age_top200[,c("Age", "AgeGenes")]
Overall2_PMI_top200 <- Overall2_PMI_top200[,c("PMI", "PMIGenes")]
Overall2_RIN_top200 <- Overall2_RIN_top200[,c("RIN", "RINGenes")]
Overall2_Sex_top200 <- Overall2_Sex_top200[,c("Sex", "SexGenes")]
Overall2_Suicide_top200 <- Overall2_Suicide_top200[,c("Suicide", "SuicideGenes")]
Overall2_SubjectGroup_top200 <- Overall2_SubjectGroup_top200[,c("Subject.Group", "Subject.GroupGenes")]

Overall2_top200 <- cbind(Overall2_pH_top200, Overall2_Age_top200, Overall2_PMI_top200, Overall2_RIN_top200, Overall2_Sex_top200, Overall2_Suicide_top200, Overall2_SubjectGroup_top200)
write.csv(Overall2_top200, "SST_top200 variance contributors (subjectGroup).csv")



###########VIP
dds0_VIP_variancePart <- DESeqDataSet(se_VIP, design = ~ 1)
dds0_VIP_variancePart <- estimateSizeFactors(dds0_VIP_variancePart)
VIP_variancePart_Counts <- counts(dds0_VIP_variancePart)
isexprVIP <- rowSums(VIP_variancePart_Counts) > 30 & rowSums(VIP_variancePart_Counts == 0) <= 60

quantVST_VIP <- vst((dds0_VIP_variancePart)[isexprVIP,])
vstCounts_VIP <- assay(quantVST_VIP)

VIP_vstplot <- meanSdPlot(assay(quantVST_VIP))
pdf("Genes_only_VIP (filtered)_vst_noscale_variance.pdf")
VIP_vstplot$gg
dev.off()

#######If normalization looks good - do variance partition, but just 1 go
VIPinfo <- as.data.frame(colData(dds0_VIP_variancePart))
#Variance partition analyses
formOverall <- ~ Age + pH + RIN + PMI + (1|Sex) + (1|Suicide) + (1|Subject.Group)

VIP_varPart_Overall <- fitExtractVarPartModel(vstCounts_VIP, formOverall, VIPinfo)
VIP_vp_Overall <- sortCols(VIP_varPart_Overall)

pdf("Variance Parition (VIP) all covariates_group_Genesonly(filtered).pdf")
plotVarPart(VIP_vp_Overall)
dev.off()

Overall_medians_VIP <- apply(VIP_varPart_Overall, 2, median)
Overall_means_VIP <- apply(VIP_varPart_Overall, 2, mean)
Overall_output_VIP <- rbind(Overall_medians_VIP, Overall_means_VIP)
write.csv(Overall_output_VIP, "VIP_Average variance explained_group.csv")

VIP_geneMeta <- as.data.frame(rowData(dds0_VIP_variancePart))

Overall2_pH_top200 <- head(VIP_varPart_Overall[order(VIP_varPart_Overall$pH, decreasing=TRUE),], n=200)
Overall2_Age_top200 <- head(VIP_varPart_Overall[order(VIP_varPart_Overall$Age, decreasing=TRUE),], n=200)
Overall2_PMI_top200 <- head(VIP_varPart_Overall[order(VIP_varPart_Overall$PMI, decreasing=TRUE),], n=200)
Overall2_RIN_top200 <- head(VIP_varPart_Overall[order(VIP_varPart_Overall$RIN, decreasing=TRUE),], n=200)
Overall2_Sex_top200 <- head(VIP_varPart_Overall[order(VIP_varPart_Overall$Sex, decreasing=TRUE),], n=200)
Overall2_Suicide_top200 <- head(VIP_varPart_Overall[order(VIP_varPart_Overall$Suicide, decreasing=TRUE),], n=200)
Overall2_SubjectGroup_top200 <- head(VIP_varPart_Overall[order(VIP_varPart_Overall$Subject.Group, decreasing=TRUE),], n=200)

Overall2_pH_top200$pHGenes <- rownames(Overall2_pH_top200)
Overall2_Age_top200$AgeGenes <- rownames(Overall2_Age_top200)
Overall2_PMI_top200$PMIGenes <- rownames(Overall2_PMI_top200)
Overall2_RIN_top200$RINGenes <- rownames(Overall2_RIN_top200)
Overall2_Sex_top200$SexGenes <- rownames(Overall2_Sex_top200)
Overall2_Suicide_top200$SuicideGenes <- rownames(Overall2_Suicide_top200)
Overall2_SubjectGroup_top200$Subject.GroupGenes <- rownames(Overall2_SubjectGroup_top200)

Overall2_pH_top200 <- Overall2_pH_top200[,c("pH", "pHGenes")]
Overall2_Age_top200 <- Overall2_Age_top200[,c("Age", "AgeGenes")]
Overall2_PMI_top200 <- Overall2_PMI_top200[,c("PMI", "PMIGenes")]
Overall2_RIN_top200 <- Overall2_RIN_top200[,c("RIN", "RINGenes")]
Overall2_Sex_top200 <- Overall2_Sex_top200[,c("Sex", "SexGenes")]
Overall2_Suicide_top200 <- Overall2_Suicide_top200[,c("Suicide", "SuicideGenes")]
Overall2_SubjectGroup_top200 <- Overall2_SubjectGroup_top200[,c("Subject.Group", "Subject.GroupGenes")]

Overall2_top200 <- cbind(Overall2_pH_top200, Overall2_Age_top200, Overall2_PMI_top200, Overall2_RIN_top200, Overall2_Sex_top200, Overall2_Suicide_top200, Overall2_SubjectGroup_top200)
write.csv(Overall2_top200, "VIP_top200 variance contributors (subjectGroup).csv")


