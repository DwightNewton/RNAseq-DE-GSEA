library(fgsea)
library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(fdrtool)
library(ggplot2)
library(reshape2)
library(dplyr)
library(GGally)
library(Vennerable)
library(gplots)
library(venn)
options(stringsAsFactors = FALSE)
options(scipen = 999)


#First, just get numbers of up/down pathways in each cell-type/contrast and do some venn diagrams for each direction (and overall)
setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/")

filelist <- list.files(pattern=".rData")
for(i in filelist){load(i)}

#MDD
MDD_pathwaycounts <- matrix(nrow=8, ncol=30)
colnames(MDD_pathwaycounts) <- as.vector(sapply(c("PV_", "PYR23_", "PYR56_","SST_", "VIP_"), paste0, c(0, 0.5, 0.9, 1, 1.5, 2), simplify = TRUE))
rownames(MDD_pathwaycounts) <- c("Up_pval(0.05)", "Down_pval(0.05)", "Up_pval(0.01)", "Down_pval(0.01)", "Up_FDR_25", "Down_FDR_25", "Up_FDR_5", "Down_FDR_5")

GSEA_datalist <- c("PV_MDD_gseaM", "PV_MDD_gseaM0p5", "PV_MDD_gseaM0p9", "PV_MDD_gseaM1", "PV_MDD_gseaM1p5", "PV_MDD_gseaM2","PYR23_MDD_gseaM", "PYR23_MDD_gseaM0p5", "PYR23_MDD_gseaM0p9", "PYR23_MDD_gseaM1", "PYR23_MDD_gseaM1p5", "PYR23_MDD_gseaM2","PYR56_MDD_gseaM", "PYR56_MDD_gseaM0p5", "PYR56_MDD_gseaM0p9", "PYR56_MDD_gseaM1", "PYR56_MDD_gseaM1p5", "PYR56_MDD_gseaM2","SST_MDD_gseaM", "SST_MDD_gseaM0p5", "SST_MDD_gseaM0p9", "SST_MDD_gseaM1", "SST_MDD_gseaM1p5", "SST_MDD_gseaM2", "VIP_MDD_gseaM", "VIP_MDD_gseaM0p5", "VIP_MDD_gseaM0p9", "VIP_MDD_gseaM1", "VIP_MDD_gseaM1p5", "VIP_MDD_gseaM2")

for(i in 1:length(colnames((MDD_pathwaycounts)))){
  workingData <- get(GSEA_datalist[i])
  MDD_pathwaycounts[1, i] <- sum(workingData$pval<0.05 & workingData$ES > 0)
  MDD_pathwaycounts[2, i] <- sum(workingData$pval<0.05 & workingData$ES < 0)
  MDD_pathwaycounts[3, i] <- sum(workingData$pval<0.01 & workingData$ES > 0)
  MDD_pathwaycounts[4, i] <- sum(workingData$pval<0.01 & workingData$ES < 0)
  MDD_pathwaycounts[5, i] <- sum(workingData$padj<0.25 & workingData$ES > 0)
  MDD_pathwaycounts[6, i] <- sum(workingData$padj<0.25 & workingData$ES < 0)
  MDD_pathwaycounts[7, i] <- sum(workingData$padj<0.05 & workingData$ES > 0)
  MDD_pathwaycounts[8, i] <- sum(workingData$padj<0.05 & workingData$ES < 0)
}
write.csv(MDD_pathwaycounts, "MDD_pathway_counts.csv")

#Bipolar
Bipolar_pathwaycounts <- matrix(nrow=8, ncol=30)
colnames(Bipolar_pathwaycounts) <- as.vector(sapply(c("PV_", "PYR23_", "PYR56_","SST_", "VIP_"), paste0, c(0, 0.5, 0.9, 1, 1.5, 2), simplify = TRUE))
rownames(Bipolar_pathwaycounts) <- c("Up_pval(0.05)", "Down_pval(0.05)", "Up_pval(0.01)", "Down_pval(0.01)", "Up_FDR_25", "Down_FDR_25", "Up_FDR_5", "Down_FDR_5")

GSEA_datalist <- c("PV_Bipolar_gseaM", "PV_Bipolar_gseaM0p5", "PV_Bipolar_gseaM0p9", "PV_Bipolar_gseaM1", "PV_Bipolar_gseaM1p5", "PV_Bipolar_gseaM2","PYR23_Bipolar_gseaM", "PYR23_Bipolar_gseaM0p5", "PYR23_Bipolar_gseaM0p9", "PYR23_Bipolar_gseaM1", "PYR23_Bipolar_gseaM1p5", "PYR23_Bipolar_gseaM2","PYR56_Bipolar_gseaM", "PYR56_Bipolar_gseaM0p5", "PYR56_Bipolar_gseaM0p9", "PYR56_Bipolar_gseaM1", "PYR56_Bipolar_gseaM1p5", "PYR56_Bipolar_gseaM2","SST_Bipolar_gseaM", "SST_Bipolar_gseaM0p5", "SST_Bipolar_gseaM0p9", "SST_Bipolar_gseaM1", "SST_Bipolar_gseaM1p5", "SST_Bipolar_gseaM2", "VIP_Bipolar_gseaM", "VIP_Bipolar_gseaM0p5", "VIP_Bipolar_gseaM0p9", "VIP_Bipolar_gseaM1", "VIP_Bipolar_gseaM1p5", "VIP_Bipolar_gseaM2")

for(i in 1:length(colnames((Bipolar_pathwaycounts)))){
  workingData <- get(GSEA_datalist[i])
  Bipolar_pathwaycounts[1, i] <- sum(workingData$pval<0.05 & workingData$ES > 0)
  Bipolar_pathwaycounts[2, i] <- sum(workingData$pval<0.05 & workingData$ES < 0)
  Bipolar_pathwaycounts[3, i] <- sum(workingData$pval<0.01 & workingData$ES > 0)
  Bipolar_pathwaycounts[4, i] <- sum(workingData$pval<0.01 & workingData$ES < 0)
  Bipolar_pathwaycounts[5, i] <- sum(workingData$padj<0.25 & workingData$ES > 0)
  Bipolar_pathwaycounts[6, i] <- sum(workingData$padj<0.25 & workingData$ES < 0)
  Bipolar_pathwaycounts[7, i] <- sum(workingData$padj<0.05 & workingData$ES > 0)
  Bipolar_pathwaycounts[8, i] <- sum(workingData$padj<0.05 & workingData$ES < 0)
}
write.csv(Bipolar_pathwaycounts, "Bipolar_pathway_counts.csv")

#SCHIZ
SCHIZ_pathwaycounts <- matrix(nrow=8, ncol=30)
colnames(SCHIZ_pathwaycounts) <- as.vector(sapply(c("PV_", "PYR23_", "PYR56_","SST_", "VIP_"), paste0, c(0, 0.5, 0.9, 1, 1.5, 2), simplify = TRUE))
rownames(SCHIZ_pathwaycounts) <- c("Up_pval(0.05)", "Down_pval(0.05)", "Up_pval(0.01)", "Down_pval(0.01)", "Up_FDR_25", "Down_FDR_25", "Up_FDR_5", "Down_FDR_5")

GSEA_datalist <- c("PV_SCHIZ_gseaM", "PV_SCHIZ_gseaM0p5", "PV_SCHIZ_gseaM0p9", "PV_SCHIZ_gseaM1", "PV_SCHIZ_gseaM1p5", "PV_SCHIZ_gseaM2","PYR23_SCHIZ_gseaM", "PYR23_SCHIZ_gseaM0p5", "PYR23_SCHIZ_gseaM0p9", "PYR23_SCHIZ_gseaM1", "PYR23_SCHIZ_gseaM1p5", "PYR23_SCHIZ_gseaM2","PYR56_SCHIZ_gseaM", "PYR56_SCHIZ_gseaM0p5", "PYR56_SCHIZ_gseaM0p9", "PYR56_SCHIZ_gseaM1", "PYR56_SCHIZ_gseaM1p5", "PYR56_SCHIZ_gseaM2","SST_SCHIZ_gseaM", "SST_SCHIZ_gseaM0p5", "SST_SCHIZ_gseaM0p9", "SST_SCHIZ_gseaM1", "SST_SCHIZ_gseaM1p5", "SST_SCHIZ_gseaM2", "VIP_SCHIZ_gseaM", "VIP_SCHIZ_gseaM0p5", "VIP_SCHIZ_gseaM0p9", "VIP_SCHIZ_gseaM1", "VIP_SCHIZ_gseaM1p5", "VIP_SCHIZ_gseaM2")

for(i in 1:length(colnames((SCHIZ_pathwaycounts)))){
  workingData <- get(GSEA_datalist[i])
  SCHIZ_pathwaycounts[1, i] <- sum(workingData$pval<0.05 & workingData$ES > 0)
  SCHIZ_pathwaycounts[2, i] <- sum(workingData$pval<0.05 & workingData$ES < 0)
  SCHIZ_pathwaycounts[3, i] <- sum(workingData$pval<0.01 & workingData$ES > 0)
  SCHIZ_pathwaycounts[4, i] <- sum(workingData$pval<0.01 & workingData$ES < 0)
  SCHIZ_pathwaycounts[5, i] <- sum(workingData$padj<0.25 & workingData$ES > 0)
  SCHIZ_pathwaycounts[6, i] <- sum(workingData$padj<0.25 & workingData$ES < 0)
  SCHIZ_pathwaycounts[7, i] <- sum(workingData$padj<0.05 & workingData$ES > 0)
  SCHIZ_pathwaycounts[8, i] <- sum(workingData$padj<0.05 & workingData$ES < 0)
}
write.csv(SCHIZ_pathwaycounts, "SCHIZ_pathway_counts.csv")

#Heatmaps/Venns of all sig hits (seperately for 25%FDR and pval<0.05)
#P<0.05 results
#Heatmap/venns across cell-types, within disorders
#MDD
PV_MDD_SigUp <- PV_MDD_gseaM1[PV_MDD_gseaM1$pval < 0.05 & PV_MDD_gseaM1$ES > 0, ]
PV_MDD_SigDown <- PV_MDD_gseaM1[PV_MDD_gseaM1$pval < 0.05 & PV_MDD_gseaM1$ES < 0, ]
PV_MDD_vennUp <- PV_MDD_SigUp$pathway
PV_MDD_vennDown <- PV_MDD_SigDown$pathway
PV_MDD_hmap <- PV_MDD_gseaM1[PV_MDD_gseaM1$pval < 0.05, c("pathway", "ES")]
PV_MDD_hmap$ES[PV_MDD_hmap$ES > 0] <- 1
PV_MDD_hmap$ES[PV_MDD_hmap$ES < 0] <- -1

PYR23_MDD_SigUp <- PYR23_MDD_gseaM1[PYR23_MDD_gseaM1$pval < 0.05 & PYR23_MDD_gseaM1$ES > 0, ]
PYR23_MDD_SigDown <- PYR23_MDD_gseaM1[PYR23_MDD_gseaM1$pval < 0.05 & PYR23_MDD_gseaM1$ES < 0, ]
PYR23_MDD_vennUp <- PYR23_MDD_SigUp$pathway
PYR23_MDD_vennDown <- PYR23_MDD_SigDown$pathway
PYR23_MDD_hmap <- PYR23_MDD_gseaM1[PYR23_MDD_gseaM1$pval < 0.05, c("pathway", "ES")]
PYR23_MDD_hmap$ES[PYR23_MDD_hmap$ES > 0] <- 1
PYR23_MDD_hmap$ES[PYR23_MDD_hmap$ES < 0] <- -1

PYR56_MDD_SigUp <- PYR56_MDD_gseaM1[PYR56_MDD_gseaM1$pval < 0.05 & PYR56_MDD_gseaM1$ES > 0, ]
PYR56_MDD_SigDown <- PYR56_MDD_gseaM1[PYR56_MDD_gseaM1$pval < 0.05 & PYR56_MDD_gseaM1$ES < 0, ]
PYR56_MDD_vennUp <- PYR56_MDD_SigUp$pathway
PYR56_MDD_vennDown <- PYR56_MDD_SigDown$pathway
PYR56_MDD_hmap <- PYR56_MDD_gseaM1[PYR56_MDD_gseaM1$pval < 0.05, c("pathway", "ES")]
PYR56_MDD_hmap$ES[PYR56_MDD_hmap$ES > 0] <- 1
PYR56_MDD_hmap$ES[PYR56_MDD_hmap$ES < 0] <- -1

SST_MDD_SigUp <- SST_MDD_gseaM1[SST_MDD_gseaM1$pval < 0.05 & SST_MDD_gseaM1$ES > 0, ]
SST_MDD_SigDown <- SST_MDD_gseaM1[SST_MDD_gseaM1$pval < 0.05 & SST_MDD_gseaM1$ES < 0, ]
SST_MDD_vennUp <- SST_MDD_SigUp$pathway
SST_MDD_vennDown <- SST_MDD_SigDown$pathway
SST_MDD_hmap <- SST_MDD_gseaM1[SST_MDD_gseaM1$pval < 0.05, c("pathway", "ES")]
SST_MDD_hmap$ES[SST_MDD_hmap$ES > 0] <- 1
SST_MDD_hmap$ES[SST_MDD_hmap$ES < 0] <- -1

VIP_MDD_SigUp <- VIP_MDD_gseaM1[VIP_MDD_gseaM1$pval < 0.05 & VIP_MDD_gseaM1$ES > 0, ]
VIP_MDD_SigDown <- VIP_MDD_gseaM1[VIP_MDD_gseaM1$pval < 0.05 & VIP_MDD_gseaM1$ES < 0, ]
VIP_MDD_vennUp <- VIP_MDD_SigUp$pathway
VIP_MDD_vennDown <- VIP_MDD_SigDown$pathway
VIP_MDD_hmap <- VIP_MDD_gseaM1[VIP_MDD_gseaM1$pval < 0.05, c("pathway", "ES")]
VIP_MDD_hmap$ES[VIP_MDD_hmap$ES > 0] <- 1
VIP_MDD_hmap$ES[VIP_MDD_hmap$ES < 0] <- -1


testlist <- c("PV_MDD_gseaM1", "PYR23_MDD_gseaM1", "PYR56_MDD_gseaM1", "SST_MDD_gseaM1", "VIP_MDD_gseaM1")

for(i in 1:length(testlist)){
  workingdata <- get(testlist[i])
  print(testlist[i])
  print(head(workingdata[workingdata$pval < 0.05,], 20))
}



#Up Venn
par(cex=2)
venn(list(PYR23_MDD_vennUp, PYR56_MDD_vennUp, SST_MDD_vennUp, PV_MDD_vennUp, VIP_MDD_vennUp), snames=c("PYR23", "PYR56", "SST", "PV", "VIP"), zcolor = c("#F8766D", "#b01005","#00BFC4","#7CAE00", "#C77CFF"), par=TRUE, box = FALSE) 

#Down Venn
par(cex=2)
venn(list(PYR23_MDD_vennDown, PYR56_MDD_vennDown, SST_MDD_vennDown, PV_MDD_vennDown, VIP_MDD_vennDown), snames=c("PYR23", "PYR56", "SST", "PV", "VIP"), zcolor = c("#F8766D", "#b01005","#00BFC4","#7CAE00", "#C77CFF"), par=TRUE, box = FALSE)



#Heatmap
names(PV_MDD_hmap) <- c("pathway", "PV_MDD_ES")
names(PYR23_MDD_hmap) <- c("pathway", "PYR23_MDD_ES")
names(PYR56_MDD_hmap) <- c("pathway", "PYR56_MDD_ES")
names(SST_MDD_hmap) <- c("pathway", "SST_MDD_ES")
names(VIP_MDD_hmap) <- c("pathway", "VIP_MDD_ES")

merged_MDD <- merge(PYR23_MDD_hmap, merge(PYR56_MDD_hmap, merge(SST_MDD_hmap, merge(PV_MDD_hmap, VIP_MDD_hmap, by="pathway", all=TRUE), by="pathway", all=TRUE), by="pathway", all=TRUE), by="pathway", all=TRUE)
merged_MDD[is.na(merged_MDD)] <- 0



Heatmap(merged_MDD[,c(2:6)], col=colorRamp2(c(-1, 0, 1),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"), border = TRUE)



#Bipolar
PV_Bipolar_SigUp <- PV_Bipolar_gseaM1[PV_Bipolar_gseaM1$pval < 0.05 & PV_Bipolar_gseaM1$ES > 0, ]
PV_Bipolar_SigDown <- PV_Bipolar_gseaM1[PV_Bipolar_gseaM1$pval < 0.05 & PV_Bipolar_gseaM1$ES < 0, ]
PV_Bipolar_vennUp <- PV_Bipolar_SigUp$pathway
PV_Bipolar_vennDown <- PV_Bipolar_SigDown$pathway
PV_Bipolar_hmap <- PV_Bipolar_gseaM1[PV_Bipolar_gseaM1$pval < 0.05, c("pathway", "ES")]
PV_Bipolar_hmap$ES[PV_Bipolar_hmap$ES > 0] <- 1
PV_Bipolar_hmap$ES[PV_Bipolar_hmap$ES < 0] <- -1

PYR23_Bipolar_SigUp <- PYR23_Bipolar_gseaM1[PYR23_Bipolar_gseaM1$pval < 0.05 & PYR23_Bipolar_gseaM1$ES > 0, ]
PYR23_Bipolar_SigDown <- PYR23_Bipolar_gseaM1[PYR23_Bipolar_gseaM1$pval < 0.05 & PYR23_Bipolar_gseaM1$ES < 0, ]
PYR23_Bipolar_vennUp <- PYR23_Bipolar_SigUp$pathway
PYR23_Bipolar_vennDown <- PYR23_Bipolar_SigDown$pathway
PYR23_Bipolar_hmap <- PYR23_Bipolar_gseaM1[PYR23_Bipolar_gseaM1$pval < 0.05, c("pathway", "ES")]
PYR23_Bipolar_hmap$ES[PYR23_Bipolar_hmap$ES > 0] <- 1
PYR23_Bipolar_hmap$ES[PYR23_Bipolar_hmap$ES < 0] <- -1

PYR56_Bipolar_SigUp <- PYR56_Bipolar_gseaM1[PYR56_Bipolar_gseaM1$pval < 0.05 & PYR56_Bipolar_gseaM1$ES > 0, ]
PYR56_Bipolar_SigDown <- PYR56_Bipolar_gseaM1[PYR56_Bipolar_gseaM1$pval < 0.05 & PYR56_Bipolar_gseaM1$ES < 0, ]
PYR56_Bipolar_vennUp <- PYR56_Bipolar_SigUp$pathway
PYR56_Bipolar_vennDown <- PYR56_Bipolar_SigDown$pathway
PYR56_Bipolar_hmap <- PYR56_Bipolar_gseaM1[PYR56_Bipolar_gseaM1$pval < 0.05, c("pathway", "ES")]
PYR56_Bipolar_hmap$ES[PYR56_Bipolar_hmap$ES > 0] <- 1
PYR56_Bipolar_hmap$ES[PYR56_Bipolar_hmap$ES < 0] <- -1

SST_Bipolar_SigUp <- SST_Bipolar_gseaM1[SST_Bipolar_gseaM1$pval < 0.05 & SST_Bipolar_gseaM1$ES > 0, ]
SST_Bipolar_SigDown <- SST_Bipolar_gseaM1[SST_Bipolar_gseaM1$pval < 0.05 & SST_Bipolar_gseaM1$ES < 0, ]
SST_Bipolar_vennUp <- SST_Bipolar_SigUp$pathway
SST_Bipolar_vennDown <- SST_Bipolar_SigDown$pathway
SST_Bipolar_hmap <- SST_Bipolar_gseaM1[SST_Bipolar_gseaM1$pval < 0.05, c("pathway", "ES")]
SST_Bipolar_hmap$ES[SST_Bipolar_hmap$ES > 0] <- 1
SST_Bipolar_hmap$ES[SST_Bipolar_hmap$ES < 0] <- -1

VIP_Bipolar_SigUp <- VIP_Bipolar_gseaM1[VIP_Bipolar_gseaM1$pval < 0.05 & VIP_Bipolar_gseaM1$ES > 0, ]
VIP_Bipolar_SigDown <- VIP_Bipolar_gseaM1[VIP_Bipolar_gseaM1$pval < 0.05 & VIP_Bipolar_gseaM1$ES < 0, ]
VIP_Bipolar_vennUp <- VIP_Bipolar_SigUp$pathway
VIP_Bipolar_vennDown <- VIP_Bipolar_SigDown$pathway
VIP_Bipolar_hmap <- VIP_Bipolar_gseaM1[VIP_Bipolar_gseaM1$pval < 0.05, c("pathway", "ES")]
VIP_Bipolar_hmap$ES[VIP_Bipolar_hmap$ES > 0] <- 1
VIP_Bipolar_hmap$ES[VIP_Bipolar_hmap$ES < 0] <- -1

#Up Venn
par(cex=2)
venn(list(PYR23_Bipolar_vennUp, PYR56_Bipolar_vennUp, SST_Bipolar_vennUp, PV_Bipolar_vennUp, VIP_Bipolar_vennUp), snames=c("PYR23", "PYR56", "SST", "PV", "VIP"), zcolor = c("#F8766D", "#b01005","#00BFC4","#7CAE00", "#C77CFF"), par=TRUE, box = FALSE) 

#Down Venn
par(cex=2)
venn(list(PYR23_Bipolar_vennDown, PYR56_Bipolar_vennDown, SST_Bipolar_vennDown, PV_Bipolar_vennDown, VIP_Bipolar_vennDown), snames=c("PYR23", "PYR56", "SST", "PV", "VIP"), zcolor = c("#F8766D", "#b01005","#00BFC4","#7CAE00", "#C77CFF"), par=TRUE, box = FALSE)



#Heatmap
names(PV_Bipolar_hmap) <- c("pathway", "PV_Bipolar_ES")
names(PYR23_Bipolar_hmap) <- c("pathway", "PYR23_Bipolar_ES")
names(PYR56_Bipolar_hmap) <- c("pathway", "PYR56_Bipolar_ES")
names(SST_Bipolar_hmap) <- c("pathway", "SST_Bipolar_ES")
names(VIP_Bipolar_hmap) <- c("pathway", "VIP_Bipolar_ES")

merged_Bipolar <- merge(PYR23_Bipolar_hmap, merge(PYR56_Bipolar_hmap, merge(SST_Bipolar_hmap, merge(PV_Bipolar_hmap, VIP_Bipolar_hmap, by="pathway", all=TRUE), by="pathway", all=TRUE), by="pathway", all=TRUE), by="pathway", all=TRUE)
merged_Bipolar[is.na(merged_Bipolar)] <- 0

#SCHIZ
PV_SCHIZ_SigUp <- PV_SCHIZ_gseaM1[PV_SCHIZ_gseaM1$pval < 0.05 & PV_SCHIZ_gseaM1$ES > 0, ]
PV_SCHIZ_SigDown <- PV_SCHIZ_gseaM1[PV_SCHIZ_gseaM1$pval < 0.05 & PV_SCHIZ_gseaM1$ES < 0, ]
PV_SCHIZ_vennUp <- PV_SCHIZ_SigUp$pathway
PV_SCHIZ_vennDown <- PV_SCHIZ_SigDown$pathway
PV_SCHIZ_hmap <- PV_SCHIZ_gseaM1[PV_SCHIZ_gseaM1$pval < 0.05, c("pathway", "ES")]
PV_SCHIZ_hmap$ES[PV_SCHIZ_hmap$ES > 0] <- 1
PV_SCHIZ_hmap$ES[PV_SCHIZ_hmap$ES < 0] <- -1

PYR23_SCHIZ_SigUp <- PYR23_SCHIZ_gseaM1[PYR23_SCHIZ_gseaM1$pval < 0.05 & PYR23_SCHIZ_gseaM1$ES > 0, ]
PYR23_SCHIZ_SigDown <- PYR23_SCHIZ_gseaM1[PYR23_SCHIZ_gseaM1$pval < 0.05 & PYR23_SCHIZ_gseaM1$ES < 0, ]
PYR23_SCHIZ_vennUp <- PYR23_SCHIZ_SigUp$pathway
PYR23_SCHIZ_vennDown <- PYR23_SCHIZ_SigDown$pathway
PYR23_SCHIZ_hmap <- PYR23_SCHIZ_gseaM1[PYR23_SCHIZ_gseaM1$pval < 0.05, c("pathway", "ES")]
PYR23_SCHIZ_hmap$ES[PYR23_SCHIZ_hmap$ES > 0] <- 1
PYR23_SCHIZ_hmap$ES[PYR23_SCHIZ_hmap$ES < 0] <- -1

PYR56_SCHIZ_SigUp <- PYR56_SCHIZ_gseaM1[PYR56_SCHIZ_gseaM1$pval < 0.05 & PYR56_SCHIZ_gseaM1$ES > 0, ]
PYR56_SCHIZ_SigDown <- PYR56_SCHIZ_gseaM1[PYR56_SCHIZ_gseaM1$pval < 0.05 & PYR56_SCHIZ_gseaM1$ES < 0, ]
PYR56_SCHIZ_vennUp <- PYR56_SCHIZ_SigUp$pathway
PYR56_SCHIZ_vennDown <- PYR56_SCHIZ_SigDown$pathway
PYR56_SCHIZ_hmap <- PYR56_SCHIZ_gseaM1[PYR56_SCHIZ_gseaM1$pval < 0.05, c("pathway", "ES")]
PYR56_SCHIZ_hmap$ES[PYR56_SCHIZ_hmap$ES > 0] <- 1
PYR56_SCHIZ_hmap$ES[PYR56_SCHIZ_hmap$ES < 0] <- -1

SST_SCHIZ_SigUp <- SST_SCHIZ_gseaM1[SST_SCHIZ_gseaM1$pval < 0.05 & SST_SCHIZ_gseaM1$ES > 0, ]
SST_SCHIZ_SigDown <- SST_SCHIZ_gseaM1[SST_SCHIZ_gseaM1$pval < 0.05 & SST_SCHIZ_gseaM1$ES < 0, ]
SST_SCHIZ_vennUp <- SST_SCHIZ_SigUp$pathway
SST_SCHIZ_vennDown <- SST_SCHIZ_SigDown$pathway
SST_SCHIZ_hmap <- SST_SCHIZ_gseaM1[SST_SCHIZ_gseaM1$pval < 0.05, c("pathway", "ES")]
SST_SCHIZ_hmap$ES[SST_SCHIZ_hmap$ES > 0] <- 1
SST_SCHIZ_hmap$ES[SST_SCHIZ_hmap$ES < 0] <- -1

VIP_SCHIZ_SigUp <- VIP_SCHIZ_gseaM1[VIP_SCHIZ_gseaM1$pval < 0.05 & VIP_SCHIZ_gseaM1$ES > 0, ]
VIP_SCHIZ_SigDown <- VIP_SCHIZ_gseaM1[VIP_SCHIZ_gseaM1$pval < 0.05 & VIP_SCHIZ_gseaM1$ES < 0, ]
VIP_SCHIZ_vennUp <- VIP_SCHIZ_SigUp$pathway
VIP_SCHIZ_vennDown <- VIP_SCHIZ_SigDown$pathway
VIP_SCHIZ_hmap <- VIP_SCHIZ_gseaM1[VIP_SCHIZ_gseaM1$pval < 0.05, c("pathway", "ES")]
VIP_SCHIZ_hmap$ES[VIP_SCHIZ_hmap$ES > 0] <- 1
VIP_SCHIZ_hmap$ES[VIP_SCHIZ_hmap$ES < 0] <- -1

#Up Venn
par(cex=2)
venn(list(PYR23_SCHIZ_vennUp, PYR56_SCHIZ_vennUp, SST_SCHIZ_vennUp, PV_SCHIZ_vennUp, VIP_SCHIZ_vennUp), snames=c("PYR23", "PYR56", "SST", "PV", "VIP"), zcolor = c("#F8766D", "#b01005","#00BFC4","#7CAE00", "#C77CFF"), par=TRUE, box = FALSE) 

#Down Venn
par(cex=2)
venn(list(PYR23_SCHIZ_vennDown, PYR56_SCHIZ_vennDown, SST_SCHIZ_vennDown, PV_SCHIZ_vennDown, VIP_SCHIZ_vennDown), snames=c("PYR23", "PYR56", "SST", "PV", "VIP"), zcolor = c("#F8766D", "#b01005","#00BFC4","#7CAE00", "#C77CFF"), par=TRUE, box = FALSE)



#Heatmap
names(PV_SCHIZ_hmap) <- c("pathway", "PV_SCHIZ_ES")
names(PYR23_SCHIZ_hmap) <- c("pathway", "PYR23_SCHIZ_ES")
names(PYR56_SCHIZ_hmap) <- c("pathway", "PYR56_SCHIZ_ES")
names(SST_SCHIZ_hmap) <- c("pathway", "SST_SCHIZ_ES")
names(VIP_SCHIZ_hmap) <- c("pathway", "VIP_SCHIZ_ES")

merged_SCHIZ <- merge(PYR23_SCHIZ_hmap, merge(PYR56_SCHIZ_hmap, merge(SST_SCHIZ_hmap, merge(PV_SCHIZ_hmap, VIP_SCHIZ_hmap, by="pathway", all=TRUE), by="pathway", all=TRUE), by="pathway", all=TRUE), by="pathway", all=TRUE)
merged_SCHIZ[is.na(merged_SCHIZ)] <- 0





merged_pval <- merge(merged_MDD, merge(merged_Bipolar, merged_SCHIZ, by="pathway", all=TRUE), by="pathway", all=TRUE)
merged_pval[is.na(merged_pval)] <- 0



Heatmap(merged_pval[,c(2:16)], column_split = c(rep("MDD", 5), rep("BD", 5), rep("SCZ", 5)),
        col=colorRamp2(c(-1, 0, 1),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = FALSE, show_row_names = FALSE, width = unit(8, "cm"), border = TRUE)

merged_pval2 <- merged_pval[,c(1,2,7,12,3,8,13,4,9,14,5,10,15,6,11,16)]

Heatmap(merged_pval2[,c(2:16)], column_split = c(rep("PYR L2/3", 3), rep("PYR L5/6", 3), rep("SST", 3), rep("PV", 3), rep("VIP", 3)),
        col=colorRamp2(c(-1, 0, 1),  c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, show_row_dend = TRUE, show_row_names = FALSE, width = unit(8, "cm"), border = TRUE)

#CT-specific heatmaps
PYR23_hmap <- merge(PYR23_MDD_hmap, merge(PYR23_Bipolar_hmap, PYR23_SCHIZ_hmap, by="pathway", all=TRUE), by="pathway", all=TRUE)
PYR56_hmap <- merge(PYR56_MDD_hmap, merge(PYR56_Bipolar_hmap, PYR56_SCHIZ_hmap, by="pathway", all=TRUE), by="pathway", all=TRUE)
SST_hmap <- merge(SST_MDD_hmap, merge(SST_Bipolar_hmap, SST_SCHIZ_hmap, by="pathway", all=TRUE), by="pathway", all=TRUE)
PV_hmap <- merge(PV_MDD_hmap, merge(PV_Bipolar_hmap, PV_SCHIZ_hmap, by="pathway", all=TRUE), by="pathway", all=TRUE)
VIP_hmap <- merge(VIP_MDD_hmap, merge(VIP_Bipolar_hmap, VIP_SCHIZ_hmap, by="pathway", all=TRUE), by="pathway", all=TRUE)

length(grep("GO.+.+", PYR23_hmap$pathway))/nrow(PYR23_hmap)
length(grep("GO.+.+", PYR56_hmap$pathway))/nrow(PYR56_hmap)
length(grep("GO.+.+", SST_hmap$pathway))/nrow(SST_hmap)
length(grep("GO.+.+", PV_hmap$pathway))/nrow(PV_hmap)
length(grep("GO.+.+", VIP_hmap$pathway))/nrow(VIP_hmap)


PYR23_hmap$pathway <- gsub("[^[:alpha:][:digit:]]+", "\\.", PYR23_hmap$pathway)
PYR56_hmap$pathway <- gsub("[^[:alpha:][:digit:]]+", "\\.", PYR56_hmap$pathway)
SST_hmap$pathway <- gsub("[^[:alpha:][:digit:]]+", "\\.", SST_hmap$pathway)
PV_hmap$pathway <- gsub("[^[:alpha:][:digit:]]+", "\\.", PV_hmap$pathway)
VIP_hmap$pathway <- gsub("[^[:alpha:][:digit:]]+", "\\.", VIP_hmap$pathway)

PYR23_hmap[is.na(PYR23_hmap)] <- 0
PYR56_hmap[is.na(PYR56_hmap)] <- 0
SST_hmap[is.na(SST_hmap)] <- 0
PV_hmap[is.na(PV_hmap)] <- 0
VIP_hmap[is.na(VIP_hmap)] <- 0


mPYR23_hmap <- as.matrix(PYR23_hmap[,-1])
mPYR56_hmap <- as.matrix(PYR56_hmap[,-1])
mSST_hmap <- as.matrix(SST_hmap[,-1])
mPV_hmap <- as.matrix(PV_hmap[,-1])
mVIP_hmap <- as.matrix(VIP_hmap[,-1])


rownames(mPYR23_hmap) <- PYR23_hmap$pathway
rownames(mPYR56_hmap) <- PYR56_hmap$pathway
rownames(mSST_hmap) <- SST_hmap$pathway
rownames(mPV_hmap) <- PV_hmap$pathway
rownames(mVIP_hmap) <- VIP_hmap$pathway


#CT-specific heatmaps


Heatmap(mPYR23_hmap, 
        column_labels = c("MDD", "BPD", "SCZ"), 
        column_title = "PYR L2/3", col=colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, 
        show_row_dend = TRUE, 
        show_row_names = FALSE, 
        width = unit(8, "cm"), 
        border = TRUE,
        column_names_rot = 0,
        column_names_side = "top")


Heatmap(mPYR56_hmap, 
        column_labels = c("MDD", "BPD", "SCZ"), 
        column_title = "PYR L5/6", col=colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, 
        show_row_dend = TRUE, 
        show_row_names = FALSE, 
        width = unit(8, "cm"), 
        border = TRUE,
        column_names_rot = 0,
        column_names_side = "top")

Heatmap(mSST_hmap, 
        column_labels = c("MDD", "BPD", "SCZ"), 
        column_title = "SST", col=colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, 
        show_row_dend = TRUE, 
        show_row_names = FALSE, 
        width = unit(8, "cm"), 
        border = TRUE,
        column_names_rot = 0,
        column_names_side = "top")

Heatmap(mPV_hmap, 
        column_labels = c("MDD", "BPD", "SCZ"), 
        column_title = "PV", col=colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, 
        show_row_dend = TRUE, 
        show_row_names = FALSE, 
        width = unit(8, "cm"), 
        border = TRUE,
        column_names_rot = 0,
        column_names_side = "top")

Heatmap(mVIP_hmap, 
        column_labels = c("MDD", "BPD", "SCZ"), 
        column_title = "VIP", col=colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, 
        show_row_dend = TRUE, 
        show_row_names = FALSE, 
        width = unit(8, "cm"), 
        border = TRUE,
        column_names_rot = 0,
        column_names_side = "top")


#Disorder heatmaps (MDD)
MDD_hmap <- merge(PYR23_MDD_hmap, merge(PYR56_MDD_hmap, merge(SST_MDD_hmap, merge(PV_MDD_hmap, VIP_MDD_hmap, by="pathway", all=TRUE), by="pathway", all=TRUE), by="pathway", all=TRUE), by="pathway", all=TRUE)


MDD_hmap$pathway <- gsub("[^[:alpha:][:digit:]]+", "\\.", MDD_hmap$pathway)
MDD_hmap[is.na(MDD_hmap)] <- 0
mMDD_hmap <- as.matrix(MDD_hmap[,-1])


rownames(mMDD_hmap) <- MDD_hmap$pathway


Heatmap(mMDD_hmap, 
        column_labels = c("L2/3 PYR", "L5/6 PYR", "SST", "PV", "VIP"), 
        column_title = "MDD", col=colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, 
        show_row_dend = TRUE, 
        show_row_names = FALSE, 
        width = unit(8, "cm"), 
        border = TRUE,
        column_names_rot = 0,
        column_names_side = "top")






#Get GO terms to use in par-child analysis, by text search
library(fgsea)
library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(fdrtool)
library(ggplot2)
library(reshape2)
library(dplyr)
library(GGally)
library(Vennerable)
library(gplots)
library(venn)
options(stringsAsFactors = FALSE)
options(scipen = 999)

setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/")
GMTfile <- gmtPathways("Human_GO_AllPathways_with_GO_iea_April_01_2020_symbol.gmt")

setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/Parent child/")

rawParentTerms <- read.csv("GO text search results.csv")

#Some how you can't get the definition field - GO downloads are just trash # test <- get_ontology("go-plus.obo") # #Need to remove obsolete terms, CHEBI/non-GO ID's, and other general clean-up

#Before analysis: convert GMT (only GO data) and raw parent terms to dataframes, integrate GMT data into raw parent term data and:
#Test different size filters and determine overlap across categories
GMT_working <- lapply(GMTfile, paste0, collapse=" | ")

GMT_df <- data.frame("Pathway" = names(GMT_working), "Genes" = unlist(GMT_working))
#Remove all non GO stuff
GMT_df <- GMT_df[(grepl("GOBP", GMT_df$Pathway) | grepl("GOMF", GMT_df$Pathway) | grepl("GOCC", GMT_df$Pathway)), ]
GMT_df$Size <- sapply(strsplit(GMT_df$Genes, split="|", fixed=TRUE), length)
GMT_df$Ontology <- sub("%.*$", "", sub("^.*%+?", "", GMT_df$Pathway))
GMT_df$GOID <- sub("^.*%", "", GMT_df$Pathway)
GMT_df$Name <- sub("%.*$", "", GMT_df$Pathway)

#Add genelists to "rawParentTerms"
ParentTerms <- merge(rawParentTerms, GMT_df, by="GOID", all.x=TRUE)
ParentTerms <- ParentTerms[complete.cases(ParentTerms),]

#Standardize parentterm function names
unique(ParentTerms$Function)

ParentTerms$Function[ParentTerms$Function == "Mitochondrial Structure and Function"] <- "mt_function"
ParentTerms$Function[ParentTerms$Function == "Channels/Transporters"] <- "channels_trans"
ParentTerms$Function[ParentTerms$Function == "ATP generation_Bioenergetics"] <- "bioenergetics"
ParentTerms$Function[ParentTerms$Function == "Transcription" ] <- "transcription"
ParentTerms$Function[ParentTerms$Function == "DNA Modification/Epigenetics" ] <- "DNA_modification"
ParentTerms$Function[ParentTerms$Function == "Golgi"] <- "golgi"
ParentTerms$Function[ParentTerms$Function == "Protein Modification"] <- "prot_modification"
ParentTerms$Function[ParentTerms$Function == "Mt-transcription"] <- "mt_transcription"
ParentTerms$Function[ParentTerms$Function == "Translation"] <- "translation"
ParentTerms$Function[ParentTerms$Function == "Synaptic Vesicles"] <- "gen_synaptic"
ParentTerms$Function[ParentTerms$Function == "ER"] <- "er"
ParentTerms$Function[ParentTerms$Function == "Synaptic Function"] <- "gen_synaptic"
ParentTerms$Function[ParentTerms$Function == "Protein Folding"] <- "prot_folding"
ParentTerms$Function[ParentTerms$Function == "Neuropeptides"] <- "neuropeptides"
ParentTerms$Function[ParentTerms$Function == "Matrix"] <- "matrix"
ParentTerms$Function[ParentTerms$Function == "Monoamines"] <- "monoamines"
ParentTerms$Function[ParentTerms$Function == "Synaptic Structure"] <- "synapse_structure"
ParentTerms$Function[ParentTerms$Function == "Inflammatory/Immune Signaling"] <- "inflammation"
ParentTerms$Function[ParentTerms$Function == "Cytoskeleton"] <- "cytoskeleton"
ParentTerms$Function[ParentTerms$Function == "GABA"] <- "gaba"
ParentTerms$Function[ParentTerms$Function == "Growth Factors"] <- "growth_fac"
ParentTerms$Function[ParentTerms$Function == "Neurotrophins"] <- "neurotroph"
ParentTerms$Function[ParentTerms$Function == "Nucleus"] <- "nucleus"
ParentTerms$Function[ParentTerms$Function == "Mt-translation"] <- "mt_translation"
ParentTerms$Function[ParentTerms$Function == "Lysosome"] <- "cytoskeleton"
ParentTerms$Function[ParentTerms$Function == "Glutamate"] <- "glutamate"


hist(ParentTerms$Size[ParentTerms$Size < 1000 & ParentTerms$Size > 100])
length(ParentTerms$Size[ParentTerms$Size < 1000 & ParentTerms$Size > 100])

#Check size distribution of sig. results
merged_pval <- merge(merged_MDD, merge(merged_Bipolar, merged_SCHIZ, by="pathway", all=TRUE), by="pathway", all=TRUE)
merged_pval[is.na(merged_pval)] <- 0
merged_pval$GOID <- sub("^.*%", "", merged_pval$pathway)
merged_pval3 <- merge(merged_pval, GMT_df, by="GOID", all.x=TRUE)
merged_pval3 <- merged_pval3[complete.cases(merged_pval3),]
#dev.off()
hist(merged_pval3$Size)

######################## To be used to search against in text-query
ParentTerms_searchable <- ParentTerms[ParentTerms$Size > 15 & ParentTerms$Size < 500, ]

#########################To be used as parent-child starting point
#Number and Overlap of GO terms across categories after filtering to size of 50-1000
ParentTerms_filt <- ParentTerms[ParentTerms$Size > 50 & ParentTerms$Size < 1000, ]
ParentTerms_filt$Function <- as.factor(ParentTerms_filt$Function)

ParentTerms_export <- ParentTerms[ParentTerms$Size > 15 & ParentTerms$Size < 1000, ]

#Distribution shows ~half of functions have <20 terms, and ~half have 20-100
hist(summary(ParentTerms_filt$Function))
summary(as.factor(ParentTerms_filt$Function))

sum(duplicated(ParentTerms_filt$GOID))

library(GeneOverlap)
gsSummary <- rep(list("A"), length(unique(ParentTerms_filt$Function)))
functions <- as.character(unique(ParentTerms_filt$Function))
names(gsSummary) <- functions
for(i in 1:length(functions)){
  gsSummary[[i]] <- ParentTerms_filt$GOID[ParentTerms_filt$Function == functions[i]]
}


gsOverlap <- newGOM(gsSummary, gsSummary)
gsOutput <- getMatrix(gsOverlap, name="intersection")
#Only appreciable overlap is synaptic vesicles and synaptic function

#Manual rearranging and filtering
merged_pval3 <- merged_pval3[,c(1,2,3,8,13,4,9,14,5,10,15,6,11,16,7,12,17,18:22)]

save(GMT_df, GMT_working, gsOutput, merged_pval3, ParentTerms, ParentTerms_export, ParentTerms_filt, ParentTerms_searchable, file="ParentChild input files.rData")


################################################### Start here for actual parent-child analysis

##############################################   
library(fgsea)
library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(fdrtool)
library(ggplot2)
library(reshape2)
library(dplyr)
library(GGally)
library(Vennerable)
library(gplots)
library(venn)
options(stringsAsFactors = FALSE)
options(scipen = 999)

setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/Parent child")
load("ParentChild input files.rData")

#First, do text-based search of results based on search terms:
#filter out matching results and put in holding list
#need to integrate GMT data 

# toMatch <- c("A1", "A9", "A6")
# matches <- unique(grep(paste(toMatch,collapse="|"), myfile$Letter, value=TRUE))

transcription <- c("transcription")
translation <- c("translat")
mt_transcription <- c("mitochondrial transcri")
mt_translation <- c("mitochondrial translat")
bioenergetics <- c("oxidative phos", "glycoly", "ATP" )
#For this, remove bioenergetics, mt translation, and mt transcription
mt_function <- c("mitochon")
prot_folding <- c("chaperone", "heat shock", "misfolded", "unfolded")
prot_modification <- c("post translational", "post-translational", "protein modification")
DNA_modification <- c("histone", "DNA modification", "DNA methyl")

####CC only
nucleus <- c("nuclear", "nucleus", "nucleolus", "nucleolar")
er <- c("endoplasmic")
golgi <- c("golgi")
lysosome <- c("lysosome")
cytoskeleton <- c("cytoskeleton", "cytoskeletal")
#Remove "mitochondria"
matrix <- c("matrix", "extracellular matrix", "neuropil", "perineuronal", "peri-neuronal")
synapse_structure <- c("synapse", "synaptic", "neuron proj", "neurite", "axon", "dendrite", "dendritic")

###BP/MF only
gen_synaptic <- c("synapse", "synaptic", "neurite", "axon", "dendrite", "dendritic", "synaptic ves", "vesicle")
# synaptic_ves <- c("synaptic ves", "vesicle") #largely redundant with gen_synaptic
gaba <- c("GABA")
glutamate <- c("glutamate synapse", "glutamatergic synapse")
channels_trans <- c("ion transport", "ion channel", "channel", "transporter")
monoamines <- c("serotonin", "serotonergic", "noradrenaline", "noradrenerguc", "adrenergic", "dopamine", "dopminergic")
neuropeptides <- c("neuropeptide", "somatostatin", "cholecystokinin", "substance P", "Enkephalin", "vasoactive intestinal peptide" )
growth_fac <- c("growth factor")
neurotroph <- c("neurotrophin", "neurotrophic")
inflammation <- c("inflamm", "immun")

functionlist <- c("transcription","translation","mt_transcription","mt_translation","bioenergetics","mt_function","prot_folding","prot_modification","DNA_modification","nucleus","er","golgi","lysosome","cytoskeleton","matrix", "synapse_structure","gen_synaptic","gaba","glutamate","channels_trans","monoamines","neuropeptides","growth_fac","neurotroph","inflammation")


ontologylist <- list(c("GOBP","GOMF","GOCC"),
                     c("GOBP","GOMF","GOCC"),
                     c("GOBP","GOMF","GOCC"),
                     c("GOBP","GOMF","GOCC"),
                     c("GOBP","GOMF","GOCC"),
                     c("GOBP","GOMF","GOCC"),
                     c("GOBP","GOMF","GOCC"),
                     c("GOBP","GOMF","GOCC"),
                     c("GOBP","GOMF","GOCC"),
                     c("GOCC"),
                     c("GOCC"),
                     c("GOCC"),
                     c("GOCC"),
                     c("GOCC"),
                     c("GOCC"),
                     c("GOCC"),
                     c("GOBP","GOMF"),
                     c("GOBP","GOMF"),
                     c("GOBP","GOMF"),
                     c("GOBP","GOMF"),
                     c("GOBP","GOMF"),
                     c("GOBP","GOMF"),
                     c("GOBP","GOMF"),
                     c("GOBP","GOMF"),
                     c("GOBP","GOMF"))
names(ontologylist) <- functionlist

exclusionlist <- list(c("mitochondria"),
                   c("mitochondria"),
                   c(""),
                   c(""),
                   c(""),
                   c("oxidative phos"," glycoly"," ATP"," mitochondrial translat"," mitochondrial transcri"),
                   c(""),
                   c(""),
                   c(""),
                   c(""),
                   c(""),
                   c(""),
                   c(""),
                   c(""),
                   c("mitochondria"),
                   c(""),
                   c(""),
                   c(""),
                   c(""),
                   c(""),
                   c(""),
                   c(""),
                   c(""),
                   c(""),
                   c(""))
names(exclusionlist) <- functionlist


#For each results column in merged_pval3:
#for each function in functionlist, grep all matches (tolower) in mergedpval3, and get all matches between ParentTerms_searchable and mergedpval3 (for either of these, only include results from ontologies belonging to "ontology list")
#Remove duplicated terms from these two searches


#Results table, containing pathway name and NES for each contrast - 18 cols: 1 for functional group, 1 for pathway, 1 for GOID, 15 for NES's
workingTextquery <- as.data.frame(matrix(nrow = 0, ncol = 18))
names(workingTextquery) <- c("FunctionalGroup", "Pathway", "GOID", names(merged_pval3)[3:17])

GMT_df$Pathway <- tolower(GMT_df$Pathway)
# write.csv(merged_pval3, "merged_pval3_test.csv")
merged_pval3 <- as.data.frame(merged_pval3)
merged_pval3$pathway <- tolower(merged_pval3$pathway)

merged_pval_leftover <- merged_pval3

#Export duplicated parent terms - manually reassign
duplicated_parents <- ParentTerms_export[duplicated(ParentTerms_export$GOID) | duplicated(ParentTerms_export$GOID, fromLast=TRUE),]
duplicated_parents <- duplicated_parents[,c(2,1,4)]
write.csv(duplicated_parents, "actually duplicated parents (export).csv")


########import the un-duplicated terms, and filter out duplicates which were not selected - also do an NA check, as parent terms not intuitive for any overlap were just removed

#Manual de-duplication was performed
unduplicated_parents <- read.csv("actually unduplicated parents (export).csv")

for(i in 1:nrow(unduplicated_parents)){
  ParentTerms_export <- ParentTerms_export[!(ParentTerms_export$Function != unduplicated_parents$Function[i] & ParentTerms_export$GOID == unduplicated_parents$GOID[i]),]
}

#don't forget NA check
ParentTerms_export <- ParentTerms_export[!(ParentTerms_export$GOID %in% unduplicated_parents$GOID[is.na(unduplicated_parents$Function)]),]
ParentTerms_export <- ParentTerms_export[!is.na(ParentTerms_export$Function),]

#After, remove duplicates within functions (only gen_synaptic)
ParentTerms_export <- ParentTerms_export[!duplicated(ParentTerms_export$GOID),]

# then re-generate the "searchable" and "filt" objects based on size
ParentTerms_filt <- ParentTerms_export[ParentTerms_export$Size > 50 & ParentTerms_export$Size < 1000, ]
ParentTerms_searchable <- ParentTerms_export[ParentTerms_export$Size > 15 & ParentTerms_export$Size < 500, ]

for(i in 1:length(functionlist)){
  for (j in 1:length(names(merged_pval3)[3:17])){
    currentcolumn <- names(merged_pval3)[j+2]
    searchterms <- tolower(get(functionlist[i]))
    currentontology <- ontologylist[[i]]
    currentexclusion <- exclusionlist[[i]]
    
    #grep-search (by selected ontology)
    grepsearch <- unique(grep(paste0(searchterms,collapse="|"), tolower(merged_pval3$Pathway[merged_pval3[,j+2] != 0]), value=TRUE))
    
    #keep only those of proper ontologies
    grepsearch <- grep(paste0(tolower(currentontology),collapse="|"), grepsearch, value=TRUE)
    
    #Get matches based on those already identified as parent-terms
    overlap_terms <- tolower(ParentTerms_searchable$Pathway[ParentTerms_searchable$Function == functionlist[i]])
    overlapsearch <- tolower(merged_pval3$pathway[merged_pval3[,j+2] != 0][tolower(merged_pval3$pathway[merged_pval3[,j+2] != 0]) %in% overlap_terms])
    
    #Merge and de-duplicate results
    mergedsearch <- unique(c(overlapsearch, grepsearch))
    
    #exclude after merge and de-duplication of results
    if(length(currentexclusion) == 1 & currentexclusion != ""){
      #remove only if necessary
      mergedsearch <- mergedsearch[!grepl(paste0(currentexclusion,collapse="|"), mergedsearch)]
    }
    if(length(mergedsearch) != 0){
      holdingDF <- data.frame("FunctionalGroup" = functionlist[i], "Pathway" = mergedsearch)
      
      #Integrate GMT and directional data
      holdingDF <- merge(holdingDF, merged_pval3[,c(1,2,j+2)], by.x="Pathway", by.y="pathway", all.x=TRUE)
      holdingDF <- holdingDF[,c(2,1,3,4)]
      
      #merge to results object
      workingTextquery <- merge(workingTextquery, holdingDF[,c(1:3)], all=TRUE)[, union(names(workingTextquery), names(holdingDF[,c(1:3)]))]
      workingTextquery[,j+3][workingTextquery$FunctionalGroup==functionlist[i] & (workingTextquery$Pathway %in% holdingDF$Pathway)] <- holdingDF[holdingDF$Pathway %in% workingTextquery$Pathway , 4]
      
      #remove positive hits from merged_pval_leftover
      merged_pval_leftover[merged_pval_leftover$pathway %in% holdingDF$Pathway,j+2] <- 0
    }
  }
}

#After everything, run "complete cases" type of check on merged_pval_leftover to remove any terms no longer pending across all contrasts

dim(merged_pval_leftover)
merged_pval_leftover <- merged_pval_leftover[rowSums(merged_pval_leftover[,c(3:17)] == 0) < ncol(merged_pval_leftover[,c(3:17)]), ]

#looks good, all accounted for after taking into account terms that are found across multiple functions
dim(merged_pval_leftover)

ParentChildInput <- merged_pval_leftover[,-1]


#Do Parent-Child analysis on rest of results (ParentChildInput)
#Before the loop - "toupper" the pathway column
#Run Parent-child seperately for each contrast - getting multiple truthtables/summaries (1/contrast): truthtable is the important file
#At the start of each loop, pull out the necessary resultsGroupsup/down from ParentChildInput


library(readr)
library(dplyr)
library(magrittr)
library(GO.db)
library(reshape2)
library(plyr)
setwd("Text and ParentChild output/")

ParentChildInput$X1 <- toupper(ParentChildInput$pathway)
ParentChildInput <- ParentChildInput[,c(22,2:16)]

OntologyMap <- list(as.list(GOBPANCESTOR),
                    as.list(GOMFANCESTOR),
                    as.list(GOCCANCESTOR))
CurrOntology <- c("BP", "MF", "CC")

#SelectedGroups in the format of: "GOID", and "Name": name being just the GO-term name
selectedGroups_func <- ParentTerms_filt[,c(2,1,4)]
selectedGroups <- selectedGroups_func[,c(2,3)]

################### This loop is based on code written by Leon French for the publication: "doi: 10.1016/j.biopsych.2018.09.019"

for(i in 1:length(names(ParentChildInput[,2:16]))){
  currentcontrast <- names(ParentChildInput[i+1])
  #Loop for each ontology - after running everything, reattach ontology information and filter accordingly
  for(j in 1:length(OntologyMap)){
    resultGroupsUp <- ParentChildInput[ParentChildInput[,i+1] > 0,c(1,i+1)]
    resultGroupsDown <- ParentChildInput[ParentChildInput[,i+1] < 0,c(1,i+1)]
    resultGroups <- bind_rows(resultGroupsUp, resultGroupsDown)
    resultGroups %<>% filter(grepl("GO:", X1)) %>% mutate(goID = gsub(".*%","" , X1))
    WorkingOntology <- OntologyMap[[j]]
    
    #############################Parent-child script as before
    # Remove GO IDs that do not have any ancestor (there are none!)
    # Match GOID from the results and filter it
    # Sort in decreasing order
    # Make a dataframe of ancestor count for all the GOIDs in our result
    # Rename the column as "occuranceInTable"
    # Add one more column of GOID (why?)
    # Filter out "all"
    WorkingOntology <- WorkingOntology[!is.na(WorkingOntology)] 
    goAncestorCounts <- table(unlist(WorkingOntology[resultGroups$goID])) 
    goAncestorCounts <- sort(goAncestorCounts, decreasing = T ) 
    df <- data.frame(matrix(unlist(goAncestorCounts), nrow=length(goAncestorCounts), byrow=T), row.names = names(goAncestorCounts))
    head(df)
    colnames(df)[1] <- "occuranceInTable"
    df$GOID <- names(goAncestorCounts)
    head(df)
    df <- tbl_df(df) %>% filter(GOID != "all") 
    head(df)
    df <- mutate(rowwise(df), ancestorCount = length(unlist(WorkingOntology[GOID])))
    head(df)
    df <- mutate(rowwise(df), name = Term(GOID)) %>% arrange(ancestorCount, desc(occuranceInTable))
    head(df)
    df$percentOfGOgroups <- df$occuranceInTable/nrow(resultGroups) * 100
    head(df)
    t(WorkingOntology[resultGroups$goID])
    
    ##################################################
    resultGroups <- as.data.frame(resultGroups)
    head(rownames(resultGroups))
    # For each selectedGOID that is in selectedGroups$GOID print the GOID and the associated GO terms
    for(selectedGOID in selectedGroups$GOID) { #for each targetted GO group
      print(selectedGOID)
      print(Term(selectedGOID))
    }
    
    # For each rowNumber that is in resultGroups (and given by commands rownames(resultGroups)) print the GOID and the associated GO terms
    
    for(rowNumber in rownames(resultGroups)) { #for each row in the result table
      print(goID <- resultGroups[rowNumber, "goID"])
      print(rowName <- Term(selectedGOID))
    }
    if(is.na(rowName)) { rowName <- selectedGOID }
    resultGroups[rowNumber, rowName]= selectedGOID %in% unlist(WorkingOntology[goID]) #test if it is one of the ancestors/parents
    #######################  
    resultGroups <- as.data.frame(resultGroups)
    head(rownames(resultGroups))
    for(selectedGOID in selectedGroups$GOID) { #for each targett  ed GO group
      print(selectedGOID)
      print(Term(selectedGOID))
      for(rowNumber in rownames(resultGroups)) { #for each row in the result table
        goID <- resultGroups[rowNumber, "goID"]
        rowName <- Term(selectedGOID)
        if(is.na(rowName)) { rowName <- selectedGOID }
        resultGroups[rowNumber, rowName]= selectedGOID %in% unlist(WorkingOntology[goID]) #test if it is one of the ancestors/parents
      }
    }
    
    write.csv(resultGroups,paste0("Sandbox_TruthTable_", currentcontrast, "_", CurrOntology[j],".csv"))
  }
}

#Afterwards, need to reduce truthtables into enriched GO-terms belonging to functions of interest
#Load each truth table, extract contrast information from name, reduce table down to columns with colsums > 0 and rowsums > 0
TTlist <- list.files(pattern = "TruthTable")
reducedTTresults <- as.data.frame(matrix(nrow=0, ncol=2))
names(reducedTTresults) <- c("Function", "Pathway")

ParentTerms_filt$Match <- gsub("[^[:alpha:][:digit:]]+", "\\.", tolower(ParentTerms_filt$Name))
#Address multiple periods by reducing them all to 1
ParentTerms_filt$Match <- gsub("([\\.])\\1+", "\\1", ParentTerms_filt$Match)

for(i in 1:length(TTlist)){
  workingTT <- read.csv(TTlist[i])
  workingTT <- workingTT[,-1]
  TTcontrast <- gsub("^.*Table_", "", TTlist[i])
  TTcontrast <- substr(TTcontrast, 1, nchar(TTcontrast)-7)
  
  #Remove anything not captured by parent-child analysis
  workingTT <- workingTT[,!c(rep(FALSE, 3), colSums((workingTT[,c(4:ncol(workingTT))])) == 0)]
  if(ncol(workingTT) > 3){
    workingTT <- workingTT[!(rowSums(as.data.frame(workingTT[,c(4:ncol(workingTT))])) == 0),]
    
    #For each column (parent-term), identify which terms are enriched, and which function the parent-term is from - create df that is merged back to "reducedTTresults"
    for (j in colnames(workingTT[4:ncol(workingTT)])){
      workingFunction <- as.character(ParentTerms_filt$Function[tolower(ParentTerms_filt$Match) %in% gsub("([\\.])\\1+", "\\1", tolower(j))])
      workingRes <- data.frame("Function" = rep(workingFunction, sum(workingTT[,j])))
      workingPathways <- workingTT$X1[(workingTT[,j])]
      workingRes$Pathway <- workingPathways
      workingRes[,TTcontrast] <- workingTT[workingTT$X1 %in% workingRes$Pathway, 2]
      reducedTTresults <- merge(reducedTTresults, workingRes, all=TRUE)
    }
  }

} 



#Afterwards, add with textquery results, and reduce duplicates within each function, assess overlap using geneOverlap package, and generate confidence/redundancy scores

names(workingTextquery)
names(workingTextquery)[1] <- "Function"
workingTextquery$Pathway <- toupper(workingTextquery$Pathway)

names(reducedTTresults)
reducedTTresults <- reducedTTresults[,c(1,2,13,14,12,10,11,9,7,8,6,16,17,15,4,5,3)]
reducedTTresults$GOID <- gsub("^.*%", "", reducedTTresults$Pathway)
reducedTTresults  <- reducedTTresults[,c(1,2,18,3:17)]

sapply(workingTextquery, class)
sapply(reducedTTresults, class)

names(workingTextquery) == names(reducedTTresults)

#Merge together
mergedOutput <- rbind(workingTextquery, reducedTTresults)

sum(duplicated(mergedOutput$Pathway))

#Normalize hits across methods (i.e. directionality not showing up for all contrasts b/c of merging issues)
FunctionUnique <- unique(mergedOutput$Function)

for(i in 1:length(FunctionUnique)){
  PathwayUnique <- unique(mergedOutput$Pathway[mergedOutput$Function == FunctionUnique[i]])
  for(j in 1:length(PathwayUnique)){
    res <- mergedOutput[mergedOutput$Function == FunctionUnique[i] & mergedOutput$Pathway == PathwayUnique[j], ]
    if(nrow(res) > 1){
      for(k in 1:ncol(mergedOutput[,c(4:18)])){
        replacement <- mergedOutput[(mergedOutput$Function == FunctionUnique[i] & mergedOutput$Pathway == PathwayUnique[j]), k+3]
        replacement <- replacement[!is.na(replacement)]
        
        if(length(replacement) > 0){
          mergedOutput[mergedOutput$Function == FunctionUnique[i] & mergedOutput$Pathway == PathwayUnique[j], k+3] <- replacement
        }
      }
    }
  }
}

test <- mergedOutput[duplicated(mergedOutput$Pathway) | duplicated(mergedOutput$Pathway, fromLast=TRUE),]
  
#Remove wrong ontologies
# ontologylist
# GMT_df$Pathway <- toupper(GMT_df$Pathway)
# mergedOutput <- merge(mergedOutput, GMT_df)[, union(names(mergedOutput), names(GMT_df))]
# 
# for(i in 1:length(ontologylist)){
#   currFunc <- functionlist[i]
#   ontologies <- ontologylist[[i]]
#   mergedOutput <- mergedOutput[!(isFALSE(mergedOutput$Ontology[mergedOutput$Function == currFunc] %in% ontologies)),]
# }

mergedOutput2 <- mergedOutput %>% group_by(Function , Pathway , GOID , PYR23_MDD_ES , PYR23_Bipolar_ES , PYR23_SCHIZ_ES , PYR56_MDD_ES , PYR56_Bipolar_ES , PYR56_SCHIZ_ES , SST_MDD_ES , SST_Bipolar_ES , SST_SCHIZ_ES , PV_MDD_ES , PV_Bipolar_ES , PV_SCHIZ_ES , VIP_MDD_ES , VIP_Bipolar_ES , VIP_SCHIZ_ES)%>%distinct(`Pathway`)

sum(duplicated(mergedOutput2$GOID))

#assess overlap using geneOverlap package
library(GeneOverlap)
gsMergedSummary <- rep(list("A"), length(unique(mergedOutput2$Function)))
functions <- as.character(unique(mergedOutput2$Function))
names(gsMergedSummary) <- functions
for(i in 1:length(functions)){
  gsMergedSummary[[i]] <- mergedOutput2$GOID[mergedOutput2$Function == functions[i]]
}


gsMergedOverlap <- newGOM(gsMergedSummary, gsMergedSummary)
gsmergedOutput <- getMatrix(gsMergedOverlap, name="intersection")

write.csv(gsmergedOutput, "Text-Parent-Child output - overlap across functions.csv")



#generate confidence/redundancy scores to address potential effects of redundancy in parent-child hits
#Get leading edges
setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/")

load("GSEA results (GSEAparam1).rData")

names(mergedOutput2)
resultsList <- c("PYR23_MDD_gseaM1", "PYR23_Bipolar_gseaM1", "PYR23_SCHIZ_gseaM1", "PYR56_MDD_gseaM1", "PYR56_Bipolar_gseaM1", "PYR56_SCHIZ_gseaM1", "SST_MDD_gseaM1", "SST_Bipolar_gseaM1", "SST_SCHIZ_gseaM1", "PV_MDD_gseaM1", "PV_Bipolar_gseaM1", "PV_SCHIZ_gseaM1", "VIP_MDD_gseaM1", "VIP_Bipolar_gseaM1", "VIP_SCHIZ_gseaM1")

setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/Parent child/Text and ParentChild output")
GMT_df$Pathway <- toupper(GMT_df$Pathway)
mergedOutput2 <- merge(mergedOutput2, GMT_df)[, union(names(mergedOutput2), names(GMT_df))]


#Generate scores
confidenceScores <- as.data.frame(matrix(nrow=25, ncol=16))
names(confidenceScores) <- c("Function", names(mergedOutput2)[c(4:18)])
confidenceScores$Function <- functionlist
LE_scores <- confidenceScores
Pathway_scores <- confidenceScores

library(EnvStats)

for(i in 1:length(resultsList)){
  LEreference <- get(resultsList[i])
  
  for(j in 1:nrow(confidenceScores)){
    currFunc <- confidenceScores$Function[j]
    pathways <- as.character(unique(mergedOutput2$Pathway[mergedOutput2$Function==currFunc & (mergedOutput2[,i+3] != 0)]))
    pathways <- pathways[!is.na(pathways)]
    if(length(pathways) > 1){
      
      #scores for pathways gene lists
      gsWorkingPathways <- rep(list("A"), length(pathways))
      names(gsWorkingPathways) <- pathways
      for(k in 1:length(pathways)){
        gsWorkingPathways[[k]] <- gsub(" ", "", unlist(strsplit(mergedOutput2$Genes[mergedOutput2$Pathway == pathways[k]], split="|", fixed=TRUE)))
      }
      
      gsOverlapPathways <- newGOM(gsWorkingPathways, gsWorkingPathways)
      gsOutputPathways <- getMatrix(gsOverlapPathways, name="Jaccard")
      for(k in 1:nrow(gsOutputPathways)){
        gsOutputPathways[k,c(k:nrow(gsOutputPathways))] <- NA
      }
      genesetScore <- mean(gsOutputPathways, na.rm=TRUE)
      #scores for leading-edges
      gsWorkingLEs <- rep(list("A"), length(pathways))
      names(gsWorkingLEs) <- pathways
      for(k in 1:length(pathways)){
        gsWorkingLEs[[k]] <- unlist(LEreference$leadingEdge[LEreference$pathway == pathways[k]])
      }
      
      gsOverlapLEs <- newGOM(gsWorkingLEs, gsWorkingLEs)
      gsOutputLEs <- getMatrix(gsOverlapLEs, name="Jaccard")
      for(k in 1:nrow(gsOutputLEs)){
        gsOutputLEs[k,c(k:nrow(gsOutputLEs))] <- NA
      }
      LEScore <- mean(gsOutputLEs, na.rm=TRUE)
      
      #Add scores to holding DFs
      Pathway_scores[j,i+1] <- genesetScore
      LE_scores[j,i+1] <- LEScore
      confidenceScores[j,i+1] <- min(genesetScore, LEScore)
    }
  }
}

write.csv(confidenceScores, "GSEA parent-child redundancy-confidence scores.csv")
#Generate heatmaps and output tables (in same order) - for each cell-type
mergedOutputOrdered <- mergedOutput2[order(match(mergedOutput2$Function, functionlist)),c(2,1,3:22)]

#Cell-specific comparisons across disorders, integrate actual NES (from GSEAresults files) into output tables

#PYR23
PYR23_Table <- mergedOutputOrdered[,c(1,22,21,3,grep("PYR23", names(mergedOutputOrdered)))]
#Remove NA rows
PYR23_Table <- PYR23_Table[!(rowSums(is.na(PYR23_Table[,c(5:7)])) == 3),]

#Remove Pathway column, Order by schiz, BD, then MDD, lastly by function
PYR23_Table <- PYR23_Table[order(PYR23_Table[,7], decreasing = TRUE),]
PYR23_Table <- PYR23_Table[order(PYR23_Table[,6], decreasing = TRUE),]
PYR23_Table <- PYR23_Table[order(PYR23_Table[,5], decreasing = TRUE),]

PYR23_Table <- PYR23_Table[order(match(PYR23_Table$Function, functionlist)),]
PYR23_Table[is.na(PYR23_Table)] <- 0

write.csv(PYR23_Table, "PYR23_Output_Table.csv")
rownames(PYR23_Table) <- paste0(PYR23_Table$Name,"|",PYR23_Table$Function)

PYR23_func_hmap <- as.matrix(PYR23_Table[,-c(1:4)])

#Lazy and inefficient - but that's how it goes sometimes!
for(i in 1:nrow(PYR23_Table)){
  PYR23_Table$PYR23_MDD_ES[i] <- PYR23_MDD_gseaM1$NES[gsub("%.*%.*$", "",PYR23_MDD_gseaM1$pathway) == PYR23_Table$Name[i]]
  PYR23_Table$PYR23_Bipolar_ES[i] <- PYR23_Bipolar_gseaM1$NES[gsub("%.*%.*$", "",PYR23_Bipolar_gseaM1$pathway) == PYR23_Table$Name[i]]
  PYR23_Table$PYR23_SCHIZ_ES[i] <- PYR23_SCHIZ_gseaM1$NES[gsub("%.*%.*$", "",PYR23_SCHIZ_gseaM1$pathway) == PYR23_Table$Name[i]]
  
  if(PYR23_func_hmap[i,1] == 0){
    PYR23_Table$PYR23_MDD_ES[i] <- NA
  }
 
  if(PYR23_func_hmap[i,2] == 0){
    PYR23_Table$PYR23_Bipolar_ES[i] <- NA
  }
  
  if(PYR23_func_hmap[i,3] == 0){
    PYR23_Table$PYR23_SCHIZ_ES[i] <- NA
  }
}

PYR23_Table[is.na(PYR23_Table)] <- 0

write.csv(PYR23_Table, "PYR23 output table_NES.csv")


Heatmap(PYR23_func_hmap, 
        column_labels = c("MDD", "BPD", "SCZ"), 
        column_title = "PYR L2/3", col=colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, 
        cluster_rows = FALSE,
        cluster_row_slices = FALSE,
        row_split = factor(gsub("^.*\\|","",rownames(PYR23_func_hmap)), levels = unique(gsub("^.*\\|","",rownames(PYR23_func_hmap)))),
        show_row_names = FALSE, 
        width = unit(8, "cm"), 
        border = TRUE,
        column_names_rot = 0,
        column_names_side = "top")

#PYR56
PYR56_Table <- mergedOutputOrdered[,c(1,22,21,3,grep("PYR56", names(mergedOutputOrdered)))]
#Remove NA rows
PYR56_Table <- PYR56_Table[!(rowSums(is.na(PYR56_Table[,c(5:7)])) == 3),]

#Remove Pathway column, Order by schiz, BD, then MDD, lastly by function
PYR56_Table <- PYR56_Table[order(PYR56_Table[,7], decreasing = TRUE),]
PYR56_Table <- PYR56_Table[order(PYR56_Table[,6], decreasing = TRUE),]
PYR56_Table <- PYR56_Table[order(PYR56_Table[,5], decreasing = TRUE),]

PYR56_Table <- PYR56_Table[order(match(PYR56_Table$Function, functionlist)),]
PYR56_Table[is.na(PYR56_Table)] <- 0

write.csv(PYR56_Table, "PYR56_Output_Table.csv")
rownames(PYR56_Table) <- paste0(PYR56_Table$Name,"|",PYR56_Table$Function)

PYR56_func_hmap <- as.matrix(PYR56_Table[,-c(1:4)])
#Lazy and inefficient - but that's how it goes sometimes!
for(i in 1:nrow(PYR56_Table)){
  PYR56_Table$PYR56_MDD_ES[i] <- PYR56_MDD_gseaM1$NES[gsub("%.*%.*$", "",PYR56_MDD_gseaM1$pathway) == PYR56_Table$Name[i]]
  PYR56_Table$PYR56_Bipolar_ES[i] <- PYR56_Bipolar_gseaM1$NES[gsub("%.*%.*$", "",PYR56_Bipolar_gseaM1$pathway) == PYR56_Table$Name[i]]
  PYR56_Table$PYR56_SCHIZ_ES[i] <- PYR56_SCHIZ_gseaM1$NES[gsub("%.*%.*$", "",PYR56_SCHIZ_gseaM1$pathway) == PYR56_Table$Name[i]]
  
  if(PYR56_func_hmap[i,1] == 0){
    PYR56_Table$PYR56_MDD_ES[i] <- NA
  }
  
  if(PYR56_func_hmap[i,2] == 0){
    PYR56_Table$PYR56_Bipolar_ES[i] <- NA
  }
  
  if(PYR56_func_hmap[i,3] == 0){
    PYR56_Table$PYR56_SCHIZ_ES[i] <- NA
  }
}

PYR56_Table[is.na(PYR56_Table)] <- 0

write.csv(PYR56_Table, "PYR56 output table_NES.csv")


Heatmap(PYR56_func_hmap, 
        column_labels = c("MDD", "BPD", "SCZ"), 
        column_title = "PYR L5/6", col=colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, 
        cluster_rows = FALSE,
        cluster_row_slices = FALSE,
        row_split = factor(gsub("^.*\\|","",rownames(PYR56_func_hmap)), levels = unique(gsub("^.*\\|","",rownames(PYR56_func_hmap)))),
        show_row_names = FALSE, 
        width = unit(8, "cm"), 
        border = TRUE,
        column_names_rot = 0,
        column_names_side = "top")

#SST
SST_Table <- mergedOutputOrdered[,c(1,22,21,3,grep("SST", names(mergedOutputOrdered)))]
#Remove NA rows
SST_Table <- SST_Table[!(rowSums(is.na(SST_Table[,c(5:7)])) == 3),]

#Remove Pathway column, Order by schiz, BD, then MDD, lastly by function
SST_Table <- SST_Table[order(SST_Table[,7], decreasing = TRUE),]
SST_Table <- SST_Table[order(SST_Table[,6], decreasing = TRUE),]
SST_Table <- SST_Table[order(SST_Table[,5], decreasing = TRUE),]

SST_Table <- SST_Table[order(match(SST_Table$Function, functionlist)),]
SST_Table[is.na(SST_Table)] <- 0

write.csv(SST_Table, "SST_Output_Table.csv")
rownames(SST_Table) <- paste0(SST_Table$Name,"|",SST_Table$Function)

SST_func_hmap <- as.matrix(SST_Table[,-c(1:4)])
#Lazy and inefficient - but that's how it goes sometimes!
for(i in 1:nrow(SST_Table)){
  SST_Table$SST_MDD_ES[i] <- SST_MDD_gseaM1$NES[gsub("%.*%.*$", "",SST_MDD_gseaM1$pathway) == SST_Table$Name[i]]
  SST_Table$SST_Bipolar_ES[i] <- SST_Bipolar_gseaM1$NES[gsub("%.*%.*$", "",SST_Bipolar_gseaM1$pathway) == SST_Table$Name[i]]
  SST_Table$SST_SCHIZ_ES[i] <- SST_SCHIZ_gseaM1$NES[gsub("%.*%.*$", "",SST_SCHIZ_gseaM1$pathway) == SST_Table$Name[i]]
  
  if(SST_func_hmap[i,1] == 0){
    SST_Table$SST_MDD_ES[i] <- NA
  }
  
  if(SST_func_hmap[i,2] == 0){
    SST_Table$SST_Bipolar_ES[i] <- NA
  }
  
  if(SST_func_hmap[i,3] == 0){
    SST_Table$SST_SCHIZ_ES[i] <- NA
  }
}

SST_Table[is.na(SST_Table)] <- 0

write.csv(SST_Table, "SST output table_NES.csv")


Heatmap(SST_func_hmap, 
        column_labels = c("MDD", "BPD", "SCZ"), 
        column_title = "SST", col=colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, 
        cluster_rows = FALSE,
        cluster_row_slices = FALSE,
        row_split = factor(gsub("^.*\\|","",rownames(SST_func_hmap)), levels = unique(gsub("^.*\\|","",rownames(SST_func_hmap)))),
        show_row_names = FALSE, 
        width = unit(8, "cm"), 
        border = TRUE,
        column_names_rot = 0,
        column_names_side = "top")

#PV
PV_Table <- mergedOutputOrdered[,c(1,22,21,3,grep("PV", names(mergedOutputOrdered)))]
#Remove NA rows
PV_Table <- PV_Table[!(rowSums(is.na(PV_Table[,c(5:7)])) == 3),]

#Remove Pathway column, Order by schiz, BD, then MDD, lastly by function
PV_Table <- PV_Table[order(PV_Table[,7], decreasing = TRUE),]
PV_Table <- PV_Table[order(PV_Table[,6], decreasing = TRUE),]
PV_Table <- PV_Table[order(PV_Table[,5], decreasing = TRUE),]

PV_Table <- PV_Table[order(match(PV_Table$Function, functionlist)),]
PV_Table[is.na(PV_Table)] <- 0

write.csv(PV_Table, "PV_Output_Table.csv")
rownames(PV_Table) <- paste0(PV_Table$Name,"|",PV_Table$Function)

PV_func_hmap <- as.matrix(PV_Table[,-c(1:4)])
#Lazy and inefficient - but that's how it goes sometimes!
for(i in 1:nrow(PV_Table)){
  PV_Table$PV_MDD_ES[i] <- PV_MDD_gseaM1$NES[gsub("%.*%.*$", "",PV_MDD_gseaM1$pathway) == PV_Table$Name[i]]
  PV_Table$PV_Bipolar_ES[i] <- PV_Bipolar_gseaM1$NES[gsub("%.*%.*$", "",PV_Bipolar_gseaM1$pathway) == PV_Table$Name[i]]
  PV_Table$PV_SCHIZ_ES[i] <- PV_SCHIZ_gseaM1$NES[gsub("%.*%.*$", "",PV_SCHIZ_gseaM1$pathway) == PV_Table$Name[i]]
  
  if(PV_func_hmap[i,1] == 0){
    PV_Table$PV_MDD_ES[i] <- NA
  }
  
  if(PV_func_hmap[i,2] == 0){
    PV_Table$PV_Bipolar_ES[i] <- NA
  }
  
  if(PV_func_hmap[i,3] == 0){
    PV_Table$PV_SCHIZ_ES[i] <- NA
  }
}

PV_Table[is.na(PV_Table)] <- 0

write.csv(PV_Table, "PV output table_NES.csv")

Heatmap(PV_func_hmap, 
        column_labels = c("MDD", "BPD", "SCZ"), 
        column_title = "PV", col=colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, 
        cluster_rows = FALSE,
        cluster_row_slices = FALSE,
        row_split = factor(gsub("^.*\\|","",rownames(PV_func_hmap)), levels = unique(gsub("^.*\\|","",rownames(PV_func_hmap)))),
        show_row_names = FALSE, 
        width = unit(8, "cm"), 
        border = TRUE,
        column_names_rot = 0,
        column_names_side = "top")

#VIP
VIP_Table <- mergedOutputOrdered[,c(1,22,21,3,grep("VIP", names(mergedOutputOrdered)))]
#Remove NA rows
VIP_Table <- VIP_Table[!(rowSums(is.na(VIP_Table[,c(5:7)])) == 3),]

#Remove Pathway column, Order by schiz, BD, then MDD, lastly by function
VIP_Table <- VIP_Table[order(VIP_Table[,7], decreasing = TRUE),]
VIP_Table <- VIP_Table[order(VIP_Table[,6], decreasing = TRUE),]
VIP_Table <- VIP_Table[order(VIP_Table[,5], decreasing = TRUE),]

VIP_Table <- VIP_Table[order(match(VIP_Table$Function, functionlist)),]
VIP_Table[is.na(VIP_Table)] <- 0

write.csv(VIP_Table, "VIP_Output_Table.csv")
rownames(VIP_Table) <- paste0(VIP_Table$Name,"|",VIP_Table$Function)

VIP_func_hmap <- as.matrix(VIP_Table[,-c(1:4)])
#Lazy and inefficient - but that's how it goes sometimes!
for(i in 1:nrow(VIP_Table)){
  VIP_Table$VIP_MDD_ES[i] <- VIP_MDD_gseaM1$NES[gsub("%.*%.*$", "",VIP_MDD_gseaM1$pathway) == VIP_Table$Name[i]]
  VIP_Table$VIP_Bipolar_ES[i] <- VIP_Bipolar_gseaM1$NES[gsub("%.*%.*$", "",VIP_Bipolar_gseaM1$pathway) == VIP_Table$Name[i]]
  VIP_Table$VIP_SCHIZ_ES[i] <- VIP_SCHIZ_gseaM1$NES[gsub("%.*%.*$", "",VIP_SCHIZ_gseaM1$pathway) == VIP_Table$Name[i]]
  
  if(VIP_func_hmap[i,1] == 0){
    VIP_Table$VIP_MDD_ES[i] <- NA
  }
  
  if(VIP_func_hmap[i,2] == 0){
    VIP_Table$VIP_Bipolar_ES[i] <- NA
  }
  
  if(VIP_func_hmap[i,3] == 0){
    VIP_Table$VIP_SCHIZ_ES[i] <- NA
  }
}

VIP_Table[is.na(VIP_Table)] <- 0

write.csv(VIP_Table, "VIP output table_NES.csv")


Heatmap(VIP_func_hmap, 
        column_labels = c("MDD", "BPD", "SCZ"), 
        column_title = "VIP", col=colorRamp2(c(-1, 0, 1), c("dodgerblue4", "white", "firebrick3")),
        cluster_columns = FALSE, 
        cluster_rows = FALSE,
        cluster_row_slices = FALSE,
        row_split = factor(gsub("^.*\\|","",rownames(VIP_func_hmap)), levels = unique(gsub("^.*\\|","",rownames(VIP_func_hmap)))),
        show_row_names = FALSE, 
        width = unit(8, "cm"), 
        border = TRUE,
        column_names_rot = 0,
        column_names_side = "top")
