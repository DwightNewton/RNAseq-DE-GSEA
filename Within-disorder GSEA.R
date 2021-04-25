library(fgsea)
library(DESeq2)
library(BiocParallel)
options(stringsAsFactors = FALSE)
options(scipen = 999)
#Run fgsea using signed -log10(pvalue)
#Will also use new "Multilevel" approach in fgsea - more accurate FDR correction and estimation of very low pvalues

setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/")

#Load required data and DESeqDataSet (dds) files to pull gene names (exported files only had gene stable IDs)
load("FDRtooled_results.rData")
filelist <- list.files(pattern="dds")
for(i in filelist){load(i)}

setwd("GSEA/")
RefGeneLists <- gmtPathways("Human_GO_AllPathways_with_GO_iea_April_01_2020_symbol.gmt")

register(MulticoreParam(workers = 10))

PVgenes <- as.data.frame(rowData(dds_PV_exons_groupCT_age_sex_group)[,c(1:8)])
PYR23genes <- as.data.frame(rowData(dds_PYR23_exons_groupCT_age_sex_group)[,c(1:8)])
PYR56genes <- as.data.frame(rowData(dds_PYR56_exons_groupCT_age_sex_group)[,c(1:8)])
SSTgenes <- as.data.frame(rowData(dds_SST_exons_groupCT_age_sex_group)[,c(1:8)])
VIPgenes <- as.data.frame(rowData(dds_VIP_exons_groupCT_age_sex_group)[,c(1:8)])

#MDD
#PV
fdrT_CT_PV_MDD <- as.data.frame(fdrT_CT_PV_MDD)
fdrT_CT_PV_MDD$rank <- -1 * log10(fdrT_CT_PV_MDD$pvalue) * sign(fdrT_CT_PV_MDD$log2FoldChange)
fdrT_CT_PV_MDD <- fdrT_CT_PV_MDD[order(fdrT_CT_PV_MDD$rank, decreasing = TRUE),]
fdrT_CT_PV_MDD$Gene.stable.ID <- rownames(fdrT_CT_PV_MDD)
fdrT_CT_PV_MDD$Gene.name <- PVgenes$Gene.name[match(fdrT_CT_PV_MDD$Gene.stable.ID, PVgenes$Gene.stable.ID)]
fdrT_CT_PV_MDD <- fdrT_CT_PV_MDD[!duplicated(fdrT_CT_PV_MDD$Gene.name),]
PV_MDD_ranklist <- fdrT_CT_PV_MDD$rank
names(PV_MDD_ranklist) <- fdrT_CT_PV_MDD$Gene.name

PV_MDD_gsea <- fgsea(pathways = RefGeneLists,
                     stats = PV_MDD_ranklist, 
                     minSize = 15,
                     maxSize = 500, 
                     nperm = 1000000, 
                     nproc = 10)

#PYR23
fdrT_CT_PYR23_MDD <- as.data.frame(fdrT_CT_PYR23_MDD)
fdrT_CT_PYR23_MDD$rank <- -1 * log10(fdrT_CT_PYR23_MDD$pvalue) * sign(fdrT_CT_PYR23_MDD$log2FoldChange)
fdrT_CT_PYR23_MDD <- fdrT_CT_PYR23_MDD[order(fdrT_CT_PYR23_MDD$rank, decreasing = TRUE),]
fdrT_CT_PYR23_MDD$Gene.stable.ID <- rownames(fdrT_CT_PYR23_MDD)
fdrT_CT_PYR23_MDD$Gene.name <- PYR23genes$Gene.name[match(fdrT_CT_PYR23_MDD$Gene.stable.ID, PYR23genes$Gene.stable.ID)]
fdrT_CT_PYR23_MDD <- fdrT_CT_PYR23_MDD[!duplicated(fdrT_CT_PYR23_MDD$Gene.name),]
PYR23_MDD_ranklist <- fdrT_CT_PYR23_MDD$rank
names(PYR23_MDD_ranklist) <- fdrT_CT_PYR23_MDD$Gene.name

PYR23_MDD_gsea <- fgsea(pathways = RefGeneLists,
                     stats = PYR23_MDD_ranklist, 
                     minSize = 15,
                     maxSize = 500, 
                     nperm = 1000000,
                     nproc = 10)

#PYR56
fdrT_CT_PYR56_MDD <- as.data.frame(fdrT_CT_PYR56_MDD)
fdrT_CT_PYR56_MDD$rank <- -1 * log10(fdrT_CT_PYR56_MDD$pvalue) * sign(fdrT_CT_PYR56_MDD$log2FoldChange)
fdrT_CT_PYR56_MDD <- fdrT_CT_PYR56_MDD[order(fdrT_CT_PYR56_MDD$rank, decreasing = TRUE),]
fdrT_CT_PYR56_MDD$Gene.stable.ID <- rownames(fdrT_CT_PYR56_MDD)
fdrT_CT_PYR56_MDD$Gene.name <- PYR56genes$Gene.name[match(fdrT_CT_PYR56_MDD$Gene.stable.ID, PYR56genes$Gene.stable.ID)]
fdrT_CT_PYR56_MDD <- fdrT_CT_PYR56_MDD[!duplicated(fdrT_CT_PYR56_MDD$Gene.name),]
PYR56_MDD_ranklist <- fdrT_CT_PYR56_MDD$rank
names(PYR56_MDD_ranklist) <- fdrT_CT_PYR56_MDD$Gene.name

PYR56_MDD_gsea <- fgsea(pathways = RefGeneLists,
                     stats = PYR56_MDD_ranklist, 
                     minSize = 15,
                     maxSize = 500, 
                     nperm = 1000000,
                     nproc = 10)

#SST
fdrT_CT_SST_MDD <- as.data.frame(fdrT_CT_SST_MDD)
fdrT_CT_SST_MDD$rank <- -1 * log10(fdrT_CT_SST_MDD$pvalue) * sign(fdrT_CT_SST_MDD$log2FoldChange)
fdrT_CT_SST_MDD <- fdrT_CT_SST_MDD[order(fdrT_CT_SST_MDD$rank, decreasing = TRUE),]
fdrT_CT_SST_MDD$Gene.stable.ID <- rownames(fdrT_CT_SST_MDD)
fdrT_CT_SST_MDD$Gene.name <- SSTgenes$Gene.name[match(fdrT_CT_SST_MDD$Gene.stable.ID, SSTgenes$Gene.stable.ID)]
fdrT_CT_SST_MDD <- fdrT_CT_SST_MDD[!duplicated(fdrT_CT_SST_MDD$Gene.name),]
SST_MDD_ranklist <- fdrT_CT_SST_MDD$rank
names(SST_MDD_ranklist) <- fdrT_CT_SST_MDD$Gene.name

SST_MDD_gsea <- fgsea(pathways = RefGeneLists,
                     stats = SST_MDD_ranklist, 
                     minSize = 15,
                     maxSize = 500, 
                     nperm = 1000000,
                     nproc = 10)

#VIP
fdrT_CT_VIP_MDD <- as.data.frame(fdrT_CT_VIP_MDD)
fdrT_CT_VIP_MDD$rank <- -1 * log10(fdrT_CT_VIP_MDD$pvalue) * sign(fdrT_CT_VIP_MDD$log2FoldChange)
fdrT_CT_VIP_MDD <- fdrT_CT_VIP_MDD[order(fdrT_CT_VIP_MDD$rank, decreasing = TRUE),]
fdrT_CT_VIP_MDD$Gene.stable.ID <- rownames(fdrT_CT_VIP_MDD)
fdrT_CT_VIP_MDD$Gene.name <- VIPgenes$Gene.name[match(fdrT_CT_VIP_MDD$Gene.stable.ID, VIPgenes$Gene.stable.ID)]
fdrT_CT_VIP_MDD <- fdrT_CT_VIP_MDD[!duplicated(fdrT_CT_VIP_MDD$Gene.name),]
VIP_MDD_ranklist <- fdrT_CT_VIP_MDD$rank
names(VIP_MDD_ranklist) <- fdrT_CT_VIP_MDD$Gene.name

VIP_MDD_gsea <- fgsea(pathways = RefGeneLists,
                     stats = VIP_MDD_ranklist, 
                     minSize = 15,
                     maxSize = 500, 
                     nperm = 1000000,
                     nproc = 10)

#Bipolar
#PV
fdrT_CT_PV_Bipolar <- as.data.frame(fdrT_CT_PV_Bipolar)
fdrT_CT_PV_Bipolar$rank <- -1 * log10(fdrT_CT_PV_Bipolar$pvalue) * sign(fdrT_CT_PV_Bipolar$log2FoldChange)
fdrT_CT_PV_Bipolar <- fdrT_CT_PV_Bipolar[order(fdrT_CT_PV_Bipolar$rank, decreasing = TRUE),]
fdrT_CT_PV_Bipolar$Gene.stable.ID <- rownames(fdrT_CT_PV_Bipolar)
fdrT_CT_PV_Bipolar$Gene.name <- PVgenes$Gene.name[match(fdrT_CT_PV_Bipolar$Gene.stable.ID, PVgenes$Gene.stable.ID)]
fdrT_CT_PV_Bipolar <- fdrT_CT_PV_Bipolar[!duplicated(fdrT_CT_PV_Bipolar$Gene.name),]
PV_Bipolar_ranklist <- fdrT_CT_PV_Bipolar$rank
names(PV_Bipolar_ranklist) <- fdrT_CT_PV_Bipolar$Gene.name

PV_Bipolar_gsea <- fgsea(pathways = RefGeneLists,
                     stats = PV_Bipolar_ranklist, 
                     minSize = 15,
                     maxSize = 500, 
                     nperm = 1000000,
                     nproc = 10)

#PYR23
fdrT_CT_PYR23_Bipolar <- as.data.frame(fdrT_CT_PYR23_Bipolar)
fdrT_CT_PYR23_Bipolar$rank <- -1 * log10(fdrT_CT_PYR23_Bipolar$pvalue) * sign(fdrT_CT_PYR23_Bipolar$log2FoldChange)
fdrT_CT_PYR23_Bipolar <- fdrT_CT_PYR23_Bipolar[order(fdrT_CT_PYR23_Bipolar$rank, decreasing = TRUE),]
fdrT_CT_PYR23_Bipolar$Gene.stable.ID <- rownames(fdrT_CT_PYR23_Bipolar)
fdrT_CT_PYR23_Bipolar$Gene.name <- PYR23genes$Gene.name[match(fdrT_CT_PYR23_Bipolar$Gene.stable.ID, PYR23genes$Gene.stable.ID)]
fdrT_CT_PYR23_Bipolar <- fdrT_CT_PYR23_Bipolar[!duplicated(fdrT_CT_PYR23_Bipolar$Gene.name),]
PYR23_Bipolar_ranklist <- fdrT_CT_PYR23_Bipolar$rank
names(PYR23_Bipolar_ranklist) <- fdrT_CT_PYR23_Bipolar$Gene.name

PYR23_Bipolar_gsea <- fgsea(pathways = RefGeneLists,
                        stats = PYR23_Bipolar_ranklist, 
                        minSize = 15,
                        maxSize = 500, 
                        nperm = 1000000,
                        nproc = 10)

#PYR56
fdrT_CT_PYR56_Bipolar <- as.data.frame(fdrT_CT_PYR56_Bipolar)
fdrT_CT_PYR56_Bipolar$rank <- -1 * log10(fdrT_CT_PYR56_Bipolar$pvalue) * sign(fdrT_CT_PYR56_Bipolar$log2FoldChange)
fdrT_CT_PYR56_Bipolar <- fdrT_CT_PYR56_Bipolar[order(fdrT_CT_PYR56_Bipolar$rank, decreasing = TRUE),]
fdrT_CT_PYR56_Bipolar$Gene.stable.ID <- rownames(fdrT_CT_PYR56_Bipolar)
fdrT_CT_PYR56_Bipolar$Gene.name <- PYR56genes$Gene.name[match(fdrT_CT_PYR56_Bipolar$Gene.stable.ID, PYR56genes$Gene.stable.ID)]
fdrT_CT_PYR56_Bipolar <- fdrT_CT_PYR56_Bipolar[!duplicated(fdrT_CT_PYR56_Bipolar$Gene.name),]
PYR56_Bipolar_ranklist <- fdrT_CT_PYR56_Bipolar$rank
names(PYR56_Bipolar_ranklist) <- fdrT_CT_PYR56_Bipolar$Gene.name

PYR56_Bipolar_gsea <- fgsea(pathways = RefGeneLists,
                        stats = PYR56_Bipolar_ranklist, 
                        minSize = 15,
                        maxSize = 500, 
                        nperm = 1000000,
                        nproc = 10)

#SST
fdrT_CT_SST_Bipolar <- as.data.frame(fdrT_CT_SST_Bipolar)
fdrT_CT_SST_Bipolar$rank <- -1 * log10(fdrT_CT_SST_Bipolar$pvalue) * sign(fdrT_CT_SST_Bipolar$log2FoldChange)
fdrT_CT_SST_Bipolar <- fdrT_CT_SST_Bipolar[order(fdrT_CT_SST_Bipolar$rank, decreasing = TRUE),]
fdrT_CT_SST_Bipolar$Gene.stable.ID <- rownames(fdrT_CT_SST_Bipolar)
fdrT_CT_SST_Bipolar$Gene.name <- SSTgenes$Gene.name[match(fdrT_CT_SST_Bipolar$Gene.stable.ID, SSTgenes$Gene.stable.ID)]
fdrT_CT_SST_Bipolar <- fdrT_CT_SST_Bipolar[!duplicated(fdrT_CT_SST_Bipolar$Gene.name),]
SST_Bipolar_ranklist <- fdrT_CT_SST_Bipolar$rank
names(SST_Bipolar_ranklist) <- fdrT_CT_SST_Bipolar$Gene.name

SST_Bipolar_gsea <- fgsea(pathways = RefGeneLists,
                      stats = SST_Bipolar_ranklist, 
                      minSize = 15,
                      maxSize = 500, 
                      nperm = 1000000,
                      nproc = 10)

#VIP
fdrT_CT_VIP_Bipolar <- as.data.frame(fdrT_CT_VIP_Bipolar)
fdrT_CT_VIP_Bipolar$rank <- -1 * log10(fdrT_CT_VIP_Bipolar$pvalue) * sign(fdrT_CT_VIP_Bipolar$log2FoldChange)
fdrT_CT_VIP_Bipolar <- fdrT_CT_VIP_Bipolar[order(fdrT_CT_VIP_Bipolar$rank, decreasing = TRUE),]
fdrT_CT_VIP_Bipolar$Gene.stable.ID <- rownames(fdrT_CT_VIP_Bipolar)
fdrT_CT_VIP_Bipolar$Gene.name <- VIPgenes$Gene.name[match(fdrT_CT_VIP_Bipolar$Gene.stable.ID, VIPgenes$Gene.stable.ID)]
fdrT_CT_VIP_Bipolar <- fdrT_CT_VIP_Bipolar[!duplicated(fdrT_CT_VIP_Bipolar$Gene.name),]
VIP_Bipolar_ranklist <- fdrT_CT_VIP_Bipolar$rank
names(VIP_Bipolar_ranklist) <- fdrT_CT_VIP_Bipolar$Gene.name

VIP_Bipolar_gsea <- fgsea(pathways = RefGeneLists,
                      stats = VIP_Bipolar_ranklist, 
                      minSize = 15,
                      maxSize = 500, 
                      nperm = 1000000,
                      nproc = 10)

#SCHIZ
#PV
fdrT_CT_PV_SCHIZ <- as.data.frame(fdrT_CT_PV_SCHIZ)
fdrT_CT_PV_SCHIZ$rank <- -1 * log10(fdrT_CT_PV_SCHIZ$pvalue) * sign(fdrT_CT_PV_SCHIZ$log2FoldChange)
fdrT_CT_PV_SCHIZ <- fdrT_CT_PV_SCHIZ[order(fdrT_CT_PV_SCHIZ$rank, decreasing = TRUE),]
fdrT_CT_PV_SCHIZ$Gene.stable.ID <- rownames(fdrT_CT_PV_SCHIZ)
fdrT_CT_PV_SCHIZ$Gene.name <- PVgenes$Gene.name[match(fdrT_CT_PV_SCHIZ$Gene.stable.ID, PVgenes$Gene.stable.ID)]
fdrT_CT_PV_SCHIZ <- fdrT_CT_PV_SCHIZ[!duplicated(fdrT_CT_PV_SCHIZ$Gene.name),]
PV_SCHIZ_ranklist <- fdrT_CT_PV_SCHIZ$rank
names(PV_SCHIZ_ranklist) <- fdrT_CT_PV_SCHIZ$Gene.name

PV_SCHIZ_gsea <- fgsea(pathways = RefGeneLists,
                     stats = PV_SCHIZ_ranklist, 
                     minSize = 15,
                     maxSize = 500, 
                     nperm = 1000000,
                     nproc = 10)

#PYR23
fdrT_CT_PYR23_SCHIZ <- as.data.frame(fdrT_CT_PYR23_SCHIZ)
fdrT_CT_PYR23_SCHIZ$rank <- -1 * log10(fdrT_CT_PYR23_SCHIZ$pvalue) * sign(fdrT_CT_PYR23_SCHIZ$log2FoldChange)
fdrT_CT_PYR23_SCHIZ <- fdrT_CT_PYR23_SCHIZ[order(fdrT_CT_PYR23_SCHIZ$rank, decreasing = TRUE),]
fdrT_CT_PYR23_SCHIZ$Gene.stable.ID <- rownames(fdrT_CT_PYR23_SCHIZ)
fdrT_CT_PYR23_SCHIZ$Gene.name <- PYR23genes$Gene.name[match(fdrT_CT_PYR23_SCHIZ$Gene.stable.ID, PYR23genes$Gene.stable.ID)]
fdrT_CT_PYR23_SCHIZ <- fdrT_CT_PYR23_SCHIZ[!duplicated(fdrT_CT_PYR23_SCHIZ$Gene.name),]
PYR23_SCHIZ_ranklist <- fdrT_CT_PYR23_SCHIZ$rank
names(PYR23_SCHIZ_ranklist) <- fdrT_CT_PYR23_SCHIZ$Gene.name

PYR23_SCHIZ_gsea <- fgsea(pathways = RefGeneLists,
                        stats = PYR23_SCHIZ_ranklist, 
                        minSize = 15,
                        maxSize = 500, 
                        nperm = 1000000,
                        nproc = 10)

#PYR56
fdrT_CT_PYR56_SCHIZ <- as.data.frame(fdrT_CT_PYR56_SCHIZ)
fdrT_CT_PYR56_SCHIZ$rank <- -1 * log10(fdrT_CT_PYR56_SCHIZ$pvalue) * sign(fdrT_CT_PYR56_SCHIZ$log2FoldChange)
fdrT_CT_PYR56_SCHIZ <- fdrT_CT_PYR56_SCHIZ[order(fdrT_CT_PYR56_SCHIZ$rank, decreasing = TRUE),]
fdrT_CT_PYR56_SCHIZ$Gene.stable.ID <- rownames(fdrT_CT_PYR56_SCHIZ)
fdrT_CT_PYR56_SCHIZ$Gene.name <- PYR56genes$Gene.name[match(fdrT_CT_PYR56_SCHIZ$Gene.stable.ID, PYR56genes$Gene.stable.ID)]
fdrT_CT_PYR56_SCHIZ <- fdrT_CT_PYR56_SCHIZ[!duplicated(fdrT_CT_PYR56_SCHIZ$Gene.name),]
PYR56_SCHIZ_ranklist <- fdrT_CT_PYR56_SCHIZ$rank
names(PYR56_SCHIZ_ranklist) <- fdrT_CT_PYR56_SCHIZ$Gene.name

PYR56_SCHIZ_gsea <- fgsea(pathways = RefGeneLists,
                        stats = PYR56_SCHIZ_ranklist, 
                        minSize = 15,
                        maxSize = 500, 
                        nperm = 1000000,
                        nproc = 10)

#SST
fdrT_CT_SST_SCHIZ <- as.data.frame(fdrT_CT_SST_SCHIZ)
fdrT_CT_SST_SCHIZ$rank <- -1 * log10(fdrT_CT_SST_SCHIZ$pvalue) * sign(fdrT_CT_SST_SCHIZ$log2FoldChange)
fdrT_CT_SST_SCHIZ <- fdrT_CT_SST_SCHIZ[order(fdrT_CT_SST_SCHIZ$rank, decreasing = TRUE),]
fdrT_CT_SST_SCHIZ$Gene.stable.ID <- rownames(fdrT_CT_SST_SCHIZ)
fdrT_CT_SST_SCHIZ$Gene.name <- SSTgenes$Gene.name[match(fdrT_CT_SST_SCHIZ$Gene.stable.ID, SSTgenes$Gene.stable.ID)]
fdrT_CT_SST_SCHIZ <- fdrT_CT_SST_SCHIZ[!duplicated(fdrT_CT_SST_SCHIZ$Gene.name),]
SST_SCHIZ_ranklist <- fdrT_CT_SST_SCHIZ$rank
names(SST_SCHIZ_ranklist) <- fdrT_CT_SST_SCHIZ$Gene.name

SST_SCHIZ_gsea <- fgsea(pathways = RefGeneLists,
                      stats = SST_SCHIZ_ranklist, 
                      minSize = 15,
                      maxSize = 500, 
                      nperm = 1000000,
                      nproc = 10)

#VIP
fdrT_CT_VIP_SCHIZ <- as.data.frame(fdrT_CT_VIP_SCHIZ)
fdrT_CT_VIP_SCHIZ$rank <- -1 * log10(fdrT_CT_VIP_SCHIZ$pvalue) * sign(fdrT_CT_VIP_SCHIZ$log2FoldChange)
fdrT_CT_VIP_SCHIZ <- fdrT_CT_VIP_SCHIZ[order(fdrT_CT_VIP_SCHIZ$rank, decreasing = TRUE),]
fdrT_CT_VIP_SCHIZ$Gene.stable.ID <- rownames(fdrT_CT_VIP_SCHIZ)
fdrT_CT_VIP_SCHIZ$Gene.name <- VIPgenes$Gene.name[match(fdrT_CT_VIP_SCHIZ$Gene.stable.ID, VIPgenes$Gene.stable.ID)]
fdrT_CT_VIP_SCHIZ <- fdrT_CT_VIP_SCHIZ[!duplicated(fdrT_CT_VIP_SCHIZ$Gene.name),]
VIP_SCHIZ_ranklist <- fdrT_CT_VIP_SCHIZ$rank
names(VIP_SCHIZ_ranklist) <- fdrT_CT_VIP_SCHIZ$Gene.name

VIP_SCHIZ_gsea <- fgsea(pathways = RefGeneLists,
                      stats = VIP_SCHIZ_ranklist, 
                      minSize = 15,
                      maxSize = 500, 
                      nperm = 1000000,
                      nproc = 10)

save(PV_MDD_gsea, PYR23_MDD_gsea, PYR56_MDD_gsea, SST_MDD_gsea, VIP_MDD_gsea, PV_Bipolar_gsea, PYR23_Bipolar_gsea, PYR56_Bipolar_gsea, SST_Bipolar_gsea, VIP_Bipolar_gsea, PV_SCHIZ_gsea, PYR23_SCHIZ_gsea, PYR56_SCHIZ_gsea, SST_SCHIZ_gsea, VIP_SCHIZ_gsea, file="GSEA results (1M permutation).rData")

save(PV_MDD_ranklist, PYR23_MDD_ranklist, PYR56_MDD_ranklist, SST_MDD_ranklist, VIP_MDD_ranklist, PV_Bipolar_ranklist, PYR23_Bipolar_ranklist, PYR56_Bipolar_ranklist, SST_Bipolar_ranklist, VIP_Bipolar_ranklist, PV_SCHIZ_ranklist, PYR23_SCHIZ_ranklist, PYR56_SCHIZ_ranklist, SST_SCHIZ_ranklist, VIP_SCHIZ_ranklist, file="GSEA rankedlists (1M permutation).rData")
