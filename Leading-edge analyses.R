#GSEA Heatmap and Enrichment Map leading-edge analysis
library(fgsea)
library(DESeq2)
library(ggplot2)
options(stringsAsFactors = FALSE)

substrLeft <- function (x, n=13){
  substr(x, 1, nchar(x)-n)}

#ParentChild/heatmap results first

#Require: GSEA results and leading edges (both in fgsea results), GSEA results assigned to clusters (from output tables), and DE stats (from "PITT_tetrad_cohort_SCT_contrasts.csv")
#GSEA results
setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output")

load("GSEA rankedlists (GSEAparam1).rData")
load("GSEA results (GSEAparam1).rData")
load("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/FDRtooled_results.rData")

PVheatmap <- read.csv("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/Parent child/Text and ParentChild output/Heatmaps/PV output table_NES.csv")
PYR23heatmap <- read.csv("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/Parent child/Text and ParentChild output/Heatmaps/PYR23 output table_NES.csv")
PYR56heatmap <- read.csv("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/Parent child/Text and ParentChild output/Heatmaps/PYR56 output table_NES.csv")
SSTheatmap <- read.csv("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/Parent child/Text and ParentChild output/Heatmaps/SST output table_NES.csv")
VIPheatmap <- read.csv("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/Parent child/Text and ParentChild output/Heatmaps/VIP output table_NES.csv")

#Genesets for each contrast
#PV
PVMDDhmap <- PVheatmap[,c(1:5,6)]
PVBipolarhmap <- PVheatmap[,c(1:5,7)]
PVSCHIZhmap <- PVheatmap[,c(1:5,8)]

#PYR23
PYR23MDDhmap <- PYR23heatmap[,c(1:5,6)]
PYR23Bipolarhmap <- PYR23heatmap[,c(1:5,7)]
PYR23SCHIZhmap <- PYR23heatmap[,c(1:5,8)]

#PYR56
PYR56MDDhmap <- PYR56heatmap[,c(1:5,6)]
PYR56Bipolarhmap <- PYR56heatmap[,c(1:5,7)]
PYR56SCHIZhmap <- PYR56heatmap[,c(1:5,8)]

#SST
SSTMDDhmap <- SSTheatmap[,c(1:5,6)]
SSTBipolarhmap <- SSTheatmap[,c(1:5,7)]
SSTSCHIZhmap <- SSTheatmap[,c(1:5,8)]

#VIP
VIPMDDhmap <- VIPheatmap[,c(1:5,6)]
VIPBipolarhmap <- VIPheatmap[,c(1:5,7)]
VIPSCHIZhmap <- VIPheatmap[,c(1:5,8)]

DEGlist <- c("fdrT_CT_PV_MDD", "fdrT_CT_PV_Bipolar", "fdrT_CT_PV_SCHIZ", "fdrT_CT_PYR23_MDD", "fdrT_CT_PYR23_Bipolar", "fdrT_CT_PYR23_SCHIZ", "fdrT_CT_PYR56_MDD", "fdrT_CT_PYR56_Bipolar", "fdrT_CT_PYR56_SCHIZ", "fdrT_CT_SST_MDD", "fdrT_CT_SST_Bipolar", "fdrT_CT_SST_SCHIZ", "fdrT_CT_VIP_MDD", "fdrT_CT_VIP_Bipolar", "fdrT_CT_VIP_SCHIZ")
GSEAreslist <- c("PV_MDD_gseaM1", "PV_Bipolar_gseaM1", "PV_SCHIZ_gseaM1", "PYR23_MDD_gseaM1", "PYR23_Bipolar_gseaM1", "PYR23_SCHIZ_gseaM1", "PYR56_MDD_gseaM1", "PYR56_Bipolar_gseaM1", "PYR56_SCHIZ_gseaM1", "SST_MDD_gseaM1", "SST_Bipolar_gseaM1", "SST_SCHIZ_gseaM1", "VIP_MDD_gseaM1", "VIP_Bipolar_gseaM1", "VIP_SCHIZ_gseaM1")
heatmaplist <- c("PVMDDhmap", "PVBipolarhmap", "PVSCHIZhmap", "PYR23MDDhmap", "PYR23Bipolarhmap", "PYR23SCHIZhmap", "PYR56MDDhmap", "PYR56Bipolarhmap", "PYR56SCHIZhmap", "SSTMDDhmap", "SSTBipolarhmap", "SSTSCHIZhmap", "VIPMDDhmap", "VIPBipolarhmap", "VIPSCHIZhmap")

workingResults <- as.data.frame(matrix(nrow=1, ncol=107))
names(workingResults) <- c("Cell_Type","Contrast", "Function", "Direction","Mean_NES", "Enriched_Terms", "Union", paste0("Gene_", 1:100))

for (i in 1:length(DEGlist)){
  contrast <- gsub("_.*$", "", gsub("^.*CT_", "", DEGlist[i]))
  celltype <- gsub("^.*_", "", gsub("^.*CT_", "", DEGlist[i]))
  workinghmap <- get(heatmaplist[i])
  workingfunctions <- unique(workinghmap$Function)
  
  workingGSEA <- get(GSEAreslist[i])
  
  for(j in 1:length(workingfunctions)){
    
    hmapterm <- workingfunctions[j]
    directions <- c("Up", "Down")
    
    for(k in 1:length(directions)){
      currentDirection <- directions[k]
      
      if(currentDirection == "Up"){termxDirectionMembers <- workinghmap$Name[(workinghmap$Function == hmapterm) & (workinghmap[,6] > 0)]}
      if(currentDirection == "Up"){directGSEA <- workingGSEA[workingGSEA$NES > 0]}
      
      if(currentDirection == "Down"){termxDirectionMembers <- workinghmap$Name[(workinghmap$Function == hmapterm) & (workinghmap[,6] < 0)]} 
      if(currentDirection == "Down"){directGSEA <- workingGSEA[workingGSEA$NES < 0]}
      if(length(termxDirectionMembers) != 0){
        directGSEA$pathway <- gsub("%+.*", "", directGSEA$pathway)
        
        subsettedGSEA <- directGSEA[directGSEA$pathway %in% termxDirectionMembers,]
        LEs <- subsettedGSEA$leadingEdge
        Intersect <- Reduce(intersect, LEs)
        
        Histo <- NULL
        for(l in 1:length(LEs)){
          Histo <- c(Histo, LEs[[l]])
        }
        Histo <- paste(Histo, sep = "", collapse = "|")
        
        Name <- paste(termxDirectionMembers, sep="", collapse="_x_")
        Magnitude <- mean(subsettedGSEA$NES)
        Magnitude <- Magnitude[!is.na(Magnitude)] 
        Magnitude <- mean(Magnitude)
        Direction <- currentDirection
        
        ROWINFO <- c(celltype, contrast, hmapterm, currentDirection, Magnitude, Name, Histo, Intersect)
        workingResults <- rbind.na(workingResults, ROWINFO)
        rownames(workingResults) <- c(1:length(rownames(workingResults)))
      }
    }
  }
}

write.csv(workingResults[-1,], file="Parent-Child leading edge OVERVIEW SUMMARY")

library(dplyr)

#plotting histograms of merged LE's
histoResults <- workingResults[is.na(workingResults[,7]) == FALSE ,]
histoResults <- histoResults[histoResults[,7] != "",]
histoResults <- histoResults[,c(1:7)]
histoResults$CT_Contrast <- paste0(histoResults$Contrast, "_", histoResults$Cell_Type)


#EnsemblIDs to GeneNames
library(DESeq2)
load("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/seFullData_5p-3bp_3p-15bp.rData")
EnsemblID_DB <- as.data.frame(rowData(seFullData))[,c(1,6)]

setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/Leading edge analyses/Parent Child leading edges/")
workingMaxcounts <- matrix(ncol = 200, nrow= 1)

for(i in 1:nrow(histoResults)){
  #De-string union of LE terms into factors
  histoVec <- unlist(strsplit(histoResults$Union[i], "\\|"))
  histoGOterms <- unlist(strsplit(histoResults$Enriched_Terms[i], "_x_"))
  
  ###Need to create a dataframe for ggplot2
  histoDF <- as.data.frame(as.factor(histoVec))
  histoDF <- histoDF %>% group_by(`as.factor(histoVec)`) %>% summarize(no_rows=length(`as.factor(histoVec)`))
  names(histoDF) <- c("Gene_Name", "Count")
  histoDF$Gene_Name <- as.character(histoDF$Gene_Name)
  histoDF$Count <- as.numeric(histoDF$Count)
  clustername <- histoResults$Function[i]
  histplot <- ggplot(histoDF, aes(x=reorder(Gene_Name,-Count), y=Count)) + geom_bar(stat="identity") 
  histplot <- histplot + ggtitle(paste(unlist(histoGOterms), collapse=" ,")) + labs(caption= paste0(histoResults$CT_Contrast[i],"_",histoResults$Function[i])) + theme(text=element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1)) + ylab(paste0("Count (nterms=",length(histoGOterms),")"))
  
  pdf(paste0(as.character(histoResults$CT_Contrast[i]),"_",as.character(histoResults$Function[i]),"_", histoResults$Direction[i],".pdf"))
  print(histplot)
  dev.off()
  
  histoDF$Gene_Name <- as.character(histoDF$Gene_Name)
  
  #subset all genes overlapping across at least half of total # of genesets 
  maxcounts <- histoDF$Gene_Name[histoDF$Count >= floor(length(histoGOterms)/2)]
  mydata <- c(histoResults$Cell_Type[i], length(histoGOterms), clustername, histoResults$Direction[i], length(maxcounts), histoResults$Mean_NES[i], maxcounts)
  workingMaxcounts <- rbind.na(workingMaxcounts, mydata)
  
  #want further output files containing gene ID, name, LFC,  pvalue for each list of genes, and number of terms each gene was enriched in
  workingContrast <- histoResults$CT_Contrast[i]
  workingDirection <- histoResults$Direction[i]
  if(workingDirection == "Up"){directSign <- 1}
  if(workingDirection == "Down"){directSign <- -1}
  
  maxcountDF <- as.data.frame(maxcounts)
  colnames(maxcountDF) <- "Gene.name"
  maxcountDF$n_Enriched <- histoDF$Count[histoDF$Gene_Name %in% maxcountDF$Gene.name]
  
  workingDEGdf <- get(DEGlist[grep(workingContrast, DEGlist)])
  workingDEGdf <- as.data.frame(workingDEGdf)
  workingDEGdf <- workingDEGdf[workingDEGdf$log2FoldChange * directSign > 0,]
  workingDEGdf$EnsemblID <- rownames(workingDEGdf)
  workingDEGdf$Gene.name <- EnsemblID_DB$Gene.name[EnsemblID_DB$Gene.stable.ID %in% workingDEGdf$EnsemblID]
  
  contrastResults <- merge(maxcountDF, workingDEGdf, by="Gene.name")
  contrastResults <- contrastResults[order(contrastResults$padj),]
  write.csv(contrastResults, file=paste0(workingContrast,"_", clustername, "_", workingDirection,".csv"), row.names = FALSE)
}

colnames(workingMaxcounts)[1:6] <- c("Cell_Type", "#_of_terms","Function", "Direction", "#_of_Genes", "Mean_NES")
workingMaxcounts <- workingMaxcounts[-1,]
write.csv(workingMaxcounts, file="SUMMARY_Leading edge (most overlapping genes).csv")



############################EnrichmentMap leading edge analysis

#Need to load node files, cluster ID files, and GSEA/DEG results as before
setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output")
load("GSEA rankedlists (GSEAparam1).rData")
load("GSEA results (GSEAparam1).rData")
load("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/FDRtooled_results.rData")

setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/Leading edge analyses/EnrichmentMap leading edges/")
clusterIDfiles <- list.files(pattern="tsv")
nodefiles <- list.files(pattern="csv")

for(i in clusterIDfiles){
  assign(paste0(substr(i, 1, attr(gregexpr("^.*_+", i)[[1]], "match.length")-1),"_clusterIDs"), read.table(i, sep = "\t", header = TRUE, quote="\""))
}
for(i in nodefiles){
  assign(paste0(substr(i, 1, attr(gregexpr("^.*_+", i)[[1]], "match.length")-1),"_nodes"), read.csv(i))
}

#Some hard-coding, not ideal but it worked!
DEGlist <- c("fdrT_CT_PV_MDD", "fdrT_CT_PV_Bipolar", "fdrT_CT_PV_SCHIZ", "fdrT_CT_PYR23_MDD", "fdrT_CT_PYR23_Bipolar", "fdrT_CT_PYR23_SCHIZ", "fdrT_CT_PYR56_MDD", "fdrT_CT_PYR56_Bipolar", "fdrT_CT_PYR56_SCHIZ", "fdrT_CT_SST_MDD", "fdrT_CT_SST_Bipolar", "fdrT_CT_SST_SCHIZ", "fdrT_CT_VIP_MDD", "fdrT_CT_VIP_Bipolar", "fdrT_CT_VIP_SCHIZ")

GSEAreslist <- c("PV_MDD_gseaM1", "PV_Bipolar_gseaM1", "PV_SCHIZ_gseaM1", "PYR23_MDD_gseaM1", "PYR23_Bipolar_gseaM1", "PYR23_SCHIZ_gseaM1", "PYR56_MDD_gseaM1", "PYR56_Bipolar_gseaM1", "PYR56_SCHIZ_gseaM1", "SST_MDD_gseaM1", "SST_Bipolar_gseaM1", "SST_SCHIZ_gseaM1", "VIP_MDD_gseaM1", "VIP_Bipolar_gseaM1", "VIP_SCHIZ_gseaM1")

clusterIDlist <- c("PV_MDD_clusterIDs", "PV_Bipolar_clusterIDs", "PV_SCHIZ_clusterIDs", "PYR23_MDD_clusterIDs", "PYR23_Bipolar_clusterIDs", "PYR23_SCHIZ_clusterIDs", "PYR56_MDD_clusterIDs", "PYR56_Bipolar_clusterIDs", "PYR56_SCHIZ_clusterIDs", "SST_MDD_clusterIDs", "SST_Bipolar_clusterIDs", "SST_SCHIZ_clusterIDs", "VIP_MDD_clusterIDs", "VIP_Bipolar_clusterIDs", "VIP_SCHIZ_clusterIDs")

nodelist <- c("PV_MDD_nodes", "PV_Bipolar_nodes", "PV_SCHIZ_nodes", "PYR23_MDD_nodes", "PYR23_Bipolar_nodes", "PYR23_SCHIZ_nodes", "PYR56_MDD_nodes", "PYR56_Bipolar_nodes", "PYR56_SCHIZ_nodes", "SST_MDD_nodes", "SST_Bipolar_nodes", "SST_SCHIZ_nodes", "VIP_MDD_nodes", "VIP_Bipolar_nodes", "VIP_SCHIZ_nodes")

workingResults <- as.data.frame(matrix(nrow=1, ncol=107))
names(workingResults) <- c("Cell_Type","Contrast", "Cluster_Name", "Direction","Mean_NES", "Enriched_Terms", "Union", paste0("Gene_", 1:100))

for (i in 1:length(DEGlist)){
  celltype <- gsub("_.*$", "", gsub("^.*CT_", "", DEGlist[i]))
  contrast <- gsub("^.*_", "", gsub("^.*CT_", "", DEGlist[i]))
  workingClusterIDdf <- get(clusterIDlist[i])
  workingnodes <- get(nodelist[i])
  workingGSEA <- get(GSEAreslist[i])
  
  for(j in 1:nrow(workingClusterIDdf)){
    workingClusterID <- eval(parse(text=workingClusterIDdf[j,3]))
    workingCluster <- workingnodes[workingnodes$X__mclCluster %in% workingClusterID,]
    workingTerms <- workingCluster$EnrichmentMap..Name
    workingGSEACluster <- workingGSEA[workingGSEA$pathway %in% workingTerms, ]
    
    #Get LEs and summarize
    LEs <- workingGSEACluster$leadingEdge
    Intersect <- Reduce(intersect, LEs)
        
    Histo <- NULL
    for(l in 1:length(LEs)){
      Histo <- c(Histo, LEs[[l]])
    }
    Histo <- paste(Histo, sep = "", collapse = "|")
        
    Name <- paste(workingTerms, sep="", collapse="_x_")
    Magnitude <- mean(workingGSEACluster$NES)
    if(sign(Magnitude) == -1){Direction <- "Down"}
    if(sign(Magnitude) == 1){Direction <- "Up"}
        
    ROWINFO <- c(celltype, contrast, workingClusterIDdf$Cluster[j], Direction, Magnitude, Name, Histo, Intersect)
    workingResults <- rbind.na(workingResults, ROWINFO)
    rownames(workingResults) <- c(1:length(rownames(workingResults)))
  }
}

write.csv(workingResults[-1,], file="EnrichmentMap leading edge OVERVIEW SUMMARY.csv")

library(dplyr)

#plotting histograms of merged LE's
histoResults <- workingResults[is.na(workingResults[,7]) == FALSE ,]
histoResults <- histoResults[histoResults[,7] != "",]
histoResults <- histoResults[,c(1:7)]
histoResults$CT_Contrast <- paste0(histoResults$Cell_Type, "_", histoResults$Contrast)

rownames(histoResults) <- c(1:length(rownames(histoResults)))
#EnsemblIDs to GeneNames
library(DESeq2)
load("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/seFullData_5p-3bp_3p-15bp.rData")
EnsemblID_DB <- as.data.frame(rowData(seFullData))[,c(1,6)]

setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/Leading edge analyses/EnrichmentMap leading edges/")
workingMaxcounts <- matrix(ncol = 200, nrow= 1)

for(i in 1:nrow(histoResults)){
  #De-string union of LE terms into factors
  histoVec <- unlist(strsplit(histoResults$Union[i], "\\|"))
  histoGOterms <- unlist(strsplit(histoResults$Enriched_Terms[i], "_x_"))
  
  
  
  ###Need to create a dataframe for ggplot2
  histoDF <- as.data.frame(as.factor(histoVec))
  histoDF <- histoDF %>% group_by(`as.factor(histoVec)`) %>% summarize(no_rows=length(`as.factor(histoVec)`))
  names(histoDF) <- c("Gene_Name", "Count")
  histoDF$Gene_Name <- as.character(histoDF$Gene_Name)
  histoDF$Count <- as.numeric(histoDF$Count)
  clustername <- histoResults$Cluster_Name[i]
  histplot <- ggplot(histoDF, aes(x=reorder(Gene_Name,-Count), y=Count)) + geom_bar(stat="identity") 
  histplot <- histplot + ggtitle(paste(unlist(histoGOterms), collapse=" ,")) + labs(caption= paste0(histoResults$CT_Contrast[i],"_",histoResults$Cluster_Name[i])) + theme(text=element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1)) + ylab(paste0("Count (nterms=",length(histoGOterms),")"))
  
  pdf(paste0(as.character(histoResults$CT_Contrast[i]),"_",as.character(histoResults$Cluster_Name[i]),".pdf"))
  print(histplot)
  dev.off()
  
  histoDF$Gene_Name <- as.character(histoDF$Gene_Name)
  
  #subset all genes overlapping across at least half of total # of genesets 
  maxcounts <- histoDF$Gene_Name[histoDF$Count >= floor(length(histoGOterms)/2)]
  mydata <- c(histoResults$Cell_Type[i], length(histoGOterms), clustername, histoResults$Direction[i], length(maxcounts), histoResults$Mean_NES[i], maxcounts)
  workingMaxcounts <- rbind.na(workingMaxcounts, mydata)
  
  #want further output files containing gene ID, name, LFC,  pvalue for each list of genes, and number of terms each gene was enriched in
  workingContrast <- histoResults$CT_Contrast[i]
  workingDirection <- histoResults$Direction[i]
  maxcountDF <- as.data.frame(maxcounts)
  colnames(maxcountDF) <- "Gene.name"
  maxcountDF$n_Enriched <- histoDF$Count[histoDF$Gene_Name %in% maxcountDF$Gene.name]
  
  workingDEGdf <- get(DEGlist[grep(workingContrast, DEGlist)])
  workingDEGdf <- as.data.frame(workingDEGdf)
  workingDEGdf$EnsemblID <- rownames(workingDEGdf)
  workingDEGdf$Gene.name <- EnsemblID_DB$Gene.name[EnsemblID_DB$Gene.stable.ID %in% workingDEGdf$EnsemblID]
  
  contrastResults <- merge(maxcountDF, workingDEGdf, by="Gene.name")
  contrastResults <- contrastResults[order(contrastResults$padj),]
  write.csv(contrastResults, file=paste0(workingContrast,"_", clustername, ".csv"), row.names = FALSE)
}

colnames(workingMaxcounts)[1:6] <- c("Cell_Type", "#_of_terms","Function", "Direction", "#_of_Genes", "Mean_NES")
workingMaxcounts <- workingMaxcounts[-1,]
write.csv(workingMaxcounts, file="SUMMARY_Leading edge (most overlapping genes).csv")








