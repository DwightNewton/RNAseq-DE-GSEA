library(fgsea)
library(DESeq2)
options(stringsAsFactors = FALSE)
options(scipen = 999)


#First, just get numbers of up/down pathways in each cell-type/contrast and do some venn diagrams for each direction (and overall)
setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/")

filelist <- list.files(pattern=".rData")
for(i in filelist){load(i)}

setwd("Parent child/")
load("ParentChild input files.rData")

setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/")
GMTfile <- gmtPathways("Human_GO_AllPathways_with_GO_iea_April_01_2020_symbol.gmt")
GMT_working <- lapply(GMTfile, paste0, collapse=" | ")
GMT_df <- data.frame("Pathway" = names(GMT_working), "Genes" = unlist(GMT_working))
GMT_df$Size <- sapply(strsplit(GMT_df$Genes, split="|", fixed=TRUE), length)
setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/Enrichment Map/")


#For each contrast, need "pos enrich", "neg enrich", GMT, "ranks", and "expressions" files

datalist <- c("PV_MDD_gseaM1", "PV_Bipolar_gseaM1", "PV_SCHIZ_gseaM1","PYR23_MDD_gseaM1", "PYR23_Bipolar_gseaM1", "PYR23_SCHIZ_gseaM1","PYR56_MDD_gseaM1","PYR56_Bipolar_gseaM1", "PYR56_SCHIZ_gseaM1","SST_MDD_gseaM1", "SST_Bipolar_gseaM1", "SST_SCHIZ_gseaM1","VIP_MDD_gseaM1", "VIP_Bipolar_gseaM1", "VIP_SCHIZ_gseaM1")
rnklist <- c("PV_MDD_ranklist", "PV_Bipolar_ranklist", "PV_SCHIZ_ranklist","PYR23_MDD_ranklist", "PYR23_Bipolar_ranklist", "PYR23_SCHIZ_ranklist","PYR56_MDD_ranklist","PYR56_Bipolar_ranklist", "PYR56_SCHIZ_ranklist","SST_MDD_ranklist", "SST_Bipolar_ranklist", "SST_SCHIZ_ranklist","VIP_MDD_ranklist", "VIP_Bipolar_ranklist", "VIP_SCHIZ_ranklist")

for(i in 1:length(datalist)){
  workingGSEA <- get(datalist[i])
  workingcontrast <- gsub("_g.*$", "", datalist[i])
  workingrnk <- get(rnklist[i])
    
  dir.create(paste0("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/Enrichment Map/", workingcontrast))
  setwd(workingcontrast)
  
  workingGSEA <- merge(workingGSEA, GMT_df, by.x = "pathway", by.y = "Pathway")
  
  #Pos enrich and Neg enrich
  Enrich <- data.frame("NAME"=workingGSEA$pathway, 
                          "GS<br> follow link to MSigDB"=workingGSEA$pathway, 
                          "GS DETAILS" = NA, 
                          "SIZE"=workingGSEA$Size,
                          "ES"=workingGSEA$ES,
                          "NES"=workingGSEA$NES,
                          "NOM p-val"=workingGSEA$pval,
                          "FDR q-val"=workingGSEA$padj,
                          "FWER p-val"=workingGSEA$padj,
                          "RANK AT MAX"=NA,
                          "LEADING EDGE"=NA)
  #Separate terms and get max-rank of leading-edge genes
  PosEnrich <- Enrich[Enrich$NES > 0,]
  PosRanks <- PosEnrich[,c("NAME", "RANK.AT.MAX")]
  PosRanks <- merge(PosRanks, workingGSEA, by.x="NAME", by.y="pathway")
  PosRanks$RANK.AT.MAX <- sapply(PosRanks$leadingEdge, match, table=names(workingrnk))
  PosRanks$ACTUAL.MAX.RANK <- sapply(PosRanks$RANK.AT.MAX, max)
  PosEnrich$RANK.AT.MAX <- PosRanks$ACTUAL.MAX.RANK[PosRanks$NAME %in% PosEnrich$NAME]
  
  NegEnrich <- Enrich[Enrich$NES < 0,]
  NegRanks <- NegEnrich[,c("NAME", "RANK.AT.MAX")]
  NegRanks <- merge(NegRanks, workingGSEA, by.x="NAME", by.y="pathway")
  NegRanks$RANK.AT.MAX <- sapply(NegRanks$leadingEdge, match, table=names(rev(workingrnk)))
  NegRanks$ACTUAL.MAX.RANK <- sapply(NegRanks$RANK.AT.MAX, max)
  NegEnrich$RANK.AT.MAX <- NegRanks$ACTUAL.MAX.RANK[NegRanks$NAME %in% NegEnrich$NAME]
  
  #specify in write.table
  names(Enrich) <- c("NAME","GS<br> follow link to MSigDB","GS DETAILS","SIZE","ES","NES","NOM p-val","FDR q-val","FWER p-val","RANK AT MAX","LEADING EDGE")
 
  #Ranks file, just gene names, "null" columns, and ranking value 
  RanksFile <- data.frame("NAME"=names(workingrnk),
                          "DESCRIPTION"="null",
                          "GENE_SYMBOL"="null",
                          "GENE_TITLE"="null",
                          "SCORE"=workingrnk)

  #Expressions file is just the ranked list again, but only the gene name and ranking value
  write.table(workingrnk, paste0("rankedlist.rnk"), quote = FALSE, col.names = FALSE, sep="\t")
  
  #Output other files(Seem to be tab-separated)
  # write.table(RanksFile, paste0("ranked_gene_list_na_pos_versus_na_neg_", workingcontrast, ".xls"), quote = FALSE, row.names = FALSE, sep="\t")
  write.table(PosEnrich, paste0("gsea_report_for_na_pos_", ".xls"), quote = FALSE, row.names = FALSE, sep="\t")
  write.table(NegEnrich, paste0("gsea_report_for_na_neg_", ".xls"), quote = FALSE, row.names = FALSE, sep="\t")

  setwd("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/Full Dataset/Differential Expression/GSEA/GSEA output/Enrichment Map/")
}





