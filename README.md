# RNAseq-DE-GSEA
This repository will cover many differential expression (DE) based downstream analyses of RNAseq data. DESeq2 was used for differential expression analysis, gene-set enrichment analysis (GSEA) for pathway enrichemnt, and a combination of parent-child and network-based visualizations were used to interpret the results. Network-based visualizations require Cytoscape to visualize output.


1. For datasets with detailed clinical/demographic data, explore contribution of relevant covariates.
   * `Covariate analysis.R`
   * In this post-mortem RNAseq dataset, RNA/tissue quality measures (pH, PMI, RIN), age, sex, and cause of death (suicide vs non-suicide) are of partiular interest.
   * If medication status and co-morbid disorders were more strongly represented, these would be additional candidates.
   * Main outputs are: plots of % gene variance explained for each covariate, lists of top *n* genes explained by covariates of interest, correlation matrices and clustering of covariates by co-correlation
