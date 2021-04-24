# RNAseq-DE-GSEA
This repository will cover many differential expression (DE) based downstream analyses of RNAseq data. DESeq2 was used for differential expression analysis, gene-set enrichment analysis (GSEA) for pathway enrichemnt, and a combination of parent-child and network-based visualizations were used to interpret the results. Network-based visualizations require Cytoscape to visualize output.

These analyses were conducted on the Sibille Lab "PITT Tetrad" dataset. Laser-dissected populations of L2/3 Pyrmaidal, L5/6 Pyramidal, SST-neurons, PV-neurons, and VIP-neurons were collected (130 cells/sample) and sequenced on the NovaSeq platform. Cells were dissected from the subgenual anterior cingulate cortex, from a post-mortem cohort containing 19 tetrads  

1. For datasets with detailed clinical/demographic data, explore contribution of relevant covariates.
   * `Covariate analysis.R`
   * In this post-mortem RNAseq dataset, RNA/tissue quality measures (pH, PMI, RIN), age, sex, and cause of death (suicide vs non-suicide) are of partiular interest.
   * If medication status and co-morbid disorders were more strongly represented, these would be additional candidates.
   * Main outputs are: plots of % gene variance explained for each covariate, lists of top *n* genes explained by covariates of interest, correlation matrices and clustering of covariates by co-correlation
