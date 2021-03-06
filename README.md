# RNAseq-DE-GSEA
This repository will cover many differential expression (DE) based downstream analyses of RNAseq data. DESeq2 was used for differential expression analysis, gene-set enrichment analysis (GSEA) for pathway enrichemnt, and a combination of parent-child and network-based visualizations were used to interpret the results. Network-based visualizations require Cytoscape to visualize output.

These analyses were conducted on the Sibille Lab "PITT Tetrad" dataset. Laser-dissected populations of L2/3 Pyramidal, L5/6 Pyramidal, SST-neurons, PV-neurons, and VIP-neurons were collected (130 cells/sample) and sequenced on the NovaSeq platform. Cells were dissected from the subgenual anterior cingulate cortex, from a post-mortem cohort containing 19 matched tetrads of major depressive disorder (MDD), bipolar disorder (BPD), schizophrenia (SCZ), and control subjects.

1. For datasets with detailed clinical/demographic data, explore contribution of relevant covariates.
   * `Covariate analysis.R`
   * In this post-mortem RNAseq dataset, RNA/tissue quality measures (pH, PMI, RIN), age, sex, and cause of death (suicide vs non-suicide) are of partiular interest.
   * If medication status and co-morbid disorders were more strongly represented, these would be additional candidates.
   * Outputs: plots of % gene variance explained for each covariate, lists of top *n* genes explained by covariates of interest, correlation matrices and clustering of covariates by co-correlation

2. Perform DE analysis with DESeq2
   * `Differential expression.R`
   * DE script can be seperated by cell-type and parallelized to save time, though this typically doesn't run overly long (~1-2 hours on CAMH cluster).
   * Saves results files as .rData files in DESeq2 format, gives most flexibility for downstream uses.
   * `Cell-type heatmaps and venns.R`: Extracts venn digrams of detected genes across cell-types and heatmaps of gene expression across cell-types (all genes and those at *n*% FDR).
   * `Contrast lists and FPKM table.R`: Exports FPKM (counts normalized for sequencing depth and gene length) matrix, and DE genelists showing fold-change and q-value for MDD-Control, BPD-Control, and SCZ-Control contrasts in each cell-type.

3. GSEA to determine biological pathways altered in each cell-type.
   * `Within-disorder GSEA.R` - performs GSEA for each contrast and cell-type.
   * A .gmt file is required, containing the universe of gene-sets/pathways you wish to query against. A good resource is: http://baderlab.org/GeneSets.

4. Parent-child analysis of GSEA results.
   * `Parent-child analysis and visualization.R` generates many outputs: tables summarizing # of GSEA hits at different significant thresholds across GSEA settings, cell-types, and contrasts. It also performs a parent-child analysis examining the degree of enrichment in biological functions of interest leveraging the Gene Ontology strucutre.
   * Gene Ontology (GO) has a heirarchal structure, allowing for analysis of parent (i.e. higher-up) and child (i.e. lower-down) relationships.
   * GO-terms representing many biological functions relevant to cell biology and psychiatric disorders were used as parent-terms.
   * Any significantly enriched GO-term in any contrast that (1) was one of the parent terms or (2) was a child-term of that parent was categorized as enriched for that function.
   * Heatmaps visualizing the results, and overlaps across constrasts, are generated.

5. EnrichmentMap Visualization of GSEA results.
   * `EnrichmentMap output script.R` is a handy script that generates "traditional" GSEA-style output files that are compatible with EnrichmentMap (Cytoscape plugin).
   * `Cytoscape import commands.txt.` imports these networks into Cytoscape using the Cytoscape command language since the GUI-based import functions can crash with large files.

6. Leading-edge analysis of parent-child and EnrichmentMap analyses.
   * `Leading-edge analyses.R` performs important downstream analyses which identify key genes of interest.
   * For a given Cytoscape cluster or Parent-Child function of interest, GSEA leading-edge genes (genes which drove enrichment of the function) are identified which overlap across a majority of pathway hits in that cluster/function.
   * Also generates histograms showing the most overlapping leading-edge genes in each cluster/function.
