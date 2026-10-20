# Scenic-results

<!-- badges: start -->

<!-- badges: end -->

This repo has the data and scripts needed for performing marker gene and marker peak enrichment on the target genes and target regions of SCENIC+ predicted regulons. SCENIC+ was run on each of the three stages of differentiation- SC, DE and S2 on the WT and KO combined dataset (unpaired scRNA and scATAC-seq). The directories contain the following files:

-   `data`: The direct eRegulons predicted by SCENIC+ for each of the three stages

    -   Differentially accessible marker peaks between KO vs WT from scATAC-seq data per stage

    -   Differntially expressed genes from bulk RNA-seq per stage

-   `functions:` All functions needed to run the scripts

-   `scripts`: Scripts for gene enrichment and peak enrichment
