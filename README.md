# DNA methylation analysis pipeline (Illumina 450K, longitudinal design)

R pipeline for a repeated-measures DNA methylation study (same individuals at
three time points: **A1**, **A2**, **A3**).

**This repository contains code only.** No data, results, or sample identifiers
are included. Supply your own input file.

## What it does
Probe filtering (detection p-value, SNPs, sex chromosomes, non-CpG, cross-reactive)
-> QC (MDS, PCA, beta densities) -> limma differential methylation (fixed-effect and
`duplicateCorrelation` models, stage trend) -> GO/KEGG (missMethyl) and Reactome ->
DMRcate -> differential variability (`varFit`) -> WGCNA co-methylation modules and
hub CpGs -> methylation state transitions -> Horvath/PhenoAge epigenetic clocks ->
integration of clocks, modules and variable CpGs.

## Input format
A CSV (`cfg$input_csv` in the script) with a `TargetID` column plus, per sample,
`<Individual>_<Stage>_Beta` and `<Individual>_<Stage>_Pval` columns, where
`<Stage>` is one of `A1`, `A2`, `A3`.

## Usage
```r
# edit the cfg list at the top of methylation_pipeline.R, then:
source("methylation_pipeline.R")
```
Outputs are written to `results/`.

## Dependencies
R >= 4.1 with Bioconductor:
```r
BiocManager::install(c("limma","minfi","IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "missMethyl","DMRcate","WGCNA","clusterProfiler","org.Hs.eg.db","ReactomePA",
  "methylclock","matrixStats"))
install.packages(c("data.table","lme4","effectsize","mclust","ggalluvial","UpSetR",
  "pheatmap","RColorBrewer","reshape2","ggplot2","tibble","purrr","tidyr","dplyr"))
remotes::install_github("markgene/maxprobes")
```

## Suggested .gitignore
```
data/
results/
*.csv
*.RData
*.rds
*.idat
.Rhistory
.Rproj.user/
```
