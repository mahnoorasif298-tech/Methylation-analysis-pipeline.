###############################################################################
# Longitudinal DNA methylation analysis pipeline (Illumina 450K arrays)
#
# Author : Mahnoor Asif
#
# Design : repeated measures on the same individuals at three time points
#            A1 = first time point (baseline)
#            A2 = second time point
#            A3 = third time point
#          Any other time point present in the input file is excluded
#          (only the stages listed in cfg$stages are analysed).
#
# Input  : one CSV file (path in cfg$input_csv) with a "TargetID" column
#          (CpG probe IDs) and, for every sample, one beta-value column and
#          one detection p-value column named
#              <Individual>_<Stage>_Beta   and   <Individual>_<Stage>_Pval
#
# Sections
#    1  Setup, configuration, helper functions
#    2  Load data and build sample metadata
#    3  Probe filtering
#    4  Quality control (PCA, beta density)
#    5  Probe-wise differential methylation (limma) and enrichment
#    6  Paired / longitudinal model (limma duplicateCorrelation) + trend
#    7  Later-stage analysis (A3 vs A2)
#    8  Region-level analysis (DMRcate)
#    9  Differential variability (missMethyl::varFit)
#   10  Co-methylation network analysis (WGCNA)
#   11  Methylation state transitions
#   12  Epigenetic clocks (Horvath, PhenoAge)
#   13  Integration of clocks, WGCNA modules and DV CpGs
#
# This repository contains CODE ONLY. No data, results, or sample identifiers.
###############################################################################


# =============================================================================
# 1. SETUP
# =============================================================================
options(stringsAsFactors = FALSE)
set.seed(2025)

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(missMethyl)
  library(DMRcate)
  library(WGCNA)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(maxprobes)
  library(methylclock)
  library(lme4)
  library(effectsize)
  library(mclust)
  library(ggalluvial)
  library(UpSetR)
  library(matrixStats)
  library(pheatmap)
  library(RColorBrewer)
  library(reshape2)
  library(ggplot2)
  library(tibble)
  library(purrr)
  library(tidyr)
  library(dplyr)   # loaded last so dplyr verbs are not masked
})

cfg <- list(
  input_csv        = "data/methylation_input.csv",  # NOT included in the repo
  out_dir          = "results",                     # git-ignored
  stages           = c("A1", "A2", "A3"),
  det_p_thresh     = 0.01,    # detection p-value cutoff
  det_p_min_prop   = 0.75,    # probe must pass in >= 75% of samples
  fdr              = 0.05,
  logfc_min        = 0.5,     # |logFC| (M-values) for later-stage CpGs
  dmr_lambda       = 1000,
  dmr_C            = 2,
  dv_k             = 4,       # k-means clusters of variance trajectories
  dv_clusters_test = c("1", "2", "3"),  # clusters tested for enrichment
  delta_var_min    = 0.05,    # increase in variance (A3 - A1) for enrichment
  top_unstable_n   = 1000,
  wgcna_top_n      = 20000,   # most variable CpGs entering WGCNA
  wgcna_soft_power = 9,       # NULL = choose from pickSoftThreshold (R2 >= 0.8)
  wgcna_min_module = 30,
  wgcna_merge_cut  = 0.25,
  wgcna_trait_cor  = 0.4,     # |module-trait correlation| to call a module
  hub_mm           = 0.8,     # module membership cutoff for hub CpGs
  hub_gs           = 0.4,     # gene significance cutoff for hub CpGs
  state_breaks     = c(-Inf, 0.3, 0.7, Inf),  # Low / Medium / High beta
  dir_delta        = 0.1,     # |delta beta| for hyper/hypo labels
  effect_delta     = 0.15,    # |delta beta| for high-confidence dynamic CpGs
  run_mclust       = TRUE     # model-based clustering (slow on large arrays)
)

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
outf <- function(f) file.path(cfg$out_dir, f)

# ---- helper functions -------------------------------------------------------

# beta -> M-value (exact 0 / 1 replaced to avoid infinite values)
beta2m <- function(b) {
  b[b <= 0] <- 1e-4
  b[b >= 1] <- 0.9999
  log2(b / (1 - b))
}

save_csv <- function(df, file, row.names = FALSE) {
  if (!is.null(df)) write.csv(df, outf(file), row.names = row.names)
}

save_gg <- function(p, file, w = 7, h = 6) {
  if (!is.null(p)) ggsave(outf(file), p, width = w, height = h, dpi = 300)
}

# gometh wrapper: skips very small CpG sets and never stops the pipeline
run_gometh <- function(sig, bg, collection) {
  if (length(sig) < 5) return(NULL)
  tryCatch(
    missMethyl::gometh(sig.cpg = sig, all.cpg = bg, collection = collection),
    error = function(e) {
      message("gometh (", collection, ") failed: ", conditionMessage(e))
      NULL
    }
  )
}

# dot plot of significant terms from a gometh result
dotplot_enrich <- function(df, term_col, title, top_n = 15) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- df[df$FDR < cfg$fdr, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df <- head(df[order(df$FDR), , drop = FALSE], top_n)
  terms <- unique(as.character(df[[term_col]]))
  df <- df[!duplicated(df[[term_col]]), , drop = FALSE]
  df$Term <- factor(df[[term_col]], levels = rev(terms))
  df$GeneRatio <- df$DE / df$N
  ggplot(df, aes(x = GeneRatio, y = Term, size = DE, color = -log10(FDR))) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = "steelblue", high = "darkred") +
    labs(title = title, x = "Gene ratio (DE / N)", y = NULL,
         size = "DE genes", color = "-log10(FDR)") +
    theme_bw(base_size = 12)
}

# GO + KEGG enrichment for a CpG set; writes CSVs and dot plots
enrich_both <- function(sig, bg, tag, title = tag) {
  res <- list(GO = run_gometh(sig, bg, "GO"), KEGG = run_gometh(sig, bg, "KEGG"))
  save_csv(res$GO,   paste0("GO_",   tag, ".csv"), row.names = TRUE)
  save_csv(res$KEGG, paste0("KEGG_", tag, ".csv"), row.names = TRUE)
  save_gg(dotplot_enrich(res$GO,   "TERM",        paste("GO:",   title)), paste0("GO_",   tag, ".pdf"))
  save_gg(dotplot_enrich(res$KEGG, "Description", paste("KEGG:", title)), paste0("KEGG_", tag, ".pdf"))
  res
}

# keep significant terms with at least min_N genes
filter_enrich <- function(df, term_col, min_N = 5) {
  if (is.null(df)) return(NULL)
  df <- df[df$FDR < cfg$fdr & df$N >= min_N, c(term_col, "N", "DE", "FDR"), drop = FALSE]
  df[order(df$FDR), , drop = FALSE]
}


# =============================================================================
# 2. LOAD DATA AND SAMPLE METADATA
# =============================================================================
raw <- as.data.frame(data.table::fread(cfg$input_csv, header = TRUE))
raw <- raw[complete.cases(raw), ]                           # drop rows with NA
raw <- raw[, colSums(is.na(raw)) == 0, drop = FALSE]        # drop columns with NA
stopifnot("TargetID" %in% colnames(raw))

beta.mat <- as.matrix(raw[, grep("Beta", colnames(raw))])
pval.mat <- as.matrix(raw[, grep("Pval", colnames(raw))])
rownames(beta.mat) <- rownames(pval.mat) <- raw$TargetID
stopifnot(ncol(beta.mat) == ncol(pval.mat))

strip_suffix <- function(x) sub("[._]?(Beta|Pval)$", "", x, ignore.case = TRUE)
if (!identical(strip_suffix(colnames(beta.mat)), strip_suffix(colnames(pval.mat)))) {
  warning("Beta and Pval columns are not in the same sample order - please check.")
}

# sample metadata parsed from column names: <Individual>_<Stage>_Beta
parse_samples <- function(nm) {
  core <- sub("[._]?Beta$", "", nm, ignore.case = TRUE)
  data.frame(Sample     = nm,
             Individual = sub("_.*", "", core),
             Stage      = sub(".*_", "", core))
}
pheno_all <- parse_samples(colnames(beta.mat))

# MDS of all samples / all stages (before any filtering)
M_all     <- beta2m(beta.mat)
stage_all <- factor(pheno_all$Stage)
pal       <- brewer.pal(max(3, nlevels(stage_all)), "Dark2")
pdf(outf("MDS_all_samples.pdf"), width = 7, height = 6)
plotMDS(M_all, top = 1000, gene.selection = "common", col = pal[stage_all],
        pch = 19, cex = 1.5, main = "MDS of DNA methylation (M-values)")
legend("topright", legend = levels(stage_all), col = pal[seq_len(nlevels(stage_all))],
       pch = 19, bty = "n")
dev.off()
rm(M_all)

# keep only the stages of interest
keep  <- pheno_all$Stage %in% cfg$stages
pheno <- pheno_all[keep, ]
pheno$Stage      <- factor(pheno$Stage, levels = cfg$stages)
pheno$Individual <- factor(make.names(pheno$Individual))
pheno$StageNum   <- as.integer(pheno$Stage)
rownames(pheno)  <- pheno$Sample
beta.mat <- beta.mat[, keep, drop = FALSE]
pval.mat <- pval.mat[, keep, drop = FALSE]


# =============================================================================
# 3. PROBE FILTERING
# =============================================================================
ann <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

filter_log <- list()
apply_filter <- function(m, keep_rows, label) {
  filter_log[[label]] <<- data.frame(Step = label, Before = nrow(m), After = sum(keep_rows))
  m[keep_rows, , drop = FALSE]
}

# (a) detection p-value
good_detp <- rowMeans(pval.mat < cfg$det_p_thresh, na.rm = TRUE) >= cfg$det_p_min_prop
b1 <- apply_filter(beta.mat, good_detp, "Detection p-value")

# (b) present in the 450K annotation
b2 <- apply_filter(b1, rownames(b1) %in% rownames(ann), "In 450K annotation")

# (c) probes overlapping known SNPs
a2 <- ann[rownames(b2), ]
b3 <- apply_filter(b2, is.na(a2$Probe_rs) | a2$Probe_rs == "", "No SNP at probe")

# (d) sex chromosomes
a3 <- ann[rownames(b3), ]
b4 <- apply_filter(b3, !(a3$chr %in% c("chrX", "chrY")), "Autosomes only")

# (e) CpG probes only (drop non-CpG "ch" probes)
b5 <- apply_filter(b4, grepl("^cg", rownames(b4)), "CpG probes only")

# (f) cross-reactive probes
xreactive  <- maxprobes::xreactive_probes(array_type = "450K")
beta.final <- apply_filter(b5, !(rownames(b5) %in% xreactive), "Not cross-reactive")

filter_summary <- do.call(rbind, filter_log)
print(filter_summary)
save_csv(filter_summary, "probe_filtering_summary.csv")
rm(b1, b2, b3, b4, b5, a2, a3)

stopifnot(identical(colnames(beta.final), pheno$Sample))
Mval <- beta2m(beta.final)

anno_slim <- ann[, c("Name", "chr", "pos", "UCSC_RefGene_Name",
                     "UCSC_RefGene_Group", "Relation_to_Island")]


# =============================================================================
# 4. QUALITY CONTROL
# =============================================================================
pca  <- prcomp(t(beta.final), scale. = TRUE)
expl <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pdf(outf("QC_PCA.pdf"), width = 7, height = 6)
barplot(expl[1:10], main = "Variance explained by top 10 PCs",
        xlab = "Principal component", ylab = "% variance")
plot(pca$x[, 1], pca$x[, 2], pch = 19, col = as.integer(pheno$Individual),
     xlab = sprintf("PC1 (%.1f%%)", expl[1]), ylab = sprintf("PC2 (%.1f%%)", expl[2]),
     main = "PCA: inter-individual variation")
legend("topright", legend = levels(pheno$Individual),
       col = seq_len(nlevels(pheno$Individual)), pch = 19, cex = 0.7)
dev.off()
print(summary(lm(pca$x[, 1] ~ pheno$Individual)))

dens <- apply(beta.final, 2, density)
pdf(outf("QC_beta_density.pdf"), width = 7, height = 6)
plot(dens[[1]], main = "Density of filtered beta values", xlab = "Beta value",
     col = 1, lwd = 2, ylim = c(0, max(sapply(dens, function(d) max(d$y)))))
for (i in seq_along(dens)[-1]) lines(dens[[i]], col = i, lwd = 2)
dev.off()
rm(dens)


# =============================================================================
# 5. PROBE-WISE DIFFERENTIAL METHYLATION (individual as fixed effect)
# =============================================================================
design_fixed <- model.matrix(~ 0 + Stage + Individual, data = pheno)
colnames(design_fixed) <- sub("^(Stage|Individual)", "", colnames(design_fixed))

cont_fixed <- makeContrasts(
  A2_vs_A1  = A2 - A1,
  A3_vs_A1  = A3 - A1,
  A3_vs_A2 = A3 - A2,
  levels = design_fixed
)

fitA  <- lmFit(Mval, design_fixed)
fitA2 <- eBayes(contrasts.fit(fitA, cont_fixed))

comparisons <- colnames(cont_fixed)
tt_A <- sapply(comparisons, function(cmp)
  topTable(fitA2, coef = cmp, number = Inf, adjust.method = "BH"), simplify = FALSE)
for (cmp in comparisons) save_csv(tt_A[[cmp]], paste0("TopTable_", cmp, ".csv"), row.names = TRUE)

sig_counts_A <- data.frame(
  Comparison = comparisons,
  Significant_CpGs = sapply(tt_A, function(x) sum(x$adj.P.Val < cfg$fdr))
)
print(sig_counts_A)

p <- ggplot(sig_counts_A, aes(Comparison, Significant_CpGs, fill = Comparison)) +
  geom_col() + theme_minimal(base_size = 14) + theme(legend.position = "none") +
  labs(title = "Significant CpGs per comparison", y = "Number of CpGs (FDR < 0.05)")
save_gg(p, "Significant_CpGs_per_comparison.pdf")

# volcano plots
volc <- do.call(rbind, lapply(comparisons, function(n)
  data.frame(logFC = tt_A[[n]]$logFC, adj.P.Val = tt_A[[n]]$adj.P.Val, Comparison = n)))
volc$Significant <- volc$adj.P.Val < cfg$fdr
p <- ggplot(volc, aes(logFC, -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey60", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(cfg$fdr), linetype = "dashed") +
  facet_wrap(~ Comparison, scales = "free") +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano plots", x = "log2 fold change (M-values)", y = "-log10(FDR)")
save_gg(p, "Volcano_all_comparisons.pdf", w = 12, h = 5)
rm(volc)

# heatmap: top 100 significant CpGs, A3 vs A2
sig_ids <- head(rownames(tt_A$A3_vs_A2)[tt_A$A3_vs_A2$adj.P.Val < cfg$fdr], 100)
if (length(sig_ids) > 1) {
  cols    <- pheno$Sample[pheno$Stage %in% c("A2", "A3")]
  ann_col <- data.frame(Stage = as.character(pheno[cols, "Stage"]), row.names = cols)
  pheatmap(beta.final[sig_ids, cols], scale = "row", show_rownames = FALSE,
           clustering_distance_cols = "correlation",
           annotation_col = ann_col,
           annotation_colors = list(Stage = c(A2 = "purple", A3 = "forestgreen")),
           main = "Top CpGs (A3 vs A2)", filename = outf("Heatmap_top100_A3_vs_A2.pdf"))
}

# missMethyl GO / KEGG for every comparison (background = all tested CpGs)
all_cpgs  <- rownames(Mval)
enrich_A  <- list()
for (cmp in comparisons) {
  sig <- rownames(tt_A[[cmp]])[tt_A[[cmp]]$adj.P.Val < cfg$fdr]
  enrich_A[[cmp]] <- enrich_both(sig, all_cpgs, paste0("missMethyl_", cmp), cmp)
}
enrich_summary <- data.frame(
  Comparison = comparisons,
  Sig_GO   = sapply(comparisons, function(c) { r <- enrich_A[[c]]$GO;   if (is.null(r)) 0 else sum(r$FDR < cfg$fdr) }),
  Sig_KEGG = sapply(comparisons, function(c) { r <- enrich_A[[c]]$KEGG; if (is.null(r)) 0 else sum(r$FDR < cfg$fdr) })
)
print(enrich_summary)
save_csv(enrich_summary, "Summary_pathways_all_comparisons.csv")

# annotate the top 100 CpGs per comparison; gene-based GO (BP) and Reactome
gene_enrich_A <- list()
for (cmp in comparisons) {
  top_ids <- head(rownames(tt_A[[cmp]]), 100)
  save_csv(data.frame(CpG = top_ids,
                      Gene = ann[top_ids, "UCSC_RefGene_Name"],
                      Relation_to_Island = ann[top_ids, "Relation_to_Island"],
                      Genomic_Region = ann[top_ids, "UCSC_RefGene_Group"]),
           paste0("Annotated_top100_", cmp, ".csv"))

  genes <- unique(unlist(strsplit(ann[top_ids, "UCSC_RefGene_Name"], ";")))
  genes <- genes[!is.na(genes) & genes != ""]
  entrez <- tryCatch(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)$ENTREZID,
                     error = function(e) character(0))
  if (length(entrez) >= 5) {
    gene_enrich_A[[cmp]] <- list(
      GO_BP = enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE),
      Reactome = enrichPathway(gene = entrez, organism = "human",
                               pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)
    )
    rx <- gene_enrich_A[[cmp]]$Reactome
    if (!is.null(rx) && nrow(as.data.frame(rx)) > 0) {
      ggsave(outf(paste0("Reactome_", cmp, ".pdf")),
             enrichplot::dotplot(rx, showCategory = 20) +
               ggtitle(paste("Reactome pathway enrichment -", cmp)),
             width = 8, height = 7)
    }
  }
}


# =============================================================================
# 6. PAIRED / LONGITUDINAL MODEL (within-individual correlation) + TREND
# =============================================================================
design_B <- model.matrix(~ 0 + Stage, data = pheno)
colnames(design_B) <- sub("^Stage", "", colnames(design_B))
block <- pheno$Individual

dupcor <- duplicateCorrelation(Mval, design_B, block = block)
message("Estimated within-individual correlation: ", round(dupcor$consensus.correlation, 4))

fitB <- lmFit(Mval, design_B, block = block, correlation = dupcor$consensus.correlation)
cont_B <- makeContrasts(
  A2_vs_A1  = A2 - A1,
  A3_vs_A1  = A3 - A1,
  A3_vs_A2 = A3 - A2,
  levels = design_B
)
fitB2 <- eBayes(contrasts.fit(fitB, cont_B))

tt_B <- list()
for (cmp in colnames(cont_B)) {
  tt <- topTable(fitB2, coef = cmp, number = Inf, sort.by = "P")
  tt$CpG <- rownames(tt)
  tt <- merge(tt, anno_slim, by.x = "CpG", by.y = "Name", all.x = TRUE, sort = FALSE)
  tt <- tt[order(tt$P.Value), ]
  tt_B[[cmp]] <- tt
  save_csv(tt, paste0("limma_longitudinal_", cmp, "_annotated.csv"))
}
sig_summary_B <- data.frame(Contrast = names(tt_B),
                            Significant_CpGs = sapply(tt_B, function(x) sum(x$adj.P.Val < cfg$fdr, na.rm = TRUE)))
print(sig_summary_B)
save_csv(sig_summary_B, "limma_longitudinal_significant_counts.csv")

# linear trend across ordered stages (A1 -> A2 -> A3 = 1, 2, 3)
design_trend <- model.matrix(~ StageNum, data = pheno)
dupcor_trend <- duplicateCorrelation(Mval, design_trend, block = block)
fit_trend <- eBayes(lmFit(Mval, design_trend, block = block,
                          correlation = dupcor_trend$consensus.correlation))
trend_tt <- topTable(fit_trend, coef = "StageNum", number = Inf)
trend_tt$CpG <- rownames(trend_tt)
save_csv(merge(trend_tt, anno_slim, by.x = "CpG", by.y = "Name", all.x = TRUE),
         "limma_trend_stage_numeric.csv")
save(fitB, fitB2, dupcor, dupcor_trend, fit_trend, file = outf("limma_fit_objects.RData"))


# =============================================================================
# 7. LATER-STAGE ANALYSIS (A3 vs A2)
# =============================================================================
a3_a2 <- tt_B$A3_vs_A2
save_csv(a3_a2, "A3_vs_A2_all_CpGs.csv")

sig_cpgs  <- a3_a2[a3_a2$adj.P.Val < cfg$fdr & abs(a3_a2$logFC) >= cfg$logfc_min, ]
hyper_ids <- sig_cpgs$CpG[sig_cpgs$logFC > 0]
hypo_ids  <- sig_cpgs$CpG[sig_cpgs$logFC < 0]
save_csv(sig_cpgs, "A3_vs_A2_significant_CpGs.csv")
save_csv(sig_cpgs[sig_cpgs$logFC > 0, ], "A3_vs_A2_hypermethylated_CpGs.csv")
save_csv(sig_cpgs[sig_cpgs$logFC < 0, ], "A3_vs_A2_hypomethylated_CpGs.csv")

post_enrich <- list(
  hyper = enrich_both(hyper_ids, all_cpgs, "stage3_hyper", "hypermethylated CpGs (A3 vs A2)"),
  hypo  = enrich_both(hypo_ids,  all_cpgs, "stage3_hypo",  "hypomethylated CpGs (A3 vs A2)")
)
for (d in names(post_enrich)) {
  save_csv(filter_enrich(post_enrich[[d]]$GO,   "TERM"),        paste0("GO_",   d, "_filtered.csv"))
  save_csv(filter_enrich(post_enrich[[d]]$KEGG, "Description"), paste0("KEGG_", d, "_filtered.csv"))
}

trend_sig <- trend_tt[trend_tt$adj.P.Val < cfg$fdr, ]
save_csv(trend_sig, "Stage_trend_significant_CpGs.csv")

# volcano
vdf <- a3_a2
vdf$Significant <- vdf$adj.P.Val < cfg$fdr & abs(vdf$logFC) >= cfg$logfc_min
p <- ggplot(vdf, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(color = Significant), alpha = 0.6) +
  theme_minimal() +
  labs(title = "Volcano plot: A3 vs A2", x = "log2 fold change (M-values)", y = "-log10(FDR)")
save_gg(p, "Volcano_A3_vs_A2.pdf")
rm(vdf)

# heatmap of the top 200 CpGs and mean trajectories of the top 10
top200 <- head(sig_cpgs$CpG[order(sig_cpgs$adj.P.Val)], 200)
if (length(top200) > 1) {
  cols <- pheno$Sample[pheno$Stage %in% c("A2", "A3")]
  pheatmap(Mval[top200, cols], scale = "row", show_rownames = FALSE,
           annotation_col = data.frame(Stage = as.character(pheno[cols, "Stage"]), row.names = cols),
           main = "Top CpGs (A3 vs A2, M-values)",
           filename = outf("Heatmap_top200_A3_vs_A2.pdf"))

  top10 <- top200[1:min(10, length(top200))]
  traj  <- reshape2::melt(Mval[top10, , drop = FALSE], varnames = c("CpG", "Sample"),
                          value.name = "Mvalue")
  traj  <- merge(traj, pheno[, c("Sample", "Stage")], by = "Sample")
  traj  <- aggregate(Mvalue ~ CpG + Stage, data = traj, FUN = mean)
  p <- ggplot(traj, aes(Stage, Mvalue, group = CpG, color = CpG)) +
    geom_line() + geom_point() + theme_minimal() +
    labs(title = "Mean methylation trajectories of top CpGs", y = "Mean M-value")
  save_gg(p, "Lineplot_top10_CpGs.pdf")
}


# =============================================================================
# 8. REGION-LEVEL ANALYSIS (DMRcate)
#    NOTE: this model has stage only (individual is not modelled here).
# =============================================================================
M_dmr <- Mval
M_dmr[!is.finite(M_dmr)] <- NA
M_dmr <- M_dmr[rownames(M_dmr) %in% rownames(ann), ]
M_dmr <- M_dmr[rowSums(is.na(M_dmr)) < ncol(M_dmr), ]

all_dmrs <- list()
for (cmp in colnames(cont_B)) {
  message("DMRcate: ", cmp)
  annot <- tryCatch(
    cpg.annotate(datatype = "array", object = M_dmr, what = "M",
                 arraytype = "450K", analysis.type = "differential",
                 design = design_B, contrasts = TRUE, cont.matrix = cont_B,
                 coef = cmp, fdr = cfg$fdr),
    error = function(e) { message("cpg.annotate failed for ", cmp, ": ", conditionMessage(e)); NULL })
  if (is.null(annot)) next

  dmrc <- tryCatch(dmrcate(annot, lambda = cfg$dmr_lambda, C = cfg$dmr_C),
                   error = function(e) { message("dmrcate failed for ", cmp, ": ", conditionMessage(e)); NULL })
  if (is.null(dmrc)) next

  all_dmrs[[cmp]] <- as.data.frame(extractRanges(dmrc, genome = "hg19"))
  save_csv(all_dmrs[[cmp]], paste0("DMRs_", cmp, ".csv"))
}
dmr_summary <- data.frame(
  Comparison = names(all_dmrs),
  Total_DMRs = sapply(all_dmrs, nrow),
  Significant_DMRs = sapply(all_dmrs, function(df)
    if ("min_smoothed_fdr" %in% names(df)) sum(df$min_smoothed_fdr < cfg$fdr) else NA_integer_)
)
print(dmr_summary)
save_csv(dmr_summary, "DMR_summary.csv")
rm(M_dmr)


# =============================================================================
# 9. DIFFERENTIAL VARIABILITY (missMethyl::varFit)
#    NOTE: topVar() uses its default coefficient here. Pairwise comparisons
#    can be obtained with contrasts.varFit() (see the missMethyl vignette).
# =============================================================================
fit_var <- varFit(Mval, design = design_B)
DV_all  <- topVar(fit_var, number = nrow(Mval))
DV_all$CpG <- rownames(DV_all)
DV_all$FDR <- p.adjust(DV_all$P.Value, method = "BH")
DV_sig <- DV_all[DV_all$FDR < cfg$fdr, ]
save_csv(DV_all, "DV_results_all.csv")
message("Differentially variable CpGs (FDR < ", cfg$fdr, "): ", nrow(DV_sig))

# per-stage beta variance for every CpG, and direction of change
stage_beta_var <- sapply(cfg$stages, function(s)
  rowVars(beta.final[, pheno$Stage == s, drop = FALSE]))
rownames(stage_beta_var) <- rownames(beta.final)

direction <- function(a, b) ifelse(a - b > 0, "Increased", ifelse(a - b < 0, "Decreased", "No change"))
var_dir <- data.frame(
  A2_vs_A1  = direction(stage_beta_var[, "A2"], stage_beta_var[, "A1"]),
  A3_vs_A1  = direction(stage_beta_var[, "A3"], stage_beta_var[, "A1"]),
  A3_vs_A2 = direction(stage_beta_var[, "A3"], stage_beta_var[, "A2"])
)
var_dir_long <- do.call(rbind, lapply(names(var_dir), function(n) {
  tb <- table(factor(var_dir[[n]], levels = c("Increased", "Decreased")))
  data.frame(Comparison = n, Direction = names(tb), Percentage = 100 * as.numeric(tb) / sum(tb))
}))
p <- ggplot(var_dir_long, aes(Comparison, Percentage, fill = Direction)) +
  geom_col() + scale_fill_manual(values = c("firebrick", "steelblue")) +
  theme_bw(base_size = 14) +
  labs(title = "CpGs with increased vs decreased variance", y = "Percentage of CpGs")
save_gg(p, "Variance_direction_percentage.pdf")

# trajectories of M-value variance for the DV CpGs
sig_dv_ids <- DV_sig$CpG
if (length(sig_dv_ids) >= cfg$dv_k * 5) {

  stage_var_M <- sapply(cfg$stages, function(s)
    rowVars(Mval[sig_dv_ids, pheno$Stage == s, drop = FALSE]))
  rownames(stage_var_M) <- sig_dv_ids
  stage_var_M <- as.data.frame(stage_var_M)

  # slope of a least-squares line through the three stage variances
  var_trends <- stage_var_M
  var_trends$Trend_Slope <- (stage_var_M[, 3] - stage_var_M[, 1]) / 2
  var_trends$Trend_Direction <- ifelse(var_trends$Trend_Slope > 0,
                                       "Increasing variability", "Decreasing variability")
  print(table(var_trends$Trend_Direction))

  set.seed(123)
  km <- kmeans(scale(stage_var_M), centers = cfg$dv_k)
  var_trends$Cluster <- factor(km$cluster)
  print(table(var_trends$Cluster))
  save_csv(var_trends, "DV_variance_trends.csv", row.names = TRUE)

  p <- ggplot(var_trends, aes(Trend_Slope)) + geom_histogram(bins = 50) +
    theme_bw() + labs(title = "Distribution of variance trends",
                      x = "Slope (change in variance per stage)")
  save_gg(p, "DV_trend_slope_histogram.pdf")

  # mean variance trajectory per cluster
  traj_long <- var_trends %>%
    tibble::rownames_to_column("CpG") %>%
    dplyr::select(CpG, dplyr::all_of(cfg$stages), Cluster) %>%
    tidyr::pivot_longer(dplyr::all_of(cfg$stages), names_to = "Stage", values_to = "Variance")
  traj_long$Stage <- factor(traj_long$Stage, levels = cfg$stages)
  cl_means <- traj_long %>%
    dplyr::group_by(Cluster, Stage) %>%
    dplyr::summarise(Mean_Variance = mean(Variance), .groups = "drop")
  p <- ggplot(cl_means, aes(Stage, Mean_Variance, group = Cluster)) +
    geom_line(linewidth = 1.2) + geom_point(size = 3) +
    facet_wrap(~ Cluster, scales = "free_y") + theme_bw(base_size = 14) +
    labs(title = "Mean variance trajectory per cluster", y = "Mean variance (M-values)")
  save_gg(p, "DV_cluster_trajectories.pdf", w = 9, h = 7)

  # GO / KEGG for the clusters of interest (background = all tested CpGs)
  cluster_enrich <- list()
  for (cl in cfg$dv_clusters_test) {
    ids <- rownames(var_trends)[var_trends$Cluster == cl]
    cluster_enrich[[cl]] <- enrich_both(ids, all_cpgs, paste0("DV_cluster_", cl),
                                        paste("variance cluster", cl))
  }

  # CpGs whose variance increases most at A3 relative to baseline
  delta_var  <- stage_var_M[, "A3"] - stage_var_M[, "A1"]
  names(delta_var) <- rownames(stage_var_M)
  top_inc    <- names(delta_var)[delta_var > cfg$delta_var_min]
  bg_dv      <- rownames(stage_var_M)  # background = DV CpGs (review if needed)
  inc_enrich <- enrich_both(top_inc, bg_dv, "DV_increased_at_A3", "variance increased at A3")

  # most unstable CpGs across all stages
  overall_var  <- apply(stage_var_M[, cfg$stages], 1, var)
  top_unstable <- names(sort(overall_var, decreasing = TRUE))[seq_len(min(cfg$top_unstable_n, length(overall_var)))]
  unst_enrich  <- enrich_both(top_unstable, bg_dv, "DV_most_unstable", "most unstable CpGs")

  unst_annot <- ann[top_unstable, c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group")]
  unst_annot$Gene_Symbol <- sapply(strsplit(unst_annot$UCSC_RefGene_Name, ";"), `[`, 1)
  unst_annot$OverallVariance <- overall_var[top_unstable]
  save_csv(unst_annot[order(-unst_annot$OverallVariance), ], "Most_unstable_CpGs_annotated.csv")
} else {
  message("Too few differentially variable CpGs for trajectory clustering; skipping.")
}


# =============================================================================
# 10. CO-METHYLATION NETWORK ANALYSIS (WGCNA)
# =============================================================================
allowWGCNAThreads()

cpg_var  <- setNames(rowVars(beta.final), rownames(beta.final))
top_cpgs <- names(sort(cpg_var, decreasing = TRUE))[seq_len(min(cfg$wgcna_top_n, length(cpg_var)))]
datExpr  <- t(beta.final[top_cpgs, ])                 # samples x CpGs

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
pheno_w <- pheno[rownames(datExpr), ]

pdf(outf("WGCNA_sample_tree.pdf"), width = 9, height = 6)
plot(hclust(dist(datExpr), method = "average"),
     main = "Sample clustering (outlier check)", sub = "", xlab = "")
dev.off()

traitMat <- model.matrix(~ 0 + Stage, data = pheno_w)
colnames(traitMat) <- sub("^Stage", "", colnames(traitMat))

# soft-threshold power
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)
sft_r2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
pdf(outf("WGCNA_soft_threshold.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], sft_r2, type = "n", xlab = "Soft threshold (power)",
     ylab = "Scale-free topology fit (signed R^2)")
text(sft$fitIndices[, 1], sft_r2, labels = powers, col = "red")
abline(h = 0.80, col = "blue")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n",
     xlab = "Soft threshold (power)", ylab = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

soft_power <- if (!is.null(cfg$wgcna_soft_power)) cfg$wgcna_soft_power else
  if (any(sft_r2 >= 0.8)) sft$fitIndices[which(sft_r2 >= 0.8)[1], 1] else 6
message("WGCNA soft-thresholding power: ", soft_power)

net <- blockwiseModules(
  datExpr, power = soft_power, networkType = "signed", TOMType = "signed",
  minModuleSize = cfg$wgcna_min_module, reassignThreshold = 0,
  mergeCutHeight = cfg$wgcna_merge_cut, numericLabels = TRUE,
  pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 0
)

# module colours and eigengenes are derived from the SAME labels, so that
# module names ("blue", "turquoise", ...) are used consistently downstream
moduleColors <- labels2colors(net$colors)
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
print(table(moduleColors))

pdf(outf("WGCNA_dendrogram.pdf"), width = 10, height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colours", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# module-trait relationships (biweight midcorrelation)
moduleTraitCor    <- bicor(MEs, traitMat, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))
rownames(moduleTraitCor) <- rownames(moduleTraitPvalue) <- sub("^ME", "", rownames(moduleTraitCor))
pheatmap(moduleTraitCor, cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = matrix(sprintf("%.2f\n(%.1g)", moduleTraitCor, moduleTraitPvalue),
                                  nrow = nrow(moduleTraitCor)),
         main = "Module-trait relationships (bicor)",
         filename = outf("WGCNA_module_trait_heatmap.pdf"), width = 7, height = 8)

# module eigengenes by stage
ME_long <- reshape2::melt(cbind(Sample = rownames(MEs), as.data.frame(MEs)), id.vars = "Sample")
ME_long$Stage <- pheno_w[as.character(ME_long$Sample), "Stage"]
p <- ggplot(ME_long, aes(Stage, value)) + geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y") + theme_bw() +
  labs(title = "Module eigengenes across stages", y = "Eigengene")
save_gg(p, "WGCNA_eigengenes_by_stage.pdf", w = 12, h = 9)

# CpG-level module membership (MM) and gene significance (GS) for A2
MM <- as.data.frame(cor(datExpr, MEs, use = "pairwise.complete.obs"))
colnames(MM) <- sub("^ME", "MM", colnames(MEs))
probeInfo <- data.frame(Probe  = colnames(datExpr),
                        Module = moduleColors,
                        GS     = as.numeric(cor(datExpr, traitMat[, "A2"], use = "pairwise.complete.obs")),
                        MM, check.names = FALSE)
save_csv(probeInfo, "WGCNA_probe_module_membership.csv")

# A2-associated modules and their hub CpGs
a2_tab <- data.frame(Module      = sub("^ME", "", colnames(MEs)),
                      Correlation = as.numeric(bicor(MEs, traitMat[, "A2"], use = "pairwise.complete.obs")))
a2_tab$Pvalue <- corPvalueStudent(a2_tab$Correlation, nSamples = nrow(datExpr))
a2_tab <- a2_tab[order(-abs(a2_tab$Correlation)), ]
save_csv(a2_tab, "WGCNA_A2_module_correlations.csv")

A2_modules <- a2_tab$Module[abs(a2_tab$Correlation) >= cfg$wgcna_trait_cor & a2_tab$Module != "grey"]
message("A2-associated modules: ", paste(A2_modules, collapse = ", "))

hub_list <- lapply(setNames(A2_modules, A2_modules), function(mod) {
  d <- probeInfo[probeInfo$Module == mod, ]
  d[abs(d[[paste0("MM", mod)]]) >= cfg$hub_mm & abs(d$GS) >= cfg$hub_gs, ]
})
hub_all <- dplyr::bind_rows(hub_list)
message("Hub CpGs identified: ", nrow(hub_all))

if (nrow(hub_all) > 0) {
  hub_annot <- merge(hub_all, anno_slim, by.x = "Probe", by.y = "Name", all.x = TRUE)
  save_csv(hub_annot, "WGCNA_A2_hub_CpGs_annotated.csv")

  # GO / KEGG per module (background = CpGs in the network)
  for (mod in names(hub_list)) {
    enrich_both(hub_list[[mod]]$Probe, colnames(datExpr), paste0("WGCNA_module_", mod),
                paste("module", mod))
  }

  # genomic context of hub CpGs
  hub_annot$Gene_Location <- "Intergenic"
  hub_annot$Gene_Location[grepl("TSS200|TSS1500", hub_annot$UCSC_RefGene_Group)] <- "Promoter"
  hub_annot$Gene_Location[grepl("Body", hub_annot$UCSC_RefGene_Group)] <- "Gene body"

  island_summary <- hub_annot %>%
    dplyr::filter(!is.na(Relation_to_Island), Relation_to_Island != "") %>%
    dplyr::group_by(Module, Relation_to_Island) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(Module) %>%
    dplyr::mutate(Percent = 100 * Count / sum(Count))
  gene_loc_summary <- hub_annot %>%
    dplyr::group_by(Module, Gene_Location) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(Module) %>%
    dplyr::mutate(Percent = 100 * Count / sum(Count))
  save_csv(gene_loc_summary, "WGCNA_hub_gene_location_by_module.csv")

  save_gg(ggplot(island_summary, aes(Module, Percent, fill = Relation_to_Island)) +
            geom_col() + theme_bw(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(title = "CpG island context of hub CpGs", y = "Percentage of CpGs", fill = "Context"),
          "WGCNA_hub_island_context.pdf", w = 8, h = 6)
  save_gg(ggplot(gene_loc_summary, aes(Module, Percent, fill = Gene_Location)) +
            geom_col() + theme_bw(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(title = "Gene-centric location of hub CpGs", y = "Percentage of CpGs", fill = "Location"),
          "WGCNA_hub_gene_location.pdf", w = 8, h = 6)
}

# bar plot of the strongest A2 module correlations
top_mod <- head(a2_tab, 8)
top_mod$Direction <- ifelse(top_mod$Correlation > 0, "Positive", "Negative")
save_gg(ggplot(top_mod, aes(reorder(Module, Correlation), Correlation, fill = Direction)) +
          geom_col(width = 0.7) + coord_flip() + theme_bw(base_size = 12) +
          labs(title = "Modules most associated with A2", x = "Module", y = "Correlation with A2"),
        "WGCNA_A2_top_modules.pdf")

save(moduleColors, MEs, moduleTraitCor, moduleTraitPvalue, probeInfo, hub_list,
     file = outf("WGCNA_results.RData"))


# =============================================================================
# 11. METHYLATION STATE TRANSITIONS (A1 -> A2 -> A3)
# =============================================================================
mean_beta <- data.frame(
  CpG = rownames(beta.final),
  sapply(cfg$stages, function(s) rowMeans(beta.final[, pheno$Stage == s, drop = FALSE])),
  check.names = FALSE
)
rownames(mean_beta) <- NULL

assign_state <- function(x) cut(x, breaks = cfg$state_breaks, labels = c("Low", "Medium", "High"))

mean_beta <- mean_beta %>%
  dplyr::mutate(
    A1_state  = assign_state(A1),
    A2_state = assign_state(A2),
    A3_state = assign_state(A3),
    delta_A2 = A2 - A1,
    delta_A3 = A3 - A2,
    Direction = dplyr::case_when(delta_A3 >  cfg$dir_delta ~ "Hypermethylation",
                                 delta_A3 < -cfg$dir_delta ~ "Hypomethylation",
                                 TRUE ~ "Stable"),
    Category = dplyr::case_when(
      A1_state != A2_state & A2_state == A3_state ~ "Persistent change",
      A1_state != A2_state & A1_state  == A3_state ~ "Reversible change",
      A1_state == A2_state & A2_state != A3_state ~ "Late change",
      TRUE ~ "Stable")
  )

# optional: model-based clustering of mean beta profiles into 3 states
if (cfg$run_mclust) {
  mc <- Mclust(mean_beta[, cfg$stages], G = 3, verbose = FALSE)
  mean_beta$Cluster <- mc$classification
  cl_mean <- rowMeans(aggregate(mean_beta[, cfg$stages],
                                by = list(Cluster = mean_beta$Cluster), FUN = mean)[, -1])
  cl_lab  <- setNames(c("Low", "Medium", "High")[rank(cl_mean)], sort(unique(mean_beta$Cluster)))
  mean_beta$Global_State <- unname(cl_lab[as.character(mean_beta$Cluster)])
}

# summaries and plots
cat_summary <- dplyr::count(mean_beta, Category, name = "n")
print(cat_summary)
save_gg(ggplot(cat_summary, aes(reorder(Category, n), n, fill = Category)) +
          geom_col(width = 0.7) + coord_flip() + theme_minimal(base_size = 13) +
          theme(legend.position = "none") +
          labs(title = "CpG methylation response categories", x = NULL, y = "Number of CpGs"),
        "State_category_counts.pdf")

dir_summary <- dplyr::count(mean_beta, Direction, name = "n")
dir_summary$Direction <- factor(dir_summary$Direction,
                                levels = c("Hypomethylation", "Stable", "Hypermethylation"))
save_gg(ggplot(dir_summary, aes(Direction, n, fill = Direction)) +
          geom_col(width = 0.7) + theme_minimal(base_size = 13) +
          theme(legend.position = "none") +
          labs(title = "Directionality of methylation change (A3 vs A2)", x = NULL, y = "Number of CpGs"),
        "State_direction_counts.pdf")

flow <- dplyr::count(mean_beta, A1_state, A2_state, A3_state, name = "n")
save_gg(ggplot(flow, aes(axis1 = A1_state, axis2 = A2_state, axis3 = A3_state, y = n)) +
          geom_alluvium(aes(fill = A1_state), alpha = 0.8) +
          geom_stratum(width = 0.25) +
          geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
          scale_x_discrete(limits = cfg$stages) +
          theme_minimal(base_size = 13) + labs(title = "CpG methylation state transitions"),
        "State_transitions_alluvial.pdf", w = 8, h = 6)

# annotation (CpG island context is collapsed from N/S shore/shelf)
anno_df <- data.frame(CpG = as.character(ann$Name),
                      Gene = as.character(ann$UCSC_RefGene_Name),
                      Relation_to_Island = sub("^[NS]_", "", as.character(ann$Relation_to_Island)),
                      stringsAsFactors = FALSE)
anno_df$Relation_to_Island <- factor(anno_df$Relation_to_Island,
                                     levels = c("Island", "Shore", "Shelf", "OpenSea"))
mean_beta_ann <- dplyr::left_join(mean_beta, anno_df, by = "CpG")

ctx <- mean_beta_ann %>%
  dplyr::filter(Category != "Stable", !is.na(Relation_to_Island)) %>%
  dplyr::count(Category, Relation_to_Island, name = "n")
save_gg(ggplot(ctx, aes(Relation_to_Island, n, fill = Category)) +
          geom_col(position = "fill") + coord_flip() + theme_minimal(base_size = 13) +
          labs(title = "Genomic context of dynamic CpGs", x = "CpG island context", y = "Proportion"),
        "State_genomic_context.pdf")

# gene-level aggregation
gene_long <- mean_beta_ann %>%
  dplyr::filter(!is.na(Gene), Gene != "") %>%
  tidyr::separate_rows(Gene, sep = ";")
gene_level <- gene_long %>%
  dplyr::filter(Category != "Stable") %>%
  dplyr::group_by(Gene, Category) %>%
  dplyr::summarise(mean_delta = mean(delta_A3), n_CpGs = dplyr::n(), .groups = "drop")
save_csv(gene_level, "State_gene_level_results.csv")

# high-confidence dynamic CpGs -> GO (BP) enrichment per category
high_conf <- gene_long %>%
  dplyr::filter(Category != "Stable",
                abs(delta_A2) >= cfg$effect_delta | abs(delta_A3) >= cfg$effect_delta)
gene_lists <- high_conf %>%
  dplyr::distinct(Gene, Category) %>%
  dplyr::group_by(Category) %>%
  dplyr::summarise(Genes = list(unique(Gene)), .groups = "drop") %>%
  tibble::deframe()

valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
run_GO_safe <- function(genes) {
  genes <- intersect(genes, valid_symbols)
  if (length(genes) < 5) return(NULL)
  enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
           pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
}
go_by_category <- purrr::map(gene_lists, run_GO_safe)
enrich_df <- purrr::imap_dfr(go_by_category, function(x, nm) {
  if (is.null(x)) return(NULL)
  df <- as.data.frame(x)
  if (nrow(df) > 0) df$Category <- nm
  df
})

# pathways significant ONLY in the "Reversible change" category
if (nrow(enrich_df) > 0) {
  path_sets <- enrich_df %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::distinct(Category, Description) %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(Categories = list(Category), .groups = "drop")
  unique_reversible <- path_sets %>%
    dplyr::filter(purrr::map_int(Categories, length) == 1,
                  purrr::map_chr(Categories, 1) == "Reversible change") %>%
    dplyr::pull(Description)
  unique_rev_tab <- enrich_df %>%
    dplyr::filter(Category == "Reversible change", Description %in% unique_reversible) %>%
    dplyr::arrange(p.adjust)
  save_csv(unique_rev_tab, "GO_unique_reversible_pathways.csv")
  if (nrow(unique_rev_tab) > 0) {
    save_gg(unique_rev_tab %>% dplyr::slice_head(n = 15) %>%
              ggplot(aes(reorder(Description, -Count), Count)) +
              geom_col(fill = "#2C7FB8") + coord_flip() + theme_minimal(base_size = 13) +
              labs(title = "Pathways unique to reversible-change CpGs", x = NULL, y = "Number of genes"),
            "GO_unique_reversible_pathways.pdf", w = 9, h = 6)
  }
}
save_csv(mean_beta_ann, "State_CpG_level_results.csv")


# =============================================================================
# 12. EPIGENETIC CLOCKS (Horvath, PhenoAge) AND WITHIN-INDIVIDUAL ACCELERATION
#     No chronological age is used: acceleration (EAA) is expressed relative
#     to each individual's own mean across the analysed stages.
#     NOTE: check the clock column names against your methylclock version.
# =============================================================================
clk <- DNAmAge(beta.final, clocks = c("Horvath", "Levine"))
stopifnot(nrow(clk) == nrow(pheno))
if ("id" %in% names(clk) && !identical(as.character(clk$id), pheno$Sample)) {
  warning("Clock output IDs differ from sample names; row order is assumed to match.")
}

clock_df <- pheno %>%
  dplyr::mutate(Horvath = clk$Horvath, PhenoAge = clk$Levine) %>%
  dplyr::group_by(Individual) %>%
  dplyr::mutate(Horvath_EAA = Horvath  - mean(Horvath,  na.rm = TRUE),
                Pheno_EAA   = PhenoAge - mean(PhenoAge, na.rm = TRUE)) %>%
  dplyr::ungroup()
save_csv(clock_df, "Clock_estimates.csv")

paired_stats <- function(df, var, s1, s2) {
  w <- df %>%
    dplyr::filter(Stage %in% c(s1, s2)) %>%
    dplyr::select(Individual, Stage, value = dplyr::all_of(var)) %>%
    tidyr::pivot_wider(names_from = Stage, values_from = value) %>%
    tidyr::drop_na()
  tst <- wilcox.test(w[[s1]], w[[s2]], paired = TRUE, exact = FALSE)
  es  <- effectsize::rank_biserial(w[[s1]], w[[s2]], paired = TRUE)
  data.frame(Variable = var, Comparison = paste(s2, "vs", s1), n_pairs = nrow(w),
             statistic = unname(tst$statistic), p_value = tst$p.value,
             rank_biserial = es[[1]])
}
clock_tests <- do.call(rbind, list(
  paired_stats(clock_df, "Pheno_EAA",   "A1",  "A2"),
  paired_stats(clock_df, "Pheno_EAA",   "A2", "A3"),
  paired_stats(clock_df, "Horvath_EAA", "A1",  "A2"),
  paired_stats(clock_df, "Horvath_EAA", "A2", "A3")
))
print(clock_tests)
save_csv(clock_tests, "Clock_paired_tests.csv")

# linear mixed-effects models (random intercept per individual)
lmer_pheno   <- lmer(Pheno_EAA   ~ Stage + (1 | Individual), data = clock_df)
lmer_horvath <- lmer(Horvath_EAA ~ Stage + (1 | Individual), data = clock_df)
print(summary(lmer_pheno))
print(summary(lmer_horvath))

p_a1_a2  <- signif(clock_tests$p_value[clock_tests$Variable == "Pheno_EAA" & clock_tests$Comparison == "A2 vs A1"], 3)
p_a2_a3 <- signif(clock_tests$p_value[clock_tests$Variable == "Pheno_EAA" & clock_tests$Comparison == "A3 vs A2"], 3)
ymax <- max(clock_df$Pheno_EAA, na.rm = TRUE)

save_gg(ggplot(clock_df, aes(Stage, Pheno_EAA, group = Individual)) +
          geom_line(alpha = 0.5) + geom_point(size = 2) + theme_minimal(base_size = 13) +
          labs(title = "PhenoAge acceleration across stages (per individual)",
               y = "Epigenetic age acceleration"),
        "Clock_PhenoAge_EAA_lines.pdf")
save_gg(ggplot(clock_df, aes(Stage, Pheno_EAA)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.6) +
          geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
          annotate("text", x = 1.5, y = ymax,       label = paste0("A2 vs A1: p = ", p_a1_a2)) +
          annotate("text", x = 2.5, y = ymax * 0.9, label = paste0("A3 vs A2: p = ", p_a2_a3)) +
          theme_bw(base_size = 13) +
          labs(title = "PhenoAge acceleration by stage", y = "PhenoAge EAA"),
        "Clock_PhenoAge_EAA_boxplot.pdf")


# =============================================================================
# 13. INTEGRATION: CLOCKS x WGCNA MODULES x DIFFERENTIAL VARIABILITY
# =============================================================================
common <- intersect(rownames(MEs), clock_df$Sample)
MEs_sub   <- MEs[common, , drop = FALSE]
clock_sub <- clock_df[match(common, clock_df$Sample), ]
stopifnot(identical(rownames(MEs_sub), clock_sub$Sample))

module_EAA <- data.frame(
  Module      = sub("^ME", "", colnames(MEs_sub)),
  Correlation = apply(MEs_sub, 2, function(x) cor(x, clock_sub$Pheno_EAA, use = "pairwise.complete.obs")),
  Pvalue      = apply(MEs_sub, 2, function(x) cor.test(x, clock_sub$Pheno_EAA)$p.value)
)
module_EAA$AdjP <- p.adjust(module_EAA$Pvalue, method = "BH")
module_EAA <- module_EAA[order(-abs(module_EAA$Correlation)), ]
save_csv(module_EAA, "WGCNA_module_vs_PhenoEAA.csv")

sig_mod <- module_EAA[module_EAA$AdjP < cfg$fdr, ]
if (nrow(sig_mod) > 0) {
  save_gg(ggplot(sig_mod, aes(reorder(Module, Correlation), Correlation, fill = Correlation > 0)) +
            geom_col(width = 0.8) + coord_flip() +
            scale_fill_manual(values = c("steelblue", "firebrick")) +
            theme_bw(base_size = 12) + theme(legend.position = "none") +
            labs(title = "Modules associated with PhenoAge acceleration", x = "Module",
                 y = "Correlation with PhenoAge EAA"),
          "WGCNA_modules_vs_PhenoEAA.pdf")
}

# eigengene heat map across stages
MEs_mat <- as.matrix(MEs_sub)
ann_state <- data.frame(State = clock_sub$Stage, row.names = clock_sub$Sample)
pheatmap(MEs_mat, annotation_row = ann_state, scale = "none", clustering_method = "average",
         annotation_colors = list(State = c(A1 = "#4DAF4A", A2 = "#377EB8", A3 = "#E41A1C")),
         main = "Module eigengenes across stages",
         filename = outf("WGCNA_eigengene_heatmap.pdf"), width = 8, height = 9)

# CpG sets: ageing-associated A2 modules, DV CpGs, stage-associated module CpGs
A2_aging_modules <- intersect(A2_modules, sig_mod$Module)
aging_hub   <- probeInfo[probeInfo$Module %in% A2_aging_modules & abs(probeInfo$GS) >= cfg$hub_gs, ]
aging_CpGs  <- unique(aging_hub$Probe)
dv_CpGs     <- unique(DV_sig$CpG)
module_CpGs <- unique(probeInfo$Probe[probeInfo$Module %in% A2_aging_modules])
save_csv(aging_hub, "A2_aging_hub_CpGs.csv")

if (length(aging_CpGs) > 0 && length(dv_CpGs) > 0 && length(module_CpGs) > 0) {
  all_sets <- unique(c(aging_CpGs, dv_CpGs, module_CpGs))
  upset_data <- data.frame(Aging     = as.integer(all_sets %in% aging_CpGs),
                           DV        = as.integer(all_sets %in% dv_CpGs),
                           Modules = as.integer(all_sets %in% module_CpGs))
  pdf(outf("UpSet_aging_DV_modules.pdf"), width = 8, height = 6)
  print(upset(upset_data, sets = c("Aging", "DV", "Modules"), order.by = "freq",
              keep.order = TRUE, mainbar.y.label = "CpG intersections",
              sets.x.label = "CpGs per set"))
  dev.off()
} else {
  message("One or more CpG sets are empty; UpSet plot skipped.")
}


# =============================================================================
# SESSION INFO (for reproducibility)
# =============================================================================
writeLines(capture.output(sessionInfo()), outf("sessionInfo.txt"))
message("Pipeline finished. Outputs written to: ", normalizePath(cfg$out_dir))
