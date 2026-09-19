#!/usr/bin/env Rscript

# ============================================================================
# edgeR Differential Gene Expression Analysis — OSCC (Tumor vs Normal)
#
# Features:
#   - TMM normalisation, filterByExpr gene filtering
#   - glmQLFit / glmQLFTest (quasi-likelihood F-test)
#   - Optional batch-covariate adjustment (boolean toggle)
#   - Exploratory PCA (PC1×PC2 and PC1×PC3) on variance-stabilised logCPM
#   - Volcano, MA, BCV, QL-dispersion, sample-distance, and library-size plots
#   - Heatmaps: top-25 up, top-25 down, and combined (FDR ≤ 0.05, ranked by logFC)
#   - Expression boxplots for top DE genes
#   - Full and significant result tables, normalised count exports
#
# References:
#   Chen Y, Lun ATL, Smyth GK (2016). "From reads to genes to pathways:
#     differential expression analysis of RNA-Seq experiments using Rsubread
#     and the edgeR quasi-likelihood pipeline." F1000Research 5:1438.
#   Robinson MD, McCarthy DJ, Smyth GK (2010). "edgeR: a Bioconductor
#     package for differential expression analysis of digital gene
#     expression data." Bioinformatics 26(1):139–140.
# ============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(matrixStats)
  library(gridExtra)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

# ============================================================================
# 1. USER-CONFIGURABLE SETTINGS
# ============================================================================

# --- Input / output paths ---------------------------------------------------
counts_file <- ifelse(
  length(commandArgs(trailingOnly = TRUE)) >= 1,
  commandArgs(trailingOnly = TRUE)[1],
  "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project/03_processed/counts/OSCC_counts_matched.tsv"
)
meta_file <- ifelse(
  length(commandArgs(trailingOnly = TRUE)) >= 2,
  commandArgs(trailingOnly = TRUE)[2],
  "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project/03_processed/metadata/OSCC_metadata_matched.tsv"
)
outdir <- ifelse(
  length(commandArgs(trailingOnly = TRUE)) >= 3,
  commandArgs(trailingOnly = TRUE)[3],
  "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project/05_results/edgeR"
)

# --- Biological design ------------------------------------------------------
group_col  <- "Tissue_Type"
normal_lvl <- "Normal"
tumor_lvl  <- "Tumor"

# --- Batch correction toggle ------------------------------------------------
# TRUE  = include Plate + TSS in design matrix  (~Plate + TSS + Tissue_Type)
# FALSE = simple design                         (~Tissue_Type)
use_batch <- FALSE
batch_cols <- c("Plate", "TSS")

# --- Filtering / significance thresholds ------------------------------------
fdr_cutoff <- 0.05
lfc_cutoff <- 1

# --- Plot controls -----------------------------------------------------------
top_n_heatmap       <- 50     # overall top-DE heatmap
top_n_directional   <- 25     # up / down heatmaps
top_n_label_volcano <- 15
top_n_boxplot       <- 12

# ============================================================================
# 2. DIRECTORY SETUP & HELPERS
# ============================================================================

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
plot_dir <- file.path(outdir, "plots")
res_dir  <- file.path(outdir, "results")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir,  showWarnings = FALSE, recursive = TRUE)

# --- Run log ----------------------------------------------------------------
log_file <- file.path(outdir, "edgeR_run_log.txt")
sink(log_file, split = TRUE)
cat("================================================================\n")
cat("edgeR OSCC DGEA — run started:", format(Sys.time()), "\n")
cat("================================================================\n")
cat("Counts file :", counts_file, "\n")
cat("Metadata file:", meta_file, "\n")
cat("Output dir   :", outdir, "\n")
cat("Batch adjust :", use_batch, "\n\n")

# --- Publication ggplot theme -----------------------------------------------
theme_pub <- theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "black"),
    legend.title  = element_text(face = "bold"),
    panel.grid    = element_blank()
  )

# --- Dual-format saver (PDF + PNG) ------------------------------------------
save_pdf_png <- function(plot_obj, filename_base, width = 8, height = 6) {
  ggsave(file.path(plot_dir, paste0(filename_base, ".pdf")),
         plot = plot_obj, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(plot_dir, paste0(filename_base, ".png")),
         plot = plot_obj, width = width, height = height, dpi = 300)
}

# --- Pheatmap dual-format saver ---------------------------------------------
save_pheatmap <- function(pheatmap_expr, filename_base, width = 10, height = 11) {
  pdf(file.path(plot_dir, paste0(filename_base, ".pdf")), width = width, height = height)
  eval(pheatmap_expr)
  dev.off()

  png(file.path(plot_dir, paste0(filename_base, ".png")),
      width = width * 300, height = height * 300, res = 300)
  eval(pheatmap_expr)
  dev.off()
}

# --- Factor helper ----------------------------------------------------------
safe_factor <- function(x, ref = NULL) {
  x <- as.factor(x)
  if (!is.null(ref) && ref %in% levels(x)) x <- relevel(x, ref = ref)
  x
}

# ============================================================================
# 3. LOAD DATA
# ============================================================================

counts_df <- read.delim(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
meta_df   <- read.delim(meta_file,   check.names = FALSE, stringsAsFactors = FALSE)

# --- Validate required columns ----------------------------------------------
stopifnot(
  "'gene_name' column missing from counts file" = "gene_name" %in% colnames(counts_df),
  "'File_ID' column missing from metadata"      = "File_ID"   %in% colnames(meta_df),
  "Group column missing from metadata"          = group_col   %in% colnames(meta_df)
)

# --- Match samples between counts and metadata ------------------------------
sample_ids <- intersect(colnames(counts_df)[-1], meta_df$File_ID)
if (length(sample_ids) < 4) stop("Too few shared samples after matching counts ↔ metadata.")

counts_df <- counts_df[, c("gene_name", sample_ids)]
meta_df   <- meta_df[match(sample_ids, meta_df$File_ID), , drop = FALSE]
stopifnot(identical(colnames(counts_df)[-1], meta_df$File_ID))

# --- Build numeric count matrix ---------------------------------------------
gene_ids  <- counts_df$gene_name
count_mat <- as.matrix(counts_df[, -1, drop = FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- make.unique(gene_ids)
colnames(count_mat) <- meta_df$File_ID

# ============================================================================
# 4. BASIC QC SUMMARIES
# ============================================================================

meta_df[[group_col]] <- safe_factor(meta_df[[group_col]], ref = normal_lvl)

cat("Samples matched :", ncol(count_mat), "\n")
cat("Genes in counts :", nrow(count_mat), "\n")
cat("\nGroup distribution:\n")
print(table(meta_df[[group_col]]))

write.table(meta_df, file.path(res_dir, "matched_metadata_used.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================================
# 5. BUILD DGEList, FILTER, NORMALISE
# ============================================================================

dge <- DGEList(
  counts = count_mat,
  genes  = data.frame(gene_name = rownames(count_mat), stringsAsFactors = FALSE)
)
dge$samples$group <- meta_df[[group_col]]

# --- filterByExpr: adaptive minimum expression filter -----------------------
keep <- filterByExpr(dge, group = meta_df[[group_col]])
dge  <- dge[keep, , keep.lib.sizes = FALSE]
cat("\nGenes retained after filterByExpr:", nrow(dge), "\n")

write.table(
  data.frame(gene_name = rownames(dge)),
  file.path(res_dir, "genes_retained_after_filterByExpr.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# --- TMM normalisation ------------------------------------------------------
dge <- calcNormFactors(dge, method = "TMM")

cat("\nNormalisation factors (first 10 samples):\n")
print(head(dge$samples[, c("group", "lib.size", "norm.factors")], 10))

# ============================================================================
# 6. EXPLORATORY QC PLOTS  (before model fitting)
# ============================================================================

# --- 6a. Library size bar chart ---------------------------------------------
lib_df <- data.frame(
  sample      = colnames(dge),
  group       = meta_df[[group_col]],
  lib_size    = dge$samples$lib.size,
  norm_factor = dge$samples$norm.factors
)

p_lib <- ggplot(lib_df, aes(x = reorder(sample, -lib_size), y = lib_size / 1e6, fill = group)) +
  geom_col(width = 0.8) +
  labs(title = "Library sizes by sample", x = NULL, y = "Library size (millions)") +
  theme_pub +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top")
save_pdf_png(p_lib, "QC_library_sizes", width = 12, height = 5)

# --- 6b. MDS plot (edgeR native) -------------------------------------------
mds_png <- file.path(plot_dir, "QC_MDS.png")
mds_pdf <- file.path(plot_dir, "QC_MDS.pdf")

group_fac   <- meta_df[[group_col]]
group_cols  <- setNames(
  brewer.pal(max(3, nlevels(group_fac)), "Set1")[seq_len(nlevels(group_fac))],
  levels(group_fac)
)
sample_cols <- group_cols[as.character(group_fac)]

png(mds_png, width = 2400, height = 1800, res = 300)
plotMDS(dge, labels = FALSE, col = sample_cols, pch = 19,
        main = "MDS plot (leading logFC dimensions)")
legend("topright", legend = names(group_cols), col = group_cols, pch = 19, bty = "n")
dev.off()

pdf(mds_pdf, width = 8, height = 6)
plotMDS(dge, labels = FALSE, col = sample_cols, pch = 19,
        main = "MDS plot (leading logFC dimensions)")
legend("topright", legend = names(group_cols), col = group_cols, pch = 19, bty = "n")
dev.off()

# --- 6c. PCA on logCPM (variance-stabilised, batch-corrected if applicable) -
#
# NOTE: For exploratory PCA we use logCPM with prior.count = 2 for better
# variance stabilisation of low-count genes (analogous to DESeq2's VST idea).
# When batch covariates are modelled, we apply limma::removeBatchEffect()
# to the logCPM *for visualisation only* — the DE model itself handles
# batch effects through the design matrix.

logcpm <- cpm(dge, log = TRUE, prior.count = 2)

# Batch-corrected logCPM for exploratory plots only
batch_cols_present <- batch_cols[batch_cols %in% colnames(meta_df)]

if (use_batch && length(batch_cols_present) > 0) {
  cat("\nBatch columns used for exploratory correction:", paste(batch_cols_present, collapse = ", "), "\n")
  for (b in batch_cols_present) {
    cat("\nBatch × Group cross-tab for", b, ":\n")
    print(table(meta_df[[b]], meta_df[[group_col]]))
  }

  # Build batch model matrix for removeBatchEffect
  batch_design <- model.matrix(
    as.formula(paste("~", paste(batch_cols_present, collapse = " + "))),
    data = meta_df
  )
  # removeBatchEffect needs the batch as covariates and preserves the group effect
  logcpm_corrected <- removeBatchEffect(
    logcpm,
    covariates = batch_design[, -1, drop = FALSE],
    design     = model.matrix(~ meta_df[[group_col]])
  )
  cat("Applied removeBatchEffect on logCPM for exploratory PCA/heatmaps.\n")
} else {
  logcpm_corrected <- logcpm
  if (use_batch) {
    cat("\nWARNING: use_batch = TRUE but no batch columns found in metadata.\n")
    cat("Proceeding without batch correction.\n")
  }
}

# Run PCA on the (optionally corrected) logCPM
pca_res <- prcomp(t(logcpm_corrected), center = TRUE, scale. = TRUE)
var_pct <- round(summary(pca_res)$importance[2, ] * 100, 1)

pca_df <- data.frame(
  sample = rownames(pca_res$x),
  PC1    = pca_res$x[, 1],
  PC2    = pca_res$x[, 2],
  PC3    = pca_res$x[, 3],
  group  = meta_df[[group_col]]
)

# PC1 vs PC2
p_pca12 <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(level = 0.95, linewidth = 0.7, show.legend = FALSE) +
  labs(
    title = "PCA — PC1 vs PC2 (logCPM)",
    x     = paste0("PC1 (", var_pct[1], "%)"),
    y     = paste0("PC2 (", var_pct[2], "%)")
  ) +
  theme_pub +
  theme(legend.position = "top")
save_pdf_png(p_pca12, "QC_PCA_PC1_vs_PC2", width = 7, height = 6)

# PC1 vs PC3
p_pca13 <- ggplot(pca_df, aes(PC1, PC3, color = group)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(level = 0.95, linewidth = 0.7, show.legend = FALSE) +
  labs(
    title = "PCA — PC1 vs PC3 (logCPM)",
    x     = paste0("PC1 (", var_pct[1], "%)"),
    y     = paste0("PC3 (", var_pct[3], "%)")
  ) +
  theme_pub +
  theme(legend.position = "top")
save_pdf_png(p_pca13, "QC_PCA_PC1_vs_PC3", width = 7, height = 6)

# Scree plot (top 10 PCs)
n_pcs_show <- min(10, ncol(pca_res$x))
scree_df <- data.frame(
  PC       = factor(paste0("PC", seq_len(n_pcs_show)), levels = paste0("PC", seq_len(n_pcs_show))),
  Variance = var_pct[seq_len(n_pcs_show)]
)

p_scree <- ggplot(scree_df, aes(PC, Variance)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_text(aes(label = paste0(Variance, "%")), vjust = -0.4, size = 3.2) +
  labs(title = "Scree plot — variance explained per PC", x = NULL, y = "% Variance explained") +
  theme_pub
save_pdf_png(p_scree, "QC_PCA_scree_plot", width = 8, height = 5)

# --- 6d. Sample distance heatmap -------------------------------------------
sample_dist <- as.dist(1 - cor(logcpm_corrected, method = "pearson"))
ann_col_qc  <- data.frame(Group = meta_df[[group_col]], row.names = colnames(logcpm_corrected))

save_pheatmap(
  quote(pheatmap(
    as.matrix(sample_dist),
    clustering_distance_rows = sample_dist,
    clustering_distance_cols = sample_dist,
    clustering_method        = "complete",
    annotation_col           = ann_col_qc,
    annotation_row           = ann_col_qc,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = "Sample distance heatmap (1 − Pearson r)"
  )),
  "QC_sample_distance_heatmap",
  width = 9, height = 8
)

# ============================================================================
# 7. DESIGN MATRIX & BATCH REVIEW
# ============================================================================

if (use_batch && length(batch_cols_present) > 0) {
  # Ensure batch columns are factors
  for (b in batch_cols_present) {
    meta_df[[b]] <- as.factor(meta_df[[b]])
  }
  design_terms   <- c(batch_cols_present, group_col)
  design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
  cat("\nDesign formula (with batch):", deparse(design_formula), "\n")
} else {
  design_formula <- as.formula(paste("~", group_col))
  cat("\nDesign formula (no batch):", deparse(design_formula), "\n")
}

design <- model.matrix(design_formula, data = meta_df)
colnames(design) <- make.names(colnames(design))

cat("Design matrix columns:", paste(colnames(design), collapse = ", "), "\n")
cat("Design matrix rank:", qr(design)$rank, "/", ncol(design), "\n")

# Check for rank deficiency (can happen with confounded batch variables)
if (qr(design)$rank < ncol(design)) {
  warning("Design matrix is rank-deficient — some batch levels may be confounded with the group. Consider dropping batch covariates.")
}

# ============================================================================
# 8. DISPERSION ESTIMATION & MODEL FITTING
# ============================================================================

dge <- estimateDisp(dge, design = design, robust = TRUE)

# --- BCV plot ---------------------------------------------------------------
pdf(file.path(plot_dir, "edgeR_BCV_plot.pdf"), width = 7, height = 6)
plotBCV(dge, main = "Biological coefficient of variation (BCV)")
dev.off()
png(file.path(plot_dir, "edgeR_BCV_plot.png"), width = 2100, height = 1800, res = 300)
plotBCV(dge, main = "Biological coefficient of variation (BCV)")
dev.off()

cat("\nCommon dispersion :", dge$common.dispersion, "\n")
cat("Common BCV        :", sqrt(dge$common.dispersion), "\n")

# --- Quasi-likelihood fit ---------------------------------------------------
fit <- glmQLFit(dge, design = design, robust = TRUE)

# --- QL dispersion plot -----------------------------------------------------
pdf(file.path(plot_dir, "edgeR_QLDisp_plot.pdf"), width = 7, height = 6)
plotQLDisp(fit, main = "Quasi-likelihood dispersion")
dev.off()
png(file.path(plot_dir, "edgeR_QLDisp_plot.png"), width = 2100, height = 1800, res = 300)
plotQLDisp(fit, main = "Quasi-likelihood dispersion")
dev.off()

# ============================================================================
# 9. DIFFERENTIAL EXPRESSION TEST — Tumor vs Normal
# ============================================================================

# Identify the Tumor coefficient in the design matrix
coef_name <- grep(
  paste0("^", make.names(group_col), make.names(tumor_lvl)),
  colnames(design), value = TRUE
)
if (length(coef_name) == 0) {
  coef_name <- setdiff(
    grep(paste0("^", make.names(group_col)), colnames(design), value = TRUE),
    "(Intercept)"
  )
}
if (length(coef_name) == 0) stop("Cannot identify the Tumor vs Normal coefficient in the design.")
coef_idx <- match(coef_name[1], colnames(design))

cat("\nTesting coefficient:", coef_name[1], "(index", coef_idx, ")\n")

qlf <- glmQLFTest(fit, coef = coef_idx)

# --- Extract and annotate results -------------------------------------------
tt <- topTags(qlf, n = Inf, sort.by = "PValue")$table

# topTags returns: gene_name (from genes slot), logFC, logCPM, F, PValue, FDR
tt$significant <- with(tt, FDR < fdr_cutoff & abs(logFC) >= lfc_cutoff)
tt$direction   <- ifelse(tt$logFC > 0, "Up in Tumor", "Down in Tumor")

# Save full results
write.table(
  tt,
  file.path(res_dir, "edgeR_all_results_tumor_vs_normal.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Save significant subset
sig_tt <- tt %>% filter(FDR < fdr_cutoff, abs(logFC) >= lfc_cutoff)
write.table(
  sig_tt,
  file.path(res_dir, "edgeR_significant_genes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Summary counts
n_up   <- sum(sig_tt$logFC > 0)
n_down <- sum(sig_tt$logFC < 0)

summary_df <- data.frame(
  metric = c(
    "Genes tested",
    "Genes retained after filtering",
    paste0("Significant (FDR<", fdr_cutoff, " & |log2FC|>=", lfc_cutoff, ")"),
    "Up in Tumor",
    "Down in Tumor"
  ),
  value = c(nrow(tt), nrow(dge), nrow(sig_tt), n_up, n_down)
)
write.table(summary_df, file.path(res_dir, "summary_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n--- DE Summary ---\n")
print(summary_df)

# ============================================================================
# 10. VOLCANO PLOT
# ============================================================================

volc_df <- tt %>%
  mutate(
    minus_log10_FDR = -log10(pmax(FDR, 1e-300)),
    colour = case_when(
      FDR < fdr_cutoff & logFC >=  lfc_cutoff ~ "Up in Tumor",
      FDR < fdr_cutoff & logFC <= -lfc_cutoff ~ "Down in Tumor",
      TRUE ~ "NS"
    )
  )

top_volc <- volc_df %>%
  filter(colour != "NS") %>%
  arrange(FDR) %>%
  slice_head(n = top_n_label_volcano)

p_volcano <- ggplot(volc_df, aes(x = logFC, y = minus_log10_FDR)) +
  geom_point(aes(color = colour), alpha = 0.7, size = 1.5) +
  scale_color_manual(
    values = c("Up in Tumor" = "#D62728", "Down in Tumor" = "#1F77B4", "NS" = "grey70"),
    name   = NULL
  ) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(
    data = top_volc, aes(label = gene_name),
    size = 3, max.overlaps = Inf, box.padding = 0.35,
    point.padding = 0.2, segment.size = 0.2, show.legend = FALSE
  ) +
  labs(
    title = "Volcano plot — Tumor vs Normal (edgeR QLF)",
    x     = expression(log[2] ~ Fold ~ Change),
    y     = expression(-log[10] ~ FDR)
  ) +
  theme_pub +
  theme(legend.position = "top")
save_pdf_png(p_volcano, "DEG_volcano_plot", width = 8, height = 7)

# ============================================================================
# 11. MA PLOT
# ============================================================================

ma_df <- tt %>%
  mutate(
    colour = case_when(
      FDR < fdr_cutoff & logFC >=  lfc_cutoff ~ "Up in Tumor",
      FDR < fdr_cutoff & logFC <= -lfc_cutoff ~ "Down in Tumor",
      TRUE ~ "NS"
    )
  )

p_ma <- ggplot(ma_df, aes(x = logCPM, y = logFC)) +
  geom_point(aes(color = colour), alpha = 0.7, size = 1.3) +
  scale_color_manual(
    values = c("Up in Tumor" = "#D62728", "Down in Tumor" = "#1F77B4", "NS" = "grey70"),
    name   = NULL
  ) +
  geom_hline(yintercept = c(-lfc_cutoff, 0, lfc_cutoff),
             linetype = c("dashed", "solid", "dashed"),
             color    = c("grey40", "black", "grey40"),
             linewidth = c(0.5, 0.3, 0.5)) +
  labs(
    title = "MA plot — Tumor vs Normal (edgeR QLF)",
    x     = expression(Average ~ log[2] ~ CPM),
    y     = expression(log[2] ~ Fold ~ Change)
  ) +
  theme_pub +
  theme(legend.position = "top")
save_pdf_png(p_ma, "DEG_MA_plot", width = 8, height = 6)

# ============================================================================
# 12. HEATMAPS — Top 25 Up, Top 25 Down, Combined
# ============================================================================

# Use the (optionally batch-corrected) logCPM for heatmap display.
# Row z-scoring is applied by pheatmap (scale = "row").
# Gene selection: FDR ≤ 0.05 (no LFC floor), ranked by logFC.

# Annotation
ann_col <- data.frame(
  Group = meta_df[[group_col]],
  row.names = colnames(logcpm_corrected)
)
ann_colors <- list(
  Group = setNames(
    brewer.pal(max(3, nlevels(ann_col$Group)), "Set1")[seq_len(nlevels(ann_col$Group))],
    levels(ann_col$Group)
  )
)
heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# Significant genes at FDR ≤ 0.05 (no LFC cutoff for heatmap selection)
sig_for_heatmap <- tt %>% filter(FDR <= fdr_cutoff)

# --- Top 25 Upregulated (highest logFC) ------------------------------------
top_up <- sig_for_heatmap %>%
  filter(logFC > 0) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = top_n_directional) %>%
  pull(gene_name)

# Intersect with available rows (safety)
top_up <- intersect(top_up, rownames(logcpm_corrected))

if (length(top_up) >= 2) {
  save_pheatmap(
    quote(pheatmap(
      logcpm_corrected[top_up, , drop = FALSE],
      color              = heat_colors,
      scale              = "row",
      cluster_rows       = TRUE,
      cluster_cols       = TRUE,
      clustering_method  = "ward.D2",
      annotation_col     = ann_col,
      annotation_colors  = ann_colors,
      show_colnames      = FALSE,
      border_color       = NA,
      fontsize_row       = 8,
      main = paste0("Top ", length(top_up), " Upregulated Genes in Tumor\n(FDR \u2264 ", fdr_cutoff, ", ranked by log2FC)")
    )),
    "Heatmap_top25_upregulated"
  )
  cat("\nHeatmap saved: top", length(top_up), "upregulated genes\n")
} else {
  cat("\nSkipping upregulated heatmap: fewer than 2 significant upregulated genes.\n")
}

# --- Top 25 Downregulated (most negative logFC) ----------------------------
top_down <- sig_for_heatmap %>%
  filter(logFC < 0) %>%
  arrange(logFC) %>%
  slice_head(n = top_n_directional) %>%
  pull(gene_name)

top_down <- intersect(top_down, rownames(logcpm_corrected))

if (length(top_down) >= 2) {
  save_pheatmap(
    quote(pheatmap(
      logcpm_corrected[top_down, , drop = FALSE],
      color              = heat_colors,
      scale              = "row",
      cluster_rows       = TRUE,
      cluster_cols       = TRUE,
      clustering_method  = "ward.D2",
      annotation_col     = ann_col,
      annotation_colors  = ann_colors,
      show_colnames      = FALSE,
      border_color       = NA,
      fontsize_row       = 8,
      main = paste0("Top ", length(top_down), " Downregulated Genes in Tumor\n(FDR \u2264 ", fdr_cutoff, ", ranked by log2FC)")
    )),
    "Heatmap_top25_downregulated"
  )
  cat("Heatmap saved: top", length(top_down), "downregulated genes\n")
} else {
  cat("Skipping downregulated heatmap: fewer than 2 significant downregulated genes.\n")
}

# --- Combined: Top 25 Up + Top 25 Down ------------------------------------
top_combined <- c(top_up, top_down)

if (length(top_combined) >= 2) {
  save_pheatmap(
    quote(pheatmap(
      logcpm_corrected[top_combined, , drop = FALSE],
      color              = heat_colors,
      scale              = "row",
      cluster_rows       = TRUE,
      cluster_cols       = TRUE,
      clustering_method  = "ward.D2",
      annotation_col     = ann_col,
      annotation_colors  = ann_colors,
      show_colnames      = FALSE,
      border_color       = NA,
      fontsize_row       = 7,
      main = paste0("Top ", length(top_up), " Up + ", length(top_down),
                     " Down Regulated Genes\n(FDR \u2264 ", fdr_cutoff, ", ranked by log2FC)")
    )),
    "Heatmap_top25_up_and_down_combined",
    width = 10, height = 13
  )
  cat("Heatmap saved: combined", length(top_combined), "genes\n")
}

# ============================================================================
# 13. EXPRESSION BOXPLOTS — Top DE genes
# ============================================================================

top_gene_list <- tt %>%
  filter(FDR < fdr_cutoff) %>%
  arrange(FDR) %>%
  slice_head(n = min(top_n_boxplot, nrow(.))) %>%
  pull(gene_name)

if (length(top_gene_list) >= 1) {
  expr_long <- as.data.frame(t(logcpm_corrected[top_gene_list, , drop = FALSE])) %>%
    rownames_to_column("sample") %>%
    mutate(Group = meta_df[[group_col]]) %>%
    pivot_longer(
      cols      = all_of(top_gene_list),
      names_to  = "gene_name",
      values_to = "logCPM"
    )

  p_expr <- ggplot(expr_long, aes(x = Group, y = logCPM, fill = Group)) +
    geom_boxplot(outlier.size = 0.4, width = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 0.6) +
    facet_wrap(~ gene_name, scales = "free_y", ncol = 3) +
    labs(title = "Top DE genes — expression distribution", x = NULL, y = "logCPM") +
    theme_pub +
    theme(legend.position = "none", strip.text = element_text(face = "bold", size = 8))
  save_pdf_png(p_expr, "DEG_top_gene_boxplots", width = 12, height = 10)
}

# ============================================================================
# 14. EXPORT NORMALISED COUNTS
# ============================================================================

norm_cpm    <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
logcpm_out  <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)

write.table(
  data.frame(gene_name = rownames(norm_cpm), norm_cpm, check.names = FALSE),
  file.path(res_dir, "normalized_counts_cpm.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(gene_name = rownames(logcpm_out), logcpm_out, check.names = FALSE),
  file.path(res_dir, "logCPM_matrix.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ============================================================================
# 15. GENE-LEVEL QC — Mean-variance relationship
# ============================================================================

gene_qc <- data.frame(
  mean_logCPM = rowMeans(logcpm_out),
  sd_logCPM   = rowSds(logcpm_out)
)

p_sd <- ggplot(gene_qc, aes(mean_logCPM, sd_logCPM)) +
  geom_point(alpha = 0.35, size = 0.7, color = "steelblue") +
  geom_smooth(method = "loess", se = FALSE, color = "firebrick", linewidth = 0.8) +
  labs(title = "Gene variability after TMM normalisation", x = "Mean logCPM", y = "SD of logCPM") +
  theme_pub
save_pdf_png(p_sd, "QC_gene_variability", width = 7, height = 6)

# ============================================================================
# 16. SESSION INFO & CLEANUP
# ============================================================================

writeLines(capture.output(sessionInfo()), con = file.path(outdir, "sessionInfo.txt"))

cat("\n================================================================\n")
cat("edgeR DGEA complete.\n")
cat("Batch adjustment:", use_batch, "\n")
cat("Results in:", outdir, "\n")
cat("================================================================\n")
sink()

message("edgeR pipeline finished. Output: ", outdir)
