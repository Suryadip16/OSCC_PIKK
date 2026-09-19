# =========================================================
# STEP 3: Batch effect diagnostics for TCGA OSCC
# =========================================================

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(readr)

project_root <- "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project"

counts_file <- file.path(project_root, "03_processed/counts/OSCC_counts_matched_noUnspecifiedSites.tsv")
meta_file   <- file.path(project_root, "03_processed/metadata/OSCC_metadata_matched_noUnspecifiedSites.tsv")

out_dir_plots <- file.path(project_root, "05_results/plots/batch_effects/_noUnspecifiedSites")
out_dir_tables <- file.path(project_root, "05_results/tables/batch_effects/_noUnspecifiedSites")
out_dir_logs   <- file.path(project_root, "09_logs/_noUnspecifiedSites")

dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_logs, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load data
# -----------------------------
counts <- read.delim(counts_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
meta   <- read.delim(meta_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# -----------------------------
# Validate expected columns
# -----------------------------
required_meta <- c("File_ID", "Condition", "Plate", "TSS")
missing_meta <- setdiff(required_meta, colnames(meta))
if (length(missing_meta) > 0) {
  stop(paste("Missing metadata columns:", paste(missing_meta, collapse = ", ")))
}

stopifnot("gene_name" %in% colnames(counts))

# -----------------------------
# Align and prepare count matrix
# -----------------------------
meta$File_ID     <- as.character(meta$File_ID)
meta$Condition <- as.factor(meta$Condition)
meta$Plate       <- as.factor(meta$Plate)
meta$TSS         <- as.factor(meta$TSS)

gene_names <- counts$gene_name
count_mat <- as.matrix(counts[, -1, drop = FALSE])
rownames(count_mat) <- gene_names
storage.mode(count_mat) <- "numeric"

# Ensure column order matches metadata
count_mat <- count_mat[, meta$File_ID, drop = FALSE]
stopifnot(identical(colnames(count_mat), meta$File_ID))

# -----------------------------
# Save batch distribution tables
# -----------------------------
tissue_plate <- table(meta$Condition, meta$Plate)
tissue_tss   <- table(meta$Condition, meta$TSS)

write.table(as.data.frame.matrix(tissue_plate),
            file = file.path(out_dir_tables, "TissueType_by_Plate.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

write.table(as.data.frame.matrix(tissue_tss),
            file = file.path(out_dir_tables, "TissueType_by_TSS.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# Also save sample-level annotation
sample_annot <- meta[, c("File_ID", "Condition", "Plate", "TSS")]
write.table(sample_annot,
            file = file.path(out_dir_tables, "OSCC_sample_annotation_for_batch_QC.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# Build DESeq2 object for VST/PCA only
# Use a biological design first, not batch-adjusted yet
# -----------------------------
dds0 <- DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData   = meta,
  design    = ~ Condition
)

# Filter out zero-count genes for transformation
dds0 <- dds0[rowSums(counts(dds0)) > 0, ]

# Variance stabilizing transformation
vsd <- vst(dds0, blind = TRUE)
vsd_mat <- assay(vsd)

# Select most variable genes for PCA
gene_var <- apply(vsd_mat, 1, var)
top_n <- min(1000, length(gene_var))
top_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(top_n)]

pca <- prcomp(t(vsd_mat[top_genes, ]), center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Condition = meta$Condition[match(rownames(pca$x), meta$File_ID)],
  Plate = meta$Plate[match(rownames(pca$x), meta$File_ID)],
  TSS = meta$TSS[match(rownames(pca$x), meta$File_ID)]
)

var_explained <- (pca$sdev^2) / sum(pca$sdev^2) * 100
pc1_lab <- paste0("PC1 (", round(var_explained[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_explained[2], 1), "%)")

# -----------------------------
# Publication-quality plotting theme
# -----------------------------
theme_pub <- theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

# PCA by Condition
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "PCA of VST-transformed counts", x = pc1_lab, y = pc2_lab, color = "Tissue Type") +
  theme_pub

# PCA by Plate
p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Plate)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "PCA colored by Plate", x = pc1_lab, y = pc2_lab, color = "Plate") +
  theme_pub

# PCA by TSS
p3 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = TSS)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "PCA colored by TSS", x = pc1_lab, y = pc2_lab, color = "TSS") +
  theme_pub

# Save plots
ggsave(file.path(out_dir_plots, "PCA_by_TissueType.pdf"), p1, width = 8, height = 6, device = cairo_pdf)
ggsave(file.path(out_dir_plots, "PCA_by_TissueType.png"), p1, width = 8, height = 6, dpi = 300)

ggsave(file.path(out_dir_plots, "PCA_by_Plate.pdf"), p2, width = 9, height = 6, device = cairo_pdf)
ggsave(file.path(out_dir_plots, "PCA_by_Plate.png"), p2, width = 9, height = 6, dpi = 300)

ggsave(file.path(out_dir_plots, "PCA_by_TSS.pdf"), p3, width = 10, height = 6, device = cairo_pdf)
ggsave(file.path(out_dir_plots, "PCA_by_TSS.png"), p3, width = 10, height = 6, dpi = 300)

# -----------------------------
# Batch summary statistics
# -----------------------------
batch_summary <- meta %>%
  group_by(Condition) %>%
  summarise(
    n_samples = n(),
    n_plate_levels = n_distinct(Plate),
    n_tss_levels = n_distinct(TSS),
    .groups = "drop"
  )

write.table(batch_summary,
            file = file.path(out_dir_tables, "batch_summary_by_tissue.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Simple diagnostic text report
report_lines <- c(
  "BATCH EFFECT DIAGNOSTIC SUMMARY",
  "--------------------------------",
  paste0("Total samples: ", nrow(meta)),
  paste0("Total genes in matrix: ", nrow(count_mat)),
  paste0("Condition levels: ", paste(levels(meta$Condition), collapse = ", ")),
  paste0("Number of Plate levels: ", nlevels(meta$Plate)),
  paste0("Number of TSS levels: ", nlevels(meta$TSS)),
  "",
  "Condition x Plate table:",
  capture.output(print(tissue_plate)),
  "",
  "Condition x TSS table:",
  capture.output(print(tissue_tss))
)

writeLines(report_lines, file.path(out_dir_logs, "batch_effect_diagnostic_report.txt"))

cat("\nBatch effect diagnostic files written to:\n")
cat(out_dir_plots, "\n")
cat(out_dir_tables, "\n")
cat(file.path(out_dir_logs, "batch_effect_diagnostic_report.txt"), "\n")