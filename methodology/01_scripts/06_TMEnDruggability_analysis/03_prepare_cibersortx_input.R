# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  CIBERSORTx Matrix Preparation Pipeline (OSCC Bulk RNA-seq)                 ║
# ║  Script: 03_prepare_cibersortx_input.R                                      ║
# ║  Purpose: Format logCPM matrix into linear-space CIBERSORTx-ready format    ║
# ║           for remote high-performance deconvolution on Stanford servers.    ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

cat("\n════════════════════════════════════════════════════════════════\n")
cat("  CIBERSORTx Matrix Preparation\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("════════════════════════════════════════════════════════════════\n\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})

base_dir   <- "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project"
meta_file  <- file.path(base_dir, "05_results/edgeR_withBatch/results/matched_metadata_used.tsv")
expr_file  <- file.path(base_dir, "05_results/edgeR_withBatch/results/logCPM_matrix.tsv")
out_dir    <- file.path(base_dir, "05_results/localisation_outputs/tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load Metadata & Filter Tumors
meta <- read_tsv(meta_file, show_col_types = FALSE)
tumor_samples <- meta %>%
  filter(Tissue_Type == "Tumor") %>%
  pull(File_ID)

cat(sprintf("  Found %d primary tumor samples with clinical follow-up.\n", length(tumor_samples)))

# 2. Load Expression Matrix
cat("  Loading logCPM expression matrix... ")
expr_df <- read_tsv(expr_file, show_col_types = FALSE)
cat("done.\n")

gene_col <- colnames(expr_df)[1]
colnames(expr_df)[1] <- "GeneSymbol"

# Ensure common samples
common_cols <- intersect(colnames(expr_df), tumor_samples)
cat(sprintf("  Matched %d tumor sample columns in expression matrix.\n", length(common_cols)))

# Subset to GeneSymbol + Tumor samples
ciber_df <- expr_df %>% select(GeneSymbol, all_of(common_cols))

# Convert logCPM (log2(CPM + 1)) to linear space: 2^logCPM - 1 (non-negative)
# CIBERSORTx LM22 assumes linear space input
cat("  Converting log2(CPM+1) expression to non-negative linear space...\n")
mat_linear <- as.matrix(ciber_df[, -1])
mat_linear <- pmax(2^mat_linear - 1, 0)
ciber_linear_df <- bind_cols(GeneSymbol = ciber_df$GeneSymbol, as.data.frame(mat_linear))

# Remove duplicate gene symbols (keep first)
ciber_linear_df <- ciber_linear_df %>%
  filter(!is.na(GeneSymbol) & GeneSymbol != "") %>%
  distinct(GeneSymbol, .keep_all = TRUE)

# Export TSV
out_file <- file.path(out_dir, "cibersortx_input_matrix.tsv")
write_tsv(ciber_linear_df, out_file)

cat("\n════════════════════════════════════════════════════════════════\n")
cat("  CIBERSORTx Input Matrix Ready!\n")
cat(sprintf("  File: %s\n", out_file))
cat(sprintf("  Dimensions: %d genes × %d tumor samples\n", nrow(ciber_linear_df), ncol(ciber_linear_df) - 1))
cat("════════════════════════════════════════════════════════════════\n\n")

cat("  ╔════════════════════════════════════════════════════════════════╗\n")
cat("  ║  NEXT STEPS FOR STANFORD CIBERSORTx (ZERO LOCAL CPU LOAD):    ║\n")
cat("  ╠════════════════════════════════════════════════════════════════╣\n")
cat("  ║  1. Go to: https://cibersortx.stanford.edu                     ║\n")
cat("  ║  2. Log in / Sign up (Free academic access).                   ║\n")
cat("  ║  3. Click 'Run CIBERSORTx' -> 'Impute Cell Fractions'.         ║\n")
cat("  ║  4. Upload Mixture file: cibersortx_input_matrix.tsv           ║\n")
cat("  ║  5. Signature Matrix: LM22 (22 human immune cell phenotypes).  ║\n")
cat("  ║  6. Permutations: 100 or 500 (runs on Stanford cloud).         ║\n")
cat("  ║  7. Download the result CSV/TSV table and save as:             ║\n")
cat("  ║     05_results/localisation_outputs/tables/cibersortx_results.tsv\n")
cat("  ╚════════════════════════════════════════════════════════════════╝\n\n")
