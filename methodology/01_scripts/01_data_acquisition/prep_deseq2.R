# =========================================================
# STEP 2: Low-count filtering + DESeq2 object preparation
# =========================================================

library(DESeq2)
library(edgeR)

project_root <- "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project"

counts_file <- file.path(project_root, "03_processed/counts/OSCC_counts_matched.tsv")
meta_file   <- file.path(project_root, "03_processed/metadata/OSCC_metadata_matched.tsv")

out_filtered_counts <- file.path(project_root, "03_processed/counts/OSCC_counts_filtered.tsv")
out_deseq_rds       <- file.path(project_root, "03_processed/counts/dds_preDESeq2.rds")
out_filter_log      <- file.path(project_root, "09_logs/gene_filtering_log.txt")

dir.create(dirname(out_filtered_counts), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_deseq_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_filter_log), recursive = TRUE, showWarnings = FALSE)

# Load matched data
counts <- read.delim(counts_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
meta   <- read.delim(meta_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# Ensure column names are as expected
stopifnot("gene_name" %in% colnames(counts))
stopifnot("File_ID" %in% colnames(meta))
stopifnot("Tissue_Type" %in% colnames(meta))

# Set row names and create matrix
gene_names <- counts$gene_name
count_mat <- as.matrix(counts[, -1, drop = FALSE])
rownames(count_mat) <- gene_names

# Make sure all count values are numeric/integer-like
storage.mode(count_mat) <- "numeric"

# Use Tissue_Type as the grouping variable
meta$Tissue_Type <- factor(meta$Tissue_Type)

# Check group balance
group_table <- table(meta$Tissue_Type)

# Low-expression filtering using edgeR
# This keeps genes with sufficient expression in at least a few samples
keep <- filterByExpr(count_mat, group = meta$Tissue_Type)

filtered_mat <- count_mat[keep, , drop = FALSE]

# Save filtered counts
filtered_df <- data.frame(gene_name = rownames(filtered_mat), filtered_mat, check.names = FALSE)
write.table(filtered_df, out_filtered_counts, sep = "\t", quote = FALSE, row.names = FALSE)

# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(filtered_mat),
  colData   = meta,
  design    = ~ Tissue_Type
)

# Optional: remove genes with all zeros after filtering
dds <- dds[rowSums(counts(dds)) > 0, ]

# Save DESeq2 object
saveRDS(dds, out_deseq_rds)

# Log summary
log_lines <- c(
  paste0("Original genes: ", nrow(count_mat)),
  paste0("Genes retained after filtering: ", nrow(filtered_mat)),
  paste0("Genes removed: ", nrow(count_mat) - nrow(filtered_mat)),
  paste0("Samples: ", ncol(count_mat)),
  paste0("Group counts: ", paste(names(group_table), group_table, sep = "=", collapse = "; "))
)

writeLines(log_lines, out_filter_log)

cat("\n--- GENE FILTERING SUMMARY ---\n")
cat(paste(log_lines, collapse = "\n"))
cat("\n")