# =========================================================
# STEP 1: Load counts + metadata, validate sample matching
# =========================================================

project_root <- "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project"

counts_file <- file.path(project_root, "01_data/raw_counts/OSCC_Coding_Counts_corrected_noUnspecifiedSites.tsv")
meta_file   <- file.path(project_root, "01_data/raw_metadata/TCGA_OSCC_Metadata_Batch_Added_noUnspecifiedSites_Cleaned.tsv")

out_counts   <- file.path(project_root, "03_processed/counts/OSCC_counts_matched_noUnspecifiedSites.tsv")
out_meta     <- file.path(project_root, "03_processed/metadata/OSCC_metadata_matched_noUnspecifiedSites.tsv")
out_log      <- file.path(project_root, "09_logs/sample_matching_log_noUnspecifiedSites.txt")
out_inventory <- file.path(project_root, "01_data/sample_inventory/OSCC_sample_inventory_noUnspecifiedSites.tsv")

# Create output folders if needed
dir.create(dirname(out_counts), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_meta), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_log), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_inventory), recursive = TRUE, showWarnings = FALSE)

# Read data
counts <- read.delim(counts_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
meta   <- read.delim(meta_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# Basic column checks
required_count_col <- "gene_name"
required_meta_cols  <- c("File_ID", "Condition")

if (!required_count_col %in% colnames(counts)) {
  stop("counts file must have a first column named 'gene_name'")
}

missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))
if (length(missing_meta_cols) > 0) {
  stop(paste("metadata file is missing columns:", paste(missing_meta_cols, collapse = ", ")))
}

# Clean up gene names and sample IDs
counts$gene_name <- as.character(counts$gene_name)
meta$`File_ID`    <- as.character(meta$`File_ID`)
meta$Condition    <- as.character(meta$Condition)

# Separate matrix and sample IDs
count_ids <- colnames(counts)[-1]
meta_ids  <- meta$`File_ID`

# Check for duplicate sample IDs in metadata
if (any(duplicated(meta_ids))) {
  dup_meta <- unique(meta_ids[duplicated(meta_ids)])
  stop(paste("Duplicate File_ID(s) found in metadata:", paste(dup_meta, collapse = ", ")))
}

# Find overlap
common_ids <- intersect(count_ids, meta_ids)
only_in_counts <- setdiff(count_ids, meta_ids)
only_in_meta    <- setdiff(meta_ids, count_ids)

# Stop if no overlap
if (length(common_ids) == 0) {
  stop("No matching sample IDs found between counts and metadata.")
}

# Reorder counts columns to match metadata order
meta_matched <- meta[match(common_ids, meta$`File_ID`), , drop = FALSE]
counts_matched <- counts[, c("gene_name", common_ids), drop = FALSE]

# Ensure metadata order matches counts order exactly
stopifnot(identical(colnames(counts_matched)[-1], meta_matched$`File_ID`))

# Convert count columns to numeric safely
counts_matched[, -1] <- lapply(counts_matched[, -1, drop = FALSE], function(x) as.numeric(x))

# Optional: check for NA introduced during numeric conversion
na_counts <- sum(is.na(as.matrix(counts_matched[, -1, drop = FALSE])))

# Make a sample inventory table
sample_inventory <- data.frame(
  File_ID = common_ids,
  Condition = meta_matched$Condition,
  stringsAsFactors = FALSE
)

# Save matched data
write.table(counts_matched, out_counts, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(meta_matched, out_meta, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sample_inventory, out_inventory, sep = "\t", quote = FALSE, row.names = FALSE)

# Write matching log
log_lines <- c(
  paste0("Total samples in counts: ", length(count_ids)),
  paste0("Total samples in metadata: ", length(meta_ids)),
  paste0("Matched samples: ", length(common_ids)),
  paste0("Samples only in counts: ", length(only_in_counts)),
  paste0("Samples only in metadata: ", length(only_in_meta)),
  paste0("NA values introduced during numeric conversion: ", na_counts)
)

writeLines(log_lines, out_log)

# Print summary to console
cat("\n--- SAMPLE MATCHING SUMMARY ---\n")
cat(paste(log_lines, collapse = "\n"))
cat("\n\nFirst few matched samples:\n")
print(head(sample_inventory))

if (length(only_in_counts) > 0) {
  cat("\nSamples present only in counts:\n")
  print(only_in_counts)
}

if (length(only_in_meta) > 0) {
  cat("\nSamples present only in metadata:\n")
  print(only_in_meta)
}

