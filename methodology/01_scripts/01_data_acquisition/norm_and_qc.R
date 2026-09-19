# =========================================================
# STEP 3: DESeq2 normalization + QC (PCA + sample distance)
# =========================================================

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

project_root <- "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project"

dds_file <- file.path(project_root, "03_processed/counts/dds_preDESeq2.rds")
meta_file <- file.path(project_root, "03_processed/metadata/OSCC_metadata_matched.tsv")

out_dir_tables <- file.path(project_root, "05_results/tables/normalization_qc")
out_dir_plots  <- file.path(project_root, "05_results/plots/normalization_qc")
out_dir_rds    <- file.path(project_root, "03_processed/normalized")

dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_plots,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_rds,    recursive = TRUE, showWarnings = FALSE)

# Load DESeq2 object and metadata
dds <- readRDS(dds_file)
meta <- read.delim(meta_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# =========================================================
# Attach metadata safely
# =========================================================

meta$File_ID <- as.character(meta$File_ID)

# Match metadata order to DESeq2 object
meta <- meta[match(colnames(dds), meta$File_ID), ]

stopifnot(
  nrow(meta) == ncol(dds),
  identical(colnames(dds), meta$File_ID)
)

# Add metadata columns individually
colData(dds)$File_ID <- meta$File_ID
colData(dds)$Tissue_Type <- factor(meta$Tissue_Type)

# Keep other metadata columns if present
for(col in setdiff(colnames(meta), c("File_ID", "Tissue_Type"))){
  colData(dds)[[col]] <- meta[[col]]
}

# Run DESeq2
dds <- DESeq(
  dds,
  test = "Wald",
  fitType = "parametric"
)

saveRDS(
  dds,
  file.path(out_dir_rds, "dds_DESeq2_fitted.rds")
)

# -----------------------------
# Normalized counts
# -----------------------------
norm_counts <- counts(dds, normalized = TRUE)
norm_df <- data.frame(gene_name = rownames(norm_counts), norm_counts, check.names = FALSE)

write.table(
  norm_df,
  file = file.path(out_dir_tables, "OSCC_normalized_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# =========================================================
# Save size factors
# =========================================================

size_factors <- data.frame(
  File_ID = as.character(colData(dds)$File_ID),
  SizeFactor = sizeFactors(dds),
  Tissue_Type = as.character(colData(dds)$Tissue_Type),
  stringsAsFactors = FALSE
)

write.table(
  size_factors,
  file = file.path(out_dir_tables, "OSCC_size_factors.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# Variance-stabilizing transformation for QC
# -----------------------------
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

saveRDS(vsd, file.path(out_dir_rds, "OSCC_vst_blind.rds"))

# -----------------------------
# PCA on top variable genes
# -----------------------------
gene_var <- apply(vsd_mat, 1, var)
top_n <- min(1000, length(gene_var))
top_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(top_n)]

pca <- prcomp(t(vsd_mat[top_genes, ]), center = TRUE, scale. = FALSE)
var_exp <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

pca_df <- data.frame(
  Sample = colData(dds)$File_ID,
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Tissue_Type = as.character(colData(dds)$Tissue_Type),
  stringsAsFactors = FALSE
)

pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = Tissue_Type)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    title = "PCA of VST-transformed OSCC RNA-seq data",
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
    color = "Tissue Type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir_plots, "PCA_VST_by_TissueType.pdf"),
  pca_plot, width = 8, height = 6, device = cairo_pdf
)
ggsave(
  file.path(out_dir_plots, "PCA_VST_by_TissueType.png"),
  pca_plot, width = 8, height = 6, dpi = 300
)

# -----------------------------
# Sample distance heatmap (clean version)
# -----------------------------
sample_dist <- dist(t(vsd_mat[top_genes, ]))
sample_dist_mat <- as.matrix(sample_dist)

anno <- data.frame(Tissue_Type = dds$Tissue_Type)
rownames(anno) <- colnames(dds)

# Save sample order for reference
sample_order <- data.frame(
  File_ID = colnames(dds),
  Tissue_Type = as.character(dds$Tissue_Type),
  stringsAsFactors = FALSE
)

write.table(
  sample_order,
  file = file.path(out_dir_tables, "OSCC_sample_order_for_heatmap.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf(file.path(out_dir_plots, "SampleDistanceHeatmap_clean.pdf"), width = 10, height = 9)
pheatmap(
  sample_dist_mat,
  annotation_col = anno,
  annotation_row = anno,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  fontsize = 8,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Sample-to-sample distances (VST, top variable genes)"
)
dev.off()

png(file.path(out_dir_plots, "SampleDistanceHeatmap_clean.png"), width = 2800, height = 2600, res = 300)
pheatmap(
  sample_dist_mat,
  annotation_col = anno,
  annotation_row = anno,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  fontsize = 8,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Sample-to-sample distances (VST, top variable genes)"
)
dev.off()

# -----------------------------
# Save PCA coordinates for record
# -----------------------------
write.table(
  pca_df,
  file = file.path(out_dir_tables, "OSCC_PCA_coordinates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------------
# Summary report
# -----------------------------
report <- c(
  "NORMALIZATION / QC SUMMARY",
  "--------------------------",
  paste0("Number of samples: ", ncol(dds)),
  paste0("Number of genes after filtering: ", nrow(dds)),
  paste0("Median size factor: ", round(median(sizeFactors(dds)), 4)),
  paste0("Size factor range: ", round(min(sizeFactors(dds)), 4), " to ", round(max(sizeFactors(dds)), 4)),
  paste0("PC1 variance explained: ", round(var_exp[1], 2), "%"),
  paste0("PC2 variance explained: ", round(var_exp[2], 2), "%")
)

writeLines(report, file.path(out_dir_tables, "normalization_qc_summary.txt"))

cat(paste(report, collapse = "\n"))