# =========================================================
# STEP 5: Intersect edgeR DEGs with PIKK gene universe
# =========================================================

library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

project_root <- "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project"

# Input files
edger_file <- file.path(project_root, "05_results/tables/deseq2/OSCC_DESeq2_Tumor_vs_Normal_shrunk.tsv")
vsd_file   <- file.path(project_root, "03_processed/normalized/OSCC_vst_blind.rds")
pikk_file  <- file.path(project_root, "10_references/PIKK_related_gene_universe.csv")

# Output folders
out_dir_tables <- file.path(project_root, "05_results/tables/pikk_intersection")
out_dir_plots  <- file.path(project_root, "05_results/plots/pikk_intersection")
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_plots,  recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------
# Load inputs
# -------------------------------------------------
deg_df <- read.delim(edger_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
pikk_df <- read.csv(pikk_file, stringsAsFactors = FALSE)

# vsd <- readRDS(vsd_file)
# vsd_mat <- assay(vsd)

# -------------------------------------------------
# Standardize gene symbols
# -------------------------------------------------
deg_df <- deg_df %>%
  mutate(
    gene_symbol = toupper(gene_name)
  )

pikk_df <- pikk_df %>%
  mutate(
    GeneSymbol = toupper(GeneSymbol)
  ) %>%
  distinct(PIKK_Group, GeneSymbol, GeneRole, .keep_all = TRUE)

# -------------------------------------------------
# Intersect DEGs with PIKK universe
# -------------------------------------------------
intersect_df <- deg_df %>%
  inner_join(
    pikk_df,
    by = c("gene_symbol" = "GeneSymbol")
  ) %>%
  mutate(
    regulation = case_when(
      !is.na(padj) & padj <= 0.05 & log2FoldChange >=  1 ~ "Up in Tumor",
      !is.na(padj) & padj <= 0.05 & log2FoldChange <= -1 ~ "Down in Tumor",
      TRUE ~ "Not significant"
    ),
    abs_log2FC = abs(log2FoldChange)
  )

# Significant dysregulated PIKK-associated genes only
pikk_sig <- intersect_df %>%
  filter(padj <= 0.05 & abs_log2FC >= 0.5) %>%
  arrange(padj, desc(abs_log2FC))

# -------------------------------------------------
# Save tables
# -------------------------------------------------
write.table(
  intersect_df,
  file = file.path(out_dir_tables, "OSCC_PIKK_extended_DEG_intersection_all.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  pikk_sig,
  file = file.path(out_dir_tables, "OSCC_PIKK_extended_DEG_intersection_significant.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

colnames(pikk_sig)

summary_by_regulation <- pikk_sig %>%
  dplyr::count(PIKK_Group, GeneRole, regulation, name = "n") %>%
  dplyr::arrange(PIKK_Group, desc(n))

write.table(
  summary_by_regulation,
  file = file.path(out_dir_tables, "OSCC_PIKK_extended_intersection_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Simple counts by direction
direction_counts <- pikk_sig %>%
  dplyr::count(regulation, name = "n")

write.table(
  direction_counts,
  file = file.path(out_dir_tables, "OSCC_PIKK_extended_direction_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------------------------------
# Console summary
# -------------------------------------------------
cat("\n========================================\n")
cat("PIKK INTERSECTION SUMMARY\n")
cat("========================================\n")
cat("Total edgeR genes loaded: ", nrow(deg_df), "\n")
cat("Genes in PIKK universe: ", nrow(pikk_df), "\n")
cat("Intersected genes: ", nrow(intersect_df), "\n")
cat("Significant PIKK-associated DEGs: ", nrow(pikk_sig), "\n")
cat("\nDirection counts:\n")
print(direction_counts)

# -------------------------------------------------
# Heatmap of intersected significant PIKK DEGs
# -------------------------------------------------
if (nrow(pikk_sig) >= 2) {
  
  heat_genes <- pikk_sig %>%
    # dplyr::filter(log2FoldChange > 0) %>%
    dplyr::pull(gene_name)
  
  heat_genes <- intersect(heat_genes, rownames(vsd_mat))
  
  heat_mat <- vsd_mat[heat_genes, , drop = FALSE]
  
  # Row-wise z-score scaling
  heat_mat_scaled <- t(scale(t(heat_mat)))
  heat_mat_scaled[is.na(heat_mat_scaled)] <- 0
  
  anno_col <- data.frame(
    Tissue_Type = factor(vsd$Tissue_Type)
  )
  rownames(anno_col) <- colnames(vsd_mat)
  
  heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
  
  pdf(
    file.path(out_dir_plots, "PIKK_extended_significant_DEG_heatmap.pdf"),
    width = 10,
    height = max(6, 0.22 * nrow(heat_mat_scaled) + 3)
  )
  pheatmap(
    heat_mat_scaled,
    color = heat_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = anno_col,
    show_colnames = FALSE,
    border_color = NA,
    fontsize_row = 8,
    clustering_method = "mcquitty",
    main = "Significant PIKK-associated DEGs in OSCC"
  )
  dev.off()
  
  png(
    file.path(out_dir_plots, "PIKK_extended_significant_DEG_heatmap.png"),
    width = 2600,
    height = max(1800, 70 * nrow(heat_mat_scaled)),
    res = 300
  )
  pheatmap(
    heat_mat_scaled,
    color = heat_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    # scale = "row",
    annotation_col = anno_col,
    show_colnames = FALSE,
    border_color = NA,
    fontsize_row = 8,
    clustering_method = "average",
    main = "Significant PIKK-associated DEGs in OSCC"
  )
  dev.off()
  
} else {
  message("Not enough significant intersected PIKK genes for a heatmap.")
}

# -------------------------------------------------
# Optional: grouped barplot for publication use
# -------------------------------------------------
if (nrow(pikk_sig) > 0) {
  bar_df <- pikk_sig %>%
    dplyr::count(PIKK_Group, regulation) %>%
    dplyr::group_by(PIKK_Group) %>%
    mutate(total = sum(n)) %>%
    ungroup()
  
  bar_plot <- ggplot(bar_df, aes(x = PIKK_Group, y = n, fill = regulation)) +
    geom_col(position = "stack", width = 0.75) +
    labs(
      title = "Dysregulated PIKK-associated genes in OSCC",
      x = "PIKK gene family / pathway group",
      y = "Number of genes",
      fill = "Direction"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(
    file.path(out_dir_plots, "PIKK_extended_dysregulated_gene_barplot.pdf"),
    bar_plot,
    width = 10,
    height = 6
  )
  
  ggsave(
    file.path(out_dir_plots, "PIKK_extended_dysregulated_gene_barplot.png"),
    bar_plot,
    width = 10,
    height = 6,
    dpi = 300
  )
}
