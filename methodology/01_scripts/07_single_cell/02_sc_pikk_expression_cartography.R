# ============================================================================
# PIKK Single-Cell Integration — Script 02: Expression Cartography
# Script:  02_sc_pikk_expression_cartography.R
# Purpose: Map all 7 Diamond Panel genes onto the GSE103322 single-cell atlas.
#          Generate UMAP feature plots, dot plots, violin plots, and
#          compartmentalisation statistics (malignant vs immune vs stromal).
#
# Input:   05_results/single_cell_outputs/rds/oscc_gse103322_seurat.rds
#
# Output:  Publication-quality figure panels (PDF + PNG, 300 dpi)
#          Summary statistics tables (TSV)
#
# Module:  1 — PIKK Diamond Panel Cellular Cartography
# ============================================================================

cat("\n================================================================\n")
cat("  Module 1: PIKK Expression Cartography (Single-Cell)\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(scales)
})

set.seed(42)

# ---- Paths ------------------------------------------------------------------
base_dir <- "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project"
out_dir  <- file.path(base_dir, "05_results/single_cell_outputs")
rds_dir  <- file.path(out_dir, "rds")
tbl_dir  <- file.path(out_dir, "tables")
fig_dir  <- file.path(out_dir, "figures")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Helper: dual save ------------------------------------------------------
save_dual <- function(plot_obj, filename, w = 10, h = 8, dpi = 300) {
  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), plot_obj,
         width = w, height = h, device = "pdf")
  ggsave(file.path(fig_dir, paste0(filename, ".png")), plot_obj,
         width = w, height = h, dpi = dpi, device = "png")
  cat("       Saved:", filename, "(PDF + PNG)\n")
}

# ---- Publication theme -------------------------------------------------------
theme_sc <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      plot.subtitle    = element_text(size = base_size, hjust = 0.5, color = "grey40"),
      axis.title       = element_text(face = "bold"),
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold", size = base_size),
      plot.margin      = margin(10, 10, 10, 10)
    )
}

# ---- 1. Load Seurat Object --------------------------------------------------
cat("[1/7] Loading Seurat object ...\n")
seu <- readRDS(file.path(rds_dir, "oscc_gse103322_seurat.rds"))
cat("       Loaded:", ncol(seu), "cells ×", nrow(seu), "genes\n")

# Standardise metadata column names to clean, syntactic R names
# (Removes spaces and parentheses which break Seurat DimPlot/LabelClusters and ggplot)
if (!"CellType_Major" %in% colnames(seu@meta.data)) {
  if ("Celltype (major-lineage)" %in% colnames(seu@meta.data)) {
    seu$CellType_Major <- seu@meta.data[["Celltype (major-lineage)"]]
  } else if ("Celltype..major.lineage." %in% colnames(seu@meta.data)) {
    seu$CellType_Major <- seu@meta.data[["Celltype..major.lineage."]]
  }
}
if (!"CellType_Malignancy" %in% colnames(seu@meta.data)) {
  if ("Celltype (malignancy)" %in% colnames(seu@meta.data)) {
    seu$CellType_Malignancy <- seu@meta.data[["Celltype (malignancy)"]]
  } else if ("Celltype..malignancy." %in% colnames(seu@meta.data)) {
    seu$CellType_Malignancy <- seu@meta.data[["Celltype..malignancy."]]
  }
}
if (!"CellType_Minor" %in% colnames(seu@meta.data)) {
  if ("Celltype (minor-lineage)" %in% colnames(seu@meta.data)) {
    seu$CellType_Minor <- seu@meta.data[["Celltype (minor-lineage)"]]
  } else if ("Celltype..minor.lineage." %in% colnames(seu@meta.data)) {
    seu$CellType_Minor <- seu@meta.data[["Celltype..minor.lineage."]]
  }
}

# Diamond Panel definition
diamond_panel <- c("PLK1", "CDK2", "TOPBP1", "RAD51", "FANCI", "KAT2B", "DEPTOR")
diamond_present <- diamond_panel[diamond_panel %in% rownames(seu)]
cat("       Diamond Panel genes found:", length(diamond_present), "/", length(diamond_panel), "\n")
cat("       Genes:", paste(diamond_present, collapse = ", "), "\n")

if (length(diamond_present) == 0) {
  stop("No Diamond Panel genes found in the Seurat object.")
}

# Gene annotations for labelling
gene_info <- data.frame(
  gene      = c("PLK1",   "CDK2",   "TOPBP1", "RAD51",  "FANCI",  "KAT2B",  "DEPTOR"),
  branch    = c("ATR",    "ATM",    "ATR",    "ATR",    "ATR",    "TRRAP",  "mTOR"),
  direction = c("Up",     "Up",     "Up",     "Up",     "Up",     "Down",   "Down"),
  role      = c("Mitotic override",
                "G1/S + HR phosphorylation",
                "ATR-activating scaffold",
                "Central HR recombinase",
                "FA/ICL repair",
                "Histone acetyltransferase",
                "mTOR endogenous brake"),
  stringsAsFactors = FALSE
) %>% filter(gene %in% diamond_present)

# ---- 2. Cell-Type Overview UMAP (Reference Panel) ---------------------------
cat("\n[2/7] Generating cell-type reference UMAP ...\n")

# Define consistent colour palette for cell types
major_lineages <- sort(unique(seu@meta.data$CellType_Major))
n_types <- length(major_lineages)

# Curated palette — distinguish malignant, immune, stromal
ct_colours <- c(
  "Malignant"      = "#E63946",   # Vivid red
  "CD8T"           = "#457B9D",   # Steel blue
  "CD8Tex"         = "#1D3557",   # Dark navy
  "CD4Tconv"       = "#A8DADC",   # Light teal
  "Mono/Macro"     = "#F4A261",   # Warm orange
  "Mast"           = "#E9C46A",   # Gold
  "Plasma"         = "#2A9D8F",   # Teal-green
  "Fibroblasts"    = "#264653",   # Dark green-grey
  "Myofibroblasts" = "#606C38",   # Olive
  "Endothelial"    = "#BC6C25",   # Brown
  "Myocyte"        = "#DDA15E"    # Tan
)
# Fallback for any unexpected cell types
for (ct in major_lineages) {
  if (!ct %in% names(ct_colours)) {
    ct_colours[ct] <- "grey60"
  }
}

p_overview <- DimPlot(
  seu,
  reduction  = "umap",
  group.by   = "CellType_Major",
  cols       = ct_colours,
  pt.size    = 0.4,
  label      = TRUE,
  label.size = 3.5,
  repel      = TRUE
) +
  labs(
    title    = "GSE103322 — OSCC Single-Cell Landscape",
    subtitle = paste0("n = ", ncol(seu), " cells | Puram et al., Cell 2017 | TISCH2 annotations"),
    color    = "Cell Type"
  ) +
  theme_sc() +
  NoAxes()

save_dual(p_overview, "01_umap_celltype_overview", w = 10, h = 8)

# ---- 3. UMAP Feature Plots — Diamond Panel ----------------------------------
cat("\n[3/7] Generating UMAP feature plots for Diamond Panel genes (individual + composite) ...\n")

feature_plots <- lapply(diamond_present, function(g) {
  info <- gene_info %>% filter(gene == g)
  direction_arrow <- ifelse(info$direction == "Up", "↑", "↓")
  
  # Individual standalone feature plot
  p_feat <- FeaturePlot(
    seu,
    features   = g,
    reduction  = "umap",
    pt.size    = 0.5,
    order      = TRUE,   # Plot expressing cells on top
    min.cutoff = "q5",
    max.cutoff = "q95",
    cols       = viridis::magma(50)
  ) +
    labs(
      title    = paste0(g, " Expression (", info$branch, " Axis ", direction_arrow, ")"),
      subtitle = paste0(info$role, " | GSE103322 OSCC atlas")
    ) +
    theme_sc(base_size = 11) +
    NoAxes() +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"))
  
  # Save individual feature plot (PDF + PNG)
  save_dual(p_feat, paste0("02_umap_feature_", g), w = 7, h = 6)
  
  # Compact version for composite multi-panel grid
  p_grid_feat <- p_feat +
    labs(
      title    = paste0(g, " (", info$branch, " ", direction_arrow, ")"),
      subtitle = info$role
    ) +
    theme_sc(base_size = 9) +
    NoAxes() +
    theme(plot.subtitle = element_text(size = 8, color = "grey50"))
  
  return(p_grid_feat)
})
names(feature_plots) <- diamond_present

# Composite multi-panel figure
n_genes <- length(diamond_present)
ncol_grid <- min(4, n_genes)
nrow_grid <- ceiling(n_genes / ncol_grid)

p_features_grid <- wrap_plots(feature_plots, ncol = ncol_grid) +
  plot_annotation(
    title    = "PIKK Diamond Panel — Single-Cell Expression Cartography",
    subtitle = "Colour: normalised expression (magma) | Cells ordered by expression level",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
    )
  )

save_dual(p_features_grid, "02_umap_pikk_feature_plots", w = 16, h = 4 * nrow_grid + 1)

# ---- 4. Dot Plot — Expression by Cell Type -----------------------------------
cat("\n[4/7] Generating dot plot (expression × % detected by cell type) ...\n")

# Set identity to major lineage for Seurat plotting functions
Idents(seu) <- "CellType_Major"

p_dotplot <- DotPlot(
  seu,
  features  = diamond_present,
  cols      = c("lightgrey", "#E63946"),
  dot.scale = 8
) +
  labs(
    title    = "PIKK Diamond Panel — Cell-Type Expression Profile",
    subtitle = "Dot size: % cells expressing | Colour: scaled mean expression",
    x        = "Gene",
    y        = "Cell Type"
  ) +
  theme_sc() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 10)
  ) +
  coord_flip()

save_dual(p_dotplot, "03_dotplot_pikk_by_celltype", w = 10, h = 7)

# ---- 5. Box Plots — Compartment Expression (Malignant vs Immune vs Stromal) --
cat("\n[5/7] Generating compartmentalised box plots (individual + composite) ...\n")

# Extract normalised expression for Diamond Panel genes + metadata
expr_data <- FetchData(seu, vars = c(diamond_present, "CellType_Malignancy", "CellType_Major"))

compartment_colours <- c(
  "Malignant cells" = "#E63946",
  "Immune cells"    = "#457B9D",
  "Stromal cells"   = "#2A9D8F"
)

box_plots <- lapply(diamond_present, function(g) {
  info <- gene_info %>% filter(gene == g)
  direction_arrow <- ifelse(info$direction == "Up", "↑", "↓")
  
  plot_df <- data.frame(
    Expression  = expr_data[[g]],
    Compartment = factor(expr_data$CellType_Malignancy,
                         levels = c("Malignant cells", "Immune cells", "Stromal cells"))
  )
  
  pct_pos <- round(sum(plot_df$Expression > 0) / nrow(plot_df) * 100, 1)
  
  # Individual box plot (standalone, easily interpretable)
  p_box <- ggplot(plot_df, aes(x = Compartment, y = Expression, fill = Compartment)) +
    geom_boxplot(
      width          = 0.45,
      outlier.size   = 0.8,
      outlier.alpha  = 0.25,
      outlier.colour = "grey35",
      alpha          = 0.85
    ) +
    stat_summary(
      fun    = mean,
      geom   = "point",
      shape  = 23,          # White diamond for mean expression
      size   = 3,
      fill   = "white",
      color  = "black",
      stroke = 0.8
    ) +
    scale_fill_manual(values = compartment_colours, guide = "none") +
    labs(
      title    = paste0(g, " Expression by Tumour Compartment"),
      subtitle = paste0(info$branch, " branch (", direction_arrow, ") | ", info$role,
                        "\nBox: median ± IQR | Diamond (◇): mean | ", pct_pos, "% cells expressing"),
      x        = NULL,
      y        = "Normalised Expression"
    ) +
    theme_sc(base_size = 11) +
    theme(
      axis.text.x   = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )
  
  # Save individual box plot (PDF + PNG)
  save_dual(p_box, paste0("04_boxplot_compartment_", g), w = 7, h = 6)
  
  # Return compact version for composite grid
  p_grid_box <- p_box +
    labs(
      title    = paste0(g, " (", info$branch, " ", direction_arrow, ")"),
      subtitle = info$role
    ) +
    theme_sc(base_size = 9) +
    theme(
      axis.text.x   = element_text(face = "bold", size = 8),
      plot.subtitle = element_text(size = 8, color = "grey50")
    )
  
  return(p_grid_box)
})
names(box_plots) <- diamond_present

# Composite multi-panel box plot figure
p_box_grid <- wrap_plots(box_plots, ncol = min(4, n_genes)) +
  plot_annotation(
    title    = "PIKK Diamond Panel — Malignant vs Immune vs Stromal Compartments",
    subtitle = "Box: median ± IQR | White diamond (◇): mean expression | Outliers: points",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
    )
  )

save_dual(p_box_grid, "04_boxplot_pikk_compartment_grid", w = 16, h = 4 * ceiling(n_genes / 4) + 1)

# ---- 6. Quantitative Compartmentalisation Analysis --------------------------
cat("\n[6/7] Running compartmentalisation statistics ...\n")

# Kruskal-Wallis + Dunn's post-hoc for malignant vs immune vs stromal
compartment_stats <- lapply(diamond_present, function(g) {
  
  cat("       ", g, ": ")
  
  # Expression per compartment
  expr_vec      <- expr_data[[g]]
  compartment   <- expr_data[["CellType_Malignancy"]]
  
  # Summary statistics
  summary_df <- data.frame(
    gene        = g,
    compartment = unique(compartment),
    stringsAsFactors = FALSE
  )
  
  summary_df <- expr_data %>%
    group_by(CellType_Malignancy) %>%
    summarise(
      mean_expr     = mean(.data[[g]], na.rm = TRUE),
      median_expr   = median(.data[[g]], na.rm = TRUE),
      pct_nonzero   = sum(.data[[g]] > 0, na.rm = TRUE) / n() * 100,
      n_cells       = n(),
      .groups       = "drop"
    ) %>%
    mutate(gene = g)
  
  # Kruskal-Wallis test
  kw_test <- kruskal.test(expr_vec ~ compartment)
  summary_df$kw_pvalue <- kw_test$p.value
  summary_df$kw_chi2   <- kw_test$statistic
  
  # Pairwise Wilcoxon (Dunn's approximation) for post-hoc
  pw_test <- pairwise.wilcox.test(expr_vec, compartment, p.adjust.method = "BH")
  pw_matrix <- pw_test$p.value
  
  # Extract specific comparisons
  comparisons <- c("Malignant_vs_Immune", "Malignant_vs_Stromal", "Immune_vs_Stromal")
  pw_pvals <- rep(NA, 3)
  
  if (!is.null(pw_matrix)) {
    rn <- rownames(pw_matrix)
    cn <- colnames(pw_matrix)
    
    # Malignant vs Immune
    if ("Malignant cells" %in% rn && "Immune cells" %in% cn) {
      pw_pvals[1] <- pw_matrix["Malignant cells", "Immune cells"]
    } else if ("Immune cells" %in% rn && "Malignant cells" %in% cn) {
      pw_pvals[1] <- pw_matrix["Immune cells", "Malignant cells"]
    }
    
    # Malignant vs Stromal
    if ("Malignant cells" %in% rn && "Stromal cells" %in% cn) {
      pw_pvals[2] <- pw_matrix["Malignant cells", "Stromal cells"]
    } else if ("Stromal cells" %in% rn && "Malignant cells" %in% cn) {
      pw_pvals[2] <- pw_matrix["Stromal cells", "Malignant cells"]
    }
    
    # Immune vs Stromal
    if ("Immune cells" %in% rn && "Stromal cells" %in% cn) {
      pw_pvals[3] <- pw_matrix["Immune cells", "Stromal cells"]
    } else if ("Stromal cells" %in% rn && "Immune cells" %in% cn) {
      pw_pvals[3] <- pw_matrix["Stromal cells", "Immune cells"]
    }
  }
  
  summary_df$pw_malig_vs_immune  <- pw_pvals[1]
  summary_df$pw_malig_vs_stromal <- pw_pvals[2]
  summary_df$pw_immune_vs_stromal <- pw_pvals[3]
  
  cat("KW p =", formatC(kw_test$p.value, format = "e", digits = 2), "\n")
  
  return(summary_df)
}) %>% bind_rows()

# Determine predominant compartment for each gene
compartment_stats <- compartment_stats %>%
  group_by(gene) %>%
  mutate(
    is_predominant = mean_expr == max(mean_expr),
    predominant_compartment = CellType_Malignancy[which.max(mean_expr)]
  ) %>%
  ungroup()

# Save
write_tsv(compartment_stats, file.path(tbl_dir, "pikk_gene_celltype_compartment_stats.tsv"))
cat("       Compartmentalisation statistics saved.\n")

# Summary table for easy interpretation
cat("\n       --- Predominant Compartment Summary ---\n")
predominant_summary <- compartment_stats %>%
  filter(is_predominant) %>%
  select(gene, predominant_compartment, mean_expr, pct_nonzero, kw_pvalue) %>%
  mutate(kw_pvalue = formatC(kw_pvalue, format = "e", digits = 2))
print(as.data.frame(predominant_summary), row.names = FALSE)

# ---- 7. Detailed Cell-Type Expression Summary --------------------------------
cat("\n[7/7] Generating detailed cell-type expression summary ...\n")

Idents(seu) <- "CellType_Major"

celltype_expr_summary <- lapply(diamond_present, function(g) {
  expr_data %>%
    group_by(CellType_Major) %>%
    summarise(
      mean_expr   = mean(.data[[g]], na.rm = TRUE),
      median_expr = median(.data[[g]], na.rm = TRUE),
      pct_nonzero = sum(.data[[g]] > 0, na.rm = TRUE) / n() * 100,
      sd_expr     = sd(.data[[g]], na.rm = TRUE),
      n_cells     = n(),
      .groups     = "drop"
    ) %>%
    mutate(gene = g) %>%
    arrange(desc(mean_expr))
}) %>% bind_rows()

write_tsv(celltype_expr_summary, file.path(tbl_dir, "pikk_gene_celltype_expression_summary.tsv"))
cat("       Detailed expression summary saved.\n")

# ---- Heatmap-style bar chart: % cells expressing per compartment -------------
cat("\n       Generating compartmentalisation bar chart ...\n")

pct_data <- compartment_stats %>%
  select(gene, CellType_Malignancy, pct_nonzero, mean_expr) %>%
  mutate(gene = factor(gene, levels = rev(diamond_present)))

p_pct_bar <- ggplot(pct_data, aes(x = gene, y = pct_nonzero,
                                   fill = CellType_Malignancy)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = compartment_colours, name = "Compartment") +
  coord_flip() +
  labs(
    title    = "PIKK Diamond Panel — % Cells Expressing per Compartment",
    subtitle = "Threshold: expression > 0 in normalised data",
    x        = NULL,
    y        = "% Cells Expressing"
  ) +
  theme_sc() +
  theme(
    axis.text.y = element_text(face = "italic", size = 11),
    legend.position = "bottom"
  )

save_dual(p_pct_bar, "05_barplot_pikk_pct_expressing", w = 9, h = 6)

# ---- Combined Dot Plot: Major Lineage (detailed) ----------------------------
cat("\n       Generating detailed major-lineage dot plot ...\n")

Idents(seu) <- "CellType_Major"

# Reorder cell types: malignant first, then immune, then stromal
lineage_order <- c("Malignant", "CD8T", "CD8Tex", "CD4Tconv", "Mono/Macro",
                   "Mast", "Plasma", "Fibroblasts", "Myofibroblasts",
                   "Endothelial", "Myocyte")
# Keep only those present in the data
lineage_order <- lineage_order[lineage_order %in% levels(Idents(seu))]
remaining <- setdiff(levels(Idents(seu)), lineage_order)
lineage_order <- c(lineage_order, remaining)

Idents(seu) <- factor(as.character(Idents(seu)), levels = lineage_order)

p_dotplot_detailed <- DotPlot(
  seu,
  features  = diamond_present,
  cols      = c("lightgrey", "#B5179E"),
  dot.scale = 8,
  dot.min   = 0.01
) +
  labs(
    title    = "PIKK Diamond Panel — Detailed Cell-Type Resolution",
    subtitle = "11 cell types | Dot size: % expressing | Colour: scaled mean expression",
    x        = "Gene",
    y        = "Cell Type (Major Lineage)"
  ) +
  theme_sc() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11),
    axis.text.y = element_text(size = 10)
  )

save_dual(p_dotplot_detailed, "06_dotplot_pikk_detailed_lineage", w = 10, h = 7)

# ==============================================================================
cat("\n================================================================\n")
cat("  Module 1 complete:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("  Figures generated: 20 (7 individual UMAP feature plots + 1 composite grid,\n")
cat("                      7 individual compartment box plots + 1 composite grid,\n")
cat("                      1 cell-type overview UMAP, 2 dot plots, 1 bar plot)\n")
cat("  Tables generated:  2\n")
cat("  Next: Run 03_sc_ddr_programme_scoring.R\n")
cat("================================================================\n\n")
