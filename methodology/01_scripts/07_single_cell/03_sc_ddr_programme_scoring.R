# ============================================================================
# PIKK Single-Cell Integration — Script 03: DDR Programme Scoring
# Script:  03_sc_ddr_programme_scoring.R
# Purpose: Use AUCell to score a custom PIKK/DDR gene set per cell, revealing
#          intra-tumour heterogeneity in DNA damage response activity.
#          Stratify malignant cells into DDR-high vs DDR-low, perform
#          differential expression, and examine patient-level variation.
#
# Input:   05_results/single_cell_outputs/rds/oscc_gse103322_seurat.rds
#
# Output:  AUCell score tables, UMAP overlays, ridge/density plots,
#          DDR-high vs DDR-low DE tables, patient-level DDR proportions.
#
# Module:  2 — DDR Programme Scoring & Malignant Cell Heterogeneity
# ============================================================================

cat("\n================================================================\n")
cat("  Module 2: DDR Programme Scoring & Malignant Heterogeneity\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(GSEABase)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(scales)
  library(ggridges)
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

# Compartment colours (consistent with Script 02)
compartment_colours <- c(
  "Malignant cells" = "#E63946",
  "Immune cells"    = "#457B9D",
  "Stromal cells"   = "#2A9D8F"
)

ct_colours <- c(
  "Malignant"      = "#E63946",
  "CD8T"           = "#457B9D",
  "CD8Tex"         = "#1D3557",
  "CD4Tconv"       = "#A8DADC",
  "Mono/Macro"     = "#F4A261",
  "Mast"           = "#E9C46A",
  "Plasma"         = "#2A9D8F",
  "Fibroblasts"    = "#264653",
  "Myofibroblasts" = "#606C38",
  "Endothelial"    = "#BC6C25",
  "Myocyte"        = "#DDA15E"
)

# ---- 1. Load Seurat Object --------------------------------------------------
cat("[1/9] Loading Seurat object ...\n")
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

# ---- 2. Define DDR Gene Sets ------------------------------------------------
cat("\n[2/9] Defining PIKK/DDR gene sets for AUCell scoring ...\n")

# === Gene Set A: PIKK-OSCC ATR-Axis Signature (10 genes) ===
# Inclusion: All members are in our 24-hub gene list AND have either
#   significant edgeR DE (|logFC| >= 1.0, FDR < 0.05) OR
#   significant KM survival association (p < 0.05).
# Rationale: This captures the ATR-centric DDR programme identified as the
#   dominant PIKK axis in our bulk RNA-seq analysis (14/22 ATR interactors
#   upregulated, 63.6% — the branch-asymmetry finding).
atr_axis_genes <- c(
  "CHEK1",   # Hub rank 4 | edgeR logFC +1.30 | ATR signal transducer
  "RAD51",   # Hub rank 5 | edgeR logFC +1.37 | KM p=0.0035 | Diamond Panel
  "PLK1",    # Hub rank 6 | edgeR logFC +1.94 | Adj HR=1.54, p=0.004 | Diamond Panel
  "TOPBP1",  # Hub rank 19 | KM p=0.0045 | Diamond Panel
  "CDK2",    # Hub rank 14 | KM p=1e-4, Adj HR=1.60 | Diamond Panel
  "EXO1",    # Hub rank 9 | edgeR logFC +1.72 | DNA end-resection nuclease
  "FANCI",   # Hub rank 16 | edgeR logFC +1.19 | KM p=0.031 | Diamond Panel
  "BRCA1",   # Hub rank 2 | edgeR logFC +1.09 | Eigenvector centrality = 1.000
  "CDC45",   # Hub rank 13 | edgeR logFC +1.61 | Replication initiation helicase
  "E2F1"     # Hub rank 8 | edgeR logFC +1.73 | S-phase transcription factor
)

# === Gene Set B: PIKK-OSCC Study-Validated Hub Interactome (17 genes) ===
# Inclusion criteria (must satisfy >= 1):
#   1. Significantly DE in edgeR (|logFC| >= 1.0, FDR < 0.05), OR
#   2. Top-10 hub gene by composite score, OR
#   3. Diamond Panel member (final biomarker shortlist), OR
#   4. KM survival p < 0.05
# EXCLUDED: ATR, ATM, MTOR, TRRAP, SMG1 (core kinases but NOT DE in our study),
#   CHEK2, H2AFX, PARP1, MSH6, POLD1, RPTOR, RICTOR, RUVBL1, KAT2A, TTI1
#   (low-ranked hubs with |logFC| < 1.0, no survival signal).
pikk_hub_genes <- c(
  # Diamond Panel (7 genes — the final biomarker shortlist)
  "PLK1",    # Hub rank 6 | edgeR +1.94 | Adj HR=1.54 | Best multivariate gene
  "CDK2",    # Hub rank 14 | KM p=1e-4 | Adj HR=1.60 | Strongest KM discriminator
  "TOPBP1",  # Hub rank 19 | KM p=0.0045 | ATR-activating scaffold
  "RAD51",   # Hub rank 5 | edgeR +1.37 | KM p=0.0035 | Central HR recombinase
  "FANCI",   # Hub rank 16 | edgeR +1.19 | KM p=0.031 | FA/ICL repair
  "KAT2B",   # Hub rank 7 | edgeR -1.85 | TRRAP axis HAT (tumour suppressor)
  "DEPTOR",  # Diamond Panel | edgeR -1.78 | mTOR endogenous brake (silenced)
  # Top-10 hub genes not in Diamond Panel
  "PRKDC",   # Hub rank 1 | edgeR +1.02 | Top composite score | Network bottleneck
  "BRCA1",   # Hub rank 2 | edgeR +1.09 | Degree=29, eigenvector=1.000
  "BRCA2",   # Hub rank 3 | edgeR +1.24 | HR co-hub with BRCA1
  "CHEK1",   # Hub rank 4 | edgeR +1.30 | ATR checkpoint kinase
  "EXO1",    # Hub rank 9 | edgeR +1.72 | DNA end-resection nuclease
  "E2F1",    # Hub rank 8 | edgeR +1.73 | Oncogenic S-phase driver
  # Significant DE genes in PIKK interactome (|logFC| >= 1.0)
  "EGFR",    # Hub rank 12 | edgeR +1.96 | KM p=0.0053 | PRKDC axis
  "AURKA",   # edgeR +2.02 | Highest logFC in PIKK interactome | Mitotic kinase
  "AURKB",   # edgeR +1.72 | Aurora kinase B | Chromosome segregation
  "CDC45"    # Hub rank 13 | edgeR +1.61 | Replication origin firing
)

# Filter to genes present in the dataset
all_genes <- rownames(seu)

atr_present  <- atr_axis_genes[atr_axis_genes %in% all_genes]
pikk_present <- pikk_hub_genes[pikk_hub_genes %in% all_genes]

cat("       ATR-Axis Signature:     ", length(atr_present), "/", length(atr_axis_genes),
    " genes present\n")
cat("       Genes:", paste(atr_present, collapse = ", "), "\n")
cat("       PIKK Hub Interactome:   ", length(pikk_present), "/", length(pikk_hub_genes),
    " genes present\n")
cat("       Genes:", paste(pikk_present, collapse = ", "), "\n")

# Construct GeneSet objects for AUCell
gs_atr  <- GeneSet(atr_present, setName = "ATR_DDR_Axis")
gs_pikk <- GeneSet(pikk_present, setName = "PIKK_Hub_Interactome")
geneSets <- GeneSetCollection(list(gs_atr, gs_pikk))

# ---- 3. Run AUCell -----------------------------------------------------------
cat("\n[3/9] Running AUCell scoring ...\n")

# Extract expression matrix (normalised counts) — robust across Seurat v4 and v5
expr_matrix <- tryCatch(
  GetAssayData(seu, layer = "data"),
  error = function(e) GetAssayData(seu, slot = "data")
)
cat("       Expression matrix:", nrow(expr_matrix), "genes ×", ncol(expr_matrix), "cells\n")

# Step 1: Build gene rankings per cell
cat("       Building AUC rankings ...\n")
cells_rankings <- AUCell_buildRankings(
  expr_matrix,
  nCores    = 1,         # Safe single-threaded
  plotStats = FALSE,
  verbose   = FALSE
)

# Step 2: Calculate AUC scores
cat("       Calculating AUC scores ...\n")
cells_AUC <- AUCell_calcAUC(
  geneSets,
  cells_rankings,
  aucMaxRank = ceiling(0.05 * nrow(cells_rankings)),  # Top 5% threshold
  verbose    = FALSE
)

# Extract scores safely by name
auc_mat <- getAUC(cells_AUC)
auc_scores <- data.frame(
  AUC_ATR_DDR  = as.numeric(auc_mat["ATR_DDR_Axis", colnames(seu)]),
  AUC_PIKK_Hub = as.numeric(auc_mat["PIKK_Hub_Interactome", colnames(seu)]),
  row.names    = colnames(seu)
)

cat("       AUC score ranges:\n")
cat("         ATR DDR Axis:          [", round(min(auc_scores$AUC_ATR_DDR), 4), ",",
    round(max(auc_scores$AUC_ATR_DDR), 4), "]\n")
cat("         PIKK Hub Interactome:  [", round(min(auc_scores$AUC_PIKK_Hub), 4), ",",
    round(max(auc_scores$AUC_PIKK_Hub), 4), "]\n")

# Add AUCell scores to Seurat metadata
seu <- AddMetaData(seu, auc_scores)

# ---- 4. UMAP Overlay of DDR Scores ------------------------------------------
cat("\n[4/9] Generating UMAP overlays of AUCell DDR scores ...\n")

# ATR DDR Axis score
p_umap_atr <- FeaturePlot(
  seu,
  features  = "AUC_ATR_DDR",
  reduction = "umap",
  pt.size   = 0.4,
  order     = TRUE,
  cols      = viridis::inferno(50)
) +
  labs(
    title    = "ATR–DDR Axis Activity Score",
    subtitle = paste0(length(atr_present), "-gene signature | AUCell (top 5% ranking threshold)")
  ) +
  theme_sc() +
  NoAxes()

# PIKK Hub Interactome score
p_umap_pikk <- FeaturePlot(
  seu,
  features  = "AUC_PIKK_Hub",
  reduction = "umap",
  pt.size   = 0.4,
  order     = TRUE,
  cols      = viridis::inferno(50)
) +
  labs(
    title    = "PIKK-OSCC Hub Interactome Activity Score",
    subtitle = paste0(length(pikk_present), "-gene study-validated signature | AUCell (top 5% ranking threshold)")
  ) +
  theme_sc() +
  NoAxes()

p_umap_combined <- p_umap_atr + p_umap_pikk +
  plot_annotation(
    title = "PIKK/DDR Programme Activity — Single-Cell Resolution",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

save_dual(p_umap_combined, "07_umap_ddr_aucell_overlay", w = 16, h = 7)

# ---- 5. Ridge Plot — DDR Scores by Cell Type ---------------------------------
cat("\n[5/9] Generating ridge plots of DDR scores by cell type ...\n")

ridge_data <- seu@meta.data %>%
  select(CellType_Major, AUC_ATR_DDR, AUC_PIKK_Hub) %>%
  rename(CellType = CellType_Major)

# Order cell types by median ATR score (malignant should be on top)
ct_order <- ridge_data %>%
  group_by(CellType) %>%
  summarise(med = median(AUC_ATR_DDR, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(CellType)
ridge_data$CellType <- factor(ridge_data$CellType, levels = ct_order)

# Ensure ct_colours covers all cell types
for (ct in unique(ridge_data$CellType)) {
  if (!ct %in% names(ct_colours)) ct_colours[ct] <- "grey60"
}

p_ridge_atr <- ggplot(ridge_data, aes(x = AUC_ATR_DDR, y = CellType, fill = CellType)) +
  geom_density_ridges(
    scale       = 1.5,
    alpha       = 0.7,
    rel_min_height = 0.01,
    quantile_lines = TRUE,
    quantiles      = 2     # Median line
  ) +
  scale_fill_manual(values = ct_colours, guide = "none") +
  labs(
    title    = "ATR–DDR Axis Activity by Cell Type",
    subtitle = "Ridge density | Vertical line: median | Higher score = stronger DDR programme",
    x        = "AUCell Score (ATR DDR Axis)",
    y        = NULL
  ) +
  theme_sc() +
  theme(axis.text.y = element_text(size = 10))

p_ridge_pikk <- ggplot(
  ridge_data %>% mutate(CellType = factor(CellType, levels = {
    ridge_data %>%
      group_by(CellType) %>%
      summarise(med = median(AUC_PIKK_Hub, na.rm = TRUE), .groups = "drop") %>%
      arrange(med) %>%
      pull(CellType)
  })),
  aes(x = AUC_PIKK_Hub, y = CellType, fill = CellType)
) +
  geom_density_ridges(
    scale       = 1.5,
    alpha       = 0.7,
    rel_min_height = 0.01,
    quantile_lines = TRUE,
    quantiles      = 2
  ) +
  scale_fill_manual(values = ct_colours, guide = "none") +
  labs(
    title    = "PIKK-OSCC Hub Interactome Activity by Cell Type",
    subtitle = "Ridge density | Vertical line: median",
    x        = "AUCell Score (PIKK Hub Interactome)",
    y        = NULL
  ) +
  theme_sc() +
  theme(axis.text.y = element_text(size = 10))

p_ridge_combined <- p_ridge_atr + p_ridge_pikk +
  plot_annotation(
    title = "DDR Programme Distribution Across the Tumour Microenvironment",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

save_dual(p_ridge_combined, "08_ridgeplot_ddr_by_celltype", w = 16, h = 8)

# ---- 6. Compartment-Level DDR Statistics -------------------------------------
cat("\n[6/9] Computing DDR score statistics by compartment ...\n")

ddr_compartment_stats <- seu@meta.data %>%
  group_by(CellType_Malignancy) %>%
  summarise(
    n_cells          = n(),
    mean_ATR_DDR     = mean(AUC_ATR_DDR, na.rm = TRUE),
    median_ATR_DDR   = median(AUC_ATR_DDR, na.rm = TRUE),
    sd_ATR_DDR       = sd(AUC_ATR_DDR, na.rm = TRUE),
    mean_PIKK_Hub    = mean(AUC_PIKK_Hub, na.rm = TRUE),
    median_PIKK_Hub  = median(AUC_PIKK_Hub, na.rm = TRUE),
    sd_PIKK_Hub      = sd(AUC_PIKK_Hub, na.rm = TRUE),
    .groups          = "drop"
  )

cat("\n       --- DDR Score by Compartment ---\n")
print(as.data.frame(ddr_compartment_stats), row.names = FALSE)

# Kruskal-Wallis for compartment differences
kw_atr  <- kruskal.test(AUC_ATR_DDR ~ CellType_Malignancy, data = seu@meta.data)
kw_pikk <- kruskal.test(AUC_PIKK_Hub ~ CellType_Malignancy, data = seu@meta.data)
cat("\n       KW test (ATR DDR by compartment):  p =", formatC(kw_atr$p.value, format = "e", digits = 2), "\n")
cat("       KW test (PIKK Hub by compartment): p =", formatC(kw_pikk$p.value, format = "e", digits = 2), "\n")

# Box/violin of DDR score by compartment
p_box_compartment <- ggplot(seu@meta.data,
       aes(x = CellType_Malignancy, y = AUC_ATR_DDR,
           fill = CellType_Malignancy)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
  scale_fill_manual(values = compartment_colours, guide = "none") +
  labs(
    title    = "ATR–DDR Programme Activity by Tumour Compartment",
    subtitle = paste0("Kruskal-Wallis p = ", formatC(kw_atr$p.value, format = "e", digits = 2)),
    x        = NULL,
    y        = "AUCell Score (ATR DDR Axis)"
  ) +
  theme_sc() +
  theme(axis.text.x = element_text(size = 11))

save_dual(p_box_compartment, "09_violin_ddr_by_compartment", w = 8, h = 6)

# ---- 7. Intra-Malignant Heterogeneity: DDR-High vs DDR-Low ------------------
cat("\n[7/9] Analysing intra-malignant DDR heterogeneity ...\n")

# Subset to malignant cells only
malignant_cells <- colnames(seu)[seu@meta.data$CellType_Malignancy == "Malignant cells"]
seu_malig <- subset(seu, cells = malignant_cells)
cat("       Malignant cell subset:", ncol(seu_malig), "cells\n")

# Stratify by ATR DDR score — upper quartile vs lower quartile
q75 <- quantile(seu_malig@meta.data$AUC_ATR_DDR, 0.75, na.rm = TRUE)
q25 <- quantile(seu_malig@meta.data$AUC_ATR_DDR, 0.25, na.rm = TRUE)

seu_malig@meta.data$DDR_status <- case_when(
  seu_malig@meta.data$AUC_ATR_DDR >= q75 ~ "DDR-High",
  seu_malig@meta.data$AUC_ATR_DDR <= q25 ~ "DDR-Low",
  TRUE ~ "DDR-Intermediate"
)

ddr_counts <- table(seu_malig@meta.data$DDR_status)
cat("       DDR-High:", ddr_counts["DDR-High"],
    "| DDR-Low:", ddr_counts["DDR-Low"],
    "| Intermediate:", ddr_counts["DDR-Intermediate"], "\n")

# Also add DDR_status back to the full Seurat object for later reference
full_ddr_status <- rep("Non-Malignant", ncol(seu))
names(full_ddr_status) <- colnames(seu)
full_ddr_status[names(seu_malig@meta.data$DDR_status)] <- seu_malig@meta.data$DDR_status
seu$DDR_status <- full_ddr_status

# UMAP showing DDR stratification in malignant cells
ddr_colours <- c(
  "DDR-High"         = "#D62828",
  "DDR-Intermediate" = "#F4A261",
  "DDR-Low"          = "#003049",
  "Non-Malignant"    = "grey85"
)

p_umap_ddr_strat <- DimPlot(
  seu,
  reduction  = "umap",
  group.by   = "DDR_status",
  cols       = ddr_colours,
  pt.size    = 0.4,
  order      = c("DDR-High", "DDR-Low", "DDR-Intermediate", "Non-Malignant")
) +
  labs(
    title    = "Malignant Cell DDR Stratification",
    subtitle = paste0("DDR-High (Q4): n=", ddr_counts["DDR-High"],
                      " | DDR-Low (Q1): n=", ddr_counts["DDR-Low"],
                      " | AUCell ATR–DDR Axis"),
    color    = "DDR Status"
  ) +
  theme_sc() +
  NoAxes()

save_dual(p_umap_ddr_strat, "10_umap_malignant_ddr_stratification", w = 10, h = 8)

# ---- 8. Differential Expression: DDR-High vs DDR-Low Malignant Cells --------
cat("\n[8/9] Running differential expression: DDR-High vs DDR-Low malignant cells ...\n")

# Subset to only DDR-High and DDR-Low cells for clean comparison
seu_de <- subset(seu_malig, subset = DDR_status %in% c("DDR-High", "DDR-Low"))
Idents(seu_de) <- "DDR_status"

cat("       DE comparison:", sum(Idents(seu_de) == "DDR-High"), "DDR-High vs",
    sum(Idents(seu_de) == "DDR-Low"), "DDR-Low malignant cells\n")

de_results <- FindMarkers(
  seu_de,
  ident.1       = "DDR-High",
  ident.2       = "DDR-Low",
  test.use      = "wilcox",
  logfc.threshold = 0.1,     # Low threshold to capture broad programme differences
  min.pct       = 0.1,
  verbose       = FALSE
)

de_results <- de_results %>%
  rownames_to_column("gene") %>%
  arrange(desc(avg_log2FC)) %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.001 ~ "***",
      p_val_adj < 0.01  ~ "**",
      p_val_adj < 0.05  ~ "*",
      TRUE              ~ "ns"
    ),
    direction = ifelse(avg_log2FC > 0, "Up_in_DDR_High", "Down_in_DDR_High")
  )

cat("       Total DE genes (adj.p < 0.05):", sum(de_results$p_val_adj < 0.05, na.rm = TRUE), "\n")
cat("         Upregulated in DDR-High:",
    sum(de_results$p_val_adj < 0.05 & de_results$avg_log2FC > 0, na.rm = TRUE), "\n")
cat("         Downregulated in DDR-High:",
    sum(de_results$p_val_adj < 0.05 & de_results$avg_log2FC < 0, na.rm = TRUE), "\n")

# Save full DE table
write_tsv(de_results, file.path(tbl_dir, "ddr_high_vs_low_malignant_DE.tsv"))

# Show top 15 up and top 15 down
cat("\n       --- Top 15 Upregulated in DDR-High Malignant Cells ---\n")
top_up <- de_results %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% head(15)
print(as.data.frame(top_up %>% select(gene, avg_log2FC, pct.1, pct.2, p_val_adj, significance)),
      row.names = FALSE)

cat("\n       --- Top 15 Downregulated in DDR-High Malignant Cells ---\n")
top_down <- de_results %>% filter(p_val_adj < 0.05, avg_log2FC < 0) %>% tail(15)
print(as.data.frame(top_down %>% select(gene, avg_log2FC, pct.1, pct.2, p_val_adj, significance)),
      row.names = FALSE)

# Check if key DDR/PIKK genes appear in the DE results
cat("\n       --- PIKK Diamond Panel in DDR-High vs DDR-Low DE ---\n")
diamond_panel <- c("PLK1", "CDK2", "TOPBP1", "RAD51", "FANCI", "KAT2B", "DEPTOR")
diamond_de <- de_results %>% filter(gene %in% diamond_panel)
if (nrow(diamond_de) > 0) {
  print(as.data.frame(diamond_de %>% select(gene, avg_log2FC, pct.1, pct.2, p_val_adj, significance)),
        row.names = FALSE)
} else {
  cat("       None of the Diamond Panel genes passed the DE threshold.\n")
}

# Volcano plot for DDR-High vs DDR-Low
de_plot_data <- de_results %>%
  mutate(
    neg_log10_padj = -log10(p_val_adj + 1e-300),
    label_gene = ifelse(
      gene %in% c(diamond_panel, head(top_up$gene, 5), tail(top_down$gene, 5)),
      gene, NA_character_
    ),
    color_group = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25  ~ "Up (DDR-High)",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down (DDR-High)",
      TRUE                                   ~ "Not significant"
    )
  )

volcano_colours <- c(
  "Up (DDR-High)"    = "#D62828",
  "Down (DDR-High)"  = "#003049",
  "Not significant"  = "grey75"
)

p_volcano <- ggplot(de_plot_data, aes(x = avg_log2FC, y = neg_log10_padj, color = color_group)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_text(aes(label = label_gene), size = 3, hjust = -0.1, vjust = 0.5,
            na.rm = TRUE, color = "black", fontface = "italic") +
  scale_color_manual(values = volcano_colours, name = "Direction") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.4) +
  labs(
    title    = "Differential Expression: DDR-High vs DDR-Low Malignant Cells",
    subtitle = paste0("Wilcoxon rank-sum test | ",
                      sum(de_plot_data$color_group == "Up (DDR-High)"), " up, ",
                      sum(de_plot_data$color_group == "Down (DDR-High)"), " down"),
    x        = "Average log₂ Fold Change (DDR-High / DDR-Low)",
    y        = "-log₁₀ Adjusted P-value"
  ) +
  theme_sc() +
  theme(legend.position = "bottom")

save_dual(p_volcano, "11_volcano_ddr_high_vs_low_malignant", w = 10, h = 8)

# ---- 9. Patient-Level DDR-High Proportions -----------------------------------
cat("\n[9/9] Computing patient-level DDR-high proportions ...\n")

patient_ddr <- seu_malig@meta.data %>%
  group_by(Patient) %>%
  summarise(
    n_malignant_cells   = n(),
    n_ddr_high          = sum(DDR_status == "DDR-High", na.rm = TRUE),
    n_ddr_low           = sum(DDR_status == "DDR-Low", na.rm = TRUE),
    pct_ddr_high        = n_ddr_high / n_malignant_cells * 100,
    mean_atr_ddr_score  = mean(AUC_ATR_DDR, na.rm = TRUE),
    mean_pikk_hub_score = mean(AUC_PIKK_Hub, na.rm = TRUE),
    stage               = first(TNMstage),
    source              = first(Source),
    .groups             = "drop"
  ) %>%
  arrange(desc(pct_ddr_high))

cat("\n       --- Patient-Level DDR-High Proportions ---\n")
print(as.data.frame(patient_ddr), row.names = FALSE)

write_tsv(patient_ddr, file.path(tbl_dir, "patient_ddr_high_proportion.tsv"))

# Bar plot: % DDR-high malignant cells per patient, coloured by stage
p_patient_bar <- ggplot(
  patient_ddr %>% mutate(Patient = reorder(Patient, pct_ddr_high)),
  aes(x = Patient, y = pct_ddr_high, fill = stage)
) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = paste0(round(pct_ddr_high, 1), "%")),
            hjust = -0.1, size = 3) +
  scale_fill_brewer(palette = "Set2", name = "TNM Stage") +
  coord_flip() +
  labs(
    title    = "Proportion of DDR-High Malignant Cells per Patient",
    subtitle = "DDR-High = upper quartile of ATR–DDR AUCell score among all malignant cells",
    x        = "Patient",
    y        = "% DDR-High Malignant Cells"
  ) +
  theme_sc() +
  theme(legend.position = "bottom") +
  ylim(0, max(patient_ddr$pct_ddr_high, na.rm = TRUE) * 1.15)

save_dual(p_patient_bar, "12_barplot_patient_ddr_high_proportion", w = 9, h = 7)

# Scatter: mean DDR score vs stage (if enough patients per stage)
if (length(unique(patient_ddr$stage[!is.na(patient_ddr$stage)])) > 1) {
  p_stage_scatter <- ggplot(
    patient_ddr %>% filter(!is.na(stage)),
    aes(x = stage, y = mean_atr_ddr_score, fill = stage)
  ) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2.5, alpha = 0.8) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(
      title    = "Mean ATR–DDR Score vs Clinical Stage",
      subtitle = "Per-patient mean AUCell score of malignant cells",
      x        = "TNM Stage",
      y        = "Mean ATR–DDR AUCell Score"
    ) +
    theme_sc()
  
  save_dual(p_stage_scatter, "13_boxplot_ddr_score_vs_stage", w = 7, h = 6)
}

# ---- Save Full AUCell Scores Table -------------------------------------------
cat("\n       Saving per-cell AUCell scores ...\n")
aucell_export <- seu@meta.data %>%
  rownames_to_column("Cell") %>%
  select(Cell, CellType_Malignancy, CellType_Major,
         Patient, Source, TNMstage,
         AUC_ATR_DDR, AUC_PIKK_Hub, DDR_status)

write_tsv(aucell_export, file.path(tbl_dir, "ddr_aucell_scores_per_cell.tsv"))

# ---- Save updated Seurat object with DDR scores -----------------------------
cat("       Saving updated Seurat object with AUCell scores ...\n")
saveRDS(seu, file.path(rds_dir, "oscc_gse103322_seurat.rds"))

# ==============================================================================
cat("\n================================================================\n")
cat("  Module 2 complete:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("  Figures generated: 7\n")
cat("  Tables generated:  3\n")
cat("  Key outputs:\n")
cat("    - AUCell DDR scores per cell\n")
cat("    - DDR-High vs DDR-Low malignant cell DE table\n")
cat("    - Patient-level DDR-high proportions\n")
cat("================================================================\n\n")
