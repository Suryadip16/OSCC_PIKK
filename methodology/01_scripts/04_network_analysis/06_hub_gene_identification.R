#!/usr/bin/env Rscript

# ============================================================================
# STEP 2.3: Hub Gene Identification
# OSCC-PIKK DDR & Checkpoint Signalling
# ============================================================================
#
# This script identifies hub genes from the PPI network (Step 2.2), correlates
# their expression with clinical variables, and classifies their pathway roles.
#
# Workflow:
#   1. Load network topology metrics (from Step 2.2)
#   2. Select hub genes by degree, betweenness, and composite score
#   3. Correlate hub gene expression with OSCC clinical variables
#   4. Classify hub genes into pathway role categories
#   5. Generate final hub gene characterisation table
#   6. Produce publication-quality visualisations
#
# Required R packages:
#   CRAN: igraph, ggplot2, ggrepel, dplyr, tidyr, readr, tibble,
#         pheatmap, RColorBrewer, scales, cowplot, reshape2, ggpubr
#
# Inputs:
#   - network_topology_metrics.tsv (from Step 2.2)
#   - OSCC_PIKK_extended_edgeR_DEG_intersection_significant.tsv
#   - logCPM_matrix.tsv (edgeR normalised expression)
#   - matched_metadata_used.tsv (clinical metadata)
#   - PIKK_related_gene_universe.csv
#   - PPI_core_network.graphml (for hub gene network highlight)
#
# Outputs:
#   - Tables: hub_genes_*.tsv, clinical correlation, characterisation table
#   - Plots : 8 publication-quality figures (PDF + PNG)

# Date   : August 2026
# ============================================================================


# ============================================================================
# 0. CONFIGURATION & PACKAGE LOADING
# ============================================================================

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(cowplot)
  library(reshape2)
})

# ggpubr for statistical comparisons on boxplots (optional)
has_ggpubr <- requireNamespace("ggpubr", quietly = TRUE)
if (has_ggpubr) {
  library(ggpubr)
} else {
  message(
    "NOTE: ggpubr not installed. Statistical annotations on boxplots ",
    "will be omitted."
  )
  message('Install with: install.packages("ggpubr")')
}

# ggraph + tidygraph for the hub network highlight
has_ggraph <- requireNamespace("ggraph", quietly = TRUE) &&
  requireNamespace("tidygraph", quietly = TRUE)
if (has_ggraph) {
  library(ggraph)
  library(tidygraph)
}

# --- Project root & paths ---------------------------------------------------
project_root <- "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project"

# Input files
topo_file <- file.path(
  project_root,
  "05_results/tables/network/network_topology_metrics.tsv"
)
pikk_deg_file <- file.path(
  project_root,
  "05_results/Deepaprabha_Intermediary_Results/OSCC_PIKK_extended_edgeR_DEG_intersection_significant.tsv"
)
logcpm_file <- file.path(
  project_root,
  "05_results/Deepaprabha_Intermediary_Results/edgeR_withBatch/results/logCPM_matrix.tsv"
)
meta_file <- file.path(
  project_root,
  "05_results/Deepaprabha_Intermediary_Results/edgeR_withBatch/results/matched_metadata_used.tsv"
)
pikk_univ_file <- file.path(
  project_root,
  "10_references/PIKK_related_gene_universe.csv"
)
graphml_file <- file.path(
  project_root,
  "05_results/network_files/PPI_core_network.graphml"
)

# Output directories
out_tables <- file.path(project_root, "05_results/tables/hub_genes")
out_plots <- file.path(project_root, "05_results/plots/hub_genes")
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)

# --- Hub gene selection parameters -------------------------------------------
top_n_degree <- 15 # top genes by degree
top_n_betweenness <- 15 # top genes by betweenness
# Composite score: z(degree) + z(betweenness) + z(|logFC|)

# --- Aesthetics ---------------------------------------------------------------
theme_pub <- theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "black"),
    legend.title  = element_text(face = "bold"),
    panel.grid    = element_blank()
  )

pikk_colors <- c(
  "ATR"   = "#E41A1C",
  "ATM"   = "#377EB8",
  "PRKDC" = "#4DAF4A",
  "MTOR"  = "#984EA3",
  "SMG1"  = "#FF7F00",
  "TRRAP" = "#A65628"
)

reg_colors <- c(
  "Up in Tumor"   = "#D73027",
  "Down in Tumor" = "#4575B4"
)

role_colors <- c(
  "Core PIKK Kinase" = "#E41A1C",
  "Checkpoint Regulator" = "#FF7F00",
  "Cell Cycle Regulator" = "#FFFF33",
  "DNA Repair Effector" = "#4DAF4A",
  "Signalling Effector" = "#377EB8",
  "Chromatin/Transcription" = "#984EA3",
  "mRNA Surveillance" = "#A65628",
  "Other" = "#999999"
)

save_pdf_png <- function(plot_obj, filename_base, width = 10, height = 8) {
  ggsave(file.path(out_plots, paste0(filename_base, ".pdf")),
    plot = plot_obj, width = width, height = height, device = cairo_pdf
  )
  ggsave(file.path(out_plots, paste0(filename_base, ".png")),
    plot = plot_obj, width = width, height = height, dpi = 300
  )
  cat("  Saved:", filename_base, "(PDF + PNG)\n")
}


# ============================================================================
# 1. LOAD INPUTS
# ============================================================================
cat("\n================================================================\n")
cat("STEP 2.3: Hub Gene Identification — OSCC-PIKK\n")
cat("================================================================\n\n")

# --- 1a. Topology metrics (from Step 2.2) ------------------------------------
topo_df <- read.delim(topo_file,
  header = TRUE, sep = "\t",
  check.names = FALSE, stringsAsFactors = FALSE
)
cat("Loaded topology metrics for", nrow(topo_df), "genes.\n")

# --- 1b. PIKK DEG intersection (for full DE info) ----------------------------
pikk_deg <- read.delim(pikk_deg_file,
  header = TRUE, sep = "\t",
  check.names = FALSE, stringsAsFactors = FALSE
)
pikk_deg$gene_symbol <- toupper(pikk_deg$gene_name)
cat("Loaded", nrow(pikk_deg), "PIKK DEGs.\n")

# --- 1c. logCPM expression matrix -------------------------------------------
cat("Loading logCPM matrix (this may take a moment)...\n")
logcpm <- read.delim(logcpm_file,
  header = TRUE, sep = "\t",
  check.names = FALSE, stringsAsFactors = FALSE,
  row.names = 1
)
cat("logCPM matrix:", nrow(logcpm), "genes x", ncol(logcpm), "samples.\n")

# --- 1d. Clinical metadata ---------------------------------------------------
meta <- read.delim(meta_file,
  header = TRUE, sep = "\t",
  check.names = FALSE, stringsAsFactors = FALSE
)
cat("Metadata:", nrow(meta), "samples.\n")

# --- 1e. PIKK gene universe --------------------------------------------------
pikk_univ <- read.csv(pikk_univ_file, stringsAsFactors = FALSE)
pikk_univ$GeneSymbol <- toupper(pikk_univ$GeneSymbol)


# ============================================================================
# 2. HUB GENE SELECTION
# ============================================================================
cat("\n--- Hub Gene Selection ---\n")

# --- 2a. Top-N by degree -----------------------------------------------------
top_degree <- topo_df %>%
  dplyr::arrange(desc(degree)) %>%
  dplyr::slice_head(n = top_n_degree)

cat("Top", top_n_degree, "by degree:\n")
cat(paste(top_degree$gene_symbol, collapse = ", "), "\n")

# --- 2b. Top-N by betweenness ------------------------------------------------
top_btwn <- topo_df %>%
  dplyr::arrange(desc(betweenness)) %>%
  dplyr::slice_head(n = top_n_betweenness)

cat("Top", top_n_betweenness, "by betweenness:\n")
cat(paste(top_btwn$gene_symbol, collapse = ", "), "\n")

# --- 2c. Union hub gene set --------------------------------------------------
hub_genes_union <- unique(c(top_degree$gene_symbol, top_btwn$gene_symbol))
cat("Union hub gene set:", length(hub_genes_union), "genes.\n")

# --- 2d. Composite hub score -------------------------------------------------
# Composite = z(degree) + z(betweenness) + z(|logFC|)
# Applied to the union set only
hub_df <- topo_df %>%
  dplyr::filter(gene_symbol %in% hub_genes_union) %>%
  dplyr::mutate(
    abs_logFC = abs(log2FoldChange),
    z_degree = as.numeric(scale(degree)),
    z_betweenness = as.numeric(scale(betweenness)),
    z_absLogFC = as.numeric(scale(abs_logFC)),
    composite_score = z_degree + z_betweenness + z_absLogFC
  ) %>%
  dplyr::arrange(desc(composite_score))

cat("\nComposite hub gene ranking (top 10):\n")
print(hub_df %>%
  dplyr::select(
    gene_symbol, degree, betweenness, log2FoldChange,
    composite_score, PIKK_Group
  ) %>%
  head(10))

# Save hub gene tables
write.table(top_degree,
  file = file.path(out_tables, "hub_genes_top_degree.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(top_btwn,
  file = file.path(out_tables, "hub_genes_top_betweenness.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(hub_df,
  file = file.path(out_tables, "hub_genes_union_composite.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)


# ============================================================================
# 3. CLINICAL CORRELATION ANALYSIS
# ============================================================================
cat("\n--- Clinical Correlation Analysis ---\n")

# Subset to tumour samples only (clinical variables like stage are tumour-specific)
meta_tumor <- meta %>%
  dplyr::filter(Tissue_Type == "Tumor")
cat("Tumour samples for clinical correlation:", nrow(meta_tumor), "\n")

# Match samples in logCPM matrix
sample_ids <- intersect(meta_tumor$File_ID, colnames(logcpm))
if (length(sample_ids) == 0) {
  # Try matching with row names or alternative ID column
  sample_ids <- intersect(meta_tumor$File_ID, colnames(logcpm))
  if (length(sample_ids) == 0) {
    message("WARNING: Cannot match sample IDs between metadata and logCPM matrix.")
    message("Metadata File_ID examples: ", paste(head(meta_tumor$File_ID), collapse = ", "))
    message("logCPM column examples: ", paste(head(colnames(logcpm)), collapse = ", "))
  }
}
cat("Matched tumour samples:", length(sample_ids), "\n")

# Extract hub gene expression for matched tumour samples
hub_genes_in_matrix <- intersect(hub_df$gene_symbol, rownames(logcpm))
if (length(hub_genes_in_matrix) == 0) {
  # Try matching with toupper
  rownames(logcpm) <- toupper(rownames(logcpm))
  hub_genes_in_matrix <- intersect(hub_df$gene_symbol, rownames(logcpm))
}
cat(
  "Hub genes found in logCPM matrix:", length(hub_genes_in_matrix),
  "of", nrow(hub_df), "\n"
)

# Build long-format expression table for plotting
hub_expr <- logcpm[hub_genes_in_matrix, sample_ids, drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_symbol") %>%
  tidyr::pivot_longer(
    cols      = -gene_symbol,
    names_to  = "File_ID",
    values_to = "logCPM"
  ) %>%
  dplyr::left_join(
    meta_tumor %>%
      dplyr::select(
        File_ID, Tissue_Type,
        diagnoses.ajcc_pathologic_stage,
        demographic.age_at_index,
        demographic.gender,
        demographic.vital_status,
        exposures.tobacco_smoking_status
      ),
    by = "File_ID"
  )

# --- Clean clinical variables ------------------------------------------------
# Simplify stage (group IVA/IVB/IVC into IV)
hub_expr <- hub_expr %>%
  dplyr::mutate(
    stage_simple = dplyr::case_when(
      grepl("Stage I$|Stage I[^V]", diagnoses.ajcc_pathologic_stage) ~ "Stage I",
      grepl("Stage II$|Stage II[^I]", diagnoses.ajcc_pathologic_stage) ~ "Stage II",
      grepl("Stage III", diagnoses.ajcc_pathologic_stage) ~ "Stage III",
      grepl("Stage IV", diagnoses.ajcc_pathologic_stage) ~ "Stage IV",
      TRUE ~ NA_character_
    ),
    stage_simple = factor(stage_simple,
      levels = c(
        "Stage I", "Stage II",
        "Stage III", "Stage IV"
      )
    ),
    vital_status = demographic.vital_status,
    smoking = exposures.tobacco_smoking_status,
    gender = demographic.gender,
    age = as.numeric(demographic.age_at_index)
  )

cat("Clinical variable levels:\n")
cat("  Stage: ", paste(levels(hub_expr$stage_simple), collapse = ", "), "\n")
cat("  Vital status: ", paste(unique(hub_expr$vital_status), collapse = ", "), "\n")

# --- 3a. Statistical tests: hub gene expression vs clinical variables --------
clinical_tests <- data.frame(
  gene_symbol = character(),
  variable = character(),
  test = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (gene in hub_genes_in_matrix) {
  gene_data <- hub_expr %>% dplyr::filter(gene_symbol == gene)

  # Stage: Kruskal-Wallis test
  stage_data <- gene_data %>% dplyr::filter(!is.na(stage_simple))
  if (length(unique(stage_data$stage_simple)) >= 2) {
    kw <- kruskal.test(logCPM ~ stage_simple, data = stage_data)
    clinical_tests <- rbind(clinical_tests, data.frame(
      gene_symbol = gene, variable = "AJCC_Stage",
      test = "Kruskal-Wallis", p_value = kw$p.value
    ))
  }

  # Vital status: Wilcoxon rank-sum test
  vital_data <- gene_data %>%
    dplyr::filter(!is.na(vital_status), vital_status %in% c("Alive", "Dead"))
  if (length(unique(vital_data$vital_status)) == 2) {
    wt <- wilcox.test(logCPM ~ vital_status, data = vital_data)
    clinical_tests <- rbind(clinical_tests, data.frame(
      gene_symbol = gene, variable = "Vital_Status",
      test = "Wilcoxon", p_value = wt$p.value
    ))
  }

  # Gender: Wilcoxon rank-sum test
  gender_data <- gene_data %>%
    dplyr::filter(!is.na(gender), gender %in% c("male", "female"))
  if (length(unique(gender_data$gender)) == 2) {
    wt_g <- wilcox.test(logCPM ~ gender, data = gender_data)
    clinical_tests <- rbind(clinical_tests, data.frame(
      gene_symbol = gene, variable = "Gender",
      test = "Wilcoxon", p_value = wt_g$p.value
    ))
  }

  # Age: Spearman correlation
  age_data <- gene_data %>% dplyr::filter(!is.na(age))
  if (nrow(age_data) >= 10) {
    cor_test <- cor.test(age_data$logCPM, age_data$age, method = "spearman")
    clinical_tests <- rbind(clinical_tests, data.frame(
      gene_symbol = gene, variable = "Age",
      test = "Spearman_correlation", p_value = cor_test$p.value
    ))
  }
}

# Adjust p-values (BH correction within each clinical variable)
# Report BOTH raw and adjusted p-values for transparency
clinical_tests <- clinical_tests %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(
    p_adjusted    = p.adjust(p_value, method = "BH"),
    sig_nominal   = p_value < 0.05, # exploratory threshold
    sig_FDR       = p_adjusted < 0.05 # conservative threshold
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(variable, p_value)

cat("\nClinical correlation tests completed:", nrow(clinical_tests), "tests.\n")
cat(
  "Nominally significant (p < 0.05):  ",
  sum(clinical_tests$sig_nominal), "\n"
)
cat(
  "FDR-significant (padj < 0.05):     ",
  sum(clinical_tests$sig_FDR), "\n"
)

write.table(clinical_tests,
  file = file.path(
    out_tables,
    "hub_gene_clinical_correlation_tests.tsv"
  ),
  sep = "\t", quote = FALSE, row.names = FALSE
)


# ============================================================================
# 4. PATHWAY ROLE CLASSIFICATION
# ============================================================================
cat("\n--- Pathway Role Classification ---\n")

# Define pathway role categories based on known DDR/checkpoint biology
# and the PIKK gene universe annotations
checkpoint_genes <- c(
  "CHEK1", "CHEK2", "TP53", "TP53BP1", "MDC1",
  "CLSPN", "TOPBP1", "RAD17", "RAD1", "RAD9A",
  "RAD9B", "HUS1", "ATRIP", "E2F1"
)

cell_cycle_genes <- c(
  "PLK1", "AURKA", "AURKB", "CCNE1", "CDK2",
  "CDC45", "CDC7", "CDKN2A", "GADD45B"
)

dna_repair_genes <- c(
  "BRCA1", "BRCA2", "RAD51", "RAD50", "RAD52",
  "EXO1", "PARP1", "XRCC4", "XRCC1", "XRCC5",
  "XRCC6", "LIG4", "FANCA", "FANCD2", "FANCI",
  "FANCG", "FANCM", "MRE11", "H2AX", "MRE11A",
  "H2AFX", "RBBP8", "POLD1", "MLH1", "MSH2",
  "MSH6", "BLM", "WRN", "DCLRE1C", "NHEJ1",
  "PAXX", "SWSAP1", "ERCC1", "ERCC4", "ERCC5",
  "ERCC8", "XPC", "RAD23B", "NBN"
)

signalling_genes <- c(
  "EGFR", "AKT1", "AKT2", "MTOR", "RPTOR",
  "RICTOR", "DEPTOR", "RHEB", "TSC1", "TSC2",
  "EIF4EBP1", "EIF4EBP2", "RPS6KB1", "RPS6KB2",
  "ULK1", "ULK2", "PRKAA1", "PRKAA2", "FKBP1A",
  "MLST8", "MAPKAP1", "PRR5", "PRR5L", "AKT1S1",
  "IRS1", "IRS2", "RRAGA", "RRAGB", "RRAGC",
  "RRAGD", "SGK1", "PRKCA", "SQSTM1", "TFEB",
  "RB1CC1", "ATG101", "LAMP1", "LAMTOR2"
)

chromatin_genes <- c(
  "KAT2A", "KAT2B", "KAT5", "KAT8", "KDM1A",
  "INO80", "TRRAP", "EP400", "DMAP1", "BRD8",
  "ACTL6A", "ING3", "MBTD1", "MORF4L1", "MORF4L2",
  "EPC1", "EPC2", "YEATS4", "RUVBL1", "RUVBL2",
  "TAF2", "TAF4", "TAF5", "TAF5L", "TAF6",
  "TAF6L", "TAF7", "TAF9", "TAF9B", "TAF10",
  "TAF12", "SUPT3H", "SUPT7L", "SUPT20H",
  "TADA1", "TADA2A", "TADA2B", "TADA3",
  "MRGBP", "MEAF6", "SGF29", "USP22",
  "ATXN7", "ATXN7L3", "VPS72", "PHF21A", "RCOR1",
  "HOXB7", "MYC", "PIAS4", "RNF144A"
)

mrna_surveillance <- c(
  "SMG1", "SMG5", "SMG6", "SMG7", "SMG8", "SMG9",
  "UPF1", "UPF2", "UPF3A", "UPF3B", "ETF1",
  "GSPT1", "GSPT2", "EIF4A3", "RBM8A", "NCBP2",
  "PABPC1"
)

core_pikk_genes <- c("ATR", "ATM", "PRKDC", "MTOR", "SMG1", "TRRAP")

# Classify hub genes
hub_df <- hub_df %>%
  dplyr::mutate(
    pathway_role = dplyr::case_when(
      gene_symbol %in% core_pikk_genes ~ "Core PIKK Kinase",
      gene_symbol %in% checkpoint_genes ~ "Checkpoint Regulator",
      gene_symbol %in% cell_cycle_genes ~ "Cell Cycle Regulator",
      gene_symbol %in% dna_repair_genes ~ "DNA Repair Effector",
      gene_symbol %in% signalling_genes ~ "Signalling Effector",
      gene_symbol %in% chromatin_genes ~ "Chromatin/Transcription",
      gene_symbol %in% mrna_surveillance ~ "mRNA Surveillance",
      TRUE ~ "Other"
    )
  )

cat("Pathway role classification:\n")
print(hub_df %>% dplyr::count(pathway_role, name = "n") %>%
  dplyr::arrange(desc(n)))


# ============================================================================
# 5. FINAL HUB GENE CHARACTERISATION TABLE
# ============================================================================
cat("\n--- Final Hub Gene Characterisation ---\n")

# Merge with clinical significance — report BOTH nominal and FDR levels
clinical_summary_nominal <- clinical_tests %>%
  dplyr::filter(sig_nominal) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(
    sig_clinical_nominal = paste(variable, collapse = "; "),
    .groups = "drop"
  )

clinical_summary_fdr <- clinical_tests %>%
  dplyr::filter(sig_FDR) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(
    sig_clinical_FDR = paste(variable, collapse = "; "),
    .groups = "drop"
  )

hub_characterisation <- hub_df %>%
  dplyr::select(
    gene_symbol, PIKK_Group, GeneRole, regulation,
    log2FoldChange, padj, degree, betweenness, closeness,
    eigenvector, composite_score, community, pathway_role
  ) %>%
  dplyr::left_join(clinical_summary_nominal, by = "gene_symbol") %>%
  dplyr::left_join(clinical_summary_fdr, by = "gene_symbol") %>%
  dplyr::mutate(
    sig_clinical_nominal = ifelse(is.na(sig_clinical_nominal), "None",
      sig_clinical_nominal
    ),
    sig_clinical_FDR = ifelse(is.na(sig_clinical_FDR), "None",
      sig_clinical_FDR
    )
  ) %>%
  dplyr::arrange(desc(composite_score))

cat("Final hub gene characterisation table:", nrow(hub_characterisation), "genes.\n")
print(hub_characterisation %>%
  dplyr::select(
    gene_symbol, PIKK_Group, pathway_role, regulation,
    composite_score, sig_clinical_nominal, sig_clinical_FDR
  ))

write.table(hub_characterisation,
  file = file.path(
    out_tables,
    "hub_gene_characterisation_table.tsv"
  ),
  sep = "\t", quote = FALSE, row.names = FALSE
)


# ============================================================================
# 6. VISUALISATIONS
# ============================================================================
cat("\n--- Hub Gene Visualisations ---\n")

# --- 6.1 Composite score lollipop chart ---------------------------------------
hub_plot_df <- hub_characterisation %>%
  dplyr::mutate(
    gene_symbol = factor(gene_symbol, levels = rev(gene_symbol))
  )

p_lollipop <- ggplot(
  hub_plot_df,
  aes(x = composite_score, y = gene_symbol)
) +
  geom_segment(aes(xend = 0, yend = gene_symbol, colour = PIKK_Group),
    linewidth = 1
  ) +
  geom_point(aes(fill = PIKK_Group, size = degree),
    shape = 21,
    colour = "grey30", stroke = 0.5
  ) +
  scale_fill_manual(values = pikk_colors, name = "PIKK Group") +
  scale_colour_manual(values = pikk_colors, guide = "none") +
  scale_size_continuous(range = c(3, 10), name = "Degree") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(
    title = "Hub Gene Composite Score Ranking",
    subtitle = "Score = z(degree) + z(betweenness) + z(|log2FC|)",
    x = "Composite Hub Score", y = NULL
  ) +
  theme_pub +
  theme(axis.text.y = element_text(face = "bold", size = 9))
save_pdf_png(p_lollipop, "hub_genes_composite_score_lollipop",
  width = 10, height = max(6, nrow(hub_plot_df) * 0.35 + 2)
)


# --- 6.2 Multi-metric dot plot (radar alternative) ----------------------------
# Show normalised degree, betweenness, |logFC| side-by-side for each hub gene
metric_long <- hub_characterisation %>%
  dplyr::select(gene_symbol,
    z_degree = degree, z_betweenness = betweenness,
    abs_logFC = log2FoldChange
  ) %>%
  dplyr::mutate(abs_logFC = abs(abs_logFC)) %>%
  # Min-max scale each metric to [0, 1]
  dplyr::mutate(
    across(
      c(z_degree, z_betweenness, abs_logFC),
      ~ (. - min(., na.rm = TRUE)) /
        (max(., na.rm = TRUE) - min(., na.rm = TRUE) + 1e-10)
    )
  ) %>%
  tidyr::pivot_longer(
    cols      = c(z_degree, z_betweenness, abs_logFC),
    names_to  = "metric",
    values_to = "scaled_value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(metric,
      "z_degree"      = "Degree",
      "z_betweenness" = "Betweenness",
      "abs_logFC"     = "|log2FC|"
    ),
    gene_symbol = factor(gene_symbol,
      levels = rev(hub_characterisation$gene_symbol)
    )
  )

p_multimetric <- ggplot(
  metric_long,
  aes(
    x = scaled_value, y = gene_symbol,
    colour = metric
  )
) +
  geom_point(size = 4, alpha = 0.85) +
  scale_colour_manual(
    values = c(
      "Degree" = "#E41A1C", "Betweenness" = "#377EB8",
      "|log2FC|" = "#4DAF4A"
    ),
    name = "Metric"
  ) +
  labs(
    title = "Hub Gene Multi-Metric Profile",
    subtitle = "Each metric min–max normalised to [0, 1]",
    x = "Normalised Score", y = NULL
  ) +
  theme_pub +
  theme(axis.text.y = element_text(face = "bold", size = 8))
save_pdf_png(p_multimetric, "hub_genes_multimetric_dotplot",
  width = 9, height = max(6, nrow(hub_characterisation) * 0.35 + 2)
)


# --- 6.3 Hub gene expression by AJCC stage (faceted boxplots) -----------------
if (length(hub_genes_in_matrix) > 0 && nrow(hub_expr) > 0) {
  # Use top hub genes (by composite score) for the faceted plot
  top_hub_for_plot <- hub_characterisation$gene_symbol[
    hub_characterisation$gene_symbol %in% hub_genes_in_matrix
  ]
  # Limit to ~12 genes for readability
  top_hub_for_plot <- head(top_hub_for_plot, 12)

  stage_plot_data <- hub_expr %>%
    dplyr::filter(
      gene_symbol %in% top_hub_for_plot,
      !is.na(stage_simple)
    ) %>%
    dplyr::mutate(
      gene_symbol = factor(gene_symbol, levels = top_hub_for_plot)
    )

  p_stage <- ggplot(
    stage_plot_data,
    aes(x = stage_simple, y = logCPM, fill = stage_simple)
  ) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8, colour = "grey30") +
    facet_wrap(~gene_symbol, scales = "free_y", ncol = 4) +
    scale_fill_brewer(palette = "YlOrRd", name = "AJCC Stage") +
    labs(
      title = "Hub Gene Expression by AJCC Pathologic Stage",
      subtitle = "OSCC tumour samples | logCPM (TMM-normalised)",
      x = "AJCC Stage", y = "logCPM"
    ) +
    theme_pub +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "none"
    )
  save_pdf_png(p_stage, "hub_gene_expression_by_stage_facet",
    width = 12, height = max(6, ceiling(length(top_hub_for_plot) / 4) * 3.5)
  )


  # --- 6.4 Hub gene expression by vital status --------------------------------
  vital_plot_data <- hub_expr %>%
    dplyr::filter(
      gene_symbol %in% top_hub_for_plot,
      !is.na(vital_status),
      vital_status %in% c("Alive", "Dead")
    ) %>%
    dplyr::mutate(
      gene_symbol = factor(gene_symbol, levels = top_hub_for_plot)
    )

  p_vital <- ggplot(
    vital_plot_data,
    aes(x = vital_status, y = logCPM, fill = vital_status)
  ) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8, colour = "grey30") +
    facet_wrap(~gene_symbol, scales = "free_y", ncol = 4) +
    scale_fill_manual(
      values = c("Alive" = "#4DAF4A", "Dead" = "#E41A1C"),
      name = "Vital Status"
    ) +
    labs(
      title = "Hub Gene Expression by Vital Status",
      subtitle = "OSCC tumour samples | logCPM (TMM-normalised)",
      x = "Vital Status", y = "logCPM"
    ) +
    theme_pub +
    theme(
      strip.text       = element_text(face = "bold", size = 9),
      legend.position  = "bottom"
    )

  # Add statistical annotations if ggpubr is available
  if (has_ggpubr) {
    p_vital <- p_vital +
      stat_compare_means(
        method = "wilcox.test", label = "p.signif",
        label.x = 1.5, label.y.npc = "top", size = 3.5
      )
  }
  save_pdf_png(p_vital, "hub_gene_expression_by_vital_status",
    width = 12, height = max(6, ceiling(length(top_hub_for_plot) / 4) * 3.5)
  )


  # --- 6.3b Stage boxplots with NOMINAL p-values (exploratory) ----------------
  # Annotated with raw Kruskal-Wallis p-values per gene facet
  if (has_ggpubr) {
    p_stage_nominal <- ggplot(
      stage_plot_data,
      aes(x = stage_simple, y = logCPM, fill = stage_simple)
    ) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.8, colour = "grey30") +
      facet_wrap(~gene_symbol, scales = "free_y", ncol = 4) +
      scale_fill_brewer(palette = "YlOrRd", name = "AJCC Stage") +
      stat_compare_means(
        method = "kruskal.test",
        label = "p.format",
        label.y.npc = "top", size = 3
      ) +
      labs(
        title = "Hub Gene Expression by AJCC Stage — Nominal p-values",
        subtitle = "Kruskal-Wallis test (unadjusted) | OSCC tumour samples",
        x = "AJCC Stage", y = "logCPM"
      ) +
      theme_pub +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "none"
      )
    save_pdf_png(p_stage_nominal, "hub_gene_expression_by_stage_nominal_pval",
      width = 12, height = max(6, ceiling(length(top_hub_for_plot) / 4) * 3.5)
    )
  }


  # --- 6.4b Vital status boxplots with NOMINAL p-values (exploratory) ---------
  # Annotated with raw Wilcoxon p-values per gene facet
  if (has_ggpubr) {
    p_vital_nominal <- ggplot(
      vital_plot_data,
      aes(x = vital_status, y = logCPM, fill = vital_status)
    ) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.8, colour = "grey30") +
      facet_wrap(~gene_symbol, scales = "free_y", ncol = 4) +
      scale_fill_manual(
        values = c("Alive" = "#4DAF4A", "Dead" = "#E41A1C"),
        name = "Vital Status"
      ) +
      stat_compare_means(
        method = "wilcox.test",
        label = "p.format",
        label.x = 1.5, label.y.npc = "top", size = 3
      ) +
      labs(
        title = "Hub Gene Expression by Vital Status — Nominal p-values",
        subtitle = "Wilcoxon test (unadjusted) | OSCC tumour samples",
        x = "Vital Status", y = "logCPM"
      ) +
      theme_pub +
      theme(
        strip.text       = element_text(face = "bold", size = 9),
        legend.position  = "bottom"
      )
    save_pdf_png(p_vital_nominal, "hub_gene_expression_by_vital_status_nominal_pval",
      width = 12, height = max(6, ceiling(length(top_hub_for_plot) / 4) * 3.5)
    )
  } else {
    message("Skipping nominal p-value annotated plots (ggpubr not installed).")
  }

  # --- 6.4c Individual plots for each nominally significant association --------
  # One standalone boxplot per gene × clinical variable (p < 0.05, unadjusted)
  cat("\n--- Individual plots for nominally significant associations ---\n")

  sig_nominal_tests <- clinical_tests %>%
    dplyr::filter(sig_nominal, gene_symbol %in% hub_genes_in_matrix)

  if (nrow(sig_nominal_tests) > 0) {
    indiv_dir <- file.path(out_plots, "individual_clinical")
    dir.create(indiv_dir, recursive = TRUE, showWarnings = FALSE)

    for (i in seq_len(nrow(sig_nominal_tests))) {
      gene <- sig_nominal_tests$gene_symbol[i]
      var <- sig_nominal_tests$variable[i]
      pval <- sig_nominal_tests$p_value[i]
      padj <- sig_nominal_tests$p_adjusted[i]
      test_name <- sig_nominal_tests$test[i]

      gene_data <- hub_expr %>% dplyr::filter(gene_symbol == gene)

      # Determine the x-axis variable and filter accordingly
      if (var == "AJCC_Stage") {
        plot_data <- gene_data %>%
          dplyr::filter(!is.na(stage_simple))
        if (nrow(plot_data) == 0) next

        p_ind <- ggplot(
          plot_data,
          aes(x = stage_simple, y = logCPM, fill = stage_simple)
        ) +
          geom_boxplot(
            alpha = 0.7, outlier.shape = 21, outlier.size = 1.5,
            colour = "grey30", width = 0.65
          ) +
          geom_jitter(width = 0.15, size = 1, alpha = 0.4, colour = "grey30") +
          scale_fill_brewer(palette = "YlOrRd", guide = "none") +
          labs(
            title = bquote(bold(.(gene)) ~ " Expression by AJCC Stage"),
            subtitle = bquote(
              "Kruskal-Wallis  " ~
                italic(p) == .(formatC(pval, format = "e", digits = 2)) ~ " | " ~
                italic(p)[adj] == .(formatC(padj, format = "e", digits = 2))
            ),
            x = "AJCC Stage", y = "logCPM (TMM-normalised)"
          ) +
          theme_pub

        if (has_ggpubr) {
          p_ind <- p_ind +
            stat_compare_means(
              method = "kruskal.test",
              label = "p.format",
              label.y.npc = "top", size = 3.5
            )
        }
      } else if (var == "Vital_Status") {
        plot_data <- gene_data %>%
          dplyr::filter(
            !is.na(vital_status),
            vital_status %in% c("Alive", "Dead")
          )
        if (nrow(plot_data) == 0) next

        p_ind <- ggplot(
          plot_data,
          aes(x = vital_status, y = logCPM, fill = vital_status)
        ) +
          geom_boxplot(
            alpha = 0.7, outlier.shape = 21, outlier.size = 1.5,
            colour = "grey30", width = 0.55
          ) +
          geom_jitter(width = 0.12, size = 1, alpha = 0.4, colour = "grey30") +
          scale_fill_manual(
            values = c("Alive" = "#4DAF4A", "Dead" = "#E41A1C"),
            guide = "none"
          ) +
          labs(
            title = bquote(bold(.(gene)) ~ " Expression by Vital Status"),
            subtitle = bquote(
              "Wilcoxon  " ~
                italic(p) == .(formatC(pval, format = "e", digits = 2)) ~ " | " ~
                italic(p)[adj] == .(formatC(padj, format = "e", digits = 2))
            ),
            x = "Vital Status", y = "logCPM (TMM-normalised)"
          ) +
          theme_pub

        if (has_ggpubr) {
          p_ind <- p_ind +
            stat_compare_means(
              method = "wilcox.test",
              label = "p.format",
              label.x = 1.5, label.y.npc = "top", size = 3.5
            )
        }
      } else if (var == "Gender") {
        plot_data <- gene_data %>%
          dplyr::filter(!is.na(gender), gender %in% c("male", "female"))
        if (nrow(plot_data) == 0) next

        p_ind <- ggplot(
          plot_data,
          aes(x = gender, y = logCPM, fill = gender)
        ) +
          geom_boxplot(
            alpha = 0.7, outlier.shape = 21, outlier.size = 1.5,
            colour = "grey30", width = 0.55
          ) +
          geom_jitter(width = 0.12, size = 1, alpha = 0.4, colour = "grey30") +
          scale_fill_manual(
            values = c("male" = "#377EB8", "female" = "#E41A1C"),
            guide = "none"
          ) +
          labs(
            title = bquote(bold(.(gene)) ~ " Expression by Gender"),
            subtitle = bquote(
              "Wilcoxon  " ~
                italic(p) == .(formatC(pval, format = "e", digits = 2)) ~ " | " ~
                italic(p)[adj] == .(formatC(padj, format = "e", digits = 2))
            ),
            x = "Gender", y = "logCPM (TMM-normalised)"
          ) +
          theme_pub

        if (has_ggpubr) {
          p_ind <- p_ind +
            stat_compare_means(
              method = "wilcox.test",
              label = "p.format",
              label.x = 1.5, label.y.npc = "top", size = 3.5
            )
        }
      } else if (var == "Age") {
        plot_data <- gene_data %>% dplyr::filter(!is.na(age))
        if (nrow(plot_data) < 10) next

        cor_res <- cor.test(plot_data$logCPM, plot_data$age, method = "spearman")

        p_ind <- ggplot(plot_data, aes(x = age, y = logCPM)) +
          geom_point(alpha = 0.5, size = 2, colour = "#377EB8") +
          geom_smooth(
            method = "lm", se = TRUE, colour = "#D73027",
            linewidth = 1, fill = "#FDD0A2"
          ) +
          labs(
            title = bquote(bold(.(gene)) ~ " Expression vs Age"),
            subtitle = bquote(
              "Spearman " ~ rho == .(round(cor_res$estimate, 3)) ~ " | " ~
                italic(p) == .(formatC(pval, format = "e", digits = 2)) ~ " | " ~
                italic(p)[adj] == .(formatC(padj, format = "e", digits = 2))
            ),
            x = "Age at diagnosis", y = "logCPM (TMM-normalised)"
          ) +
          theme_pub
      } else {
        next
      }

      # Save individual plot
      safe_gene <- gsub("[^A-Za-z0-9_]", "_", gene)
      safe_var <- gsub("[^A-Za-z0-9_]", "_", var)
      fname <- paste0(safe_gene, "_vs_", safe_var)

      ggsave(file.path(indiv_dir, paste0(fname, ".pdf")),
        plot = p_ind, width = 6, height = 5, device = cairo_pdf
      )
      ggsave(file.path(indiv_dir, paste0(fname, ".png")),
        plot = p_ind, width = 6, height = 5, dpi = 300
      )
      cat("  Saved:", fname, "(PDF + PNG)\n")
    }

    cat("Individual clinical plots saved to:", indiv_dir, "\n")
    cat("Total individual plots:", nrow(sig_nominal_tests), "\n")
  } else {
    cat("No nominally significant hub gene–clinical associations found.\n")
  }
} else {
  message("Skipping clinical expression plots (no matching expression data).")
}

# --- 6.5 Hub gene expression heatmap with clinical annotations ----------------
if (length(hub_genes_in_matrix) > 0 && length(sample_ids) > 0) {
  # Expression matrix for hub genes across tumour samples
  heat_genes <- intersect(hub_characterisation$gene_symbol, hub_genes_in_matrix)
  heat_mat <- as.matrix(logcpm[heat_genes, sample_ids, drop = FALSE])

  # Row-wise z-score scaling
  heat_mat_z <- t(scale(t(heat_mat)))
  heat_mat_z[is.na(heat_mat_z)] <- 0

  # Column annotation (clinical metadata)
  anno_col <- meta_tumor %>%
    dplyr::filter(File_ID %in% sample_ids) %>%
    dplyr::mutate(
      Stage = dplyr::case_when(
        grepl("Stage I$|Stage I[^V]", diagnoses.ajcc_pathologic_stage) ~ "I",
        grepl("Stage II$|Stage II[^I]", diagnoses.ajcc_pathologic_stage) ~ "II",
        grepl("Stage III", diagnoses.ajcc_pathologic_stage) ~ "III",
        grepl("Stage IV", diagnoses.ajcc_pathologic_stage) ~ "IV",
        TRUE ~ "NA"
      ),
      Vital_Status = demographic.vital_status,
      Gender = demographic.gender
    ) %>%
    dplyr::select(File_ID, Stage, Vital_Status, Gender) %>%
    tibble::column_to_rownames("File_ID")

  # Ensure column order matches
  anno_col <- anno_col[colnames(heat_mat_z), , drop = FALSE]

  # Row annotation (PIKK group + pathway role)
  anno_row <- hub_characterisation %>%
    dplyr::filter(gene_symbol %in% heat_genes) %>%
    dplyr::select(gene_symbol, PIKK_Group, pathway_role) %>%
    tibble::column_to_rownames("gene_symbol")
  anno_row <- anno_row[rownames(heat_mat_z), , drop = FALSE]

  # Annotation colours
  anno_colors <- list(
    PIKK_Group = pikk_colors,
    pathway_role = role_colors[unique(anno_row$pathway_role)],
    Stage = c(
      "I" = "#FFFFB2", "II" = "#FECC5C",
      "III" = "#FD8D3C", "IV" = "#E31A1C", "NA" = "grey80"
    ),
    Vital_Status = c("Alive" = "#4DAF4A", "Dead" = "#E41A1C"),
    Gender = c("male" = "#377EB8", "female" = "#E41A1C")
  )

  heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

  pdf(file.path(out_plots, "hub_gene_expression_heatmap_clinical.pdf"),
    width = 14, height = max(6, nrow(heat_mat_z) * 0.3 + 3)
  )
  pheatmap(heat_mat_z,
    color = heat_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = anno_col,
    annotation_row = anno_row,
    annotation_colors = anno_colors,
    show_colnames = FALSE,
    fontsize_row = 8,
    border_color = NA,
    clustering_method = "ward.D2",
    main = "Hub Gene Expression — OSCC Tumour Samples (z-scored logCPM)"
  )
  dev.off()

  png(file.path(out_plots, "hub_gene_expression_heatmap_clinical.png"),
    width = 4200,
    height = max(2400, nrow(heat_mat_z) * 100),
    res = 300
  )
  pheatmap(heat_mat_z,
    color = heat_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = anno_col,
    annotation_row = anno_row,
    annotation_colors = anno_colors,
    show_colnames = FALSE,
    fontsize_row = 8,
    border_color = NA,
    clustering_method = "ward.D2",
    main = "Hub Gene Expression — OSCC Tumour Samples (z-scored logCPM)"
  )
  dev.off()
  cat("  Saved: hub_gene_expression_heatmap_clinical (PDF + PNG)\n")
} else {
  message("Skipping clinical heatmap (no matching expression data).")
}


# --- 6.6 |logFC| vs degree bubble chart (colour = PIKK group, size = betweenness)
p_bubble <- ggplot(
  hub_characterisation,
  aes(x = degree, y = abs(log2FoldChange))
) +
  geom_point(aes(size = betweenness, fill = PIKK_Group),
    shape = 21, colour = "grey30", stroke = 0.5, alpha = 0.85
  ) +
  geom_text_repel(aes(label = gene_symbol),
    size = 3.2, fontface = "bold",
    max.overlaps = 25
  ) +
  scale_fill_manual(values = pikk_colors, name = "PIKK Group") +
  scale_size_continuous(
    range = c(3, 15),
    name = "Betweenness\nCentrality"
  ) +
  labs(
    title = "Hub Gene Multidimensional Profile",
    subtitle = "Degree × |log2FC| × Betweenness × PIKK Group",
    x = "Degree Centrality",
    y = expression("|" * log[2] * "FC|")
  ) +
  theme_pub
save_pdf_png(p_bubble, "hub_gene_logFC_vs_centrality_bubble",
  width = 11, height = 8
)


# --- 6.7 Pathway role classification barplot ----------------------------------
role_summary <- hub_characterisation %>%
  dplyr::count(pathway_role, regulation, name = "n") %>%
  dplyr::mutate(
    pathway_role = factor(pathway_role,
      levels = names(sort(table(hub_characterisation$pathway_role),
        decreasing = TRUE
      ))
    )
  )

p_role <- ggplot(role_summary, aes(x = pathway_role, y = n, fill = regulation)) +
  geom_col(colour = "grey30", width = 0.75) +
  scale_fill_manual(values = reg_colors, name = "Regulation") +
  labs(
    title = "Hub Gene Pathway Role Classification",
    x = NULL, y = "Number of Hub Genes"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"))
save_pdf_png(p_role, "hub_gene_role_classification_barplot",
  width = 10, height = 6
)


# --- 6.8 Hub gene network highlight ------------------------------------------
if (has_ggraph && file.exists(graphml_file)) {
  cat("Loading core network for hub gene highlight plot...\n")
  g_loaded <- read_graph(graphml_file, format = "graphml")

  # Check if the graph loaded correctly
  if (vcount(g_loaded) > 0) {
    # Mark hub genes
    V(g_loaded)$is_hub <- V(g_loaded)$name %in% hub_characterisation$gene_symbol

    # Compute layout
    set.seed(42)
    layout_coords <- layout_with_fr(g_loaded)
    V(g_loaded)$x <- layout_coords[, 1]
    V(g_loaded)$y <- layout_coords[, 2]

    # Get PIKK group for hub genes
    hub_pikk <- hub_characterisation$PIKK_Group
    names(hub_pikk) <- hub_characterisation$gene_symbol

    V(g_loaded)$hub_group <- ifelse(
      V(g_loaded)$is_hub,
      hub_pikk[V(g_loaded)$name],
      "Non-hub"
    )

    tg_loaded <- as_tbl_graph(g_loaded)

    hub_node_colors <- c(pikk_colors, "Non-hub" = "grey85")

    p_hub_net <- ggraph(tg_loaded, layout = "manual", x = x, y = y) +
      geom_edge_link(alpha = 0.15, colour = "grey60") +
      # Non-hub nodes (smaller, greyed out)
      geom_node_point(
        data = . %>% filter(!is_hub),
        aes(size = 3),
        fill = "grey85", shape = 21, colour = "grey60", stroke = 0.3
      ) +
      # Hub nodes (larger, coloured, with border)
      geom_node_point(
        data = . %>% filter(is_hub),
        aes(size = 8, fill = hub_group),
        shape = 21, colour = "black", stroke = 1.2
      ) +
      # Hub gene labels
      geom_node_text(
        data = . %>% filter(is_hub),
        aes(label = name),
        repel = TRUE, size = 3.5, fontface = "bold",
        max.overlaps = 30
      ) +
      scale_fill_manual(values = hub_node_colors, name = "PIKK Group") +
      scale_size_identity() +
      labs(
        title = paste0(
          "PPI Network — Hub Genes Highlighted (",
          sum(V(g_loaded)$is_hub), " hubs)"
        ),
        subtitle = paste0(
          vcount(g_loaded), " nodes, ",
          ecount(g_loaded), " edges"
        )
      ) +
      theme_graph(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right"
      ) +
      guides(fill = guide_legend(
        override.aes = list(size = 5),
        title.position = "top"
      ))
    save_pdf_png(p_hub_net, "hub_gene_network_highlight",
      width = 13, height = 10
    )
  } else {
    message("GraphML file loaded but has no vertices. Skipping hub network plot.")
  }
} else {
  if (!has_ggraph) {
    message("Skipping hub network plot (ggraph not installed).")
  } else {
    message(
      "Skipping hub network plot (GraphML file not found at: ",
      graphml_file, ")"
    )
    message("Run 05_ppi_network_analysis.R first to generate the network file.")
  }
}


# ============================================================================
# 7. SESSION SUMMARY
# ============================================================================
cat("\n================================================================\n")
cat("STEP 2.3 — HUB GENE IDENTIFICATION COMPLETE\n")
cat("================================================================\n")
cat("Hub genes identified (union set): ", nrow(hub_characterisation), "\n")
cat("  - Top by degree:       ", top_n_degree, "\n")
cat("  - Top by betweenness:  ", top_n_betweenness, "\n")
cat("  - Union (unique):      ", length(hub_genes_union), "\n")
cat("Clinical tests performed:           ", nrow(clinical_tests), "\n")
cat("Nominally significant (p < 0.05):   ", sum(clinical_tests$sig_nominal), "\n")
cat("FDR-significant (padj < 0.05):      ", sum(clinical_tests$sig_FDR), "\n")
cat("\nPathway role breakdown:\n")
print(table(hub_characterisation$pathway_role))
cat("\nTop-5 hub genes by composite score:\n")
print(hub_characterisation %>%
  dplyr::select(
    gene_symbol, PIKK_Group, pathway_role,
    composite_score, regulation
  ) %>%
  head(5))
cat("\nOutputs written to:\n")
cat("  Tables:", out_tables, "\n")
cat("  Plots: ", out_plots, "\n")

# Save session info
writeLines(
  capture.output(sessionInfo()),
  file.path(out_tables, "session_info_hub_gene_analysis.txt")
)

cat("\nSession info saved.\n")
cat("================================================================\n")
