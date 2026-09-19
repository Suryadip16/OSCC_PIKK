# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  PIKK Pathway Druggability Visualisation Suite in OSCC                       ║
# ║  Script: 02_pikk_druggability_plots.R                                       ║
# ║  Purpose: Generate publication-ready figures for DGIdb druggability results  ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

cat("\n════════════════════════════════════════════════════════════════\n")
cat("  PIKK Druggability Visualisation Pipeline\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("════════════════════════════════════════════════════════════════\n\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(RColorBrewer)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
base_dir <- "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project"
in_dir   <- file.path(base_dir, "05_results/localisation_outputs/tables")
fig_dir  <- file.path(base_dir, "05_results/localisation_outputs/figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

score_file <- file.path(in_dir, "pikk_druggability_composite_scores.tsv")
inter_file <- file.path(in_dir, "dgidb_raw_gene_interactions.tsv")

if (!file.exists(score_file)) {
  stop("Input file not found: ", score_file, "\nPlease run 01_pikk_druggability_analysis.py first!")
}

df_scores <- read_tsv(score_file, show_col_types = FALSE)

# Helper: Save Dual PDF + PNG
save_dual <- function(plot_obj, filename, w = 10, h = 8, dpi = 300) {
  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), plot_obj, width = w, height = h, device = "pdf")
  ggsave(file.path(fig_dir, paste0(filename, ".png")), plot_obj, width = w, height = h, dpi = dpi, device = "png")
  cat("    Saved:", filename, "(PDF + PNG)\n")
}

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 1: COMPOSITE DRUGGABILITY SCORE BREAKDOWN                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
cat("── Generating Figure 1: Composite Druggability Score Breakdown ────────\n")

# Prepare stacked data for Diamond Panel + Top Hubs
diamond_df <- df_scores %>%
  filter(is_diamond_panel == "Yes" | composite_druggability_score > 0.40) %>%
  arrange(composite_druggability_score)

diamond_df$gene_symbol <- factor(diamond_df$gene_symbol, levels = diamond_df$gene_symbol)

# Sub-signals weighted contributions
# Weights: Approval 0.35, Interaction 0.30, Diversity 0.20, Drugs 0.15
df_sub <- diamond_df %>%
  mutate(
    `Approval Status (35%)`     = approval_tier_score * 0.35,
    `Interaction Score (30%)`   = interaction_subscore * 0.30,
    `Source Diversity (20%)`    = source_diversity_subscore * 0.20,
    `Drug Count (15%)`          = drug_count_subscore * 0.15
  ) %>%
  select(gene_symbol, clinical_maturity, is_diamond_panel, composite_druggability_score,
         `Approval Status (35%)`, `Interaction Score (30%)`, `Source Diversity (20%)`, `Drug Count (15%)`) %>%
  pivot_longer(cols = c(`Approval Status (35%)`, `Interaction Score (30%)`, `Source Diversity (20%)`, `Drug Count (15%)`),
               names_to = "SubSignal", values_to = "Contribution")

# Color palette for sub-signals
signal_colors <- c(
  "Approval Status (35%)"   = "#E64B35",
  "Interaction Score (30%)" = "#4DBBD5",
  "Source Diversity (20%)"  = "#00A087",
  "Drug Count (15%)"        = "#F39B7F"
)

p1 <- ggplot(df_sub, aes(x = gene_symbol, y = Contribution, fill = SubSignal)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = signal_colors, name = "Druggability Sub-Signals") +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1.0, 0.2)) +
  geom_text(data = diamond_df,
            aes(x = gene_symbol, y = composite_druggability_score + 0.03,
                label = sprintf("%.2f (%s)", composite_druggability_score, clinical_maturity)),
            inherit.aes = FALSE, hjust = 0, size = 3.6, fontface = "bold") +
  labs(
    title = "PIKK Target Druggability & Therapeutic Actionability",
    subtitle = "Multi-Signal DGIdb v5.0 Decomposition (Approval, Evidence Confidence, Database Breadth)",
    x = NULL,
    y = "Composite Druggability Score (0.0 – 1.0)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey30", size = 11),
    axis.text.y = element_text(face = "bold", size = 12)
  )

save_dual(p1, "01_druggability_composite_score_breakdown", w = 11, h = 7)

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 2: TARGET-DRUG INTERACTION NETWORK (CLEAN BIPARTITE PLOT)          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
cat("── Generating Figure 2: Target-Drug Interaction Network ─────────────────\n")

# Helper function to clean and filter drug names
clean_drug_name <- function(dname) {
  dname <- trimws(dname)
  # Remove trailing salt phrases
  dname <- gsub(" (HYDROCHLORIDE|DIHYDROCHLORIDE|DISODIUM|SODIUM|MALEATE|MESYLATE|SULFATE|PROTEIN-BOUND|ANHYDROUS|TOSYLATE|CAMSYLATE)$", "", dname, ignore.case = TRUE)
  return(dname)
}

# Curate top informative oncology & investigational drugs per Diamond target
# (Prioritising approved antineoplastics, clinical trial inhibitors, and clean chemical probes)
diamond_drugs_list <- list()

for (i in 1:nrow(df_scores)) {
  row <- df_scores[i, ]
  gene <- row$gene_symbol
  
  if (row$is_diamond_panel != "Yes" || row$n_drugs == 0) next
  
  app_drugs <- unlist(strsplit(row$approved_drugs, "; "))
  cand_drugs <- unlist(strsplit(row$top_candidate_drugs, "; "))
  
  # Combine & clean
  all_d <- unique(c(app_drugs, cand_drugs))
  all_d <- all_d[all_d != "" & !is.na(all_d)]
  
  # Filter out raw CHEMBL IDs and overly long IUPAC chemical strings (>30 chars or complex brackets)
  clean_d <- all_d[!grepl("^CHEMBL:", all_d, ignore.case = TRUE)]
  clean_d <- clean_d[!grepl("^\\([0-9A-Z]{2,}", clean_d)] # filters long (7S)-... IUPAC
  clean_d <- clean_d[nchar(clean_d) <= 30]
  clean_d <- sapply(clean_d, clean_drug_name)
  clean_d <- unique(clean_d[clean_d != ""])
  
  # Take top 4-6 most relevant drugs
  top_d <- head(clean_d, 6)
  
  for (d in top_d) {
    is_app <- ifelse(any(sapply(app_drugs, function(x) grepl(d, x, fixed = TRUE))),
                     "FDA Approved", "Investigational / Trial")
    diamond_drugs_list[[length(diamond_drugs_list) + 1]] <- data.frame(
      gene_symbol = gene,
      drug_name = d,
      status = is_app,
      stringsAsFactors = FALSE
    )
  }
}

drug_targets_df <- bind_rows(diamond_drugs_list)

# If any target has 0 drugs (e.g. TOPBP1, FANCI), add a placeholder so the panel is comprehensive
diamond_genes_all <- c("PLK1", "CDK2", "RAD51", "KAT2B", "DEPTOR", "TOPBP1", "FANCI")
missing_genes <- setdiff(diamond_genes_all, drug_targets_df$gene_symbol)
for (mg in missing_genes) {
  drug_targets_df <- rbind(drug_targets_df, data.frame(
    gene_symbol = mg,
    drug_name = "No Direct Small Molecule (Novel Target)",
    status = "Undruggable / Novel",
    stringsAsFactors = FALSE
  ))
}

drug_targets_df$gene_symbol <- factor(drug_targets_df$gene_symbol, levels = diamond_genes_all)

# Publication-grade faceted lollipop plot for target-drug interactions
p2 <- ggplot(drug_targets_df, aes(y = drug_name, x = 1, color = status)) +
  geom_segment(aes(y = drug_name, yend = drug_name, x = 0, xend = 1),
               linetype = "solid", color = "grey82", linewidth = 0.8) +
  geom_point(size = 4.2, shape = 19) +
  geom_text(aes(label = status, x = 1.08), hjust = 0, size = 3.1, fontface = "bold") +
  facet_wrap(~ gene_symbol, scales = "free_y", ncol = 3) +
  scale_x_continuous(limits = c(0, 2.4), breaks = NULL) +
  scale_color_manual(
    values = c(
      "FDA Approved"              = "#E64B35",
      "Investigational / Trial"   = "#3C5488",
      "Undruggable / Novel"       = "#8491B4"
    ),
    name = "Clinical Regulatory Status"
  ) +
  labs(
    title = "Therapeutic Tractability: Lead Small-Molecule & Repurposing Candidates",
    subtitle = "DGIdb v5.0 Curated Interactions for the PIKK Diamond Panel Targets in OSCC",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "#F4F6F9", color = "grey75"),
    strip.text = element_text(face = "bold", size = 12, color = "#1B263B"),
    axis.text.y = element_text(face = "bold", size = 9.5, color = "grey20"),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = 5)),
    plot.subtitle = element_text(color = "grey35", size = 10.5, margin = margin(b = 15)),
    plot.margin = margin(t = 20, r = 25, b = 20, l = 25)
  )

save_dual(p2, "02_drug_target_interaction_network_bipartite", w = 12, h = 9)

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 3: THERAPEUTIC ACTIONABILITY MATRIX                               ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
cat("── Generating Figure 3: Therapeutic Actionability Heatmap ──────────────\n")

matrix_df <- df_scores %>%
  filter(is_diamond_panel == "Yes") %>%
  select(gene_symbol, clinical_maturity, composite_druggability_score, n_drugs, n_approved_drugs) %>%
  mutate(
    `Direct Small Molecule Available` = ifelse(n_drugs > 0, "Yes", "No"),
    `FDA Approved Drug Exists`       = ifelse(n_approved_drugs > 0, "Yes", "No"),
    `Clinical Phase Tier`            = clinical_maturity
  ) %>%
  pivot_longer(cols = c(`Direct Small Molecule Available`, `FDA Approved Drug Exists`, `Clinical Phase Tier`),
               names_to = "Feature", values_to = "Status")

p3 <- ggplot(matrix_df, aes(x = Feature, y = gene_symbol, fill = Status)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = Status), color = "black", fontface = "bold", size = 3.8) +
  scale_fill_brewer(palette = "Pastel1", name = "Clinical Evidence") +
  labs(
    title = "Translational Actionability Matrix: PIKK Diamond Panel",
    subtitle = "Assessment of Direct Drug Tractability in Oncology",
    x = NULL,
    y = "Gene Target"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold")
  )

save_dual(p3, "03_therapeutic_actionability_heatmap", w = 9, h = 6)

cat("\n════════════════════════════════════════════════════════════════\n")
cat("  All 3 Druggability Figures Generated Successfully!\n")
cat("════════════════════════════════════════════════════════════════\n\n")
