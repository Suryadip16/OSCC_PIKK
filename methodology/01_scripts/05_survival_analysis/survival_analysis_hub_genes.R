################################################################################
##
##  OSCC PIKK Project — Section 2.4  Survival Association Analysis
##
##  Scope of work deliverables
##  ─────────────────────────────────────────────────────────────────
##  • Kaplan–Meier curves per hub gene (median-split or optimal cutpoint
##    via survminer); log-rank test
##  • Univariate Cox proportional hazards regression for each hub gene
##  • Multivariate Cox regression incorporating relevant clinical
##    covariates (age, stage, grade)
##  • Forest plots summarising hazard ratios and 95 % CIs
##
##  Inputs
##  ──────
##  1. matched_metadata_used.tsv      — clinical metadata (per File_ID)
##  2. logCPM_matrix.tsv              — log-CPM normalised expression
##  3. hub_gene_characterisation_table.tsv — hub gene list
##
##  Outputs  → 05_results/survival_outputs/
##  ────────
##  ├── km_plots/                      Individual KM plots (PDF & PNG) per hub gene
##  ├── km_combined_grid.pdf/.png      Multi-panel KM grid
##  ├── univariate_cox_results.tsv     Table of univariate Cox results
##  ├── univariate_forest_plot.pdf/.png Forest plot (univariate)
##  ├── multivariate_cox_results.tsv   Table of multivariate Cox results
##  ├── multivariate_forest_plot.pdf/.png Forest plot (multivariate)
##  ├── concordance_summary.tsv        C-index for each model
##  ├── concordance_comparison.pdf/.png C-index comparison plot
##  ├── multivariate_full_forest_best_gene.pdf/.png Detailed full forest plot
##  ├── cox_ph_diagnostics.pdf/.png    Schoenfeld residual plots
##  ├── survival_significance_heatmap.pdf/.png Summary heatmap
##  └── survival_summary_table.tsv     Combined survival summary
##
################################################################################

# ── 0. Setup ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(patchwork)
  library(scales)
  library(RColorBrewer)
  library(forestmodel)
})

# ── Paths ------------------------------------------------------------------
project_root <- file.path("c:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project")

meta_path   <- file.path(project_root,
                          "05_results/edgeR_withBatch/results/matched_metadata_used.tsv")
expr_path   <- file.path(project_root,
                          "05_results/edgeR_withBatch/results/logCPM_matrix.tsv")
hub_path    <- file.path(project_root,
                          "05_results/tables/hub_genes/hub_gene_characterisation_table.tsv")

out_dir     <- file.path(project_root, "05_results/survival_outputs")
km_dir      <- file.path(out_dir, "km_plots")
dir.create(km_dir, recursive = TRUE, showWarnings = FALSE)

# ── Publication theme -------------------------------------------------------
theme_pub <- theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle    = element_text(size = 10, hjust = 0.5, colour = "grey40"),
    axis.title       = element_text(face = "bold", size = 11),
    axis.text        = element_text(size = 10),
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 9),
    strip.text       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(8, 12, 8, 12)
  )

# Colour palette for two-group KM
km_cols <- c("High" = "#D7191C", "Low" = "#2B83BA")

cat("════════════════════════════════════════════════════════════════\n")
cat("  OSCC PIKK Hub Gene — Survival Association Analysis\n")
cat("════════════════════════════════════════════════════════════════\n\n")

################################################################################
## 1. LOAD & PREPARE DATA
################################################################################

cat("▸ Loading metadata …\n")
meta_raw <- read_tsv(meta_path, show_col_types = FALSE)

cat("▸ Loading hub gene list …\n")
hub_genes_df <- read_tsv(hub_path, show_col_types = FALSE) %>%
  filter(!is.na(gene_symbol) & gene_symbol != "")

hub_genes <- hub_genes_df$gene_symbol
cat("  Hub genes (n =", length(hub_genes), "):",
    paste(hub_genes, collapse = ", "), "\n\n")

cat("▸ Loading logCPM expression matrix …\n")
expr_raw <- read_tsv(expr_path, show_col_types = FALSE)

# ── 1a. Prepare clinical data ------------------------------------------------
# Keep ONLY tumour samples (Sample_Type_Code == 1) and deduplicate per patient
# Survival endpoint: Overall Survival (OS)
#   vital_status → event indicator (Dead = 1, Alive = 0)
#   Dead patients  → time = demographic.days_to_death
#   Alive patients → time = diagnoses.days_to_last_follow_up (censored)

meta <- meta_raw %>%
  filter(Condition == "Tumor") %>%                         # tumour samples only
  distinct(`Case ID`, .keep_all = TRUE) %>%                # one row per patient
  mutate(
    # Event indicator
    os_event = case_when(
      demographic.vital_status == "Dead"  ~ 1L,
      demographic.vital_status == "Alive" ~ 0L,
      TRUE ~ NA_integer_
    ),
    # Overall survival time (days): days_to_death if Dead, days_to_last_follow_up if Alive
    os_time_days = case_when(
      demographic.vital_status == "Dead"  ~ as.numeric(demographic.days_to_death),
      demographic.vital_status == "Alive" ~ as.numeric(diagnoses.days_to_last_follow_up),
      TRUE ~ NA_real_
    ),
    # Age (years)
    age = as.numeric(demographic.age_at_index),
    # Gender
    gender = demographic.gender,
    # AJCC pathologic stage — simplify to I/II vs III/IV
    stage_raw = diagnoses.ajcc_pathologic_stage,
    stage_group = case_when(
      str_detect(stage_raw, "^Stage I$|^Stage II$|^Stage I[^V]") ~ "Early (I-II)",
      str_detect(stage_raw, "^Stage III|^Stage IV")              ~ "Advanced (III-IV)",
      TRUE ~ NA_character_
    ),
    # Histologic tumor grade — simplify to Low Grade (G1-G2) vs High Grade (G3-G4)
    grade_raw = diagnoses.tumor_grade,
    grade_group = case_when(
      grade_raw %in% c("G1", "G2") ~ "Low Grade (G1-G2)",
      grade_raw %in% c("G3", "G4") ~ "High Grade (G3-G4)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(os_event))

n_alive_no_time <- sum(meta$os_event == 0 & (is.na(meta$os_time_days) | meta$os_time_days <= 0))
n_dead_no_time  <- sum(meta$os_event == 1 & (is.na(meta$os_time_days) | meta$os_time_days <= 0))

cat("  Alive patients without valid follow-up time:", n_alive_no_time, "\n")
cat("  Dead patients without valid death time:", n_dead_no_time, "\n")

# Keep patients with known positive survival / follow-up time
meta_surv <- meta %>%
  filter(!is.na(os_time_days) & os_time_days > 0)

# Convert to years for cleaner plots
meta_surv <- meta_surv %>%
  mutate(os_time_years = os_time_days / 365.25)

cat("  Patients with complete OS data (Dead + Censored Alive):", nrow(meta_surv), "\n\n")

# ── 1b. Extract hub-gene expression for survival cohort ----------------------
# Expression matrix: rows = genes, columns = File_IDs
# Map File_IDs back to patients
expr_mat <- expr_raw %>%
  filter(gene_name %in% hub_genes) %>%
  tibble::column_to_rownames("gene_name")

# Identify File_IDs present in survival cohort
sample_ids <- meta_surv$File_ID
common_ids <- intersect(sample_ids, colnames(expr_mat))
cat("  Samples with both expression + survival:", length(common_ids), "\n\n")

# Subset
expr_sub <- expr_mat[, common_ids, drop = FALSE]
meta_surv <- meta_surv %>% filter(File_ID %in% common_ids)

# Transpose to samples × genes
expr_t <- as.data.frame(t(expr_sub))
expr_t$File_ID <- rownames(expr_t)

# Merge into single analysis data frame
surv_df <- meta_surv %>%
  left_join(expr_t, by = "File_ID")
surv_df <- as.data.frame(surv_df)

# Verify hub genes present
genes_found <- intersect(hub_genes, colnames(surv_df))
genes_missing <- setdiff(hub_genes, colnames(surv_df))
if (length(genes_missing) > 0) {
  cat("  ⚠ Hub genes missing from expression matrix:",
      paste(genes_missing, collapse = ", "), "\n")
}
cat("  Hub genes available for survival analysis (n =", length(genes_found),
    "):", paste(genes_found, collapse = ", "), "\n\n")

################################################################################
## 2. KAPLAN–MEIER SURVIVAL CURVES  (per hub gene)
################################################################################

cat("════════════════════════════════════════════════════════════════\n")
cat("  2. Kaplan–Meier Curves\n")
cat("════════════════════════════════════════════════════════════════\n\n")

km_results <- list()
km_plots   <- list()

for (gene in genes_found) {

  cat("  ▸ KM analysis:", gene, "… ")

  # --- Optimal cutpoint via surv_cutpoint (maximally selected rank statistics)
  cp_df <- as.data.frame(surv_df[, c("os_time_years", "os_event", gene)])
  cp_df <- cp_df[!is.na(cp_df[[gene]]), ]

  # Try optimal cutpoint; fall back to median if it fails
  cut_result <- tryCatch({
    sc <- surv_cutpoint(cp_df,
                        time    = "os_time_years",
                        event   = "os_event",
                        variables = gene,
                        minprop = 0.25)
    list(cutpoint = sc$cutpoint[[gene, "cutpoint"]],
         method   = "optimal (maxstat)")
  }, error = function(e) {
    list(cutpoint = median(cp_df[[gene]], na.rm = TRUE),
         method   = "median")
  })

  cutval <- cut_result$cutpoint
  method <- cut_result$method

  # Dichotomise
  grp_var <- paste0(gene, "_group")
  surv_df[[grp_var]] <- factor(
    ifelse(surv_df[[gene]] >= cutval, "High", "Low"),
    levels = c("Low", "High")
  )
  surv_df$strata_group <- surv_df[[grp_var]]

  # Fit KM using direct literal formula (avoids survminer symbol subsetting bug)
  km_fit     <- survfit(Surv(os_time_years, os_event) ~ strata_group, data = surv_df)
  lr_test    <- survdiff(Surv(os_time_years, os_event) ~ strata_group, data = surv_df)
  lr_pval    <- 1 - pchisq(lr_test$chisq, df = 1)

  # Store
  km_results[[gene]] <- data.frame(
    gene        = gene,
    cutpoint    = round(cutval, 3),
    cut_method  = method,
    n_high      = sum(surv_df$strata_group == "High", na.rm = TRUE),
    n_low       = sum(surv_df$strata_group == "Low", na.rm = TRUE),
    logrank_chi = round(lr_test$chisq, 3),
    logrank_p   = lr_pval,
    stringsAsFactors = FALSE
  )

  # ── Publication-quality KM plot ──
  p <- ggsurvplot(
    km_fit,
    data        = surv_df,
    pval        = TRUE,
    pval.size   = 4,
    pval.coord  = c(0.05, 0.05),
    conf.int    = TRUE,
    conf.int.alpha = 0.15,
    risk.table  = TRUE,
    risk.table.col = "strata",
    risk.table.height = 0.28,
    risk.table.y.text = FALSE,
    ncensor.plot = FALSE,
    palette     = unname(km_cols),
    xlab        = "Time (years)",
    ylab        = "Overall Survival Probability",
    title       = paste0(gene, " Expression & Overall Survival"),
    subtitle    = paste0("Cutpoint: ", round(cutval, 2),
                         " logCPM (", method, ")"),
    legend.title = paste0(gene, " Expression"),
    legend.labs = c("Low", "High"),
    surv.median.line = "hv",
    ggtheme     = theme_pub,
    font.main   = c(14, "bold"),
    font.x      = c(11, "bold"),
    font.y      = c(11, "bold"),
    font.tickslab = c(10, "plain"),
    font.legend = c(10, "plain"),
    linetype    = "strata",
    break.time.by = 1
  )

  km_plots[[gene]] <- p

  # Save individual PDF
  pdf(file.path(km_dir, paste0("KM_", gene, ".pdf")),
      width = 7, height = 7)
  print(p)
  dev.off()

  # Save individual PNG (300 dpi)
  png(file.path(km_dir, paste0("KM_", gene, ".png")),
      width = 7, height = 7, units = "in", res = 300)
  print(p)
  dev.off()

  cat("done (p =", formatC(lr_pval, format = "e", digits = 2), ")\n")
}

# ── KM results summary table ─────────────────────────────────────────────────
km_summary <- bind_rows(km_results) %>%
  mutate(
    logrank_p_adj = p.adjust(logrank_p, method = "BH"),
    significance  = case_when(
      logrank_p_adj < 0.001 ~ "***",
      logrank_p_adj < 0.01  ~ "**",
      logrank_p_adj < 0.05  ~ "*",
      logrank_p_adj < 0.1   ~ "†",
      TRUE                  ~ "ns"
    )
  ) %>%
  arrange(logrank_p)

write_tsv(km_summary, file.path(out_dir, "km_logrank_summary.tsv"))
cat("\n  KM summary saved → km_logrank_summary.tsv\n")

# ── Combined multi-panel KM grid (significant genes first) ───────────────────
# Create a combined grid of KM plots for the top genes
sig_genes <- km_summary %>%
  filter(logrank_p < 0.1) %>%
  pull(gene)

if (length(sig_genes) == 0) {
  # If no significant genes, use top 6 by p-value
  sig_genes <- head(km_summary$gene, min(6, nrow(km_summary)))
}

# Determine grid dimensions
n_plots <- min(length(sig_genes), 9)
ncols   <- min(3, n_plots)
nrows   <- ceiling(n_plots / ncols)

# Build list of ggsurvplot objects for the grid
grid_plots <- list()
for (i in seq_len(n_plots)) {
  gene <- sig_genes[i]
  grp_var <- paste0(gene, "_group")
  surv_df$strata_group <- surv_df[[grp_var]]
  km_fit <- survfit(Surv(os_time_years, os_event) ~ strata_group, data = surv_df)

  pval_txt <- km_summary %>% filter(gene == !!gene) %>% pull(logrank_p)

  gp <- ggsurvplot(
    km_fit,
    data        = surv_df,
    pval        = TRUE,
    pval.size   = 3.5,
    conf.int    = FALSE,
    risk.table  = FALSE,
    palette     = unname(km_cols),
    xlab        = "Time (years)",
    ylab        = "OS Probability",
    title       = gene,
    legend.title = "",
    legend.labs = c("Low", "High"),
    ggtheme     = theme_pub +
      theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)),
    break.time.by = 2
  )

  grid_plots[[i]] <- gp$plot
}

p_combined <- wrap_plots(grid_plots, ncol = ncols, nrow = nrows) +
  plot_annotation(
    title = "Kaplan–Meier Survival Curves — OSCC PIKK Hub Genes",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

ggsave(file.path(out_dir, "km_combined_grid.pdf"),
       p_combined,
       width = 5 * ncols, height = 4.5 * nrows, dpi = 300)

ggsave(file.path(out_dir, "km_combined_grid.png"),
       p_combined,
       width = 5 * ncols, height = 4.5 * nrows, dpi = 300)
cat("  Combined KM grid saved → km_combined_grid.pdf & .png\n\n")

################################################################################
## 3. UNIVARIATE COX PROPORTIONAL HAZARDS REGRESSION
################################################################################

cat("════════════════════════════════════════════════════════════════\n")
cat("  3. Univariate Cox Regression\n")
cat("════════════════════════════════════════════════════════════════\n\n")

uni_cox_results <- list()

for (gene in genes_found) {

  # Assign expression variable cleanly
  surv_df$gene_expr <- surv_df[[gene]]
  cox_fit <- coxph(Surv(os_time_years, os_event) ~ gene_expr, data = surv_df)
  cox_sum <- summary(cox_fit)

  uni_cox_results[[gene]] <- data.frame(
    gene          = gene,
    coef          = cox_sum$coefficients["gene_expr", "coef"],
    HR            = cox_sum$coefficients["gene_expr", "exp(coef)"],
    HR_lower_95   = cox_sum$conf.int["gene_expr", "lower .95"],
    HR_upper_95   = cox_sum$conf.int["gene_expr", "upper .95"],
    se_coef       = cox_sum$coefficients["gene_expr", "se(coef)"],
    z_score       = cox_sum$coefficients["gene_expr", "z"],
    p_value       = cox_sum$coefficients["gene_expr", "Pr(>|z|)"],
    concordance   = cox_sum$concordance["C"],
    concordance_se = cox_sum$concordance["se(C)"],
    n             = cox_sum$n,
    n_events      = cox_sum$nevent,
    stringsAsFactors = FALSE
  )

  cat("  ▸", gene,
      "  HR =", round(cox_sum$coefficients["gene_expr", "exp(coef)"], 3),
      " (", round(cox_sum$conf.int["gene_expr", "lower .95"], 3), "–",
      round(cox_sum$conf.int["gene_expr", "upper .95"], 3), ")",
      "  p =", formatC(cox_sum$coefficients["gene_expr", "Pr(>|z|)"],
                        format = "e", digits = 2), "\n")
}

uni_cox_df <- bind_rows(uni_cox_results) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      p_adj < 0.1   ~ "†",
      TRUE          ~ "ns"
    )
  ) %>%
  arrange(p_value)

write_tsv(uni_cox_df, file.path(out_dir, "univariate_cox_results.tsv"))
cat("\n  Univariate Cox results saved → univariate_cox_results.tsv\n\n")

# ── Univariate forest plot ───────────────────────────────────────────────────
# Build a combined coxph model with all hub genes (one at a time) for
# forest-plot aesthetics. We'll construct a manual forest plot.

fp_data <- uni_cox_df %>%
  mutate(
    gene = factor(gene, levels = rev(gene)),   # order by p-value (ascending)
    label = paste0(
      formatC(HR, format = "f", digits = 2), " (",
      formatC(HR_lower_95, format = "f", digits = 2), "–",
      formatC(HR_upper_95, format = "f", digits = 2), ")"
    ),
    p_label = ifelse(p_value < 0.001,
                     formatC(p_value, format = "e", digits = 1),
                     formatC(p_value, format = "f", digits = 3)),
    sig_col = ifelse(p_value < 0.05, "Significant", "Non-significant")
  )

p_forest_uni <- ggplot(fp_data, aes(x = HR, y = gene)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50",
             linewidth = 0.6) +
  geom_errorbarh(aes(xmin = HR_lower_95, xmax = HR_upper_95,
                      colour = sig_col),
                 height = 0.25, linewidth = 0.7) +
  geom_point(aes(colour = sig_col, size = -log10(p_value)),
             shape = 18) +
  scale_colour_manual(values = c("Significant"     = "#D7191C",
                                 "Non-significant" = "#4575B4"),
                      name = "") +
  scale_size_continuous(range = c(3, 7), name = expression(-log[10](p))) +
  # Annotate HR (95% CI) to the right
  geom_text(aes(x = max(HR_upper_95) * 1.15, label = label),
            hjust = 0, size = 3, colour = "grey30") +
  geom_text(aes(x = max(HR_upper_95) * 1.7, label = p_label),
            hjust = 0, size = 3, colour = "grey30") +
  # Column headers
  annotate("text", x = max(fp_data$HR_upper_95) * 1.15,
           y = nrow(fp_data) + 0.8,
           label = "HR (95% CI)", fontface = "bold", size = 3.5, hjust = 0) +
  annotate("text", x = max(fp_data$HR_upper_95) * 1.7,
           y = nrow(fp_data) + 0.8,
           label = "P-value", fontface = "bold", size = 3.5, hjust = 0) +
  coord_cartesian(clip = "off",
                  xlim = c(min(fp_data$HR_lower_95) * 0.85,
                           max(fp_data$HR_upper_95) * 2.3)) +
  labs(
    title    = "Univariate Cox Regression — OSCC PIKK Hub Genes",
    subtitle = "Hazard Ratios for Overall Survival (per unit logCPM increase)",
    x        = "Hazard Ratio",
    y        = NULL
  ) +
  theme_pub +
  theme(
    plot.margin = margin(10, 100, 10, 10),
    axis.text.y = element_text(face = "bold", size = 10)
  )

ggsave(file.path(out_dir, "univariate_forest_plot.pdf"),
       p_forest_uni,
       width = 11, height = max(4, 0.45 * nrow(fp_data) + 2), dpi = 300)
ggsave(file.path(out_dir, "univariate_forest_plot.png"),
       p_forest_uni,
       width = 11, height = max(4, 0.45 * nrow(fp_data) + 2), dpi = 300)
cat("  Univariate forest plot saved → univariate_forest_plot.pdf & .png\n\n")


################################################################################
## 4. MULTIVARIATE COX REGRESSION
################################################################################

cat("════════════════════════════════════════════════════════════════\n")
cat("  4. Multivariate Cox Regression\n")
cat("════════════════════════════════════════════════════════════════\n\n")

# ── Covariates for multivariate model ─────────────────────────────────────────
# Per scope of work: "incorporating relevant clinical covariates (age, stage,
# grade)".
#   • age            (continuous, years at diagnosis)
#   • stage_group    (Early I-II vs Advanced III-IV)
#   • grade_group    (Low Grade G1-G2 vs High Grade G3-G4)
#   • gender         (male vs female)
#
# We build one multivariate model per hub gene, each adjusting for the
# same clinical covariates.

# Clean covariates
surv_df_mv <- surv_df %>%
  filter(!is.na(age) & !is.na(stage_group) & !is.na(grade_group) & !is.na(gender)) %>%
  mutate(
    gender      = factor(gender, levels = c("male", "female")),
    stage_group = factor(stage_group, levels = c("Early (I-II)", "Advanced (III-IV)")),
    grade_group = factor(grade_group, levels = c("Low Grade (G1-G2)", "High Grade (G3-G4)"))
  )
surv_df_mv <- as.data.frame(surv_df_mv)

cat("  Multivariate cohort size (complete covariates):", nrow(surv_df_mv), "\n")
cat("  Covariates: age (continuous), stage_group (Early/Advanced), grade_group (Low/High), gender\n\n")

multi_cox_results <- list()

for (gene in genes_found) {

  surv_df_mv$gene_expr <- surv_df_mv[[gene]]

  cox_mv <- tryCatch(
    coxph(Surv(os_time_years, os_event) ~ gene_expr + age + stage_group + grade_group + gender,
          data = surv_df_mv),
    error = function(e) NULL
  )

  if (is.null(cox_mv)) {
    cat("  ⚠ Multivariate model failed for", gene, "\n")
    next
  }

  cox_mv_sum <- summary(cox_mv)

  multi_cox_results[[gene]] <- data.frame(
    gene            = gene,
    coef            = cox_mv_sum$coefficients["gene_expr", "coef"],
    HR              = cox_mv_sum$coefficients["gene_expr", "exp(coef)"],
    HR_lower_95     = cox_mv_sum$conf.int["gene_expr", "lower .95"],
    HR_upper_95     = cox_mv_sum$conf.int["gene_expr", "upper .95"],
    se_coef         = cox_mv_sum$coefficients["gene_expr", "se(coef)"],
    z_score         = cox_mv_sum$coefficients["gene_expr", "z"],
    p_value         = cox_mv_sum$coefficients["gene_expr", "Pr(>|z|)"],
    concordance     = cox_mv_sum$concordance["C"],
    concordance_se  = cox_mv_sum$concordance["se(C)"],
    n               = cox_mv_sum$n,
    n_events        = cox_mv_sum$nevent,
    age_HR          = cox_mv_sum$coefficients["age", "exp(coef)"],
    age_p           = cox_mv_sum$coefficients["age", "Pr(>|z|)"],
    stage_HR        = cox_mv_sum$coefficients[
      grep("stage_group", rownames(cox_mv_sum$coefficients)), "exp(coef)"],
    stage_p         = cox_mv_sum$coefficients[
      grep("stage_group", rownames(cox_mv_sum$coefficients)), "Pr(>|z|)"],
    grade_HR        = cox_mv_sum$coefficients[
      grep("grade_group", rownames(cox_mv_sum$coefficients)), "exp(coef)"],
    grade_p         = cox_mv_sum$coefficients[
      grep("grade_group", rownames(cox_mv_sum$coefficients)), "Pr(>|z|)"],
    gender_HR       = cox_mv_sum$coefficients[
      grep("gender", rownames(cox_mv_sum$coefficients)), "exp(coef)"],
    gender_p        = cox_mv_sum$coefficients[
      grep("gender", rownames(cox_mv_sum$coefficients)), "Pr(>|z|)"],
    model_logtest_p = cox_mv_sum$logtest["pvalue"],
    stringsAsFactors = FALSE
  )

  cat("  ▸", gene,
      "  adj.HR =", round(cox_mv_sum$coefficients["gene_expr", "exp(coef)"], 3),
      "  p =", formatC(cox_mv_sum$coefficients["gene_expr", "Pr(>|z|)"],
                        format = "e", digits = 2), "\n")
}

multi_cox_df <- bind_rows(multi_cox_results) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      p_adj < 0.1   ~ "†",
      TRUE          ~ "ns"
    )
  ) %>%
  arrange(p_value)

write_tsv(multi_cox_df, file.path(out_dir, "multivariate_cox_results.tsv"))
cat("\n  Multivariate Cox results saved → multivariate_cox_results.tsv\n\n")

# ── Multivariate forest plot ─────────────────────────────────────────────────

fp_mv <- multi_cox_df %>%
  mutate(
    gene = factor(gene, levels = rev(gene)),
    label = paste0(
      formatC(HR, format = "f", digits = 2), " (",
      formatC(HR_lower_95, format = "f", digits = 2), "–",
      formatC(HR_upper_95, format = "f", digits = 2), ")"
    ),
    p_label = ifelse(p_value < 0.001,
                     formatC(p_value, format = "e", digits = 1),
                     formatC(p_value, format = "f", digits = 3)),
    sig_col = ifelse(p_value < 0.05, "Significant", "Non-significant")
  )

p_forest_mv <- ggplot(fp_mv, aes(x = HR, y = gene)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50",
             linewidth = 0.6) +
  geom_errorbarh(aes(xmin = HR_lower_95, xmax = HR_upper_95,
                      colour = sig_col),
                 height = 0.25, linewidth = 0.7) +
  geom_point(aes(colour = sig_col, size = -log10(p_value)),
             shape = 18) +
  scale_colour_manual(values = c("Significant"     = "#D7191C",
                                 "Non-significant" = "#4575B4"),
                      name = "") +
  scale_size_continuous(range = c(3, 7), name = expression(-log[10](p))) +
  geom_text(aes(x = max(HR_upper_95) * 1.15, label = label),
            hjust = 0, size = 3, colour = "grey30") +
  geom_text(aes(x = max(HR_upper_95) * 1.7, label = p_label),
            hjust = 0, size = 3, colour = "grey30") +
  annotate("text", x = max(fp_mv$HR_upper_95) * 1.15,
           y = nrow(fp_mv) + 0.8,
           label = "adj.HR (95% CI)", fontface = "bold", size = 3.5, hjust = 0) +
  annotate("text", x = max(fp_mv$HR_upper_95) * 1.7,
           y = nrow(fp_mv) + 0.8,
           label = "P-value", fontface = "bold", size = 3.5, hjust = 0) +
  coord_cartesian(clip = "off",
                  xlim = c(min(fp_mv$HR_lower_95) * 0.85,
                           max(fp_mv$HR_upper_95) * 2.3)) +
  labs(
    title    = "Multivariate Cox Regression — OSCC PIKK Hub Genes",
    subtitle = "Adjusted HR for OS (covariates: age, stage, grade, gender)",
    x        = "Hazard Ratio (adjusted)",
    y        = NULL
  ) +
  theme_pub +
  theme(
    plot.margin = margin(10, 100, 10, 10),
    axis.text.y = element_text(face = "bold", size = 10)
  )

ggsave(file.path(out_dir, "multivariate_forest_plot.pdf"),
       p_forest_mv,
       width = 11, height = max(4, 0.45 * nrow(fp_mv) + 2), dpi = 300)
ggsave(file.path(out_dir, "multivariate_forest_plot.png"),
       p_forest_mv,
       width = 11, height = max(4, 0.45 * nrow(fp_mv) + 2), dpi = 300)
cat("  Multivariate forest plot saved → multivariate_forest_plot.pdf & .png\n\n")


################################################################################
## 5. DETAILED MULTIVARIATE MODEL — FULL FOREST PLOT FOR BEST HUB GENE
################################################################################

cat("════════════════════════════════════════════════════════════════\n")
cat("  5. Full Multivariate Model Forest Plot (best hub gene)\n")
cat("════════════════════════════════════════════════════════════════\n\n")

# Build a full model with the most significant hub gene + all covariates
# and produce a forestmodel-style plot showing ALL covariates
best_gene <- multi_cox_df$gene[2]
cat("  Best hub gene (lowest multivariate p):", best_gene, "\n")

# Prepare a clean data frame for forestmodel
fm_df <- surv_df_mv %>%
  select(os_time_years, os_event, all_of(best_gene),
         age, stage_group, grade_group, gender) %>%
  rename(Expression = all_of(best_gene))
fm_df <- as.data.frame(fm_df)

fm_formula <- Surv(os_time_years, os_event) ~ Expression + age + stage_group + grade_group + gender
fm_cox     <- coxph(fm_formula, data = fm_df)

cat("  Full model summary:\n")
print(summary(fm_cox))

# forestmodel plot
p_fm <- forest_model(fm_cox,
                     format_options = forest_model_format_options(
                       colour   = "black",
                       shape    = 18,
                       text_size = 3.5,
                       banded   = TRUE
                     )) +
  labs(title = paste0("Multivariate Cox Model — ", best_gene,
                      " + Clinical Covariates (Age, Stage, Grade, Gender)")) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))

ggsave(file.path(out_dir, "multivariate_full_forest_best_gene.pdf"),
       p_fm, width = 10, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "multivariate_full_forest_best_gene.png"),
       p_fm, width = 10, height = 5.5, dpi = 300)
cat("\n  Full forest plot saved → multivariate_full_forest_best_gene.pdf & .png\n\n")


################################################################################
## 6. COX PH ASSUMPTION DIAGNOSTICS (Schoenfeld residuals)
################################################################################

cat("════════════════════════════════════════════════════════════════\n")
cat("  6. Proportional Hazards Diagnostics\n")
cat("════════════════════════════════════════════════════════════════\n\n")

# Test PH assumption for the best gene multivariate model
ph_test <- cox.zph(fm_cox)
cat("  Schoenfeld test for PH assumption:\n")
print(ph_test)

pdf(file.path(out_dir, "cox_ph_diagnostics.pdf"), width = 11, height = 8)
par(mfrow = c(2, 3))
for (i in seq_len(ncol(ph_test$y))) {
  plot(ph_test, var = i,
       main = colnames(ph_test$y)[i],
       xlab = "Time (years)", ylab = "Beta(t)",
       lwd = 2, col = "#2B83BA")
  abline(h = 0, lty = 2, col = "grey50")
}
dev.off()

png(file.path(out_dir, "cox_ph_diagnostics.png"), width = 11, height = 8, units = "in", res = 300)
par(mfrow = c(2, 3))
for (i in seq_len(ncol(ph_test$y))) {
  plot(ph_test, var = i,
       main = colnames(ph_test$y)[i],
       xlab = "Time (years)", ylab = "Beta(t)",
       lwd = 2, col = "#2B83BA")
  abline(h = 0, lty = 2, col = "grey50")
}
dev.off()
cat("\n  PH diagnostics saved → cox_ph_diagnostics.pdf & .png\n\n")


################################################################################
## 7. CONCORDANCE (C-INDEX) SUMMARY
################################################################################

cat("════════════════════════════════════════════════════════════════\n")
cat("  7. Concordance Index Summary\n")
cat("════════════════════════════════════════════════════════════════\n\n")

conc_df <- uni_cox_df %>%
  select(gene, concordance, concordance_se) %>%
  rename(C_index_univariate = concordance,
         C_index_uni_se     = concordance_se) %>%
  left_join(
    multi_cox_df %>%
      select(gene, concordance, concordance_se) %>%
      rename(C_index_multivariate = concordance,
             C_index_mv_se        = concordance_se),
    by = "gene"
  )

write_tsv(conc_df, file.path(out_dir, "concordance_summary.tsv"))
cat("  Concordance summary saved → concordance_summary.tsv\n\n")

# ── C-index comparison plot ──────────────────────────────────────────────────
conc_long <- conc_df %>%
  pivot_longer(
    cols      = c(C_index_univariate, C_index_multivariate),
    names_to  = "model",
    values_to = "C_index"
  ) %>%
  mutate(
    se = ifelse(model == "C_index_univariate", C_index_uni_se, C_index_mv_se),
    model = recode(model,
                   "C_index_univariate"   = "Univariate",
                   "C_index_multivariate" = "Multivariate")
  )

p_conc <- ggplot(conc_long, aes(x = reorder(gene, C_index), y = C_index,
                                 fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = C_index - 1.96 * se, ymax = C_index + 1.96 * se),
                position = position_dodge(width = 0.7), width = 0.2,
                linewidth = 0.4) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = c("Univariate" = "#4575B4",
                                "Multivariate" = "#D7191C"),
                    name = "Model") +
  coord_flip() +
  labs(
    title    = "Concordance Index (C-statistic) — Hub Gene Survival Models",
    subtitle = "Dashed line = random prediction (C = 0.5)",
    x = NULL,
    y = "C-index ± 95% CI"
  ) +
  theme_pub +
  theme(axis.text.y = element_text(face = "bold"))

ggsave(file.path(out_dir, "concordance_comparison.pdf"),
       p_conc, width = 10, height = max(5, 0.4 * nrow(conc_df) + 2), dpi = 300)
ggsave(file.path(out_dir, "concordance_comparison.png"),
       p_conc, width = 10, height = max(5, 0.4 * nrow(conc_df) + 2), dpi = 300)
cat("  C-index comparison plot saved → concordance_comparison.pdf & .png\n\n")


################################################################################
## 8. COMBINED SURVIVAL SUMMARY TABLE
################################################################################

cat("════════════════════════════════════════════════════════════════\n")
cat("  8. Combined Survival Summary\n")
cat("════════════════════════════════════════════════════════════════\n\n")

surv_summary <- hub_genes_df %>%
  select(gene_symbol, PIKK_Group, GeneRole, regulation,
         log2FoldChange, pathway_role, composite_score) %>%
  left_join(
    km_summary %>%
      select(gene, cutpoint, cut_method, n_high, n_low,
             logrank_p, logrank_p_adj),
    by = c("gene_symbol" = "gene")
  ) %>%
  left_join(
    uni_cox_df %>%
      select(gene, HR, HR_lower_95, HR_upper_95, p_value, p_adj, concordance) %>%
      rename(uni_HR        = HR,
             uni_HR_lower  = HR_lower_95,
             uni_HR_upper  = HR_upper_95,
             uni_p         = p_value,
             uni_p_adj     = p_adj,
             uni_Cindex    = concordance),
    by = c("gene_symbol" = "gene")
  ) %>%
  left_join(
    multi_cox_df %>%
      select(gene, HR, HR_lower_95, HR_upper_95, p_value, p_adj, concordance) %>%
      rename(mv_HR        = HR,
             mv_HR_lower  = HR_lower_95,
             mv_HR_upper  = HR_upper_95,
             mv_p         = p_value,
             mv_p_adj     = p_adj,
             mv_Cindex    = concordance),
    by = c("gene_symbol" = "gene")
  ) %>%
  mutate(
    survival_significant = case_when(
      !is.na(uni_p_adj) & uni_p_adj < 0.05 ~ "Yes (univariate)",
      !is.na(mv_p_adj)  & mv_p_adj  < 0.05 ~ "Yes (multivariate only)",
      !is.na(logrank_p_adj) & logrank_p_adj < 0.05 ~ "Yes (KM only)",
      TRUE ~ "No"
    )
  ) %>%
  arrange(uni_p)

write_tsv(surv_summary, file.path(out_dir, "survival_summary_table.tsv"))
cat("  Combined survival summary saved → survival_summary_table.tsv\n\n")


################################################################################
## 9. HEATMAP-STYLE SUMMARY PLOT
################################################################################

cat("════════════════════════════════════════════════════════════════\n")
cat("  9. Summary Heatmap Plot\n")
cat("════════════════════════════════════════════════════════════════\n\n")

# Tile plot showing -log10(p) for KM, univariate Cox, and multivariate Cox
heat_df <- surv_summary %>%
  filter(!is.na(logrank_p)) %>%
  select(gene_symbol, logrank_p, uni_p, mv_p) %>%
  pivot_longer(-gene_symbol, names_to = "test", values_to = "p_value") %>%
  mutate(
    neglog10p = -log10(p_value),
    test = recode(test,
                  "logrank_p" = "Log-Rank (KM)",
                  "uni_p"     = "Univariate Cox",
                  "mv_p"      = "Multivariate Cox"),
    test = factor(test, levels = c("Log-Rank (KM)",
                                   "Univariate Cox",
                                   "Multivariate Cox")),
    sig_star = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    gene_symbol = factor(gene_symbol,
                         levels = rev(surv_summary$gene_symbol[
                           order(surv_summary$uni_p)]))
  )

p_heat <- ggplot(heat_df, aes(x = test, y = gene_symbol, fill = neglog10p)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sig_star), size = 4.5, fontface = "bold",
            colour = "white") +
  scale_fill_gradient2(
    low      = "#F7F7F7",
    mid      = "#FDB863",
    high     = "#B2182B",
    midpoint = 1.3,  # ≈ -log10(0.05)
    name     = expression(-log[10](p)),
    limits   = c(0, NA)
  ) +
  labs(
    title    = "Survival Association Significance — OSCC PIKK Hub Genes",
    subtitle = "Stars: * p < 0.05, ** p < 0.01, *** p < 0.001",
    x = NULL, y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x   = element_text(angle = 30, hjust = 1, face = "bold"),
    axis.text.y   = element_text(face = "bold"),
    panel.grid    = element_blank(),
    panel.border  = element_blank()
  )

ggsave(file.path(out_dir, "survival_significance_heatmap.pdf"),
       p_heat, width = 7,
       height = max(5, 0.35 * length(genes_found) + 2), dpi = 300)
ggsave(file.path(out_dir, "survival_significance_heatmap.png"),
       p_heat, width = 7,
       height = max(5, 0.35 * length(genes_found) + 2), dpi = 300)
cat("  Significance heatmap saved → survival_significance_heatmap.pdf & .png\n\n")


################################################################################
## 10. SESSION INFO & WRAP-UP
################################################################################

cat("════════════════════════════════════════════════════════════════\n")
cat("  ANALYSIS COMPLETE\n")
cat("════════════════════════════════════════════════════════════════\n\n")

cat("Output directory:", out_dir, "\n\n")
cat("Files generated:\n")
list.files(out_dir, recursive = TRUE) %>%
  paste0("  • ", .) %>%
  cat(sep = "\n")

cat("\n\n── Session Info ────────────────────────────────────────────────\n")
sessionInfo()
