# =========================================================
# STEP 6: Functional / pathway enrichment analysis (edgeR)
# Primary foreground: extended PIKK-response intersection set
# Background / GSEA ranking: edgeR DGEA all-gene results
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(ggplot2)
  library(stringr)
})

project_root <- "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project"

# -----------------------------
# Input files (edgeR-based)
# -----------------------------
ext_sig_file  <- file.path(project_root, "05_results/Deepaprabha_Intermediary_Results/OSCC_PIKK_extended_edgeR_DEG_intersection_significant.tsv")
edgeR_sig_file <- file.path(project_root, "05_results/Deepaprabha_Intermediary_Results/edgeR_withBatch/results/edgeR_significant_genes.tsv")
edgeR_all_file <- file.path(project_root, "05_results/Deepaprabha_Intermediary_Results/edgeR_withBatch/results/edgeR_all_results_tumor_vs_normal.tsv")

# -----------------------------
# Output folders (edgeR sub-directory)
# -----------------------------
out_dir_tables <- file.path(project_root, "05_results/tables/enrichment/edgeR")
out_dir_plots  <- file.path(project_root, "05_results/plots/enrichment/edgeR")
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helper: save enrichment results and plots
# -----------------------------
save_enrichment <- function(enrich_obj, prefix, out_dir_tables, out_dir_plots) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    message("No enrichment results for: ", prefix)
    return(NULL)
  }
  
  res_df <- as.data.frame(enrich_obj)
  write.table(
    res_df,
    file = file.path(out_dir_tables, paste0(prefix, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  p_dot <- dotplot(enrich_obj, showCategory = min(15, nrow(res_df))) +
    ggtitle(prefix) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold")
    )
  
  ggsave(
    file.path(out_dir_plots, paste0(prefix, "_dotplot.pdf")),
    p_dot, width = 9, height = 6, device = cairo_pdf
  )
  ggsave(
    file.path(out_dir_plots, paste0(prefix, "_dotplot.png")),
    p_dot, width = 9, height = 6, dpi = 300
  )
  
  p_bar <- barplot(enrich_obj, showCategory = min(15, nrow(res_df))) +
    ggtitle(prefix) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold")
    )
  
  ggsave(
    file.path(out_dir_plots, paste0(prefix, "_barplot.pdf")),
    p_bar, width = 9, height = 6, device = cairo_pdf
  )
  ggsave(
    file.path(out_dir_plots, paste0(prefix, "_barplot.png")),
    p_bar, width = 9, height = 6, dpi = 300
  )
  
  return(res_df)
}

# -----------------------------
# Helper: gene symbol -> Entrez
# -----------------------------
symbol_to_entrez <- function(symbols) {
  symbols <- unique(na.omit(toupper(symbols)))
  mapping <- bitr(
    symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  mapping <- mapping[!duplicated(mapping$SYMBOL), ]
  return(mapping)
}

# =========================================================
# 1) Load inputs
# =========================================================
ext_sig   <- read.delim(ext_sig_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
edgeR_sig <- read.delim(edgeR_sig_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
deg_all   <- read.delim(edgeR_all_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# Standardize gene symbol column names
if (!"gene_name" %in% names(ext_sig)) stop("Extended intersection file must contain gene_name")
if (!"gene_name" %in% names(deg_all))  stop("edgeR all-results file must contain gene_name")

ext_sig$gene_symbol <- toupper(ext_sig$gene_name)
deg_all$gene_symbol <- toupper(deg_all$gene_name)

# =========================================================
# 2) Prepare foreground gene sets
# =========================================================
# Primary set: extended PIKK-response intersection (edgeR)
ext_symbols <- unique(ext_sig$gene_symbol)

# Universe/background: all tested edgeR genes
# Keep only genes that can be mapped to Entrez IDs
universe_map <- symbol_to_entrez(deg_all$gene_symbol)

# Foreground mapping
ext_map <- symbol_to_entrez(ext_symbols)

# Save mapping summary
mapping_summary <- data.frame(
  set = c("extended_edgeR_intersection", "universe_edgeR"),
  input_genes = c(length(ext_symbols), length(unique(deg_all$gene_symbol))),
  mapped_entrez = c(nrow(ext_map), nrow(universe_map))
)

write.table(
  mapping_summary,
  file = file.path(out_dir_tables, "gene_mapping_summary_edgeR.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\nGene mapping summary (edgeR):\n")
print(mapping_summary)

# =========================================================
# 3) ORA on the extended edgeR intersection set
# =========================================================
ext_entrez <- unique(ext_map$ENTREZID)
universe_entrez <- unique(universe_map$ENTREZID)

ora_go_bp_ext <- enrichGO(
  gene          = ext_entrez,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE
)

ora_go_mf_ext <- enrichGO(
  gene          = ext_entrez,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE
)

ora_go_cc_ext <- enrichGO(
  gene          = ext_entrez,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE
)

ora_kegg_ext <- enrichKEGG(
  gene          = ext_entrez,
  universe      = universe_entrez,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

ora_reactome_ext <- enrichPathway(
  gene          = ext_entrez,
  universe      = universe_entrez,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

# Save ORA results and plots
save_enrichment(ora_go_bp_ext,    "ORA_GO_BP_edgeR_extended",    out_dir_tables, out_dir_plots)
save_enrichment(ora_go_mf_ext,    "ORA_GO_MF_edgeR_extended",    out_dir_tables, out_dir_plots)
save_enrichment(ora_go_cc_ext,    "ORA_GO_CC_edgeR_extended",    out_dir_tables, out_dir_plots)
save_enrichment(ora_kegg_ext,     "ORA_KEGG_edgeR_extended",     out_dir_tables, out_dir_plots)
save_enrichment(ora_reactome_ext, "ORA_Reactome_edgeR_extended", out_dir_tables, out_dir_plots)

# =========================================================
# ORA visualization helpers
# =========================================================

plot_ora_bubble <- function(df, prefix, out_dir_plots, top_n = 20) {
  cat("Entered function:", prefix, "\n")
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df2 <- df %>%
    mutate(
      GeneCount = as.numeric(stringr::str_extract(GeneRatio, "^[0-9]+")),
      BgCount = as.numeric(stringr::str_extract(BgRatio, "^[0-9]+")),
      negLog10Padj = -log10(p.adjust),
      Description = stringr::str_wrap(Description, width = 45)
    ) %>%
    arrange(p.adjust) %>%
    slice_head(n = min(top_n, nrow(.))) %>%
    mutate(Description = factor(Description, levels = rev(unique(Description))))
  
  cat("\n", prefix, "\n")
  cat("Rows available:", nrow(df), "\n")
  cat("Rows plotted:", nrow(df2), "\n")
  
  p <- ggplot(df2, aes(x = RichFactor, y = Description)) +
    geom_point(aes(size = GeneCount, color = negLog10Padj), alpha = 0.9) +
    scale_color_viridis_c(option = "magma") +
    labs(
      title = prefix,
      x = "Rich factor",
      y = NULL,
      size = "Gene count",
      color = expression(-log[10](p.adjust))
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold")
    )
  
  ggsave(file.path(out_dir_plots, paste0(prefix, "_bubble.pdf")), p, width = 10, height = 9)
  ggsave(file.path(out_dir_plots, paste0(prefix, "_bubble.png")), p, width = 10, height = 9, dpi = 300)
  
  return(df2)
}

plot_ora_bar <- function(df, prefix, out_dir_plots, top_n = 20) {
  cat("Entered function:", prefix, "\n")
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df2 <- df %>%
    mutate(
      negLog10Padj = -log10(p.adjust),
      Description = stringr::str_wrap(Description, width = 45)
    ) %>%
    arrange(p.adjust) %>%
    slice_head(n = min(top_n, nrow(.))) %>%
    mutate(Description = factor(Description, levels = rev(unique(Description))))
  
  cat("\n", prefix, "\n")
  cat("Rows available:", nrow(df), "\n")
  cat("Rows plotted:", nrow(df2), "\n")
  
  p <- ggplot(df2, aes(x = FoldEnrichment, y = Description)) +
    geom_col(aes(fill = negLog10Padj), width = 0.75) +
    scale_fill_viridis_c(option = "plasma") +
    labs(
      title = prefix,
      x = "Fold enrichment",
      y = NULL,
      fill = expression(-log[10](p.adjust))
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold")
    )
  
  ggsave(file.path(out_dir_plots, paste0(prefix, "_bar.pdf")), p, width = 10, height = 9)
  ggsave(file.path(out_dir_plots, paste0(prefix, "_bar.png")), p, width = 10, height = 9, dpi = 300)
  
  return(df2)
}

ora_go_bp_ext_df <- as.data.frame(ora_go_bp_ext)
ora_go_mf_ext_df <- as.data.frame(ora_go_mf_ext)
ora_go_cc_ext_df <- as.data.frame(ora_go_cc_ext)
ora_kegg_ext_df  <- as.data.frame(ora_kegg_ext)
ora_reactome_ext_df <- as.data.frame(ora_reactome_ext)

write.table(ora_go_bp_ext_df, file = file.path(out_dir_tables, "ORA_GO_BP_edgeR_extended.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ora_go_mf_ext_df, file = file.path(out_dir_tables, "ORA_GO_MF_edgeR_extended.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ora_go_cc_ext_df, file = file.path(out_dir_tables, "ORA_GO_CC_edgeR_extended.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ora_kegg_ext_df, file = file.path(out_dir_tables, "ORA_KEGG_edgeR_extended.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ora_reactome_ext_df, file = file.path(out_dir_tables, "ORA_Reactome_edgeR_extended.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_ora_bubble(ora_go_bp_ext_df, "Top20_ORA_GO_BP_edgeR_extended", out_dir_plots)
plot_ora_bar(ora_go_bp_ext_df, "Top20_ORA_GO_BP_edgeR_extended", out_dir_plots)

plot_ora_bubble(ora_go_mf_ext_df, "Top20_ORA_GO_MF_edgeR_extended", out_dir_plots)
plot_ora_bar(ora_go_mf_ext_df, "Top20_ORA_GO_MF_edgeR_extended", out_dir_plots)

plot_ora_bubble(ora_go_cc_ext_df, "Top20_ORA_GO_CC_edgeR_extended", out_dir_plots)
plot_ora_bar(ora_go_cc_ext_df, "Top20_ORA_GO_CC_edgeR_extended", out_dir_plots)

plot_ora_bubble(ora_kegg_ext_df, "Top20_ORA_KEGG_edgeR_extended", out_dir_plots)
plot_ora_bar(ora_kegg_ext_df, "Top20_ORA_KEGG_edgeR_extended", out_dir_plots)

plot_ora_bubble(ora_reactome_ext_df, "Top20_ORA_Reactome_edgeR_extended", out_dir_plots)
plot_ora_bar(ora_reactome_ext_df, "Top20_ORA_Reactome_edgeR_extended", out_dir_plots)


# =========================================================
# 5) Prepare ranked gene list for GSEA (from edgeR results)
# =========================================================
# edgeR QLF results have logFC and F columns. We construct a
# signed ranking metric: sign(logFC) * sqrt(F).
# This is analogous to DESeq2's Wald stat.

deg_rank <- deg_all %>%
  filter(!is.na(F), !is.na(logFC)) %>%
  mutate(
    gene_symbol = toupper(gene_name),
    signed_stat = sign(logFC) * sqrt(F),
    abs_stat    = abs(signed_stat)
  ) %>%
  arrange(desc(abs_stat)) %>%
  distinct(gene_symbol, .keep_all = TRUE)

rank_map <- symbol_to_entrez(deg_rank$gene_symbol)

deg_rank <- deg_rank %>%
  inner_join(rank_map, by = c("gene_symbol" = "SYMBOL")) %>%
  distinct(ENTREZID, .keep_all = TRUE)

geneList <- deg_rank$signed_stat
names(geneList) <- deg_rank$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

gsea_rank_df <- data.frame(
  ENTREZID = names(geneList),
  signed_stat = as.numeric(geneList),
  stringsAsFactors = FALSE
)

write.table(
  gsea_rank_df,
  file = file.path(out_dir_tables, "GSEA_ranked_gene_list_edgeR.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# =========================================================
# 6) GSEA
# =========================================================
gsea_go_bp <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

gsea_kegg <- gseKEGG(
  geneList      = geneList,
  organism      = "hsa",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

gsea_reactome <- gsePathway(
  geneList      = geneList,
  organism      = "human",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

# =========================================================
# GSEA visualization helper: combined + up/down split
# =========================================================

plot_gsea_outputs <- function(gsea_obj, prefix, geneList, out_dir_plots, top_n = 15) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) return(NULL)
  
  gsea_df <- as.data.frame(gsea_obj) %>%
    arrange(p.adjust)
  
  # Find signed enrichment column
  score_col <- intersect(c("NES", "enrichmentScore", "ES"), colnames(gsea_df))[1]
  if (is.na(score_col)) {
    stop("No signed enrichment score column found (NES, enrichmentScore, or ES).")
  }
  
  gsea_df <- gsea_df %>%
    mutate(score = .data[[score_col]])
  
  # top_n is total terms; split equally by direction
  top_each <- max(1, floor(top_n / 2))
  
  make_subset_obj <- function(df_sub) {
    if (is.null(df_sub) || nrow(df_sub) == 0) return(NULL)
    obj <- gsea_obj
    obj@result <- df_sub
    obj
  }
  
  save_basic_plots <- function(obj, df_use, plot_prefix, fill_title) {
    if (is.null(obj) || nrow(df_use) == 0) return(NULL)
    
    k <- min(nrow(df_use), top_each)
    
    p_dot <- dotplot(obj, showCategory = k) +
      ggtitle(fill_title) +
      theme_classic(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold")
      )
    
    ggsave(
      file.path(out_dir_plots, paste0(plot_prefix, "_dotplot.pdf")),
      p_dot, width = 9, height = 6
    )
    ggsave(
      file.path(out_dir_plots, paste0(plot_prefix, "_dotplot.png")),
      p_dot, width = 9, height = 6, dpi = 300
    )
    
    p_ridge <- ridgeplot(obj, showCategory = k) +
      ggtitle(fill_title) +
      theme_classic(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold")
      )
    
    ggsave(
      file.path(out_dir_plots, paste0(plot_prefix, "_ridgeplot.pdf")),
      p_ridge, width = 9, height = 6
    )
    ggsave(
      file.path(out_dir_plots, paste0(plot_prefix, "_ridgeplot.png")),
      p_ridge, width = 9, height = 6, dpi = 300
    )
    
    top_ids <- df_use$ID[seq_len(min(3, nrow(df_use)))]
    
    for (term_id in top_ids) {
      term_label <- df_use$Description[df_use$ID == term_id][1]
      safe_name <- stringr::str_replace_all(term_label, "[^A-Za-z0-9_]+", "_")
      
      p_enrich <- gseaplot2(
        gsea_obj,
        geneSetID = term_id,
        title = term_label,
        color = ifelse(grepl("down", plot_prefix, ignore.case = TRUE), "blue", "red")
      )
      
      ggsave(
        file.path(out_dir_plots, paste0(plot_prefix, "_", safe_name, "_gseaplot2.pdf")),
        p_enrich, width = 9, height = 6
      )
      ggsave(
        file.path(out_dir_plots, paste0(plot_prefix, "_", safe_name, "_gseaplot2.png")),
        p_enrich, width = 9, height = 6, dpi = 300
      )
    }
  }
  
  # -------------------------------------------------
  # Combined plots: top_n total
  # -------------------------------------------------
  p_dot_all <- dotplot(gsea_obj, showCategory = min(top_n, nrow(gsea_df))) +
    ggtitle(prefix) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold")
    )
  
  ggsave(file.path(out_dir_plots, paste0(prefix, "_dotplot.pdf")), p_dot_all, width = 9, height = 6)
  ggsave(file.path(out_dir_plots, paste0(prefix, "_dotplot.png")), p_dot_all, width = 9, height = 6, dpi = 300)
  
  p_ridge_all <- ridgeplot(gsea_obj, showCategory = min(top_n, nrow(gsea_df))) +
    ggtitle(prefix) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold")
    )
  
  ggsave(file.path(out_dir_plots, paste0(prefix, "_ridgeplot.pdf")), p_ridge_all, width = 9, height = 6)
  ggsave(file.path(out_dir_plots, paste0(prefix, "_ridgeplot.png")), p_ridge_all, width = 9, height = 6, dpi = 300)
  
  top_ids_all <- gsea_df$ID[seq_len(min(3, nrow(gsea_df)))]
  for (term_id in top_ids_all) {
    term_label <- gsea_df$Description[gsea_df$ID == term_id][1]
    safe_name <- stringr::str_replace_all(term_label, "[^A-Za-z0-9_]+", "_")
    
    p_enrich <- gseaplot2(gsea_obj, geneSetID = term_id, title = term_label, color = "red")
    
    ggsave(
      file.path(out_dir_plots, paste0(prefix, "_", safe_name, "_gseaplot2.pdf")),
      p_enrich, width = 9, height = 6
    )
    ggsave(
      file.path(out_dir_plots, paste0(prefix, "_", safe_name, "_gseaplot2.png")),
      p_enrich, width = 9, height = 6, dpi = 300
    )
  }
  
  # -------------------------------------------------
  # Direction-specific subsets
  # -------------------------------------------------
  gsea_up <- gsea_df %>%
    filter(!is.na(p.adjust), p.adjust < 0.05, score > 0) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_each)
  
  gsea_down <- gsea_df %>%
    filter(!is.na(p.adjust), p.adjust < 0.05, score < 0) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_each)
  
  write.table(
    gsea_up,
    file = file.path(out_dir_plots, paste0(prefix, "_upregulated_terms.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    gsea_down,
    file = file.path(out_dir_plots, paste0(prefix, "_downregulated_terms.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  # Upregulated plots
  if (nrow(gsea_up) > 0) {
    obj_up <- make_subset_obj(gsea_up)
    save_basic_plots(
      obj = obj_up,
      df_use = gsea_up,
      plot_prefix = paste0(prefix, "_upregulated"),
      fill_title = paste0(prefix, " - Upregulated pathways")
    )
  }
  
  # Downregulated plots
  if (nrow(gsea_down) > 0) {
    obj_down <- make_subset_obj(gsea_down)
    save_basic_plots(
      obj = obj_down,
      df_use = gsea_down,
      plot_prefix = paste0(prefix, "_downregulated"),
      fill_title = paste0(prefix, " - Downregulated pathways")
    )
  }
  
  return(gsea_df)
}

gsea_go_bp_df <- as.data.frame(gsea_go_bp)
gsea_kegg_df <- as.data.frame(gsea_kegg)
gsea_reactome_df <- as.data.frame(gsea_reactome)

write.table(gsea_go_bp_df, file = file.path(out_dir_tables, "GSEA_GO_BP_edgeR_allgenes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gsea_kegg_df, file = file.path(out_dir_tables, "GSEA_KEGG_edgeR_allgenes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gsea_reactome_df, file = file.path(out_dir_tables, "GSEA_Reactome_edgeR_allgenes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_gsea_outputs(gsea_go_bp, "GSEA_GO_BP_edgeR_allgenes", geneList, out_dir_plots, top_n = 20)
plot_gsea_outputs(gsea_kegg, "GSEA_KEGG_edgeR_allgenes", geneList, out_dir_plots, top_n = 20)
plot_gsea_outputs(gsea_reactome, "GSEA_Reactome_edgeR_allgenes", geneList, out_dir_plots, top_n = 20)


# =========================================================
# 8) Session summary
# =========================================================
cat("\n==============================\n")
cat("edgeR ENRICHMENT ANALYSIS COMPLETE\n")
cat("==============================\n")
cat("Extended foreground genes (input): ", length(ext_symbols), "\n")
cat("Universe mapped genes (edgeR):     ", nrow(universe_map), "\n")
cat("GSEA ranked genes (edgeR):         ", length(geneList), "\n")
cat("Ranking metric: sign(logFC) * sqrt(F)\n")
cat("Outputs written to:\n")
cat(out_dir_tables, "\n")
cat(out_dir_plots, "\n")
