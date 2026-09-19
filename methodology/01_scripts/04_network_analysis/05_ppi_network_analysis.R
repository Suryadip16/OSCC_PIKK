#!/usr/bin/env Rscript

# ============================================================================
# STEP 2.2: Gene Interaction Network Analysis
# OSCC-PIKK DDR & Checkpoint Signalling
# ============================================================================
#
# This script constructs and analyses the protein-protein interaction (PPI)
# network for the significant PIKK-associated DEGs identified from edgeR DGEA.
# Interactions are retrieved from the STRING database (v12.0, combined
# score >= 0.7, high confidence).
#
# Workflow:
#   1. Load PIKK DEG intersection results & full gene universe
#   2. Retrieve PPI data from STRING (direct API or manual TSV fallback)
#   3. Construct igraph network & annotate nodes
#   4. Compute centrality metrics (degree, betweenness, closeness, eigenvector)
#   5. Extract sub-network (DEGs + PIKK-universe first-order interactors)
#   6. Detect communities (Louvain algorithm)
#   7. Export tables, Cytoscape-compatible files, and GraphML
#   8. Generate publication-quality network & centrality visualisations
#
# Required R packages:
#   CRAN : igraph, ggraph, tidygraph, ggplot2, ggrepel, dplyr, tidyr,
#          readr, tibble, pheatmap, RColorBrewer, scales, cowplot, reshape2
#   (Optional) Bioconductor: STRINGdb
#
# Inputs:
#   - OSCC_PIKK_extended_edgeR_DEG_intersection_significant.tsv (65 DEGs)
#   - PIKK_related_gene_universe.csv (204 genes)
#   - STRING PPI data (from API or manual download)
#
# Outputs:
#   - Tables: network_topology_metrics.tsv, STRING_interactions_*.tsv,
#             community_assignments.tsv, subnetwork edges
#   - Plots : 11 publication-quality figures (PDF + PNG)
#   - Files : GraphML, Cytoscape-compatible edge/node attribute TSVs
#
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

# ggraph + tidygraph for publication-quality network layouts
has_ggraph <- requireNamespace("ggraph", quietly = TRUE) &&
  requireNamespace("tidygraph", quietly = TRUE)
if (has_ggraph) {
  library(ggraph)
  library(tidygraph)
  cat("ggraph + tidygraph loaded for network visualisation.\n")
} else {
  message("NOTE: ggraph / tidygraph not installed. Using base igraph plotting.")
  message("For publication-quality network figures, install with:")
  message('  install.packages(c("ggraph", "tidygraph"))')
}

# --- Project root & paths ---------------------------------------------------
project_root <- "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project"

# Input files
pikk_deg_file <- file.path(
  project_root,
  "05_results/Deepaprabha_Intermediary_Results/OSCC_PIKK_extended_edgeR_DEG_intersection_significant.tsv"
)
pikk_univ_file <- file.path(
  project_root,
  "10_references/PIKK_related_gene_universe.csv"
)

# Output directories
out_tables <- file.path(project_root, "05_results/tables/network")
out_plots <- file.path(project_root, "05_results/plots/network")
out_network <- file.path(project_root, "05_results/network_files")
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(out_network, recursive = TRUE, showWarnings = FALSE)

# --- STRING configuration ---------------------------------------------------
# Primary: use STRING REST API (no extra R packages needed)
# Fallback: read a manually-downloaded STRING TSV file
use_string_api <- TRUE
string_confidence <- 700 # 0–1000 scale; 700 = high confidence
string_species <- 9606 # Homo sapiens

# Path for manually-downloaded STRING export (used if API fails or is disabled)
manual_string_file <- file.path(
  project_root,
  "10_references/STRING_interactions_manual.tsv"
)

# --- Aesthetics --------------------------------------------------------------
theme_pub <- theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "black"),
    legend.title  = element_text(face = "bold"),
    panel.grid    = element_blank()
  )

# PIKK family colour palette (6 groups)
pikk_colors <- c(
  "ATR"   = "#E41A1C",
  "ATM"   = "#377EB8",
  "PRKDC" = "#4DAF4A",
  "MTOR"  = "#984EA3",
  "SMG1"  = "#FF7F00",
  "TRRAP" = "#A65628"
)

# Regulation direction colours
reg_colors <- c(
  "Up in Tumor"   = "#D73027",
  "Down in Tumor" = "#4575B4"
)

# Dual-format saver (PDF + PNG)
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
cat("STEP 2.2: PPI Network Analysis — OSCC-PIKK\n")
cat("================================================================\n\n")

# --- 1a. Significant PIKK DEGs (65 genes) -----------------------------------
pikk_deg <- read.delim(pikk_deg_file,
  header = TRUE, sep = "\t",
  check.names = FALSE, stringsAsFactors = FALSE
)
pikk_deg$gene_symbol <- toupper(pikk_deg$gene_name)
cat("Loaded", nrow(pikk_deg), "significant PIKK-associated DEGs.\n")

# --- 1b. Full PIKK gene universe (204 genes) --------------------------------
pikk_univ <- read.csv(pikk_univ_file, stringsAsFactors = FALSE)
pikk_univ$GeneSymbol <- toupper(pikk_univ$GeneSymbol)
pikk_univ <- pikk_univ %>%
  distinct(PIKK_Group, GeneSymbol, GeneRole, .keep_all = TRUE)
cat("Loaded", nrow(pikk_univ), "genes in the PIKK gene universe.\n")

# --- 1c. Gene symbol alias mapping for STRING --------------------------------
# Some gene symbols in our universe differ from STRING's preferred names
alias_map <- c(
  "H2AX"    = "H2AFX",
  "MRE11"   = "MRE11A"
)

# Create a lookup from query name → original name
original_names <- unique(pikk_univ$GeneSymbol)
query_names <- sapply(original_names, function(g) {
  if (g %in% names(alias_map)) alias_map[[g]] else g
})
alias_lookup <- data.frame(
  query_name = unname(query_names),
  original_name = original_names,
  stringsAsFactors = FALSE
) %>% distinct(query_name, .keep_all = TRUE)

cat("Gene aliases applied:", sum(original_names != query_names), "genes remapped.\n")

# DEG gene list (with aliases applied)
deg_symbols <- unique(pikk_deg$gene_symbol)
deg_query <- sapply(deg_symbols, function(g) {
  if (g %in% names(alias_map)) alias_map[[g]] else g
})

# Reverse alias lookup function
resolve_alias <- function(names_vec) {
  idx <- match(toupper(names_vec), toupper(alias_lookup$query_name))
  ifelse(!is.na(idx), alias_lookup$original_name[idx], names_vec)
}


# ============================================================================
# 2. STRING PPI DATA RETRIEVAL
# ============================================================================
cat("\n--- STRING PPI Data Retrieval ---\n")

string_interactions <- NULL

# --- 2a. Primary: STRING REST API (no extra packages) -----------------------
if (use_string_api) {
  cat("Querying STRING API for", length(query_names), "PIKK universe genes...\n")
  cat(
    "Species:", string_species, "| Confidence threshold:",
    string_confidence, "\n"
  )

  string_url <- paste0(
    "https://string-db.org/api/tsv/network?",
    "identifiers=", paste(unique(query_names), collapse = "%0d"),
    "&species=", string_species,
    "&required_score=", string_confidence,
    "&caller_identity=OSCC_PIKK_DDR_analysis"
  )

  string_interactions <- tryCatch(
    {
      raw <- read.delim(string_url, stringsAsFactors = FALSE, header = TRUE)
      cat("STRING API returned", nrow(raw), "interactions.\n")
      raw
    },
    error = function(e) {
      message("STRING API query failed: ", e$message)
      message("Falling back to manual STRING file...\n")
      NULL
    }
  )
}

# --- 2b. Fallback: manually-downloaded STRING TSV ----------------------------
if (is.null(string_interactions)) {
  if (file.exists(manual_string_file)) {
    cat("Reading manual STRING file:", manual_string_file, "\n")
    string_interactions <- read.delim(manual_string_file,
      stringsAsFactors = FALSE,
      header = TRUE,
      comment.char = "#"
    )
    cat("Loaded", nrow(string_interactions), "interactions from manual file.\n")
  } else {
    # =========================================================================
    # MANUAL STRING DOWNLOAD INSTRUCTIONS
    # =========================================================================
    # If the STRING API is unavailable, download the data manually:
    #   1. Go to https://string-db.org/
    #   2. Click "Multiple proteins" in the left menu
    #   3. Paste the gene list printed below
    #   4. Select "Homo sapiens" as the organism
    #   5. Click "Search" → then "Continue" on the next page
    #   6. Set minimum required interaction score to 0.700 (high confidence)
    #      in the Settings (bottom of page)
    #   7. Click the "Exports" tab at the top
    #   8. Click "as short tabular text output (TSV)" to download
    #   9. Save the file as:
    #      10_references/STRING_interactions_manual.tsv
    #  10. Re-run this script with  use_string_api <- FALSE
    # =========================================================================
    cat("\n======================================================\n")
    cat("STRING DATA NOT AVAILABLE\n")
    cat("======================================================\n")
    cat("Please download the STRING network manually.\n")
    cat("Paste the following gene list into STRING:\n\n")
    cat(paste(unique(query_names), collapse = "\n"), "\n\n")
    cat("Save the TSV export to:\n  ", manual_string_file, "\n")
    cat("Then re-run this script.\n")
    stop("STRING interaction data required. See instructions above.")
  }
}

# --- 2c. Standardise STRING output columns -----------------------------------
# STRING API returns: stringId_A, stringId_B, preferredName_A, preferredName_B,
#                     ncbiTaxonId, score, nscore, fscore, pscore, ascore,
#                     escore, dscore, tscore
# Manual exports may have: #node1, node2, ... , combined_score

# Detect column names and standardise
col_names <- colnames(string_interactions)

if ("preferredName_A" %in% col_names && "preferredName_B" %in% col_names) {
  # Standard STRING API format
  string_edges <- string_interactions %>%
    dplyr::rename(
      gene_A         = preferredName_A,
      gene_B         = preferredName_B,
      combined_score = score
    ) %>%
    dplyr::select(gene_A, gene_B, combined_score, everything())
} else if ("node1" %in% col_names || "X.node1" %in% col_names) {
  # Manual web export format (may have #node1 parsed as X.node1)
  n1_col <- ifelse("X.node1" %in% col_names, "X.node1", "node1")
  string_edges <- string_interactions %>%
    dplyr::rename(
      gene_A         = !!sym(n1_col),
      gene_B         = node2,
      combined_score = combined_score
    ) %>%
    dplyr::select(gene_A, gene_B, combined_score, everything())
} else {
  # Try to infer — assume first two text columns are gene names
  message("WARNING: Unexpected STRING column format. Attempting to parse...")
  message("Columns found: ", paste(col_names, collapse = ", "))
  string_edges <- string_interactions
  colnames(string_edges)[1:2] <- c("gene_A", "gene_B")
  if (!"combined_score" %in% colnames(string_edges)) {
    # Look for a 'score' column
    score_col <- grep("score", colnames(string_edges), value = TRUE)[1]
    if (!is.na(score_col)) {
      string_edges$combined_score <- string_edges[[score_col]]
    } else {
      string_edges$combined_score <- 0.7
      message("WARNING: No score column found. Setting all scores to 0.7.")
    }
  }
}

# Standardise gene names (resolve aliases back to original names)
string_edges$gene_A <- toupper(string_edges$gene_A)
string_edges$gene_B <- toupper(string_edges$gene_B)
string_edges$gene_A <- resolve_alias(string_edges$gene_A)
string_edges$gene_B <- resolve_alias(string_edges$gene_B)

# Remove self-loops and duplicate edges
string_edges <- string_edges %>%
  dplyr::filter(gene_A != gene_B) %>%
  dplyr::mutate(
    edge_key = ifelse(gene_A < gene_B,
      paste(gene_A, gene_B, sep = "::"),
      paste(gene_B, gene_A, sep = "::")
    )
  ) %>%
  dplyr::distinct(edge_key, .keep_all = TRUE) %>%
  dplyr::select(-edge_key)

cat("After cleaning:", nrow(string_edges), "unique interactions.\n")

# Save all filtered STRING interactions
write.table(string_edges,
  file = file.path(out_tables, "STRING_interactions_all_PIKK_universe.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)


# ============================================================================
# 3. NETWORK CONSTRUCTION — CORE NETWORK (DEGs only)
# ============================================================================
cat("\n--- Network Construction ---\n")

# --- 3a. Core network: edges where BOTH endpoints are significant DEGs -------
core_edges <- string_edges %>%
  dplyr::filter(gene_A %in% deg_symbols & gene_B %in% deg_symbols)

cat("Core network (DEG-DEG interactions):", nrow(core_edges), "edges.\n")

# Build node attribute table for DEGs present in core network
core_genes <- unique(c(core_edges$gene_A, core_edges$gene_B))

# Also include isolated DEGs (no STRING interactions) as disconnected nodes
all_deg_nodes <- unique(pikk_deg$gene_symbol)
isolated_degs <- setdiff(all_deg_nodes, core_genes)
cat(
  "Connected DEGs:", length(core_genes),
  "| Isolated DEGs:", length(isolated_degs), "\n"
)

# Node attribute table
node_attrs <- pikk_deg %>%
  dplyr::select(
    gene_symbol, log2FoldChange, logCPM, F, PValue, padj,
    regulation, PIKK_Group, GeneRole
  ) %>%
  dplyr::distinct(gene_symbol, .keep_all = TRUE)

# Create igraph object
g_core <- graph_from_data_frame(
  d        = core_edges[, c("gene_A", "gene_B")],
  directed = FALSE,
  vertices = node_attrs$gene_symbol
)

# Add edge attribute: combined score
E(g_core)$combined_score <- core_edges$combined_score

# Add node attributes
for (i in seq_len(nrow(node_attrs))) {
  v_name <- node_attrs$gene_symbol[i]
  if (v_name %in% V(g_core)$name) {
    V(g_core)[v_name]$logFC <- node_attrs$log2FoldChange[i]
    V(g_core)[v_name]$padj <- node_attrs$padj[i]
    V(g_core)[v_name]$regulation <- node_attrs$regulation[i]
    V(g_core)[v_name]$pikk_group <- node_attrs$PIKK_Group[i]
    V(g_core)[v_name]$gene_role <- node_attrs$GeneRole[i]
    V(g_core)[v_name]$logCPM <- node_attrs$logCPM[i]
  }
}

cat("Core igraph: ", vcount(g_core), "nodes,", ecount(g_core), "edges.\n")


# ============================================================================
# 4. TOPOLOGY METRICS
# ============================================================================
cat("\n--- Topology Metrics ---\n")

# Compute centrality metrics on the core network
topo_df <- data.frame(
  gene_symbol        = V(g_core)$name,
  degree             = degree(g_core),
  betweenness        = betweenness(g_core, normalized = TRUE),
  closeness          = closeness(g_core, normalized = TRUE),
  eigenvector        = eigen_centrality(g_core)$vector,
  clustering_coeff   = transitivity(g_core, type = "local"),
  stringsAsFactors   = FALSE
)

# Replace NaN clustering coefficients (isolated or degree-1 nodes) with 0
topo_df$clustering_coeff[is.nan(topo_df$clustering_coeff)] <- 0

# Merge with DE + PIKK annotation
topo_df <- topo_df %>%
  dplyr::left_join(node_attrs, by = "gene_symbol") %>%
  dplyr::arrange(desc(degree), desc(betweenness))

cat("Topology metrics computed for", nrow(topo_df), "nodes.\n")
cat("Top-5 by degree:\n")
print(head(topo_df %>% dplyr::select(
  gene_symbol, degree, betweenness,
  PIKK_Group, regulation
), 5))

# Save topology metrics table
write.table(topo_df,
  file = file.path(out_tables, "network_topology_metrics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("Topology metrics saved.\n")


# ============================================================================
# 5. SUB-NETWORK EXTRACTION
# ============================================================================
cat("\n--- Sub-network Extraction ---\n")

# Sub-network: edges where at least ONE endpoint is a significant DEG
# (includes non-DEG PIKK universe genes as first-order interactors)
subnetwork_edges <- string_edges %>%
  dplyr::filter(gene_A %in% deg_symbols | gene_B %in% deg_symbols)

# Identify first-order interactors (PIKK universe genes that are NOT DEGs
# but interact with DEGs)
subnetwork_genes <- unique(c(subnetwork_edges$gene_A, subnetwork_edges$gene_B))
first_order_interactors <- setdiff(subnetwork_genes, deg_symbols)
cat(
  "Sub-network genes:", length(subnetwork_genes),
  "(DEGs:", sum(subnetwork_genes %in% deg_symbols),
  ", first-order interactors:", length(first_order_interactors), ")\n"
)
cat("Sub-network edges:", nrow(subnetwork_edges), "\n")

# Build sub-network node attributes
subnetwork_node_attrs <- pikk_univ %>%
  dplyr::filter(GeneSymbol %in% subnetwork_genes) %>%
  dplyr::distinct(GeneSymbol, .keep_all = TRUE) %>%
  dplyr::left_join(
    pikk_deg %>%
      dplyr::select(gene_symbol, log2FoldChange, padj, regulation) %>%
      dplyr::distinct(gene_symbol, .keep_all = TRUE),
    by = c("GeneSymbol" = "gene_symbol")
  ) %>%
  dplyr::mutate(
    node_type  = ifelse(GeneSymbol %in% deg_symbols, "DEG", "Interactor"),
    regulation = ifelse(is.na(regulation), "Not DE", regulation)
  )

# Build sub-network igraph
g_sub <- graph_from_data_frame(
  d        = subnetwork_edges[, c("gene_A", "gene_B")],
  directed = FALSE,
  vertices = subnetwork_node_attrs$GeneSymbol
)
E(g_sub)$combined_score <- subnetwork_edges$combined_score

# Add node attributes to sub-network
for (i in seq_len(nrow(subnetwork_node_attrs))) {
  v_name <- subnetwork_node_attrs$GeneSymbol[i]
  if (v_name %in% V(g_sub)$name) {
    V(g_sub)[v_name]$pikk_group <- subnetwork_node_attrs$PIKK_Group[i]
    V(g_sub)[v_name]$gene_role <- subnetwork_node_attrs$GeneRole[i]
    V(g_sub)[v_name]$node_type <- subnetwork_node_attrs$node_type[i]
    V(g_sub)[v_name]$regulation <- subnetwork_node_attrs$regulation[i]
    V(g_sub)[v_name]$logFC <- ifelse(
      is.na(subnetwork_node_attrs$log2FoldChange[i]), 0,
      subnetwork_node_attrs$log2FoldChange[i]
    )
  }
}

# Save sub-network edges
write.table(subnetwork_edges,
  file = file.path(out_tables, "first_order_subnetwork_edges.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# --- DDR/Checkpoint-focused sub-network (ATR, ATM, PRKDC groups only) --------
ddr_groups <- c("ATR", "ATM", "PRKDC")
ddr_genes <- pikk_univ %>%
  dplyr::filter(PIKK_Group %in% ddr_groups) %>%
  dplyr::pull(GeneSymbol) %>%
  unique()

ddr_edges <- string_edges %>%
  dplyr::filter(gene_A %in% ddr_genes & gene_B %in% ddr_genes)

if (nrow(ddr_edges) > 0) {
  ddr_node_genes <- unique(c(ddr_edges$gene_A, ddr_edges$gene_B))
  ddr_node_attrs <- pikk_univ %>%
    dplyr::filter(GeneSymbol %in% ddr_node_genes) %>%
    dplyr::distinct(GeneSymbol, .keep_all = TRUE) %>%
    dplyr::left_join(
      pikk_deg %>%
        dplyr::select(gene_symbol, log2FoldChange, padj, regulation) %>%
        dplyr::distinct(gene_symbol, .keep_all = TRUE),
      by = c("GeneSymbol" = "gene_symbol")
    ) %>%
    dplyr::mutate(
      node_type  = ifelse(GeneSymbol %in% deg_symbols, "DEG", "Non-DE"),
      regulation = ifelse(is.na(regulation), "Not DE", regulation)
    )

  g_ddr <- graph_from_data_frame(
    d        = ddr_edges[, c("gene_A", "gene_B")],
    directed = FALSE,
    vertices = ddr_node_attrs$GeneSymbol
  )
  E(g_ddr)$combined_score <- ddr_edges$combined_score

  for (i in seq_len(nrow(ddr_node_attrs))) {
    v_name <- ddr_node_attrs$GeneSymbol[i]
    if (v_name %in% V(g_ddr)$name) {
      V(g_ddr)[v_name]$pikk_group <- ddr_node_attrs$PIKK_Group[i]
      V(g_ddr)[v_name]$node_type <- ddr_node_attrs$node_type[i]
      V(g_ddr)[v_name]$regulation <- ddr_node_attrs$regulation[i]
      V(g_ddr)[v_name]$logFC <- ifelse(
        is.na(ddr_node_attrs$log2FoldChange[i]), 0,
        ddr_node_attrs$log2FoldChange[i]
      )
    }
  }
  cat(
    "DDR sub-network (ATR+ATM+PRKDC):", vcount(g_ddr), "nodes,",
    ecount(g_ddr), "edges.\n"
  )
} else {
  g_ddr <- NULL
  cat("No DDR sub-network edges found at this confidence threshold.\n")
}


# ============================================================================
# 6. COMMUNITY DETECTION
# ============================================================================
cat("\n--- Community Detection (Louvain) ---\n")

set.seed(42)
louvain <- cluster_louvain(g_core)

community_df <- data.frame(
  gene_symbol = V(g_core)$name,
  community = membership(louvain),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    topo_df %>% dplyr::select(gene_symbol, degree, PIKK_Group, regulation),
    by = "gene_symbol"
  )

# Identify the dominant PIKK group per community for labelling
community_labels <- community_df %>%
  dplyr::group_by(community) %>%
  dplyr::summarise(
    n_genes = n(),
    dominant_group = names(sort(table(PIKK_Group), decreasing = TRUE))[1],
    top_gene = gene_symbol[which.max(degree)],
    .groups = "drop"
  )

cat("Detected", max(community_df$community), "communities.\n")
print(community_labels)

# Add community membership to topology table
topo_df <- topo_df %>%
  dplyr::left_join(
    community_df %>% dplyr::select(gene_symbol, community),
    by = "gene_symbol"
  )

# Save community assignments
write.table(community_df,
  file = file.path(out_tables, "community_assignments.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Update saved topology metrics (now includes community column)
write.table(topo_df,
  file = file.path(out_tables, "network_topology_metrics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Add community as vertex attribute
V(g_core)$community <- community_df$community[
  match(V(g_core)$name, community_df$gene_symbol)
]


# ============================================================================
# 7. EXPORT NETWORK FILES (Cytoscape + GraphML)
# ============================================================================
cat("\n--- Exporting Network Files ---\n")

# --- 7a. Cytoscape edge list ------------------------------------------------
cytoscape_edges <- core_edges %>%
  dplyr::transmute(
    source          = gene_A,
    target          = gene_B,
    interaction     = "pp",
    combined_score  = combined_score
  )
write.table(cytoscape_edges,
  file = file.path(out_network, "cytoscape_edge_list.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# --- 7b. Cytoscape node attributes ------------------------------------------
cytoscape_nodes <- topo_df %>%
  dplyr::select(
    gene_symbol, PIKK_Group, GeneRole, regulation,
    log2FoldChange, padj, degree, betweenness, closeness,
    eigenvector, clustering_coeff, community
  )
write.table(cytoscape_nodes,
  file = file.path(out_network, "cytoscape_node_attributes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# --- 7c. GraphML export -----------------------------------------------------
write_graph(g_core, file.path(out_network, "PPI_core_network.graphml"),
  format = "graphml"
)
cat("Cytoscape edge list, node attributes, and GraphML exported.\n")


# ============================================================================
# 8. NETWORK VISUALISATIONS
# ============================================================================
cat("\n--- Network Visualisations ---\n")

# Pre-compute layout (Fruchterman–Reingold) and store coordinates
# so all network views share the same spatial arrangement
set.seed(42)
layout_coords <- layout_with_fr(g_core)
V(g_core)$x <- layout_coords[, 1]
V(g_core)$y <- layout_coords[, 2]

# Decide which nodes get labels (top ~20% by degree to avoid clutter)
deg_threshold <- quantile(degree(g_core), 0.75)
label_genes <- V(g_core)$name[degree(g_core) >= deg_threshold]

if (has_ggraph) {
  tg_core <- as_tbl_graph(g_core)

  # --- 8.1 Full network coloured by PIKK group --------------------------------
  p_pikk <- ggraph(tg_core, layout = "manual", x = x, y = y) +
    geom_edge_link(aes(width = combined_score),
      alpha = 0.2,
      colour = "grey60", show.legend = FALSE
    ) +
    scale_edge_width(range = c(0.3, 1.5)) +
    geom_node_point(aes(size = degree(g_core), fill = pikk_group),
      shape = 21, colour = "grey30", stroke = 0.4
    ) +
    geom_node_text(aes(label = ifelse(name %in% label_genes, name, "")),
      repel = TRUE, size = 3, max.overlaps = 20,
      fontface = "bold"
    ) +
    scale_fill_manual(values = pikk_colors, name = "PIKK Group") +
    scale_size_continuous(range = c(3, 14), name = "Degree") +
    labs(
      title = "PIKK-associated DEG PPI Network (STRING ≥ 0.7)",
      subtitle = paste0(
        vcount(g_core), " nodes, ",
        ecount(g_core), " edges"
      )
    ) +
    theme_graph(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
  save_pdf_png(p_pikk, "PPI_network_by_PIKK_group", width = 13, height = 10)


  # --- 8.2 Network coloured by regulation direction ----------------------------
  p_reg <- ggraph(tg_core, layout = "manual", x = x, y = y) +
    geom_edge_link(aes(width = combined_score),
      alpha = 0.2,
      colour = "grey60", show.legend = FALSE
    ) +
    scale_edge_width(range = c(0.3, 1.5)) +
    geom_node_point(aes(size = degree(g_core), fill = regulation),
      shape = 21, colour = "grey30", stroke = 0.4
    ) +
    geom_node_text(aes(label = ifelse(name %in% label_genes, name, "")),
      repel = TRUE, size = 3, max.overlaps = 20,
      fontface = "bold"
    ) +
    scale_fill_manual(values = reg_colors, name = "Regulation") +
    scale_size_continuous(range = c(3, 14), name = "Degree") +
    labs(
      title = "PIKK DEG PPI Network — Regulation Direction",
      subtitle = paste0(
        vcount(g_core), " nodes, ",
        ecount(g_core), " edges"
      )
    ) +
    theme_graph(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
  save_pdf_png(p_reg, "PPI_network_by_regulation", width = 13, height = 10)


  # --- 8.3 Network with logFC colour gradient ----------------------------------
  logfc_vals <- V(g_core)$logFC
  logfc_lim <- max(abs(logfc_vals), na.rm = TRUE)

  p_logfc <- ggraph(tg_core, layout = "manual", x = x, y = y) +
    geom_edge_link(aes(width = combined_score),
      alpha = 0.2,
      colour = "grey60", show.legend = FALSE
    ) +
    scale_edge_width(range = c(0.3, 1.5)) +
    geom_node_point(aes(size = degree(g_core), fill = logFC),
      shape = 21, colour = "grey30", stroke = 0.4
    ) +
    geom_node_text(aes(label = ifelse(name %in% label_genes, name, "")),
      repel = TRUE, size = 3, max.overlaps = 20,
      fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = "#4575B4", mid = "white", high = "#D73027",
      midpoint = 0,
      limits = c(-logfc_lim, logfc_lim),
      name = expression(log[2] * FC)
    ) +
    scale_size_continuous(range = c(3, 14), name = "Degree") +
    labs(
      title = expression(bold("PIKK DEG PPI Network —" ~ log[2] * "FC Gradient")),
      subtitle = "Blue = down in tumour, Red = up in tumour"
    ) +
    theme_graph(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
  save_pdf_png(p_logfc, "PPI_network_logFC_gradient", width = 13, height = 10)


  # --- 8.4 DDR/Checkpoint sub-network -----------------------------------------
  if (!is.null(g_ddr) && vcount(g_ddr) > 0) {
    set.seed(42)
    layout_ddr <- layout_with_fr(g_ddr)
    V(g_ddr)$x <- layout_ddr[, 1]
    V(g_ddr)$y <- layout_ddr[, 2]

    tg_ddr <- as_tbl_graph(g_ddr)

    # Label all nodes in DDR sub-network (smaller, so labels won't clutter)
    p_ddr <- ggraph(tg_ddr, layout = "manual", x = x, y = y) +
      geom_edge_link(aes(width = combined_score),
        alpha = 0.25,
        colour = "grey50", show.legend = FALSE
      ) +
      scale_edge_width(range = c(0.5, 2)) +
      geom_node_point(
        aes(
          size = degree(g_ddr), fill = pikk_group,
          shape = node_type
        ),
        colour = "grey30", stroke = 0.5
      ) +
      geom_node_text(aes(label = name),
        repel = TRUE, size = 3.2,
        max.overlaps = 30, fontface = "bold"
      ) +
      scale_fill_manual(values = pikk_colors[ddr_groups], name = "PIKK Group") +
      scale_shape_manual(
        values = c("DEG" = 21, "Non-DE" = 23),
        name = "Node Type"
      ) +
      scale_size_continuous(range = c(4, 15), name = "Degree") +
      guides(fill = guide_legend(override.aes = list(shape = 21, size = 5))) +
      labs(
        title = "DDR & Checkpoint Signalling Sub-network (ATR + ATM + PRKDC)",
        subtitle = paste0(
          vcount(g_ddr), " nodes, ",
          ecount(g_ddr), " edges | STRING ≥ 0.7"
        )
      ) +
      theme_graph(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right"
      )
    save_pdf_png(p_ddr, "PPI_subnetwork_DDR_checkpoint", width = 14, height = 10)
  } else {
    message("Skipping DDR sub-network plot (no edges at this threshold).")
  }


  # --- 8.5 Community-coloured network ------------------------------------------
  n_communities <- max(V(g_core)$community, na.rm = TRUE)
  comm_palette <- colorRampPalette(brewer.pal(
    min(n_communities, 12), "Set3"
  ))(n_communities)

  p_comm <- ggraph(tg_core, layout = "manual", x = x, y = y) +
    geom_edge_link(aes(width = combined_score),
      alpha = 0.15,
      colour = "grey60", show.legend = FALSE
    ) +
    scale_edge_width(range = c(0.3, 1.5)) +
    geom_node_point(
      aes(
        size = degree(g_core),
        fill = factor(community)
      ),
      shape = 21, colour = "grey30", stroke = 0.4
    ) +
    geom_node_text(aes(label = ifelse(name %in% label_genes, name, "")),
      repel = TRUE, size = 3, max.overlaps = 20,
      fontface = "bold"
    ) +
    scale_fill_manual(values = comm_palette, name = "Community") +
    scale_size_continuous(range = c(3, 14), name = "Degree") +
    labs(
      title = "PIKK DEG PPI Network — Louvain Communities",
      subtitle = paste0(n_communities, " communities detected")
    ) +
    theme_graph(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
  save_pdf_png(p_comm, "PPI_communities_network", width = 13, height = 10)
} else {
  # ---- Fallback: base igraph plots ------------------------------------------
  cat("Using base igraph plotting (install ggraph for better figures).\n")

  node_cols_pikk <- pikk_colors[V(g_core)$pikk_group]
  node_sizes <- rescale(degree(g_core), to = c(5, 20))

  # 8.1 PIKK group network
  pdf(file.path(out_plots, "PPI_network_by_PIKK_group.pdf"), width = 12, height = 10)
  plot(g_core,
    layout = layout_coords,
    vertex.color = node_cols_pikk, vertex.size = node_sizes,
    vertex.label = ifelse(V(g_core)$name %in% label_genes, V(g_core)$name, NA),
    vertex.label.cex = 0.7, vertex.label.color = "black",
    vertex.frame.color = "grey50",
    edge.width = rescale(E(g_core)$combined_score, to = c(0.5, 3)),
    edge.color = adjustcolor("grey60", alpha.f = 0.4),
    main = "PIKK DEG PPI Network — by PIKK Group"
  )
  legend("bottomright",
    legend = names(pikk_colors), fill = pikk_colors,
    cex = 0.8, title = "PIKK Group", bty = "n"
  )
  dev.off()

  png(file.path(out_plots, "PPI_network_by_PIKK_group.png"),
    width = 3600, height = 3000, res = 300
  )
  plot(g_core,
    layout = layout_coords,
    vertex.color = node_cols_pikk, vertex.size = node_sizes,
    vertex.label = ifelse(V(g_core)$name %in% label_genes, V(g_core)$name, NA),
    vertex.label.cex = 0.7, vertex.label.color = "black",
    vertex.frame.color = "grey50",
    edge.width = rescale(E(g_core)$combined_score, to = c(0.5, 3)),
    edge.color = adjustcolor("grey60", alpha.f = 0.4),
    main = "PIKK DEG PPI Network — by PIKK Group"
  )
  legend("bottomright",
    legend = names(pikk_colors), fill = pikk_colors,
    cex = 0.8, title = "PIKK Group", bty = "n"
  )
  dev.off()

  cat("  Saved: PPI_network_by_PIKK_group (PDF + PNG) [base igraph]\n")
  cat("  (Install ggraph for additional network plots)\n")
}


# ============================================================================
# 9. CENTRALITY VISUALISATIONS
# ============================================================================
cat("\n--- Centrality Visualisations ---\n")

# --- 9.1 Degree distribution barplot -----------------------------------------
degree_dist <- data.frame(degree = degree(g_core)) %>%
  dplyr::count(degree, name = "count")

p_degdist <- ggplot(degree_dist, aes(x = degree, y = count)) +
  geom_col(fill = "#4575B4", colour = "grey30", width = 0.8) +
  labs(
    title = "Degree Distribution — PIKK DEG PPI Network",
    x = "Degree (number of interactions)", y = "Number of genes"
  ) +
  theme_pub
save_pdf_png(p_degdist, "degree_distribution_barplot", width = 8, height = 5)


# --- 9.2 Degree vs betweenness scatter ----------------------------------------
p_scatter <- ggplot(topo_df, aes(x = degree, y = betweenness)) +
  geom_point(aes(fill = PIKK_Group, size = abs(log2FoldChange)),
    shape = 21, colour = "grey30", stroke = 0.5, alpha = 0.85
  ) +
  geom_text_repel(
    data = topo_df %>%
      dplyr::filter(degree >= quantile(degree, 0.85) |
        betweenness >= quantile(betweenness, 0.85, na.rm = TRUE)),
    aes(label = gene_symbol),
    size = 3.2, fontface = "bold", max.overlaps = 20
  ) +
  scale_fill_manual(values = pikk_colors, name = "PIKK Group") +
  scale_size_continuous(
    range = c(2, 10),
    name = expression("|" * log[2] * "FC|")
  ) +
  labs(
    title = "Degree vs Betweenness Centrality",
    subtitle = "Hub gene candidates (upper-right quadrant)",
    x = "Degree centrality", y = "Betweenness centrality (normalised)"
  ) +
  theme_pub +
  # Add quadrant lines at median
  geom_hline(
    yintercept = median(topo_df$betweenness, na.rm = TRUE),
    linetype = "dashed", colour = "grey50", linewidth = 0.4
  ) +
  geom_vline(
    xintercept = median(topo_df$degree, na.rm = TRUE),
    linetype = "dashed", colour = "grey50", linewidth = 0.4
  )
save_pdf_png(p_scatter, "centrality_scatter_degree_vs_betweenness",
  width = 10, height = 8
)


# --- 9.3 Centrality heatmap (all metrics) ------------------------------------
heat_data <- topo_df %>%
  dplyr::select(
    gene_symbol, degree, betweenness, closeness,
    eigenvector, clustering_coeff
  ) %>%
  tibble::column_to_rownames("gene_symbol")

# Scale each metric to [0, 1]
heat_scaled <- as.data.frame(lapply(heat_data, function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) {
    return(rep(0.5, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}))
rownames(heat_scaled) <- rownames(heat_data)

# Row annotation: PIKK group
anno_row <- data.frame(
  PIKK_Group = topo_df$PIKK_Group[match(
    rownames(heat_scaled),
    topo_df$gene_symbol
  )],
  row.names = rownames(heat_scaled)
)
anno_colors <- list(PIKK_Group = pikk_colors)

heat_palette <- colorRampPalette(c("#F7FBFF", "#4292C6", "#08306B"))(100)

pdf(file.path(out_plots, "centrality_heatmap.pdf"),
  width = 8,
  height = max(6, 0.22 * nrow(heat_scaled) + 2)
)
pheatmap(as.matrix(heat_scaled),
  color = heat_palette,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = anno_row,
  annotation_colors = anno_colors,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 10,
  border_color = NA,
  main = "Centrality Metrics — PIKK DEG PPI Network (min–max scaled)"
)
dev.off()

png(file.path(out_plots, "centrality_heatmap.png"),
  width = 2400,
  height = max(1800, 70 * nrow(heat_scaled)),
  res = 300
)
pheatmap(as.matrix(heat_scaled),
  color = heat_palette,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = anno_row,
  annotation_colors = anno_colors,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 10,
  border_color = NA,
  main = "Centrality Metrics — PIKK DEG PPI Network (min–max scaled)"
)
dev.off()
cat("  Saved: centrality_heatmap (PDF + PNG)\n")


# --- 9.4 Top-20 genes by degree centrality -----------------------------------
top20_degree <- topo_df %>%
  dplyr::arrange(desc(degree)) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::mutate(gene_symbol = factor(gene_symbol,
    levels = rev(gene_symbol)
  ))

p_top_deg <- ggplot(top20_degree, aes(x = degree, y = gene_symbol)) +
  geom_col(aes(fill = PIKK_Group), colour = "grey30", width = 0.75) +
  geom_text(aes(label = degree), hjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = pikk_colors, name = "PIKK Group") +
  labs(
    title = "Top-20 Genes by Degree Centrality",
    x = "Degree (number of interactions)", y = NULL
  ) +
  theme_pub +
  theme(axis.text.y = element_text(face = "bold"))
save_pdf_png(p_top_deg, "top20_degree_barplot", width = 9, height = 7)


# --- 9.5 Top-20 genes by betweenness centrality ------------------------------
top20_between <- topo_df %>%
  dplyr::arrange(desc(betweenness)) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::mutate(gene_symbol = factor(gene_symbol,
    levels = rev(gene_symbol)
  ))

p_top_btwn <- ggplot(
  top20_between,
  aes(x = betweenness, y = gene_symbol)
) +
  geom_col(aes(fill = PIKK_Group), colour = "grey30", width = 0.75) +
  geom_text(aes(label = round(betweenness, 3)),
    hjust = -0.1,
    size = 3.5, fontface = "bold"
  ) +
  scale_fill_manual(values = pikk_colors, name = "PIKK Group") +
  labs(
    title = "Top-20 Genes by Betweenness Centrality",
    x = "Betweenness centrality (normalised)", y = NULL
  ) +
  theme_pub +
  theme(axis.text.y = element_text(face = "bold"))
save_pdf_png(p_top_btwn, "top20_betweenness_barplot", width = 9, height = 7)


# --- 9.6 PIKK group connectivity boxplot -------------------------------------
p_group_box <- ggplot(topo_df, aes(x = PIKK_Group, y = degree)) +
  geom_boxplot(aes(fill = PIKK_Group),
    alpha = 0.7, outlier.shape = NA,
    colour = "grey30"
  ) +
  geom_jitter(aes(colour = regulation),
    width = 0.2, size = 2.5,
    alpha = 0.8
  ) +
  scale_fill_manual(values = pikk_colors, guide = "none") +
  scale_colour_manual(values = reg_colors, name = "Regulation") +
  labs(
    title = "Degree Distribution by PIKK Family Group",
    x = "PIKK Group", y = "Degree (number of interactions)"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(face = "bold", angle = 0))
save_pdf_png(p_group_box, "PIKK_group_connectivity_boxplot",
  width = 9, height = 6
)


# ============================================================================
# 10. SESSION SUMMARY
# ============================================================================
cat("\n================================================================\n")
cat("STEP 2.2 — PPI NETWORK ANALYSIS COMPLETE\n")
cat("================================================================\n")
cat("Core network nodes:        ", vcount(g_core), "\n")
cat("Core network edges:        ", ecount(g_core), "\n")
cat("Sub-network nodes:         ", vcount(g_sub), "\n")
cat("Sub-network edges:         ", ecount(g_sub), "\n")
cat("Louvain communities:       ", max(community_df$community), "\n")
cat("STRING confidence threshold:", string_confidence / 1000, "\n")
cat("\nOutputs written to:\n")
cat("  Tables: ", out_tables, "\n")
cat("  Plots:  ", out_plots, "\n")
cat("  Network:", out_network, "\n")

# Save session info
writeLines(
  capture.output(sessionInfo()),
  file.path(out_tables, "session_info_network_analysis.txt")
)

cat("\nSession info saved.\n")
cat("================================================================\n")
