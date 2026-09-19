# ============================================================================
# PIKK Single-Cell Integration — Script 01: Seurat Object Construction
# Script:  01_sc_seurat_setup.R
# Purpose: Load TISCH2 GSE103322 h5 expression matrix and cell metadata,
#          construct a fully annotated Seurat object with pre-computed UMAP
#          coordinates, and save for downstream analysis.
#
# Input:   HNSC_GSE103322_expression.h5       (~5,900 cells × 18,242 genes)
#          HNSC_GSE103322_CellMetainfo_table.tsv
#
# Output:  05_results/single_cell_outputs/rds/oscc_gse103322_seurat.rds
#          05_results/single_cell_outputs/tables/sc_qc_cell_counts.tsv
#
# Dataset: Puram et al., Cell 2017 — 18 treatment-naïve HNSCC patients
#          Re-annotated by TISCH2 (major/minor lineage, malignancy)
# ============================================================================

cat("\n================================================================\n")
cat("  PIKK scRNA-seq Integration — Seurat Object Construction\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(hdf5r)
  library(dplyr)
  library(readr)
  library(tibble)
})

# ---- Paths ------------------------------------------------------------------
base_dir <- "C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project"
sc_dir   <- file.path(base_dir, "01_data/single_cell")
out_dir  <- file.path(base_dir, "05_results/single_cell_outputs")
rds_dir  <- file.path(out_dir, "rds")
tbl_dir  <- file.path(out_dir, "tables")
fig_dir  <- file.path(out_dir, "figures")

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

h5_file   <- file.path(sc_dir, "HNSC_GSE103322_expression.h5")
meta_file <- file.path(sc_dir, "HNSC_GSE103322_CellMetainfo_table.tsv")

# ---- 1. Load Cell Metadata --------------------------------------------------
cat("[1/6] Loading cell metadata ...\n")
cell_meta <- read_tsv(meta_file, show_col_types = FALSE)
cat("       Loaded", nrow(cell_meta), "cells with", ncol(cell_meta), "metadata columns\n")
cat("       Columns:", paste(colnames(cell_meta), collapse = ", "), "\n")

# ---- 2. Load h5 Expression Matrix -------------------------------------------
cat("\n[2/6] Inspecting and loading h5 expression matrix ...\n")

h5 <- H5File$new(h5_file, mode = "r")

# Inspect h5 structure to determine format
cat("       h5 top-level contents:", paste(h5$names, collapse = ", "), "\n")

# TISCH2 h5 files typically store data in one of these formats:
#   Format A: /X (dense matrix), /obs (cell names), /var (gene names)
#   Format B: /matrix/data, /matrix/indices, /matrix/indptr (sparse CSC)
#   Format C: Direct dense matrix with gene_names and cell_names datasets

# Attempt to read based on available structure
h5_names <- h5$names

if ("matrix" %in% h5_names) {
  # 10X-style sparse matrix format
  cat("       Detected: 10X/sparse matrix format\n")
  
  # Check sub-structure
  matrix_group <- h5[["matrix"]]
  cat("       matrix/ contents:", paste(matrix_group$names, collapse = ", "), "\n")
  
  data_vec    <- matrix_group[["data"]]$read()
  indices_vec <- matrix_group[["indices"]]$read()
  indptr_vec  <- matrix_group[["indptr"]]$read()
  shape       <- matrix_group[["shape"]]$read()
  
  # Gene names
  if ("features" %in% matrix_group$names) {
    feat_group <- matrix_group[["features"]]
    if ("name" %in% feat_group$names) {
      gene_names <- feat_group[["name"]]$read()
    } else if ("id" %in% feat_group$names) {
      gene_names <- feat_group[["id"]]$read()
    }
  } else if ("gene_names" %in% matrix_group$names) {
    gene_names <- matrix_group[["gene_names"]]$read()
  }
  
  # Cell barcodes
  if ("barcodes" %in% matrix_group$names) {
    cell_names <- matrix_group[["barcodes"]]$read()
  }
  
  # Reconstruct sparse matrix (CSC format: genes × cells)
  sparse_mat <- Matrix::sparseMatrix(
    i        = indices_vec + 1L,   # Convert 0-indexed to 1-indexed
    p        = indptr_vec,
    x        = data_vec,
    dims     = shape,
    repr     = "C"
  )
  rownames(sparse_mat) <- gene_names
  colnames(sparse_mat) <- cell_names
  
  expr_matrix <- sparse_mat
  cat("       Loaded sparse matrix:", nrow(expr_matrix), "genes ×", ncol(expr_matrix), "cells\n")
  
} else if ("X" %in% h5_names) {
  # AnnData-style format
  cat("       Detected: AnnData-style h5 format\n")
  
  X_dataset <- h5[["X"]]
  
  # Check if X is a group (sparse) or dataset (dense)
  if (inherits(X_dataset, "H5Group")) {
    cat("       X is sparse (group) — reading components ...\n")
    data_vec    <- X_dataset[["data"]]$read()
    indices_vec <- X_dataset[["indices"]]$read()
    indptr_vec  <- X_dataset[["indptr"]]$read()
    
    # Get gene and cell names
    if ("var" %in% h5_names) {
      var_group <- h5[["var"]]
      if ("_index" %in% var_group$names) {
        gene_names <- var_group[["_index"]]$read()
      } else if ("index" %in% var_group$names) {
        gene_names <- var_group[["index"]]$read()
      } else {
        gene_names <- var_group[[var_group$names[1]]]$read()
      }
    }
    if ("obs" %in% h5_names) {
      obs_group <- h5[["obs"]]
      if ("_index" %in% obs_group$names) {
        cell_names <- obs_group[["_index"]]$read()
      } else if ("index" %in% obs_group$names) {
        cell_names <- obs_group[["index"]]$read()
      } else {
        cell_names <- obs_group[[obs_group$names[1]]]$read()
      }
    }
    
    n_genes <- length(gene_names)
    n_cells <- length(cell_names)
    
    sparse_mat <- Matrix::sparseMatrix(
      i    = indices_vec + 1L,
      p    = indptr_vec,
      x    = data_vec,
      dims = c(n_genes, n_cells),
      repr = "C"
    )
    rownames(sparse_mat) <- gene_names
    colnames(sparse_mat) <- cell_names
    expr_matrix <- sparse_mat
    
  } else {
    cat("       X is dense — reading full matrix ...\n")
    dense_mat <- X_dataset$read()
    
    # Get gene and cell names
    if ("var" %in% h5_names) {
      var_group <- h5[["var"]]
      if ("_index" %in% var_group$names) {
        gene_names <- var_group[["_index"]]$read()
      } else if ("index" %in% var_group$names) {
        gene_names <- var_group[["index"]]$read()
      } else {
        gene_names <- var_group[[var_group$names[1]]]$read()
      }
    }
    if ("obs" %in% h5_names) {
      obs_group <- h5[["obs"]]
      if ("_index" %in% obs_group$names) {
        cell_names <- obs_group[["_index"]]$read()
      } else if ("index" %in% obs_group$names) {
        cell_names <- obs_group[["index"]]$read()
      } else {
        cell_names <- obs_group[[obs_group$names[1]]]$read()
      }
    }
    
    # dense_mat is cells × genes (AnnData convention) — transpose to genes × cells
    if (nrow(dense_mat) == length(cell_names) && ncol(dense_mat) == length(gene_names)) {
      expr_matrix <- t(dense_mat)
    } else {
      expr_matrix <- dense_mat
    }
    rownames(expr_matrix) <- gene_names
    colnames(expr_matrix) <- cell_names
  }
  
  cat("       Loaded matrix:", nrow(expr_matrix), "genes ×", ncol(expr_matrix), "cells\n")
  
} else {
  # Fallback: try to read first available dataset
  cat("       h5 structure:", paste(h5_names, collapse = ", "), "\n")
  stop("Unrecognised h5 format. Please inspect the file manually with h5$ls(recursive = TRUE)")
}

h5$close_all()
cat("       h5 file closed.\n")

# ---- 3. Align Metadata with Expression Matrix --------------------------------
cat("\n[3/6] Aligning metadata with expression matrix ...\n")

# Ensure cell names match between expression and metadata
meta_cells <- cell_meta$Cell
expr_cells <- colnames(expr_matrix)

common_cells <- intersect(meta_cells, expr_cells)
cat("       Cells in metadata:", length(meta_cells), "\n")
cat("       Cells in expression:", length(expr_cells), "\n")
cat("       Common cells:", length(common_cells), "\n")

if (length(common_cells) == 0) {
  # Try matching without the metadata Cell column directly — sometimes formatting differs
  cat("       WARNING: No direct match. Attempting fuzzy alignment ...\n")
  # Check if cell names in expression need trimming or reformatting
  cat("       First 3 metadata cells:  ", paste(head(meta_cells, 3), collapse = " | "), "\n")
  cat("       First 3 expression cells:", paste(head(expr_cells, 3), collapse = " | "), "\n")
  stop("Cell name mismatch between metadata and expression matrix. Manual inspection required.")
}

# Subset and align
expr_matrix <- expr_matrix[, common_cells]
cell_meta   <- cell_meta %>% filter(Cell %in% common_cells)
cell_meta   <- cell_meta[match(common_cells, cell_meta$Cell), ]

cat("       Final aligned dataset:", nrow(expr_matrix), "genes ×", ncol(expr_matrix), "cells\n")

# ---- 4. Construct Seurat Object ----------------------------------------------
cat("\n[4/6] Constructing Seurat object ...\n")

# Add clean syntactic column names alongside original TISCH2 columns
# (Prevents Seurat LabelClusters, DimPlot, and ggplot parsing errors)
cell_meta <- cell_meta %>%
  mutate(
    CellType_Major      = `Celltype (major-lineage)`,
    CellType_Malignancy = `Celltype (malignancy)`,
    CellType_Minor      = `Celltype (minor-lineage)`
  )

# Create Seurat object with the expression matrix
seu <- CreateSeuratObject(
  counts    = expr_matrix,
  project   = "OSCC_PIKK_GSE103322",
  meta.data = cell_meta %>% column_to_rownames("Cell")
)

cat("       Seurat object created:", ncol(seu), "cells,", nrow(seu), "genes\n")

# Inject TISCH2-provided UMAP coordinates as a dimensional reduction
umap_coords <- cell_meta %>%
  select(Cell, UMAP_1, UMAP_2) %>%
  column_to_rownames("Cell") %>%
  as.matrix()

# Ensure UMAP coordinate order matches Seurat cell order
umap_coords <- umap_coords[colnames(seu), ]
colnames(umap_coords) <- c("umap_1", "umap_2")

seu[["umap"]] <- CreateDimReducObject(
  embeddings = umap_coords,
  key        = "umap_",
  assay      = DefaultAssay(seu)
)

cat("       TISCH2 UMAP coordinates injected into Seurat object\n")

# Normalise data (log-normalisation for visualisation and AUCell downstream)
cat("       Running NormalizeData (LogNormalize, scale.factor = 10000) ...\n")
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

# ---- 5. QC Summary -----------------------------------------------------------
cat("\n[5/6] Generating QC summary tables ...\n")

# --- Cell counts by cell type (major lineage) ---
ct_major <- seu@meta.data %>%
  group_by(CellType_Major) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(desc(n_cells))

cat("\n       --- Cell Counts by Major Lineage ---\n")
print(as.data.frame(ct_major), row.names = FALSE)

# --- Cell counts by malignancy status ---
ct_malig <- seu@meta.data %>%
  group_by(CellType_Malignancy) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(desc(n_cells))

cat("\n       --- Cell Counts by Malignancy Status ---\n")
print(as.data.frame(ct_malig), row.names = FALSE)

# --- Cell counts by patient ---
ct_patient <- seu@meta.data %>%
  group_by(Patient) %>%
  summarise(
    n_cells    = n(),
    n_tumour   = sum(Source == "Tumor", na.rm = TRUE),
    n_lymph    = sum(Source == "tLN", na.rm = TRUE),
    n_normal   = sum(Source == "Normal", na.rm = TRUE),
    stage      = first(TNMstage),
    .groups    = "drop"
  ) %>%
  arrange(desc(n_cells))

cat("\n       --- Cell Counts by Patient ---\n")
print(as.data.frame(ct_patient), row.names = FALSE)

# --- Diamond Panel gene presence check ---
diamond_panel <- c("PLK1", "CDK2", "TOPBP1", "RAD51", "FANCI", "KAT2B", "DEPTOR")
present <- diamond_panel %in% rownames(seu)
cat("\n       --- Diamond Panel Gene Presence ---\n")
for (i in seq_along(diamond_panel)) {
  status <- ifelse(present[i], "FOUND", "MISSING")
  cat("       ", diamond_panel[i], ":", status, "\n")
}

if (!all(present)) {
  missing_genes <- diamond_panel[!present]
  cat("\n       WARNING: Missing genes:", paste(missing_genes, collapse = ", "), "\n")
  cat("       These genes will be excluded from downstream analyses.\n")
}

# --- Full PIKK interactome gene presence check ---
pikk_universe <- c(
  # ATR axis
  "ATR", "CHEK1", "RAD51", "PLK1", "EXO1", "CDC45", "E2F1", "FANCI",
  "BRCA1", "AURKA", "AURKB", "TOPBP1", "CHEK2", "H2AFX",
  # ATM axis
  "ATM", "BRCA2", "CDK2", "MSH6", "POLD1", "RRAGD",
  # PRKDC axis
  "PRKDC", "EGFR", "HOXB7", "PARP1",
  # TRRAP axis
  "TRRAP", "KAT2B", "RUVBL1", "KAT2A", "SUPT7L", "TAF2",
  # mTOR axis
  "MTOR", "DEPTOR", "ACTL6A", "RPTOR", "RICTOR", "TTI1",
  # SMG1 axis
  "SMG1", "NCBP2"
)

pikk_present <- pikk_universe[pikk_universe %in% rownames(seu)]
pikk_missing <- pikk_universe[!pikk_universe %in% rownames(seu)]
cat("\n       --- PIKK Universe Gene Check ---\n")
cat("       Present:", length(pikk_present), "of", length(pikk_universe), "genes\n")
if (length(pikk_missing) > 0) {
  cat("       Missing:", paste(pikk_missing, collapse = ", "), "\n")
}

# Save QC tables
qc_combined <- bind_rows(
  ct_major %>% mutate(category = "major_lineage") %>% rename(group = CellType_Major),
  ct_malig %>% mutate(category = "malignancy")    %>% rename(group = CellType_Malignancy)
)
write_tsv(qc_combined, file.path(tbl_dir, "sc_qc_cell_counts.tsv"))
write_tsv(ct_patient,  file.path(tbl_dir, "sc_qc_patient_summary.tsv"))
cat("\n       QC tables saved.\n")

# ---- 6. Save Seurat Object ---------------------------------------------------
cat("\n[6/6] Saving Seurat object ...\n")
rds_path <- file.path(rds_dir, "oscc_gse103322_seurat.rds")
saveRDS(seu, rds_path)
cat("       Saved to:", rds_path, "\n")
cat("       Object size:", format(object.size(seu), units = "auto"), "\n")

cat("\n================================================================\n")
cat("  Script 01 complete:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("  Seurat object:", ncol(seu), "cells ×", nrow(seu), "genes\n")
cat("  Diamond Panel genes found:", sum(present), "of", length(diamond_panel), "\n")
cat("  PIKK universe genes found:", length(pikk_present), "of", length(pikk_universe), "\n")
cat("  UMAP coordinates: TISCH2 pre-computed (injected)\n")
cat("  Next: Run 02_sc_pikk_expression_cartography.R\n")
cat("================================================================\n\n")
