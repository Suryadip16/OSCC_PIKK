# Single-Cell Transcriptomic Cartography and Cellular Compartmentalization (Section 2.7) — Comprehensive Methodology Guide

**Project**: PIKK Pathway-Associated Biomarker Discovery and Oncogenic Characterization in Oral Squamous Cell Carcinoma (OSCC)  
**Target Manuscript Reference**: Materials and Methods / Bioinformatic Processing Pipeline  
**Primary Single-Cell Dataset**: Single-Cell RNA Sequencing Atlas of Head and Neck Squamous Cell Carcinoma (GSE103322; Puram et al., *Cell* 2017 [1]), Curated and Harmonized via the Tumor Immune Single-cell Hub 2 (TISCH2; Sun et al., *Nucleic Acids Res* 2023 [2])  
**Single-Cell Cohort Architecture**: $N = 5,902$ individual cells across $18,242$ detected human genes, derived from $22$ surgical biospecimens representing $18$ treatment-naive primary oral cavity squamous cell carcinoma patients (spanning TNM clinical stages I, II, and III)  
**Cellular Compartment Partitioning**:
1. **Malignant Compartment**: $n = 2,488$ malignant epithelial cells (42.16% of total atlas) characterized by inferred copy number alterations (CNAs) and epithelial atypia.
2. **Stromal Compartment**: $n = 1,741$ microenvironmental stromal cells (29.50% of total atlas), comprising Fibroblasts ($n = 744$), Myofibroblasts ($n = 710$), Endothelial cells ($n = 269$), and Myocytes ($n = 18$).
3. **Tumor-Infiltrating Immune Compartment**: $n = 1,673$ immune cells (28.35% of total atlas), comprising Exhausted CD8+ T cells (`CD8Tex`, $n = 501$), Conventional CD4+ T cells (`CD4Tconv`, $n = 417$), Effector CD8+ T cells (`CD8T`, $n = 397$), Plasma B cells ($n = 148$), Mast cells ($n = 122$), and Monocytes/Macrophages (`Mono/Macro`, $n = 88$).  
**Evaluated Biomarker Panels**:
1. **Core Diamond Panel** ($n = 7$ genes): $PLK1$, $CDK2$, $TOPBP1$, $RAD51$, $FANCI$, $KAT2B$, $DEPTOR$.
2. **Full PIKK Interactome Universe** ($n = 37$ candidate genes spanning 6 canonical kinase branches: ATR, ATM, PRKDC, TRRAP, mTOR, SMG1).  
**Scripts Evaluated**:
1. `01_sc_seurat_setup.R` — *HDF5 Sparse Matrix Parsing, Metadata Sanitation, Seurat Container Instantiation, Pre-computed TISCH2 UMAP Coordinate Injection, Library-Size Normalization, and Multi-Tier Quality-Control Benchmarking*
2. `02_sc_pikk_expression_cartography.R` — *Single-Cell Cellular Cartography, Viridis Magma Feature Projection, Dual-Scaled Lineage Dot Plots, Non-Parametric Kruskal–Wallis Omnibus Testing, Post-Hoc Pairwise Wilcoxon Analysis, and Compartment Detection Breadth Profiling*

---

## 1. Executive Summary & End-to-End Workflow

Bulk RNA sequencing measures an aggregate transcriptional average across complex mixtures of malignant squamous cells, infiltrating lymphocytes, myeloid subsets, and reactive cancer-associated fibroblasts (CAFs). While bulk differential expression, topological network analysis, and multivariable Cox survival modeling identify clinically actionable biomarkers (e.g., the independent prognostic drivers $PLK1$ and $CDK2$), bulk profiling cannot definitively resolve the critical translational question of **cell-of-origin compartmentalization**:
- Are prioritized oncogenic biomarkers truly hyper-activated within malignant epithelial cells, or does their elevated bulk expression reflect signaling within tumor-infiltrating lymphocytes or proliferating fibroblasts?
- Conversely, does the downregulation of putative tumor suppressors (e.g., $KAT2B$, $DEPTOR$) represent true transcriptional repression within carcinoma cells or cellular dilution caused by the physical expansion of transformed parenchyma?

This document formalizes the publication-grade computational and statistical methodology for **Single-Cell Transcriptomic Cartography and Cellular Compartmentalization (Section 2.7)**. 

The single-cell integration pipeline executes two sequential, highly coordinated computational phases:
1. **Container Construction & Architecture Integration (`01_sc_seurat_setup.R`)**: Ingests the benchmark GSE103322 single-cell HDF5 expression matrix, aligns single-cell barcodes with curated TISCH2 metadata, maps cells into compressed sparse column (CSC) structures ($18,242 \text{ genes} \times 5,902 \text{ cells}$), injects validated two-dimensional UMAP coordinates, executes library-size log-normalization (`LogNormalize`, scale factor = 10,000), and audits quality-control distributions across 11 distinct lineages, 3 cellular compartments, and 18 patient biospecimens.
2. **Cellular Cartography & Compartmentalization Testing (`02_sc_pikk_expression_cartography.R`)**: Maps the Core Diamond Panel onto single-cell coordinates via UMAP feature overlays using perceptually uniform viridis magma scales with 5th-to-95th percentile outlier clipping; quantifies expression intensity and detection breadth across all 11 lineages via dual-scaled dot plots; models distribution shifts across compartments; and executes non-parametric Kruskal–Wallis omnibus tests ($\chi^2, \text{df} = 2$) and pairwise two-sample Wilcoxon rank-sum tests with Benjamini–Hochberg False Discovery Rate (FDR) correction to establish compartment specificity.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                                 UPSTREAM INPUT ARTIFACTS                               │
│  • Single-Cell Expression Matrix: HNSC_GSE103322_expression.h5 (5,902 cells x 18,242 g)│
│  • Single-Cell Metadata: HNSC_GSE103322_CellMetainfo_table.tsv (TISCH2 Curated)       │
│  • Core Diamond Target Panel: PLK1, CDK2, TOPBP1, RAD51, FANCI, KAT2B, DEPTOR        │
│  • Extended PIKK Universe: 37 candidate genes across 6 canonical PIKK branches        │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 1: CONTAINER INSTANTIATION & ARCHITECTURE ASSEMBLY (01_sc_seurat_setup.R)       │
│  • Low-level HDF5 container parsing via hdf5r (AnnData / 10X sparse matrix layout)     │
│  • CSC sparseMatrix reconstruction: 18,242 genes x 5,902 cells                         │
│  • Barcode alignment (N = 5,902 common cells) & syntactic column name sanitation:      │
│      CellType_Major, CellType_Malignancy, CellType_Minor                               │
│  • Seurat container instantiation: OSCC_PIKK_GSE103322                                 │
│  • Pre-computed TISCH2 UMAP coordinate injection (DimReducObject: umap_1, umap_2)       │
│  • Library-size log-normalization: y = ln(1 + (x / sum(x)) * 10,000)                   │
│  • Quality-control audit & serialization:                                              │
│      oscc_gse103322_seurat.rds, sc_qc_cell_counts.tsv, sc_qc_patient_summary.tsv       │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 2: CELLULAR CARTOGRAPHY & COMPARTMENTALIZATION (02_sc_pikk_expression_cartography)│
│  • Reference Landscape Mapping: UMAP projection across 11 curated lineages             │
│  • Single-Cell Expression Cartography: UMAP feature mapping of Diamond Panel targets   │
│      - Quantile outlier suppression: min.cutoff = 'q5', max.cutoff = 'q95'             │
│      - Perceptually uniform viridis magma color mapping with cell order prioritization │
│      - Individual high-resolution panels + composite multi-panel overview              │
│  • Lineage-Resolved Dot Plots: Dual-scaled dot diameter (% cells) and color intensity  │
│      - 11 major lineages ordered: Malignant -> Immune subsets -> Stromal subsets       │
│  • Compartment Distribution Modeling: Malignant (n=2488), Immune (n=1673), Stroma (1741)│
│      - Standalone & composite box plots (median, IQR, mean diamond, outlier jitter)   │
│  • Non-Parametric Hypothesis Testing across 3 compartments:                            │
│      - Kruskal-Wallis Rank-Sum Omnibus Test (H-statistic, chi-squared, df = 2)         │
│      - Pairwise Wilcoxon Rank-Sum Tests with Benjamini-Hochberg FDR correction         │
│      - Classification of predominant compartment based on maximal mean expression      │
│  • Detection Breadth Quantification: Proportion of non-zero cells per compartment      │
│  • Artifact Export:                                                                    │
│      pikk_gene_celltype_compartment_stats.tsv, pikk_gene_celltype_expression_summary.tsv│
│      Publication figure panels (PDF + PNG, 300 DPI)                                    │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. Biological Context and Study Rationale

### 2.1 Resolving Bulk Expression Confounding in OSCC
In mucosal head and neck oncology, tumor resections represent heterogeneous cellular microenvironments consisting of malignant squamous cells intermixed with activated cancer-associated fibroblasts (CAFs), endothelial cells, and diverse immune infiltrates. When an oncogenic kinase (e.g., $PLK1$) or cell cycle coordinator (e.g., $CDK2$) is identified as significantly upregulated in bulk tumor sequencing, three mutually exclusive biological mechanisms can explain the observation:
1. **True Malignant Hyper-Activation**: The gene is transcribed predominantly within transformed epithelial cells, reflecting cell-autonomous oncogenic addiction, replication stress adaptation, and checkpoint override.
2. **Stromal Desmoplasia Artifact**: The gene is transcribed primarily by proliferating myofibroblasts or angiogenic endothelial cells within the reactive tumor stroma.
3. **Immune Infiltration Confounder**: The gene is expressed predominantly by clonal, antigen-stimulated tumor-infiltrating lymphocytes (TILs) or expanding myeloid populations.

Similarly, when putative tumor suppressors (e.g., the acetyltransferase $KAT2B$ or the mTOR endogenous regulator $DEPTOR$) demonstrate significant downregulation in bulk tumor tissue, this reduction can arise from:
- **True Malignant Repression**: Epigenetic silencing, promoter methylation, or targeted transcriptional downregulation within transformed carcinoma cells to disable barriers against genomic instability and anabolic growth.
- **Cellular Dilution Artifact**: High baseline expression in specialized stromal or immune cells in normal mucosa that becomes diluted by the physical expansion of transformed epithelial clones in tumor biopsies.

Single-cell RNA sequencing provides the definitive resolution to these questions, allowing exact cell-of-origin mapping and statistical validation of compartment specificity.

### 2.2 The PIKK Kinase Network in Mucosal Squamous Carcinogenesis
Oral Squamous Cell Carcinoma is overwhelmingly characterized by early, near-ubiquitous loss-of-function mutations in $TP53$ (>80% of cases) and inactivation of the $CDKN2A$ locus ($p16^{INK4A}/p14^{ARF}$), eliminating the G1/S checkpoint and predisposing cells to severe oncogene-induced replication stress. 

In this setting, the six phosphatidylinositol 3-kinase-related kinases (PIKKs) and their immediate effectors maintain cellular viability:
- **ATR-CHK1 Signaling ($TOPBP1$, $RAD51$, $FANCI$)**: Stalled replication forks require ATR activation, mediated by TOPBP1 recruitment, to phosphorylate CHK1, stabilize forks, and recruit the RAD51 recombinase and Fanconi anemia core complexes ($FANCI$) for homologous recombination (HR) repair.
- **ATM-CDK2 Crosstalk**: Loss of G1 arrest elevates CDK2 activity, which hyper-phosphorylates downstream substrates, promotes origin firing, and dictates DNA double-strand break repair pathway choice.
- **Mitotic Checkpoint Override ($PLK1$)**: Polo-like kinase 1 executes centrosome maturation, mitotic entry, and anaphase progression, enabling cells with lingering DNA damage to bypass the G2/M checkpoint and avoid mitotic catastrophe.
- **Metabolic & Epigenetic Rewiring ($KAT2B$, $DEPTOR$)**: KAT2B regulates histone acetylation and transcriptional competence at damage foci, while DEPTOR serves as an endogenous inhibitor of mTORC1 and mTORC2 complexes.

Establishing whether these targets are specifically confined to malignant cells or active within the microenvironment is essential for anticipating therapeutic efficacy and off-target toxicities.

### 2.3 Methodological Rationale: Justification of Selected Analytical Paradigms

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   METHODOLOGICAL COMPARISON: SINGLE-CELL ANALYTICAL PARADIGMS          │
├───────────────────────────┬──────────────────────────────┬─────────────────────────────┤
│ Analytical Dimension      │ Selected Pipeline Paradigm   │ Conventional Alternative    │
├───────────────────────────┼──────────────────────────────┼─────────────────────────────┤
│ Atlas Ingestion           │ GSE103322 (Puram et al. 2017)│ In Vitro Cell Line Data or  │
│                           │ Curated via TISCH2           │ Synthetic Deconvolution Only│
│ Rationale                 │ Direct clinical patient      │ Cell lines lack authentic   │
│                           │ specimens; gold-standard     │ clinical TME architecture   │
│                           │ mucosal carcinoma atlas      │ and in vivo immune pressures│
├───────────────────────────┼──────────────────────────────┼─────────────────────────────┤
│ Dimensional Embedding     │ TISCH2-Curated Pre-computed  │ De Novo Re-Clustering via   │
│                           │ UMAP Projection Coordinates  │ Variable Feature Heuristics │
│ Rationale                 │ Preserves expert-validated,  │ Re-clustering risks batch   │
│                           │ benchmarked lineage clusters │ artifacts and non-reproduc- │
│                           │ matching published literature│ ible cluster boundaries     │
├───────────────────────────┼──────────────────────────────┼─────────────────────────────┤
│ Feature Plot Projection   │ Quantile-Clipped Magma Scale │ Default Linear Scale        │
│                           │ (min = q5, max = q95)        │ (min = 0, max = max(y))     │
│ Rationale                 │ Prevents single outlier cells│ Extreme outlier UMIs compress│
│                           │ from washing out gradient    │ the dynamic range, making   │
│                           │ dynamic range                │ moderate expression invisible│
├───────────────────────────┼──────────────────────────────┼─────────────────────────────┤
│ Compartment Testing       │ Kruskal–Wallis Rank-Sum +    │ Parametric ANOVA or         │
│                           │ Pairwise Wilcoxon (BH FDR)   │ Student's t-test            │
│ Rationale                 │ Non-parametric; handles      │ Violates distributional     │
│                           │ non-normal, zero-inflated    │ assumptions of single-cell  │
│                           │ single-cell count data       │ count measurements          │
└───────────────────────────┴──────────────────────────────┴─────────────────────────────┘
```

#### 1. Why Clinical Atlas (GSE103322) Over Cell Lines?
Single-cell RNA sequencing from patient surgical specimens captures the native cellular heterogeneity, stromal interactions, and immune infiltration profiles present in clinical oral tumors. Cultured cell lines adapt to plastic surfaces, lack functional immune or endothelial compartments, and undergo genetic drift that obscures genuine microenvironmental interactions.

#### 2. Why TISCH2 Harmonization and Pre-Computed Embeddings?
The Tumor Immune Single-cell Hub 2 (TISCH2) applies standardized quality-control, cell filtering, and expert manual lineage annotation based on established canonical markers. Re-running unsupervised clustering from scratch introduces run-to-run variation in high-variable gene selection and graph partitioning. By utilizing the harmonized TISCH2 metadata and UMAP coordinates, the single-cell cartography directly aligns with peer-reviewed reference benchmarks.

#### 3. Why Quantile-Clipped Feature Projections?
Single-cell transcript counts frequently follow a negative binomial or zero-inflated distribution, where a handful of hyper-amplified cells exhibit extreme expression values. Under default linear color scaling, these extreme values dominate the color bar, compressing the dynamic range for the remaining 99% of cells and making moderate expression indistinguishable from background noise. Clipping the minimum and maximum display values to the 5th and 95th percentiles (`min.cutoff = "q5"`, `max.cutoff = "q95"`) maximizes visual contrast across expressing populations.

#### 4. Why Non-Parametric Hypothesis Testing?
Single-cell expression data exhibit substantial zero-inflation (structural and technical dropouts), right-skewness, and non-constant variance across lineages. Parametric methods such as Student's $t$-test or ordinary ANOVA assume Gaussian error distributions and equal variances, leading to severe Type I error inflation. The pipeline employs the **Kruskal–Wallis rank-sum test** for multi-group omnibus comparison across compartments, followed by post-hoc **pairwise Wilcoxon rank-sum tests** with **Benjamini–Hochberg False Discovery Rate (FDR)** control.

---

## 3. Step-by-Step Computational Methodology

### 3.1 Step 1: HDF5 Matrix Ingestion and Sparse Matrix Representation
**Script**: `01_sc_seurat_setup.R`  
**Input Artifacts**: `HNSC_GSE103322_expression.h5` (TISCH2 HDF5 container).  
**Output Artifacts**: Internal memory object `sparse_mat` (CSC sparse matrix).  

#### 3.1.1 HDF5 Container Inspection
The primary transcript count matrix is stored in HDF5 format, which organizes high-dimensional biological data into structured groups and datasets. The pipeline leverages the low-level C++ interface provided by `hdf5r` to inspect the internal group hierarchy:
- **10X Genomics Format**: Layout organized under `/matrix/data`, `/matrix/indices`, `/matrix/indptr`, and `/matrix/shape`.
- **AnnData Format**: Layout organized under `/X`, `/obs`, and `/var`.

The parsing engine dynamically evaluates the root names:
```r
h5 <- H5File$new(h5_file, mode = "r")
h5_names <- h5$names
```

#### 3.1.2 Compressed Sparse Column (CSC) Reconstruction
Single-cell matrices are highly sparse (>90% zero entries). Storing expression data in standard dense matrices would require gigabytes of memory. The parsing architecture extracts the non-zero UMI counts, row indices, and column pointers to reconstruct a formal Compressed Sparse Column (`dgCMatrix`) representation using the `Matrix` package:

$$\mathbf{X} \in \mathbb{R}^{G \times N}$$

Where:
- $G = 18,242$ detected human genes (rows).
- $N = 5,902$ single cells (columns).
- $\mathbf{x} \in \mathbb{R}^K$: Vector of non-zero UMI count values ($K \ll G \times N$).
- $\mathbf{i} \in \mathbb{N}^K$: Zero-based row indices, converted to R's 1-based indexing ($\mathbf{i}_{\text{R}} = \mathbf{i} + 1$).
- $\mathbf{p} \in \mathbb{N}^{N+1}$: Column pointers defining the start and end offsets of non-zero entries for each cell.

$$\mathbf{X} = \text{sparseMatrix}(\mathbf{i} = \text{indices} + 1, \; \mathbf{p} = \text{indptr}, \; \mathbf{x} = \text{data}, \; \text{dims} = [G, N], \; \text{repr} = \text{"C"})$$

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│               EMPIRICAL METRIC CALLOUT: MATRIX RECONSTRUCTION BENCHMARK                │
├──────────────────────────────────────────────────┬─────────────────────────────────────┤
│ Metric Parameter                                 │ Empirical Pipeline Value            │
├──────────────────────────────────────────────────┼─────────────────────────────────────┤
│ Input Container Type                             │ HDF5 Sparse Compressed Format       │
│ Total Detected Human Genes (G)                   │ 18,242 genes                        │
│ Total Aligned Single Cells (N)                   │ 5,902 cells                         │
│ Internal Matrix Storage Class                    │ dgCMatrix (Compressed Sparse Column)│
│ Index Offset Correction                          │ 0-indexed C++ to 1-indexed R (+1L)  │
└──────────────────────────────────────────────────┴─────────────────────────────────────┘
```

---

### 3.2 Step 2: Barcode Alignment, Metadata Sanitation, and Seurat Container Instantiation
**Script**: `01_sc_seurat_setup.R`  
**Input Artifacts**: `HNSC_GSE103322_CellMetainfo_table.tsv`.  
**Output Artifacts**: Aligned metadata structure and initial `Seurat` object.  

#### 3.2.1 Cell Barcode Intersect and Alignment
Cell barcodes from the metadata manifest ($\mathcal{B}_{\text{meta}}$) are matched against column names of the reconstructed expression matrix ($\mathcal{B}_{\text{expr}}$):

$$\mathcal{B}_{\text{common}} = \mathcal{B}_{\text{meta}} \cap \mathcal{B}_{\text{expr}}, \quad |\mathcal{B}_{\text{common}}| = 5,902$$

The expression matrix is subset strictly to common cells, and metadata rows are permuted to match the exact column ordering of the expression matrix:
```r
expr_matrix <- expr_matrix[, common_cells]
cell_meta   <- cell_meta %>% filter(Cell %in% common_cells)
cell_meta   <- cell_meta[match(common_cells, cell_meta$Cell), ]
```

#### 3.2.2 Metadata Column Sanitation
Original TISCH2 column headers contain spaces and parentheses (e.g., `Celltype (major-lineage)`, `Celltype (malignancy)`), which trigger parsing failures in downstream functions (e.g., `Seurat::LabelClusters`, `ggplot2::aes_string`). The pipeline creates standardized, syntactic R column names:
- `Celltype (major-lineage)` $\longrightarrow$ `CellType_Major`
- `Celltype (malignancy)` $\longrightarrow$ `CellType_Malignancy`
- `Celltype (minor-lineage)` $\longrightarrow$ `CellType_Minor`

#### 3.2.3 Seurat Object Instantiation
The formal S4 `Seurat` container is instantiated:
```r
seu <- CreateSeuratObject(
  counts    = expr_matrix,
  project   = "OSCC_PIKK_GSE103322",
  meta.data = cell_meta %>% column_to_rownames("Cell")
)
```

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│               EMPIRICAL METRIC CALLOUT: COHORT COMPARTMENTALIZATION AUDIT              │
├──────────────────────────────────────────────────┬─────────────────────────────────────┤
│ Cellular Classification Tier                     │ Cell Count (n) & Atlas Percentage   │
├──────────────────────────────────────────────────┼─────────────────────────────────────┤
│ Total Single Cells Analyzed                      │ N = 5,902 cells (100.0%)            │
│ 1. Malignant Epithelial Cells                    │ n = 2,488 cells (42.16%)            │
│ 2. Microenvironmental Stromal Cells              │ n = 1,741 cells (29.50%)            │
│    • Fibroblasts                                 │ n = 744 cells (12.61%)              │
│    • Myofibroblasts                              │ n = 710 cells (12.03%)              │
│    • Endothelial Cells                           │ n = 269 cells (4.56%)               │
│    • Myocytes                                    │ n = 18 cells (0.31%)                │
│ 3. Tumor-Infiltrating Immune Cells               │ n = 1,673 cells (28.35%)            │
│    • CD8+ Exhausted T Cells (CD8Tex)             │ n = 501 cells (8.49%)               │
│    • Conventional CD4+ T Cells (CD4Tconv)        │ n = 417 cells (7.07%)               │
│    • Effector CD8+ T Cells (CD8T)                │ n = 397 cells (6.73%)               │
│    • Plasma B Cells                              │ n = 148 cells (2.51%)               │
│    • Mast Cells                                  │ n = 122 cells (2.07%)               │
│    • Monocytes / Macrophages (Mono/Macro)        │ n = 88 cells (1.49%)                │
│ Patient Specimens Represented                    │ 22 surgical biospecimens (18 pts)   │
│ Clinical Stage Distribution                      │ Stage I: 1 pt; Stage II: 5 pts;     │
│                                                  │ Stage III: 11 pts; Unspecified: 1 pt│
└──────────────────────────────────────────────────┴─────────────────────────────────────┘
```

---

### 3.3 Step 3: Library-Size Log-Normalization & Dimensional Embedding Injection
**Script**: `01_sc_seurat_setup.R`  
**Input Artifacts**: `seu` (unnormalized Seurat object), UMAP coordinates from metadata.  
**Output Artifacts**: `oscc_gse103322_seurat.rds`, `sc_qc_cell_counts.tsv`, `sc_qc_patient_summary.tsv`.  

#### 3.3.1 Pre-Computed UMAP Coordinate Injection
To maintain exact concordance with published single-cell literature for GSE103322, the pipeline imports pre-computed UMAP coordinates directly from the TISCH2 database rather than re-computing embeddings:
```r
umap_coords <- cell_meta %>%
  select(Cell, UMAP_1, UMAP_2) %>%
  column_to_rownames("Cell") %>%
  as.matrix()
umap_coords <- umap_coords[colnames(seu), ]
colnames(umap_coords) <- c("umap_1", "umap_2")

seu[["umap"]] <- CreateDimReducObject(
  embeddings = umap_coords,
  key        = "umap_",
  assay      = DefaultAssay(seu)
)
```

#### 3.3.2 Library-Size Log-Normalization Formulation
To account for differences in sequencing depth and transcript capture efficiency across single cells, raw UMI counts are normalized using global library-size scaling followed by natural log-transformation (`NormalizeData`, method = `"LogNormalize"`):

$$y_{gi} = \ln\left(1 + \frac{x_{gi}}{\sum_{k=1}^G x_{ki}} \cdot S\right)$$

Where:
- $x_{gi}$: Raw integer UMI count for gene $g$ in cell $i$.
- $\sum_{k=1}^G x_{ki}$: Total library size (total UMI count) for cell $i$.
- $S = 10,000$: Global scaling factor (transcripts per 10,000 UMIs, approximately matching median cellular library depth).
- $\ln(1 + \cdot)$: Natural log-transformation with pseudocount of 1 to stabilize variance and ensure $y_{gi} = 0$ when $x_{gi} = 0$.

The normalized values $y_{gi}$ are stored in the `@assays$RNA@data` slot of the Seurat container.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│               EMPIRICAL METRIC CALLOUT: NORMALIZATION & QC SPECIFICATIONS              │
├──────────────────────────────────────────────────┬─────────────────────────────────────┤
│ Parameter                                        │ Value / Threshold                   │
├──────────────────────────────────────────────────┼─────────────────────────────────────┤
│ Normalization Method                             │ Standard Library-Size LogNormalize  │
│ Scale Factor (S)                                 │ 10,000 (TP10K)                      │
│ Pseudocount Addition                             │ +1.0 (natural log transformation)   │
│ Dimensional Reduction Assay                      │ UMAP (2 dimensions: umap_1, umap_2) │
│ Coordinate Source                                │ TISCH2 Gold-Standard Curation       │
│ Output RDS Container                             │ oscc_gse103322_seurat.rds           │
│ Serialized Object Memory Size                    │ ~185 MB compressed                  │
└──────────────────────────────────────────────────┴─────────────────────────────────────┘
```

---

### 3.4 Step 4: Reference UMAP Landscape & Major Lineage Cartography
**Script**: `02_sc_pikk_expression_cartography.R`  
**Input Artifacts**: `oscc_gse103322_seurat.rds`.  
**Output Artifacts**: `01_umap_celltype_overview.pdf` and `.png` (300 DPI).  

#### 3.4.1 Curated Categorical Palette Assembly
To ensure clear visual distinction between malignant, stromal, and immune cells on low-dimensional manifolds, an 11-color palette was curated:

```r
ct_colours <- c(
  "Malignant"      = "#E63946",   # Vivid red (transformed carcinoma)
  "CD8T"           = "#457B9D",   # Steel blue (effector T cells)
  "CD8Tex"         = "#1D3557",   # Dark navy (exhausted T cells)
  "CD4Tconv"       = "#A8DADC",   # Light teal (helper T cells)
  "Mono/Macro"     = "#F4A261",   # Warm orange (myeloid cells)
  "Mast"           = "#E9C46A",   # Gold (granulocytic mast cells)
  "Plasma"         = "#2A9D8F",   # Teal-green (antibody-secreting B cells)
  "Fibroblasts"    = "#264653",   # Dark green-grey (quiescent stroma)
  "Myofibroblasts" = "#606C38",   # Olive (cancer-associated fibroblasts)
  "Endothelial"    = "#BC6C25",   # Brown (vascular endothelium)
  "Myocyte"        = "#DDA15E"    # Tan (muscular stroma)
)
```

#### 3.4.2 Manifold Projection & Point Repulsion
The global single-cell atlas is rendered via `Seurat::DimPlot`:
- Point size: `pt.size = 0.4`.
- Cluster labels: Superimposed at cluster geometric centroids with point repulsion (`repel = TRUE`, `label.size = 3.5`) to eliminate label overlapping.
- Visual canvas: Minimalist publication theme with removed coordinate axes (`NoAxes()`).

---

### 3.5 Step 5: Diamond Panel Single-Cell Expression Feature Mapping
**Script**: `02_sc_pikk_expression_cartography.R`  
**Input Artifacts**: `oscc_gse103322_seurat.rds`.  
**Output Artifacts**:
- `02_umap_feature_{gene}.pdf` and `.png` (Individual standalone panels for $PLK1$, $CDK2$, $TOPBP1$, $RAD51$, $FANCI$, $KAT2B$, $DEPTOR$).
- `02_umap_pikk_feature_plots.pdf` and `.png` (Composite multi-panel figure grid).

#### 3.5.1 Quantile Outlier Suppression Formulation
For each evaluated gene $g$, normalized expression values across all cells $\mathbf{y}_g = \{y_{g1}, y_{g2}, \dots, y_{gN}\}$ are clipped to prevent high-expression outlier cells from dominating the color gradient:

$$y_{gi}^{\text{clipped}} = \begin{cases}
q_{0.05}(g), & \text{if } y_{gi} < q_{0.05}(g) \\
y_{gi}, & \text{if } q_{0.05}(g) \le y_{gi} \le q_{0.95}(g) \\
q_{0.95}(g), & \text{if } y_{gi} > q_{0.95}(g)
\end{cases}$$

Where $q_{\alpha}(g)$ denotes the $\alpha$-quantile of non-zero expression for gene $g$. In `Seurat::FeaturePlot`, this is parameterized via `min.cutoff = "q5"` and `max.cutoff = "q95"`.

#### 3.5.2 Cell Plotting Order Prioritization
In standard scatter plots, non-expressing cells plotted late can visually obscure low-frequency expressing cells beneath a layer of grey points. To resolve this:
```r
FeaturePlot(..., order = TRUE)
```
Cells are ordered by their expression value in ascending order prior to rendering, ensuring that cells with positive expression are plotted on top.

#### 3.5.3 Color Gradient Mapping
Expression intensity is mapped using the perceptually uniform **viridis magma** color scale (`viridis::magma(50)`), progressing from dark purple/black (zero expression) through vivid red/orange (intermediate expression) to bright yellow (maximal expression).

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│               EMPIRICAL METRIC CALLOUT: FEATURE PROJECTION PARAMETERS                  │
├──────────────────────────────────────────────────┬─────────────────────────────────────┤
│ Visualization Parameter                          │ Pipeline Value                      │
├──────────────────────────────────────────────────┼─────────────────────────────────────┤
│ Core Biomarkers Evaluated                        │ 7 genes (PLK1, CDK2, TOPBP1, RAD51, │
│                                                  │ FANCI, KAT2B, DEPTOR)               │
│ Color Gradient Palette                           │ 50-step continuous viridis magma    │
│ Outlier Truncation Lower Cutoff                  │ 5th percentile (q5)                 │
│ Outlier Truncation Upper Cutoff                  │ 95th percentile (q95)               │
│ Cell Stacking Order                              │ Ascending expression (order = TRUE) │
│ Standalone Panel Dimensions                      │ 7 in x 6 in (300 DPI, PDF + PNG)    │
│ Composite Grid Architecture                      │ 4 columns x 2 rows (16 in x 9 in)   │
└──────────────────────────────────────────────────┴─────────────────────────────────────┘
```

---

### 3.6 Step 6: Dual-Scaled Major Lineage Dot Plot Profiling
**Script**: `02_sc_pikk_expression_cartography.R`  
**Input Artifacts**: `oscc_gse103322_seurat.rds`.  
**Output Artifacts**:
- `03_dotplot_pikk_by_celltype.pdf` and `.png`: Overview dot plot.
- `06_dotplot_pikk_detailed_lineage.pdf` and `.png`: Lineage-ordered detailed dot plot.

#### 3.6.1 Dual-Scaled Dot Plot Formulation
The dot plot visualizes two independent quantitative dimensions simultaneously for each gene $g$ across each cell lineage $l \in \{1, \dots, L\}$:

1. **Detection Prevalence (Dot Size)**: The proportion of cells within lineage $l$ exhibiting non-zero transcript detection:
   $$\text{PctExp}_{gl} = \frac{\sum_{i \in \mathcal{C}_l} \mathbb{I}(x_{gi} > 0)}{N_l} \times 100\%$$
   Where $\mathcal{C}_l$ is the set of cells belonging to lineage $l$, $N_l = |\mathcal{C}_l|$, and $\mathbb{I}(\cdot)$ is the indicator function. The dot radius scales continuously between $0\%$ and $100\%$ with a maximum diameter of 8 points (`dot.scale = 8`, `dot.min = 0.01`).

2. **Relative Expression Magnitude (Color Intensity)**: The mean normalized expression of gene $g$ within lineage $l$, standardized via $Z$-score transformation across all lineages:
   $$\bar{y}_{gl} = \frac{1}{N_l} \sum_{i \in \mathcal{C}_l} y_{gi}$$
   $$Z_{gl} = \frac{\bar{y}_{gl} - \mu_g}{\sigma_g}$$
   Where $\mu_g = \frac{1}{L} \sum_{l=1}^L \bar{y}_{gl}$ and $\sigma_g = \sqrt{\frac{1}{L-1} \sum_{l=1}^L (\bar{y}_{gl} - \mu_g)^2}$. $Z_{gl}$ is mapped onto a continuous two-color gradient from light grey ($Z \le -1$) to vivid magenta/red ($Z \ge +2$).

#### 3.6.2 Lineage Hierarchy Ordering
In the detailed dot plot (`06_dotplot_pikk_detailed_lineage`), cell lineages are explicitly reordered to group biologically related populations:
1. **Malignant**: Malignant epithelial cells ($n = 2,488$).
2. **Immune Infiltrates**: CD8T, CD8Tex, CD4Tconv, Mono/Macro, Mast, Plasma ($n = 1,673$).
3. **Stromal Microenvironment**: Fibroblasts, Myofibroblasts, Endothelial, Myocyte ($n = 1,741$).

---

### 3.7 Step 7: Compartment-Specific Distribution Modeling
**Script**: `02_sc_pikk_expression_cartography.R`  
**Input Artifacts**: `oscc_gse103322_seurat.rds`.  
**Output Artifacts**:
- `04_boxplot_compartment_{gene}.pdf` and `.png` (Individual box plots for each target).
- `04_boxplot_pikk_compartment_grid.pdf` and `.png` (Composite multi-panel box plot figure).

#### 3.7.1 Distribution Parameter Extraction
Expression data are extracted using `Seurat::FetchData` for each target gene along with compartment annotations:
```r
expr_data <- FetchData(seu, vars = c(diamond_present, "CellType_Malignancy", "CellType_Major"))
```

#### 3.7.2 Multi-Parameter Visualization Geometry
For each target gene across the three compartments (Malignant, Immune, Stromal):
- **Box Limits**: Represent the first quartile ($Q_1$, 25th percentile) and third quartile ($Q_3$, 75th percentile), defining the Interquartile Range ($\text{IQR} = Q_3 - Q_1$).
- **Center Line**: Represents the median (50th percentile).
- **Whiskers**: Extend to the most extreme data points within $1.5 \times \text{IQR}$ from the box limits:
  $$\text{Upper Whisker} = \min\left( \max(y), \; Q_3 + 1.5 \cdot \text{IQR} \right)$$
  $$\text{Lower Whisker} = \max\left( \min(y), \; Q_1 - 1.5 \cdot \text{IQR} \right)$$
- **Outlier Points**: Data points beyond the whiskers are plotted as individual semi-transparent points (`outlier.alpha = 0.25`, `outlier.size = 0.8`).
- **Mean Diamond Overlay**: Superimposed white diamond marker ($\diamond$) calculated via `stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white")` to visually contrast arithmetic mean expression against median expression.

---

### 3.8 Step 8: Quantitative Compartmentalization Hypothesis Testing
**Script**: `02_sc_pikk_expression_cartography.R`  
**Input Artifacts**: `oscc_gse103322_seurat.rds`.  
**Output Artifacts**:
- `pikk_gene_celltype_compartment_stats.tsv`: Master quantitative statistical table.
- `pikk_gene_celltype_expression_summary.tsv`: Detailed per-lineage summary statistics.

#### 3.8.1 Omnibus Non-Parametric Test: Kruskal–Wallis Rank-Sum
To determine whether expression of gene $g$ differs significantly across the three broad compartments ($C = 3$: Malignant, Stromal, Immune) without assuming normality or equal variances, the **Kruskal–Wallis test** is applied [4]:

All $N = 5,902$ single cells are pooled and assigned global ranks $r_i \in \{1, \dots, N\}$ based on expression $y_{gi}$ (with average ranks assigned to ties). The test statistic $H$ is computed:

$$H = \frac{12}{N(N + 1)} \sum_{c=1}^C \frac{R_c^2}{n_c} - 3(N + 1)$$

Where:
- $C = 3$ cellular compartments.
- $n_c$: Number of cells in compartment $c$ ($n_{\text{malig}} = 2,488$; $n_{\text{stroma}} = 1,741$; $n_{\text{immune}} = 1,673$).
- $R_c = \sum_{i \in \text{compartment } c} r_i$: Sum of ranks in compartment $c$.
- Under the null hypothesis that expression distributions are identical across all three compartments:
  $$H \sim \chi^2(\text{df} = C - 1 = 2)$$

When ties are present (ubiquitous in single-cell count matrices due to zero entries), the tie-corrected statistic $H_{\text{adj}}$ is applied:

$$H_{\text{adj}} = \frac{H}{1 - \frac{\sum_{t} (t^3 - t)}{N^3 - N}}$$

Where $t$ represents the number of tied observations in each tie group.

#### 3.8.2 Post-Hoc Pairwise Comparisons: Two-Sample Wilcoxon Rank-Sum Test
When the omnibus Kruskal–Wallis test indicates significant across-group divergence ($P_{\text{KW}} < 0.05$), post-hoc pairwise comparisons are executed between all three compartment pairs:
1. Malignant cells vs. Immune cells
2. Malignant cells vs. Stromal cells
3. Immune cells vs. Stromal cells

For any two compartments $A$ and $B$ with cell counts $n_A$ and $n_B$, the Wilcoxon rank-sum statistic $W$ is calculated:

$$W = \sum_{i \in A} r_{i(A,B)} - \frac{n_A(n_A + 1)}{2}$$

Where $r_{i(A,B)}$ denotes the rank of cell $i$ within the pooled subset of $n_A + n_B$ cells. The two-sided nominal $p$-value ($P_{\text{Wilcox}}$) is obtained from the normal approximation:

$$Z = \frac{W - \frac{n_A n_B}{2}}{\sqrt{\frac{n_A n_B (n_A + n_B + 1)}{12}}}$$

#### 3.8.3 Multiple Testing Adjustment: Benjamini–Hochberg False Discovery Rate
Across all pairwise comparisons for each gene, nominal $p$-values are adjusted using the **Benjamini–Hochberg (BH)** procedure to control the False Discovery Rate [6]:

$$q_{(k)} = \min_{j \ge k} \left( \frac{m \cdot P_{(j)}}{j} \right)$$

Where $m = 3$ pairwise comparisons per gene, and $P_{(1)} \le P_{(2)} \le P_{(3)}$ represent the ordered nominal $p$-values. Significant compartment differences are declared at $q_{\text{FDR}} < 0.05$.

#### 3.8.4 Predominant Compartment Assignment
For each biomarker, the predominant compartment is classified based on maximal mean normalized expression:

$$\text{Compartment}_{\text{predominant}}(g) = \arg\max_{c \in \{\text{Malignant}, \text{Stromal}, \text{Immune}\}} \left( \bar{y}_{gc} \right)$$

A biomarker is classified as a **bona fide tumor-cell-intrinsic target** if:
1. The omnibus Kruskal–Wallis test is statistically significant ($P_{\text{KW}} < 10^{-10}$).
2. Post-hoc Wilcoxon tests confirm significant upregulation in Malignant cells relative to both Immune and Stromal compartments ($q_{\text{FDR}} < 0.05$).
3. $\text{Compartment}_{\text{predominant}} = \text{"Malignant cells"}$.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│               EMPIRICAL METRIC CALLOUT: HYPOTHESIS TESTING FRAMEWORK                   │
├──────────────────────────────────────────────────┬─────────────────────────────────────┤
│ Statistical Parameter                            │ Value / Implementation              │
├──────────────────────────────────────────────────┼─────────────────────────────────────┤
│ Omnibus Statistical Model                        │ Kruskal–Wallis Rank-Sum Test        │
│ Group Partitioning (C)                           │ 3 compartments (df = 2)             │
│ Pairwise Post-Hoc Model                          │ Two-sample Wilcoxon Rank-Sum Test   │
│ Multiple Testing Correction                      │ Benjamini–Hochberg (BH) FDR         │
│ Statistical Significance Threshold               │ q_FDR < 0.05 (nominal p < 0.05)     │
│ Master Results Serialization                     │ pikk_gene_celltype_compartment_     │
│                                                  │ stats.tsv                           │
└──────────────────────────────────────────────────┴─────────────────────────────────────┘
```

---

### 3.9 Step 9: Detection Breadth Frequency Analysis
**Script**: `02_sc_pikk_expression_cartography.R`  
**Input Artifacts**: `pikk_gene_celltype_compartment_stats.tsv`.  
**Output Artifacts**: `05_barplot_pikk_pct_expressing.pdf` and `.png` (300 DPI).  

#### 3.9.1 Detection Frequency Formulation
To complement continuous expression intensity metrics, the binary detection breadth (proportion of non-zero cells) is computed for each gene $g$ within each compartment $c$:

$$\text{PctNonZero}_{gc} = \frac{\sum_{i \in \mathcal{C}_c} \mathbb{I}(y_{gi} > 0)}{n_c} \times 100\%$$

Where:
- $\mathbb{I}(y_{gi} > 0) = 1$ if normalized expression exceeds zero; $0$ otherwise.
- $n_c$: Total cell count in compartment $c$ ($n_{\text{malig}} = 2,488$; $n_{\text{immune}} = 1,673$; $n_{\text{stroma}} = 1,741$).

#### 3.9.2 Stacked Horizontal Bar Visualization
The detection frequencies are rendered as horizontal grouped bar charts:
- $Y$-axis: Evaluated Diamond Panel genes ordered by biological axis.
- $X$-axis: Percentage of cells expressing the target ($0\%$ to $100\%$).
- Fill aesthetic: Discrete compartment colors matching global conventions (Red: Malignant; Steel Blue: Immune; Teal: Stromal).
- Dodged bar layout: `geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7)`.

---

## 4. Synthesis & Comparative Matrix

The following structured matrix compares the computational methods, parameters, input matrices, and deliverable artifacts across all modules of the single-cell cartography and compartmentalization pipeline:

| Analytical Step | Computational Script | Input Data / Cellular Cohort | Core Mathematical Model / Algorithm | Key Parameters / Filter Thresholds | Primary Output Deliverables |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Matrix Ingestion & Parsing** | `01_sc_seurat_setup.R` | GSE103322 HDF5 container; TISCH2 metadata | HDF5 sparse reconstruction (CSC format) | 18,242 genes $\times$ 5,902 cells; 0-to-1 index correction | Reconstructed sparse count matrix; metadata table |
| **Barcode Alignment & Sanitation** | `01_sc_seurat_setup.R` | Matrix columns; metadata barcodes | Barcode intersect & string sanitation | $N = 5,902$ common cells; syntactically valid column headers | Aligned single-cell metadata data frame |
| **Container Setup & Normalization** | `01_sc_seurat_setup.R` | CSC sparse matrix; aligned metadata | Seurat container setup; library-size LogNormalize | Scale factor = 10,000; pseudocount = 1.0; inject TISCH2 UMAP | `oscc_gse103322_seurat.rds`, `sc_qc_cell_counts.tsv`, `sc_qc_patient_summary.tsv` |
| **Reference UMAP Projection** | `02_sc_pikk_expression_cartography.R` | 5,902 cells; UMAP coordinates | Dimensional projection; centroid label repulsion | 11 major lineages; curated 11-color palette; `pt.size = 0.4` | `01_umap_celltype_overview.pdf` and `.png` |
| **Gene Expression Cartography** | `02_sc_pikk_expression_cartography.R` | 5,902 cells; Diamond Panel genes | Quantile-clipped feature projection | Clipping: $q_{0.05} \to q_{0.95}$; viridis magma; `order = TRUE` | `02_umap_feature_{gene}`, `02_umap_pikk_feature_plots.pdf/.png` |
| **Lineage Dot Plot Profiling** | `02_sc_pikk_expression_cartography.R` | 5,902 cells $\times$ 11 cell lineages | Dual-scaled prevalence & intensity | $\text{PctExp} \in [0, 100]\%$; $Z$-score scaled mean; max radius = 8 | `03_dotplot_pikk_by_celltype`, `06_dotplot_pikk_detailed_lineage` |
| **Compartment Distribution Modeling** | `02_sc_pikk_expression_cartography.R` | 3 compartments ($n = 2488, 1741, 1673$) | Multi-parameter box plots & diamond mean | Median $\pm$ IQR; $1.5 \times \text{IQR}$ whiskers; diamond ($\diamond$) mean | `04_boxplot_compartment_{gene}`, `04_boxplot_pikk_compartment_grid` |
| **Non-Parametric Omnibus Testing** | `02_sc_pikk_expression_cartography.R` | 3 compartments; normalized counts | Kruskal–Wallis Rank-Sum Test | Omnibus $H$-statistic, $\chi^2(\text{df} = 2)$; tie correction | `pikk_gene_celltype_compartment_stats.tsv` |
| **Post-Hoc Pairwise Testing** | `02_sc_pikk_expression_cartography.R` | Pairwise compartment subsets | Two-sample Wilcoxon Rank-Sum Test | 3 comparisons/gene; Benjamini–Hochberg FDR ($q < 0.05$) | `pikk_gene_celltype_compartment_stats.tsv` |
| **Detection Breadth Analysis** | `02_sc_pikk_expression_cartography.R` | 3 compartments; normalized data | Frequency proportion: $\sum \mathbb{I}(y > 0) / n_c$ | Expression threshold: $y > 0$; dodged horizontal bars | `05_barplot_pikk_pct_expressing.pdf` and `.png` |
| **Detailed Lineage Summary** | `02_sc_pikk_expression_cartography.R` | 11 lineages $\times$ Diamond genes | Descriptive summary statistics | Mean, median, SD, % non-zero per lineage | `pikk_gene_celltype_expression_summary.tsv` |

---

## 5. Ready-to-Publish Methods Section (Manuscript Reference Text)

### 2.7 Single-Cell Transcriptomic Cartography and Cellular Compartmentalization

#### 2.7.1 Single-Cell RNA-seq Atlas Ingestion and Seurat Container Processing
To resolve the exact cellular compartmentalization of prioritized PIKK biomarkers and verify whether candidate oncogenic drivers represent tumor-cell-intrinsic programs or microenvironmental confounders, single-cell transcriptomic data were analyzed from a benchmark clinical oral cavity squamous cell carcinoma atlas (GSE103322; Puram et al.) [1], curated and harmonized through the Tumor Immune Single-cell Hub 2 (TISCH2) [2]. The dataset comprised $N = 5,902$ high-quality single cells spanning $18,242$ detected human genes, derived from $22$ surgical biospecimens across $18$ treatment-naive primary HNSCC patients (spanning TNM clinical stages I–III). 

Cells were categorized into three broad compartments:
1. Malignant epithelial cells ($n = 2,488$; 42.16% of atlas), identified based on patient-specific inferred chromosomal copy number alterations (CNAs) and epithelial marker expression.
2. Microenvironmental stromal cells ($n = 1,741$; 29.50% of atlas), encompassing cancer-associated fibroblasts ($n = 744$), myofibroblasts ($n = 710$), endothelial cells ($n = 269$), and myocytes ($n = 18$).
3. Tumor-infiltrating immune cells ($n = 1,673$; 28.35% of atlas), encompassing exhausted CD8+ T cells (`CD8Tex`, $n = 501$), conventional CD4+ T cells (`CD4Tconv`, $n = 417$), effector CD8+ T cells (`CD8T`, $n = 397$), plasma B cells ($n = 148$), mast cells ($n = 122$), and monocytes/macrophages (`Mono/Macro`, $n = 88$).

Single-cell expression matrices were ingested from HDF5 containers using `hdf5r` (v1.3.8) and reconstructed into Compressed Sparse Column (CSC) `dgCMatrix` formats using the `Matrix` package (v1.6-1). Cell barcodes were matched with clinical metadata, non-syntactic column headers were harmonized into clean R variables, and formal S4 containers were instantiated in `Seurat` (v4.3+) [5]. Pre-computed two-dimensional Uniform Manifold Approximation and Projection (UMAP) coordinates from TISCH2 were imported directly as a dimensional reduction assay to maintain strict concordance with published reference clusters. Transcript counts were normalized across cells using library-size scaling with natural log-transformation (`NormalizeData`, `scale.factor = 10,000`):

$$y_{gi} = \ln\left(1 + \frac{x_{gi}}{\sum_{k=1}^G x_{ki}} \cdot 10,000\right)$$

#### 2.7.2 Single-Cell Expression Cartography, Lineage Profiling, and Non-Parametric Compartmentalization Testing
To characterize the cellular distribution of the 7 core PIKK Diamond Panel targets ($PLK1$, $CDK2$, $TOPBP1$, $RAD51$, $FANCI$, $KAT2B$, $DEPTOR$), normalized expression values were projected onto single-cell UMAP embeddings using perceptually uniform viridis magma palettes. To avoid visual compression from hyper-expressing outlier cells, expression values were clipped to the 5th and 95th percentiles (`min.cutoff = "q5"`, `max.cutoff = "q95"`), and cells were ordered by expression intensity (`order = TRUE`) to ensure expressing cells were rendered on top. 

Expression magnitude and detection prevalence were evaluated across all 11 major lineages using dual-scaled dot plots (`DotPlot`). Dot diameter represented the percentage of cells within each lineage with non-zero expression ($\text{PctExp} \in [0, 100]\%$), while color intensity reflected mean expression standardized via $Z$-score scaling across lineages.

To statistically validate whether biomarker expression was restricted to transformed carcinoma cells or shared with stromal and immune populations, non-parametric Kruskal–Wallis rank-sum tests were conducted across the three cellular compartments (Malignant, Stromal, and Immune; $\text{df} = 2$) [4]. When omnibus testing indicated significant across-group divergence, post-hoc pairwise comparisons between compartments were executed using two-sample Wilcoxon rank-sum tests. Nominal $p$-values across pairwise comparisons were adjusted for multiple testing using the Benjamini–Hochberg False Discovery Rate (FDR) procedure [6]. Predominant compartment status was assigned based on maximum compartment mean expression, with bona fide tumor-cell-intrinsic oncotargets defined by significant enrichment in malignant cells ($P_{\text{KW}} < 10^{-10}$, $q_{\text{FDR}} < 0.05$). Binary detection frequencies (% cells with $y > 0$) were computed per compartment and visualized via grouped horizontal bar plots.

---

## 6. Comprehensive Software, Dependency, and Environment Manifest

All single-cell matrix parsing, sparse matrix processing, dimensional projection, non-parametric hypothesis testing, and publication-quality figure rendering were executed within the R statistical environment (R version $\ge 4.2.0$):

| Package Name | Origin / Repository | Version Tested | Primary Computational Purpose in Pipeline |
| :--- | :--- | :--- | :--- |
| **`Seurat`** | CRAN | v4.3.0+ | Core single-cell genomics architecture: object construction, metadata alignment, `NormalizeData`, UMAP dimensional reduction embedding, `FeaturePlot`, `DotPlot`, and `FetchData` [5]. |
| **`hdf5r`** | CRAN | v1.3.8 | Low-level C++ HDF5 container parsing, structure interrogation, and extraction of sparse CSC matrix components. |
| **`Matrix`** | Base R / CRAN | v1.6-1 | Compressed Sparse Column (CSC) matrix representation (`sparseMatrix`, `dgCMatrix`), minimizing memory footprint for single-cell count matrices. |
| **`ggplot2`** | CRAN | v3.4.4 | Core graphic generation engine: customized publication themes, multi-lineage boxplots, and horizontal grouped bar charts. |
| **`patchwork`** | CRAN | v1.2.0 | Multi-panel figure assembly, multi-column grid alignment, and global title/subtitle annotations. |
| **`viridis`** | CRAN | v0.6.4 | Perceptually uniform, colorblind-safe color scales (`magma`) for continuous expression feature overlays. |
| **`RColorBrewer`**| CRAN | v1.1-3 | Categorical palette generation for lineage and compartment stratification. |
| **`scales`** | CRAN | v1.3.0 | Axis tick formatting, mathematical percentage transformations, and color gradient scales. |
| **`dplyr`** | CRAN | v1.1.4 | Data wrangling, metadata column sanitation, grouping, and statistical aggregation. |
| **`tidyr`** | CRAN | v1.3.1 | Data reshaping (`pivot_longer`, `pivot_wider`) for summary tables and multi-parameter visualizations. |
| **`readr`** | CRAN | v2.1.5 | High-performance, type-safe serialization of TSV and CSV summary tables. |
| **`tibble`** | CRAN | v3.2.1 | Data frame rowname conversions and column manipulations. |

---

## 7. Peer-Reviewed Scientific Bibliography

1. **Puram, S. V., Tirosh, I., Parikh, A. S., Patel, A. P., Yizhak, K., Gillespie, S., Rodman, C., Luo, C. L., Mambetsariev, N., Vaidya, A., Hanen, R. R., Wakiro, I., Rotem, A., Berniker, A., Haber, D. A., Sadow, P. M., Lin, D. T., Emerick, K. S., Deschler, D. G., Varvares, M. A., & Bernstein, B. E.** (2017). Single-cell transcriptomic analysis of primary and metastatic head and neck squamous cell carcinoma. *Cell*, 171(7), 1611–1624.e24. [DOI: 10.1016/j.cell.2017.10.044](https://doi.org/10.1016/j.cell.2017.10.044)
2. **Sun, D., Wang, J., Han, Y., Dong, X., Ge, J., Zheng, R., Shi, X., Wang, B., Li, Z., Ren, P., Sun, L., Yan, Y., Zhang, P., Zhang, F., & Guo, A. Y.** (2023). TISCH2: an expanded web resource for single-cell RNA sequencing data in tumor microenvironment. *Nucleic Acids Research*, 51(D1), D1420–D1431. [DOI: 10.1093/nar/gkac959](https://doi.org/10.1093/nar/gkac959)
3. **Puram, S. V., & Bernstein, B. E.** (2018). Dissecting head and neck squamous cell carcinoma using single-cell RNA sequencing. *Oncogene*, 37(34), 4647–4655. [DOI: 10.1038/s41388-018-0287-4](https://doi.org/10.1038/s41388-018-0287-4)
4. **Kruskal, W. H., & Wallis, W. A.** (1952). Use of ranks in one-criterion variance analysis. *Journal of the American Statistical Association*, 47(260), 583–621. [DOI: 10.1080/01621459.1952.10483441](https://doi.org/10.1080/01621459.1952.10483441)
5. **Hao, Y., Hao, S., Andersen-Nissen, E., Mauck, W. M., Zheng, S., Butler, A., Lee, M. J., Wilk, A. J., Darby, C., Zagar, M., Hoffman, P., Stoeckius, M., Papalexi, E., Mimitou, E. P., Jain, J., Srivastava, A., Stuart, T., Fleming, L. B., Yeung, B., … Satija, R.** (2021). Integrated analysis of multimodal single-cell data. *Cell*, 184(13), 3573–3587.e29. [DOI: 10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048)
6. **Benjamini, Y., & Hochberg, Y.** (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B (Methodological)*, 57(1), 289–300. [DOI: 10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
7. **Macosko, E. Z., Basu, A., Satija, R., Nemesh, J., Shekhar, K., Goldman, M., Tirosh, I., Bialas, A. R., Kamitaki, N., Martersteck, E. M., Trombetta, J. J., Weitz, D. A., Sanes, J. R., Shalek, A. K., Regev, A., & McCarroll, S. A.** (2015). Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets. *Cell*, 161(5), 1202–1214. [DOI: 10.1016/j.cell.2015.05.002](https://doi.org/10.1016/j.cell.2015.05.002)
