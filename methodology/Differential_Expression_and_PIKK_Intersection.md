# Differential Gene Expression Analysis and PIKK Intersection (Section 2.2) — Comprehensive Methodology Guide

**Project**: PIKK Pathway-Associated Biomarker Discovery and Oncogenic Characterization in Oral Squamous Cell Carcinoma (OSCC)  
**Target Manuscript Reference**: Materials and Methods / Bioinformatic Processing Pipeline  
**Primary Dataset**: The Cancer Genome Atlas Head and Neck Squamous Cell Carcinoma (TCGA-HNSC), Curated Oral Cavity Cohort ($n = 240$ primary samples; $n = 224$ Primary Solid Tumor, $n = 16$ Solid Tissue Normal)  
**Feature Space**: 16,905 Robustly Expressed Protein-Coding Genes (pre-filtered via `filterByExpr`)  
**Target Biological Universe**: Curated PIKK Gene Universe ($n = 204$ gene-pathway functional annotations spanning $ATM$, $ATR$, $PRKDC$, $MTOR$, $TRRAP$, and $SMG1$ core kinases and their associated interactomes)  
**Scripts Evaluated**:
1. `master_edgeR_OSCC.R` — *edgeR Quasi-Likelihood Differential Expression Pipeline, TMM Normalization, Dispersion Estimation, and Exploratory Diagnostics*
2. `intersect_PIKK.R` — *Targeted Curation, Harmonization, and Two-Tier Statistical Intersection with the PIKK Pathway Universe*

---

## 1. Executive Summary & End-to-End Workflow

Differential Gene Expression Analysis (DGEA) constitutes the primary inferential engine of cancer functional genomics, allowing researchers to discern transcriptional remodeling between malignant and healthy tissues. In bulk RNA-sequencing (RNA-seq) datasets characterized by high sample dimensionality and biological heterogeneity, robust statistical inference demands modeling frameworks that explicitly capture count-based sampling error, library composition disparities, and replicate variance overdispersion.

This document outlines the rigorous, publication-grade methodology governing the **differential expression modeling and PIKK (Phosphatidylinositol 3-kinase-related kinase) family pathway intersection** for the OSCC biomarker discovery project. 

The analytical framework combines the state-of-the-art **edgeR Quasi-Likelihood (QL) F-test pipeline** (`master_edgeR_OSCC.R`) with a **targeted molecular intersection algorithm** (`intersect_PIKK.R`). This pipeline bridges the upstream audit and filtering procedures established in Section 2.1 with downstream functional enrichment, protein-protein interaction (PPI) topology, and clinical survival modeling.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                                 UPSTREAM INPUT ARTIFACTS                               │
│  • Matched Count Matrix: OSCC_counts_matched.tsv (n = 240 samples; 16,905 genes)      │
│  • Curated Metadata: OSCC_metadata_matched.tsv (Tumor = 224, Normal = 16)              │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 1: DGEList INITIALIZATION & TMM NORMALIZATION (master_edgeR_OSCC.R)             │
│  • Factor level anchoring: reference group set to "Normal"                             │
│  • Trimmed Mean of M-values (TMM) library size scaling via calcNormFactors()           │
│  • Calculation of effective library sizes (N_g* = N_g × f_g)                           │
│  • Generation of variance-stabilized log2(CPM) matrix with offset (c = 2)              │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 2: EXPLORATORY MULTI-DIMENSIONAL QUALITY CONTROL (master_edgeR_OSCC.R)         │
│  • Multidimensional Scaling (MDS) on leading logFC dimensions                          │
│  • Principal Component Analysis (PCA) & Scree eigen-decomposition                      │
│  • Pairwise sample distance evaluation using Pearson correlation dissimilarity         │
│  • Assessment of batch covariates & design matrix rank-deficiency verification         │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 3: NEGATIVE BINOMIAL DISPERSION & GLM FITTING (master_edgeR_OSCC.R)             │
│  • Empirical Bayes dispersion estimation: Common, Trended, and Tagwise (estimateDisp)  │
│  • Robust empirical Bayes shrinkage of gene-wise dispersions toward the trend curve    │
│  • Log-linear Negative Binomial GLM parameter estimation                               │
│  • Quasi-Likelihood dispersion fitting (glmQLFit) with squeezed error variance         │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 4: HYPOTHESIS TESTING & STATISTICAL STRATIFICATION (master_edgeR_OSCC.R)        │
│  • Quasi-Likelihood F-test (glmQLFTest) on Tumor vs Normal contrast vector             │
│  • Multiple testing correction via Benjamini-Hochberg False Discovery Rate (FDR)       │
│  • Dual-threshold significance filtering: FDR ≤ 0.05 and |log2FC| ≥ 1.0 (2-fold change)│
│  • Publication visualization generation: Volcano, MA, and Ward.D2 Heatmaps             │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 5: PIKK PATHWAY INTERSECTION & TARGET IDENTIFICATION (intersect_PIKK.R)        │
│  • Ingestion of curated PIKK reference universe (204 gene-pathway functional records)  │
│  • Case-insensitive HGNC gene symbol harmonization and set-theoretic inner join        │
│  • Two-tier statistical thresholding:                                                  │
│      - Strict Canonical Threshold: FDR ≤ 0.05 and |log2FC| ≥ 1.0                       │
│      - Extended Complex Threshold: FDR ≤ 0.05 and |log2FC| ≥ 0.5                       │
│  • Stratification by PIKK Group (ATM, ATR, PRKDC, MTOR, TRRAP, SMG1) & Gene Role       │
│  • Heatmap clustering using McQuitty / Average linkage on row-standardized Z-scores    │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. Biological Context and Study Rationale

### 2.1 Molecular Pathogenesis of OSCC and the PIKK Kinase Signaling Axis
Oral Squamous Cell Carcinoma (OSCC) represents the predominant malignancy of the head and neck, arising from the stratified squamous epithelium of the oral cavity proper (oral tongue, floor of mouth, alveolar ridge, hard palate, and buccal mucosa). Unlike HPV-positive oropharyngeal carcinomas that are driven by viral oncoproteins E6 and E7, oral cavity carcinomas are overwhelmingly HPV-negative and driven by prolonged exposure to exogenous carcinogens (tobacco, betel nut, alcohol).

At the molecular level, OSCC is defined by:
1. **Pervasive Inactivation of Cell Cycle Checkpoints**: High-frequency loss-of-function mutations or deletions in $TP53$ ($>80\%$) and $CDKN2A$ ($p16^{\text{INK4a}}$), disabling the classical $G_1/S$ restriction checkpoint.
2. **Replication Stress and Genomic Instability**: Oncogene-induced replication stress, unscheduled firing of replication origins, and persistent structural DNA damage.
3. **Heightened Dependence on the PIKK Kinase Superfamily**:
   In the absence of functional $p53$-mediated apoptosis and $G_1$ arrest, oral carcinoma cells become strictly reliant on the **Phosphatidylinositol 3-kinase-related kinase (PIKK)** family of atypical serine/threonine protein kinases to mitigate replication collapse, manage DNA double-strand breaks (DSBs), and coordinate metabolic rewiring.

The PIKK family comprises six giant, structurally conserved catalytic kinases that govern essential stress-response and maintenance programs:
- **ATM (Ataxia Telangiectasia Mutated)**: The primary transducer of double-strand break (DSB) signaling. Activated by the MRN complex ($MRE11$-$RAD50$-$NBN$), ATM autophosphorylates and phosphorylates downstream targets ($CHEK2$, $H2AX$, $BRCA1$) to enforce $G_2/M$ checkpoint arrest and initiate homologous recombination repair (HRR).
- **ATR (ATM and Rad3-Related)**: The master coordinator of the replication stress response. Recruited to persistent single-stranded DNA (ssDNA) coated by Replication Protein A (RPA) via its obligatory partner $ATRIP$, ATR activates $CHEK1$, stabilizing stalled replication forks and preventing catastrophic fork collapse.
- **PRKDC (DNA-Dependent Protein Kinase Catalytic Subunit, DNA-PKcs)**: The central catalytic component of non-homologous end joining (NHEJ). Recruited to DSB ends by the $XRCC5/XRCC6$ (Ku70/Ku80) heterodimer, PRKDC governs rapid, template-independent DNA ligation.
- **MTOR (Mechanistic Target of Rapamycin)**: Operates within two structurally and functionally distinct multi-protein complexes ($mTORC1$ and $mTORC2$), integrating nutrient availability, growth factor signaling, and energy status to regulate ribosome biogenesis, protein translation, autophagy, and anabolic metabolism.
- **TRRAP (Transformation/Transcription Domain-Associated Protein)**: A pseudokinase member lacking intrinsic phosphotransferase activity that functions as an indispensable scaffolding subunit in multi-protein histone acetyltransferase (HAT) complexes (e.g., SAGA, TIP60/NuA4), coupling chromatin remodeling to DNA repair and oncogenic transcription.
- **SMG1 (Suppressor with Morphogenetic Effect on Genitalia 1)**: The effector kinase of the nonsense-mediated mRNA decay (NMD) pathway. SMG1 phosphorylates $UPF1$ upon recognizing premature termination codons (PTCs), surveilling transcript fidelity and managing cellular responses to oxidative and genotoxic stress.

Systematically dissecting which PIKK family members, downstream effectors, or upstream regulatory complexes undergo transcriptional dysregulation is essential to identify vulnerabilities in the OSCC DNA damage response (DDR) machinery.

### 2.2 Methodological Rationale: Why the edgeR Quasi-Likelihood (QL) Pipeline?
While various parametric and non-parametric tools exist for differential RNA-seq analysis, the **edgeR Quasi-Likelihood (QL) framework** (Chen, Lun, & Smyth, 2016 [1]; Lund et al., 2012 [2]) was selected over conventional alternatives (such as the standard negative binomial Likelihood Ratio Test [LRT] or Wald test) due to several statistical properties:

1. **Unbalanced Cohort Architecture**:
   Our clinical cohort consists of 224 primary oral tumors and 16 matched normal mucosal controls. Standard asymptotic tests (e.g., asymptotic Wald tests or standard chi-square likelihood ratio tests) assume large, balanced sample sizes across all conditions. When evaluating 16 normal controls against 224 tumors, asymptotic approximations underestimate the uncertainty in dispersion estimation, leading to inflated type-I error rates (false positives).
2. **Gene-Specific Uncertainty via the Quasi-Likelihood $F$-Test**:
   Standard negative binomial GLMs assume that the true dispersion parameter for each gene is known without error once estimated from the trend curve. In contrast, the QL framework introduces an additional quasi-likelihood dispersion parameter ($s_i^2$) for each gene. The resulting test statistic follows an **$F$-distribution** rather than a $\chi^2$-distribution:
   $$F_i \sim F(k, d_{\text{post}})$$
   The denominator degrees of freedom ($d_{\text{post}}$) reflect the finite sample size and uncertainty in variance estimation, ensuring exact, conservative, and robust error rate control across the genome (Lund et al., 2012 [2]).
3. **Robust Empirical Bayes Dispersion Shrinkage**:
   Biological replicates in human patient cohorts exhibit variable, gene-specific dispersion. The robust empirical Bayes procedure in edgeR (`robust = TRUE`) shrinks gene-wise dispersions toward a global mean-dispersion trend while simultaneously identifying and accommodating outlier genes (e.g., hypervariable immune or mucosal differentiation genes) without allowing them to distort the global trend curve (Phipson et al., 2016 [3]).

---

## 3. Step-by-Step Computational Methodology

### 3.1 Step 1: Data Ingestion, Factor Releveling, and DGEList Construction
**Script**: `master_edgeR_OSCC.R`  
**Input Artifacts**:
- `counts_file`: Matched raw count matrix (`OSCC_counts_matched.tsv`) containing 16,905 protein-coding genes across 240 samples.
- `meta_file`: Matched clinical metadata table (`OSCC_metadata_matched.tsv`) containing sample identifiers (`File_ID`) and biological condition (`Tissue_Type`).
**Output Artifacts**:
- `matched_metadata_used.tsv`: Serialized record of the exact sample annotations passed to the modeling engine.
- Initialized `DGEList` object.

#### 3.1.1 Algorithmic Procedure
1. **Sample Intersection & Alignment Verification**:
   The script extracts the sample identifiers from the counts matrix headers (excluding the first column, `gene_name`) and metadata records:
   $$\mathcal{S}_{\text{counts}} = \text{colnames}(K)[-1], \quad \mathcal{S}_{\text{meta}} = \text{meta}\$\text{File\_ID}$$
   $$\mathcal{S}_{\text{shared}} = \mathcal{S}_{\text{counts}} \cap \mathcal{S}_{\text{meta}}$$
   A strict assertion requires $|\mathcal{S}_{\text{shared}}| \ge 4$. The columns of the count matrix and rows of the metadata table are synchronized:
   $$\text{stopifnot}(\text{identical}(\text{colnames}(K_{\text{matched}})[-1], \text{meta}_{\text{matched}}\$\text{File\_ID}))$$
2. **Matrix Casting & Duplicate Symbol Disambiguation**:
   Counts are coerced to a numeric matrix storage mode. Gene symbols are sanitized using `make.unique()` to resolve any duplicate HGNC annotations:
   $$\text{storage.mode}(K) \leftarrow \text{"numeric"}$$
3. **Reference Group Anchoring**:
   In linear modeling, differential coefficients quantify change relative to an unperturbed baseline. The biological condition variable (`Tissue_Type`) is cast to a categorical factor with the non-malignant tissue designated as the explicit reference level:
   $$\text{meta}\$\text{Tissue\_Type} \leftarrow \text{relevel}(\text{as.factor}(\text{meta}\$\text{Tissue\_Type}), \text{ref} = \text{"Normal"})$$
   This guarantees that all downstream model contrasts quantify **Tumor relative to Normal** ($+\log_2\text{FC} \implies$ upregulated in tumor; $-\log_2\text{FC} \implies$ downregulated in tumor).
4. **`DGEList` Instantiation**:
   The count matrix and sample metadata are encapsulated within an edgeR `DGEList` container:
   $$\mathcal{D} = \text{DGEList}(\text{counts} = K, \text{genes} = \text{gene\_names}, \text{samples} = \text{meta})$$

---

### 3.2 Step 2: Library Size Normalization via Trimmed Mean of M-values (TMM)
**Script**: `master_edgeR_OSCC.R`  
**Input Artifacts**: `DGEList` container $\mathcal{D}$  
**Output Artifacts**: Updated `DGEList` containing sample-specific normalization factors ($f_j$); `QC_library_sizes.pdf/.png`.

#### 3.2.1 Biological and Statistical Rationale
A pervasive challenge in transcriptomic profiling of solid tumors is **composition bias**. Malignant transformation frequently induces global transcriptional amplification (e.g., driven by $MYC$ deregulation) or massive overexpression of tissue-specific transcripts (e.g., keratins, matrix metalloproteinases). 

If a small subset of genes is expressed at extraordinarily high levels in tumor tissue, they consume a disproportionate fraction of the total sequencing reads. Consequently, all other genes in that tumor library appear artificially downregulated when scaled purely by total read count (Library Size Scaling or standard CPM):
$$\text{Naive CPM}_{ij} = \frac{K_{ij}}{\sum_i K_{ij}} \times 10^6$$
To eliminate composition artifacts, the **Trimmed Mean of M-values (TMM)** method (Robinson & Oshlack, 2010 [4]) estimates scaling factors that force the expression ratio of unperturbed housekeeping genes to remain invariant across libraries.

#### 3.2.2 Mathematical Formulation of TMM
Let $Y_{gk}$ denote the observed read count for gene $g$ in sample $k$, and let $N_k = \sum_g Y_{gk}$ denote the unadjusted library size.

1. **Selection of the Reference Library**:
   A reference sample $r$ is chosen whose 75th percentile of expression is closest to the mean 75th percentile across all samples in the dataset. All other libraries $k$ are compared against reference $r$.
2. **Calculation of Log-Fold Changes ($M$) and Absolute Abundance ($A$)**:
   For each gene $g$, the log-ratio of expression ($M_{gk}$) and the average absolute log-concentration ($A_{gk}$) are computed:
   $$M_{gk} = \log_2 \left( \frac{Y_{gk} / N_k}{Y_{gr} / N_r} \right)$$
   $$A_{gk} = \frac{1}{2} \log_2 \left( \frac{Y_{gk}}{N_k} \cdot \frac{Y_{gr}}{N_r} \right) = \frac{1}{2} \left[ \log_2 \left( \frac{Y_{gk}}{N_k} \right) + \log_2 \left( \frac{Y_{gr}}{N_r} \right) \right]$$
3. **Double Asymmetric Trimming**:
   To prevent outliers from distorting the normalization factor:
   - The top and bottom $30\%$ of genes ranked by log-fold change ($M_{gk}$) are trimmed ($\text{trim}_M = 0.30$).
   - The top and bottom $5\%$ of genes ranked by absolute expression ($A_{gk}$) are trimmed ($\text{trim}_A = 0.05$).
   Let $G^*$ denote the retained subset of doubly-trimmed genes.
4. **Precision-Weighted Average of Retained Log-Ratios**:
   The log-normalization factor $\log_2(f_k)$ is computed as the weighted average of $M_{gk}$ values across $g \in G^*$:
   $$\log_2(f_k) = \frac{\sum_{g \in G^*} w_{gk} M_{gk}}{\sum_{g \in G^*} w_{gk}}$$
   The statistical weights $w_{gk}$ are calculated using the delta method to estimate the inverse asymptotic variance of $M_{gk}$:
   $$w_{gk} = \left[ \text{Var}(M_{gk}) \right]^{-1} = \left( \frac{N_k - Y_{gk}}{N_k Y_{gk}} + \frac{N_r - Y_{gr}}{N_r Y_{gr}} \right)^{-1}$$
5. **Effective Library Size Calculation**:
   The scaling factors are centered such that their product equals 1: $\prod_{k=1}^m f_k = 1$. The **effective library size** ($N_k^*$) for each sample is:
   $$N_k^* = N_k \times f_k$$
   All downstream model fitting and CPM transformations in edgeR incorporate $N_k^*$ rather than raw library size $N_k$.

---

### 3.3 Step 3: Exploratory Multi-Dimensional Quality Control & Batch Exploration
**Script**: `master_edgeR_OSCC.R`  
**Input Artifacts**: Normalized `DGEList` object  
**Output Artifacts**:
- `QC_MDS.pdf/.png`: Multidimensional scaling plot based on leading logFC.
- `QC_PCA_PC1_vs_PC2.pdf/.png`, `QC_PCA_PC1_vs_PC3.pdf/.png`: Principal Component Analysis score plots with $95\%$ confidence ellipses.
- `QC_PCA_scree_plot.pdf/.png`: Percentage of total variance explained by the top 10 principal components.
- `QC_sample_distance_heatmap.pdf/.png`: Sample-to-sample distance matrix based on Pearson correlation dissimilarity ($1 - r$).

#### 3.3.1 Multidimensional Scaling (MDS)
MDS produces an ordination where the distances between sample pairs correspond to the **root-mean-square average of the largest $\log_2$-fold changes** between them:
$$d_{\text{MDS}}(j, k) = \sqrt{\frac{1}{|\mathcal{G}_{\text{lead}}|} \sum_{g \in \mathcal{G}_{\text{lead}}} \left( \log_2 \text{CPM}_{gj} - \log_2 \text{CPM}_{gk} \right)^2}$$
Where $\mathcal{G}_{\text{lead}}$ is the set of top 500 genes with the largest pairwise log-fold changes between samples $j$ and $k$. This metric visualizes sample separation driven by biological differences while ignoring invariant genes.

#### 3.3.2 Variance-Stabilized $\log_2(\text{CPM})$ Transformation for PCA
Linear ordination via PCA requires variance-stabilized data where the variance of a gene is decoupled from its mean abundance. The script computes $\log_2(\text{CPM})$ with an offset parameter ($c = 2$):
$$y_{ij} = \log_2 \left( \frac{K_{ij} + c}{N_j^* + 2c} \times 10^6 \right) = \log_2 \left( K_{ij} + 2 \right) - \log_2 \left( N_j^* + 4 \right) + 6 \log_2(10)$$
Adding a prior count of $c = 2$ compresses the artificially inflated variance of near-zero counts, stabilizing low-expression features without distorting high-abundance genes.

#### 3.3.3 Principal Component Analysis (PCA) & Scree Analysis
PCA is performed via Singular Value Decomposition (SVD) on the centered and standardized expression matrix $\tilde{Y}$:
$$\tilde{Y} = U \Sigma V^T$$
Where:
- The columns of $V$ define the principal component eigenvectors (loadings).
- The matrix product $U \Sigma$ defines the sample projections (principal component scores: PC1, PC2, PC3).
- The percentage of total variance explained by principal component $l$ is:
  $$\text{Var\%}(PC_l) = \frac{\sigma_l^2}{\sum_{m=1}^p \sigma_m^2} \times 100$$
Parametric bivariate normal confidence ellipses ($95\%$ coverage) are computed for each biological group:
$$(x - \bar{x})^T S^{-1} (x - \bar{x}) \le \chi_2^2(0.95)$$

#### 3.3.4 Sample Distance Metric via Pearson Dissimilarity
To assess global sample relatedness, pairwise Pearson correlation coefficients ($r_{jk}$) are computed across the normalized expression space:
$$r_{jk} = \frac{\sum_{i=1}^G (y_{ij} - \bar{y}_{\cdot j})(y_{ik} - \bar{y}_{\cdot k})}{\sqrt{\sum_{i=1}^G (y_{ij} - \bar{y}_{\cdot j})^2 \sum_{i=1}^G (y_{ik} - \bar{y}_{\cdot k})^2}}$$
The sample distance matrix is defined by the metric dissimilarity:
$$D_{jk} = 1 - r_{jk}$$
Hierarchical clustering is executed using **complete linkage**, wherein the distance between two clusters $\mathcal{C}_A$ and $\mathcal{C}_B$ is the maximum distance between any single pair of constituent samples:
$$D_{\text{complete}}(\mathcal{C}_A, \mathcal{C}_B) = \max_{j \in \mathcal{C}_A, k \in \mathcal{C}_B} D_{jk}$$

---

### 3.4 Step 4: Negative Binomial Dispersion Modeling via Empirical Bayes
**Script**: `master_edgeR_OSCC.R`  
**Input Artifacts**: `DGEList` object and design matrix $\mathbf{X}$  
**Output Artifacts**: Estimated dispersion parameters; `edgeR_BCV_plot.pdf/.png`.

#### 3.4.1 Mathematical Foundations of Negative Binomial Overdispersion
In bulk RNA-seq, the count $Y_{gi}$ for gene $g$ in library $i$ represents a discrete random variable modeled under a Negative Binomial distribution:
$$Y_{gi} \sim \text{NB}(\mu_{gi}, \phi_g)$$
With expectation and variance defined as:
$$\mathbb{E}(Y_{gi}) = \mu_{gi} = N_i^* \exp(x_i^T \beta_g)$$
$$\text{Var}(Y_{gi}) = \mu_{gi} + \phi_g \mu_{gi}^2$$
Where:
- $\mu_{gi}$ is the expected count (technical Poisson shot-noise component).
- $\phi_g$ is the **gene-specific biological dispersion** parameter.
- The **Biological Coefficient of Variation (BCV)** is the square root of the dispersion parameter:
  $$\text{BCV}_g = \sqrt{\phi_g} = \frac{\sqrt{\text{Var}_{\text{bio}}(Y_{gi})}}{\mu_{gi}}$$
  $\text{BCV}_g$ reflects the coefficient of variation with which the true abundance of gene $g$ varies between biological replicates.

#### 3.4.2 Three-Tier Empirical Bayes Dispersion Estimation (`estimateDisp`)
Because individual gene sample sizes are insufficient to estimate $\phi_g$ with high precision, edgeR executes a three-tier empirical Bayes shrinkage protocol:
1. **Common Dispersion ($\phi_{\text{common}}$)**:
   A single global dispersion parameter is calculated by maximizing the conditional log-likelihood across all genes simultaneously:
   $$\hat{\phi}_{\text{common}} = \arg\max_\phi \sum_{g=1}^G \ell_{\text{cond}}(\phi; Y_g)$$
2. **Trended Dispersion ($\phi_{\text{trend}}(\mu)$)**:
   A continuous mean-dispersion trend is fitted across all genes using local likelihood smoothing (locfit), capturing the systematic decrease in biological dispersion as mean expression increases.
3. **Tagwise Dispersions with Robust Empirical Bayes Shrinkage**:
   Individual gene dispersions $\phi_g$ are estimated by maximizing a weighted combination of the individual conditional log-likelihood and the trended likelihood:
   $$\ell_{\text{post}}(\phi_g) = \ell_{\text{cond}}(\phi_g; Y_g) + d_0 \ell_{\text{trend}}(\phi_{\text{trend}}(\mu_g))$$
   The parameter $d_0$ represents the prior degrees of freedom, governing the extent to which individual gene dispersions are squeezed toward the global trend line. Setting `robust = TRUE` identifies outliers whose sample variance deviates excessively from the trend and moderates their influence, preventing hypervariable genes from biasing the global prior (Phipson et al., 2016 [3]).

---

### 3.5 Step 5: Generalized Linear Model (GLM) Fitting and Quasi-Likelihood F-Testing
**Script**: `master_edgeR_OSCC.R`  
**Input Artifacts**: Dispersed `DGEList` and design matrix $\mathbf{X}$  
**Output Artifacts**: Fitted `DGEGLM` object; `edgeR_QLDisp_plot.pdf/.png`; full table of differential expression statistics (`edgeR_all_results_tumor_vs_normal.tsv`).

#### 3.5.1 The Quasi-Likelihood Negative Binomial Framework (`glmQLFit`)
In standard negative binomial GLM testing, the dispersion parameter $\phi_g$ is assumed to be known without error. However, biological variation fluctuates from gene to gene beyond what can be captured by $\phi_g$ alone.

The Quasi-Likelihood (QL) framework (Lund et al., 2012 [2]) introduces a gene-specific quasi-likelihood dispersion parameter ($\sigma_g^2$), such that:
$$\text{Var}(Y_{gi}) = \sigma_g^2 \left( \mu_{gi} + \phi_g \mu_{gi}^2 \right)$$
1. **Raw QL Dispersion Estimation**:
   For each gene $g$, the raw quasi-likelihood dispersion estimator $s_g^2$ is derived from the deviance $D_g$ of the fitted GLM:
   $$s_g^2 = \frac{D_g}{m - p}$$
   Where:
   - $D_g = 2 \sum_i \left[ y_{gi} \ln \left( \frac{y_{gi}}{\hat{\mu}_{gi}} \right) - \left( y_{gi} + \frac{1}{\phi_g} \right) \ln \left( \frac{y_{gi} + 1/\phi_g}{\hat{\mu}_{gi} + 1/\phi_g} \right) \right]$ is the negative binomial deviance.
   - $m$ is the total number of samples ($m = 240$).
   - $p$ is the rank of the design matrix ($p = 2$ for the unadjusted model `~ Tissue_Type`).
   - $m - p$ is the residual degrees of freedom ($240 - 2 = 238$).
2. **Empirical Bayes Squeezing of QL Dispersions**:
   Individual estimates $s_g^2$ are squeezed toward a global trended prior $s_0^2(\mu_g)$ with prior degrees of freedom $d_0$:
   $$\tilde{s}_g^2 = \frac{d_0 s_0^2(\mu_g) + (m - p) s_g^2}{d_0 + (m - p)}$$
   This shrinkage stabilizes variance estimation for every individual gene while preserving true biological hyper-variability.

#### 3.5.2 Hypothesis Testing via the Quasi-Likelihood F-Test (`glmQLFTest`)
To evaluate whether a gene is significantly differentially expressed between Primary Solid Tumor ($T$) and Normal mucosal tissue ($N$), we formulate the linear contrast hypothesis:
$$H_0: \beta_{g,\text{Tumor}} = 0 \quad \text{versus} \quad H_1: \beta_{g,\text{Tumor}} \ne 0$$
Let $\ell(\hat{\beta}_g)$ denote the unconstrained log-likelihood under $H_1$, and let $\ell(\tilde{\beta}_g)$ denote the constrained log-likelihood under $H_0$.

The Quasi-Likelihood $F$-statistic is computed as:
$$F_g = \frac{2 \left[ \ell(\hat{\beta}_g) - \ell(\tilde{\beta}_g) \right]}{k \cdot \tilde{s}_g^2} = \frac{\Delta D_g}{k \cdot \tilde{s}_g^2}$$
Where:
- $\Delta D_g$ is the change in deviance resulting from dropping the `Tissue_TypeTumor` term ($k = 1$ degree of freedom).
- $\tilde{s}_g^2$ is the squeezed quasi-likelihood dispersion.
Under $H_0$, the test statistic $F_g$ follows an exact $F$-distribution:
$$F_g \sim F(k, d_{\text{post}}), \quad d_{\text{post}} = d_0 + (m - p)$$
The exact two-sided nominal $p$-value is:
$$p_g = \mathbb{P}\left( F(1, d_{\text{post}}) \ge F_g \right)$$
Because $d_{\text{post}}$ is finite, the QL $F$-test exhibits heavier tails than the asymptotic $\chi^2$ distribution, preventing false-positive discoveries in small or unbalanced reference cohorts.

#### 3.5.3 Multiple Testing Correction via False Discovery Rate (FDR)
To control the expected proportion of false discoveries among rejected null hypotheses across all $G = 16,905$ simultaneous tests, nominal $p$-values are adjusted using the Benjamini–Hochberg (BH) step-up procedure (Benjamini & Hochberg, 1995 [5]):
1. Nominal $p$-values are ordered ascendingly: $p_{(1)} \le p_{(2)} \le \dots \le p_{(G)}$.
2. The False Discovery Rate $q$-value (adjusted $p$-value) for gene ranked $i$ is:
   $$q_{(i)} = \min_{k \ge i} \left( \frac{G \cdot p_{(k)}}{k} \right)$$
3. Statistical significance is enforced at:
   $$\text{FDR}_g \le 0.05$$

---

### 3.6 Step 6: Significance Stratification & Publication Visualizations
**Script**: `master_edgeR_OSCC.R`  
**Input Artifacts**: Annotated statistics table (`tt`)  
**Output Artifacts**:
- `edgeR_significant_genes.tsv`: Filtered table containing strictly significant DEGs.
- `summary_counts.tsv`: High-level summary of tested, retained, and directional DEG counts.
- Publication figures:
  - `DEG_volcano_plot.pdf/.png`: Volcano plot depicting statistical significance versus effect magnitude.
  - `DEG_MA_plot.pdf/.png`: MA plot depicting mean abundance versus fold change.
  - `Heatmap_top25_upregulated.pdf/.png`, `Heatmap_top25_downregulated.pdf/.png`, `Heatmap_top25_up_and_down_combined.pdf/.png`: Supervised expression heatmaps.
  - `DEG_top_gene_boxplots.pdf/.png`: Expression distribution boxplots for leading candidate genes.
  - `QC_gene_variability.pdf/.png`: LOESS curve of mean logCPM versus standard deviation.

#### 3.6.1 Dual-Threshold Significance Stratification
To ensure biological relevance alongside statistical reliability, candidate differentially expressed genes (DEGs) must satisfy two simultaneous criteria:
1. **Statistical Significance**: $\text{FDR} < 0.05$ (Benjamini–Hochberg adjusted).
2. **Biological Effect Size**: $|\log_2\text{FC}| \ge 1.0$ (representing at least a 2-fold change in expression between malignant and healthy tissues):
   $$\text{Significant Up} = \{g \mid \text{FDR}_g < 0.05 \land \log_2\text{FC}_g \ge 1.0\}$$
   $$\text{Significant Down} = \{g \mid \text{FDR}_g < 0.05 \land \log_2\text{FC}_g \le -1.0\}$$

#### 3.6.2 Volcano & MA Plot Formulations
- **Volcano Plot**: Models statistical evidence as a function of biological effect size:
  $$x_g = \log_2\text{FC}_g, \quad y_g = -\log_{10}(\text{FDR}_g)$$
  Horizontal dashed threshold is positioned at $y = -\log_{10}(0.05) \approx 1.301$, and vertical dashed thresholds at $x = \pm 1.0$. The top 15 most significant genes are dynamically labeled using `ggrepel` with point-repulsion bounding boxes.
- **MA Plot**: Models biological effect size as a function of average expression abundance:
  $$x_g = \text{Average } \log_2\text{CPM}_g, \quad y_g = \log_2\text{FC}_g$$
  Verifies that differential calls are not restricted to low-abundance artifacts and confirms that fold-change symmetry is maintained across the expression spectrum.

#### 3.6.3 Supervised Hierarchical Heatmap Clustering
For publication heatmaps:
1. Normalized expression counts are transformed to $\log_2(\text{CPM})$ with $c = 2$.
2. Row-wise standard normalization ($Z$-score scaling) is performed for each gene across all 240 samples:
   $$Z_{ij} = \frac{y_{ij} - \bar{y}_i}{\text{SD}(y_i)}$$
3. Hierarchical agglomerative clustering is executed using **Ward's minimum variance method (`ward.D2`)** based on Euclidean distances:
   $$D_{\text{Ward}}(\mathcal{C}_A, \mathcal{C}_B) = \frac{|\mathcal{C}_A| |\mathcal{C}_B|}{|\mathcal{C}_A| + |\mathcal{C}_B|} \|\bar{y}_{\mathcal{C}_A} - \bar{y}_{\mathcal{C}_B}\|^2$$
   Ward's method minimizes total within-cluster variance, generating compact, biologically interpretable dendrograms.

---

### 3.7 Step 7: Curation and Two-Tier Intersection with the PIKK Gene Universe
**Script**: `intersect_PIKK.R`  
**Input Artifacts**:
- `edger_file`: Full differential expression statistics (`OSCC_DESeq2_Tumor_vs_Normal_shrunk.tsv` or `edgeR_all_results_tumor_vs_normal.tsv`).
- `pikk_file`: Curated PIKK reference universe (`PIKK_related_gene_universe.csv`).
- `vsd_file`: Variance-stabilized expression matrix container (`OSCC_vst_blind.rds`).
**Output Artifacts**:
- `OSCC_PIKK_extended_DEG_intersection_all.tsv`: Master intersection table containing all mapped PIKK-associated genes with differential statistics.
- `OSCC_PIKK_extended_DEG_intersection_significant.tsv`: Subset of significantly dysregulated PIKK-associated candidates.
- `OSCC_PIKK_extended_intersection_summary.tsv`: Contingency summary by PIKK kinase family and gene functional role.
- `OSCC_PIKK_extended_direction_counts.tsv`: Frequency table of upregulation versus downregulation.
- Publication figures:
  - `PIKK_extended_significant_DEG_heatmap.pdf/.png`: Clustered heatmap of significant PIKK candidates.
  - `PIKK_extended_dysregulated_gene_barplot.pdf/.png`: Stacked frequency bar plot of dysregulated genes per PIKK kinase group.

#### 3.7.1 Curation Architecture of the PIKK Gene Universe
The PIKK pathway reference universe ($n = 204$ functional annotations) was manually curated from primary literature, Reactome, and KEGG pathway repositories to capture the comprehensive signaling machinery of the six human PIKK kinases.

The universe is structured into two functional tiers across all six PIKK families:
1. **Core Catalytic Kinases (`GeneRole == "Core"`)**: The fundamental catalytic enzymes ($ATM$, $ATR$, $PRKDC$, $MTOR$, $TRRAP$, $SMG1$).
2. **Associated Functional Interactors (`GeneRole == "Associated"`)**:
   - Upstream DNA damage sensors and recruitment adapters (e.g., $ATRIP$, $NBN$, $MRE11$, $RAD50$, $XRCC5$, $XRCC6$).
   - Direct checkpoint kinases and downstream phosphorylation substrates (e.g., $CHEK1$, $CHEK2$, $TP53$, $BRCA1$, $H2AX$, $MDC1$, $TOPBP1$).
   - Multi-protein complex stoichiometric subunits (e.g., $RPTOR$, $RICTOR$, $MLST8$ for MTOR; $KAT5$/TIP60, $TRRAP$, $EP400$ for chromatin remodeling; $UPF1$, $UPF2$, $UPF3A$, $SMG7$, $SMG9$ for NMD surveillance).
   - Fanconi Anemia pathway partners and homologous recombination effectors (e.g., $FANCA$, $FANCD2$, $FANCI$, $FANCM$, $BLM$, $EXO1$).

#### 3.7.2 Harmonization and Set-Theoretic Inner Join
To eliminate annotation artifacts arising from capitalization discrepancies:
1. Gene symbols across both datasets are standardized to uppercase:
   $$\text{deg\_df}\$\text{gene\_symbol} \leftarrow \text{toupper}(\text{gene\_name})$$
   $$\text{pikk\_df}\$\text{GeneSymbol} \leftarrow \text{toupper}(\text{GeneSymbol})$$
2. Redundant duplicate entries in the reference universe are filtered to preserve unique group-symbol-role combinations:
   $$\text{pikk\_curated} = \text{distinct}(\text{pikk\_df}, \text{PIKK\_Group}, \text{GeneSymbol}, \text{GeneRole})$$
3. A set-theoretic inner join maps statistical metrics onto the curated biological universe:
   $$\mathcal{I} = \text{deg\_df} \bowtie_{(\text{gene\_symbol} = \text{GeneSymbol})} \text{pikk\_curated}$$

#### 3.7.3 Two-Tier Statistical Stratification Scheme
Because core regulatory kinases and epigenetic scaffold proteins often exert potent biological effects through modest, stoichiometric shifts in expression, evaluating candidates under a single, rigid cutoff risks discarding critical biological drivers. Thus, the methodology defines two explicit tiers:

```
Table 1: Two-Tier Statistical Selection Framework for PIKK Pathway Intersection
┌───────────────────────────────┬─────────────────┬─────────────────┬──────────────────────────────────────────┐
│ Stratification Tier           │ FDR Threshold   │ log2FC Floor    │ Primary Biological Rationale             │
├───────────────────────────────┼─────────────────┼─────────────────┼──────────────────────────────────────────┤
│ Tier 1: Canonical Strict      │ FDR ≤ 0.05      │ |log2FC| ≥ 1.0  │ High-magnitude transcriptional switches; │
│                               │                 │ (≥ 2.0-fold)    │ robust primary biomarker candidates.      │
├───────────────────────────────┼─────────────────┼─────────────────┼──────────────────────────────────────────┤
│ Tier 2: Extended Complex      │ FDR ≤ 0.05      │ |log2FC| ≥ 0.5  │ Stoichiometric shifts in multi-protein   │
│                               │                 │ (≥ 1.41-fold)   │ complexes (e.g., mTORC1, SAGA, NuA4).    │
└───────────────────────────────┴─────────────────┴─────────────────┴──────────────────────────────────────────┘
```

#### 3.7.4 Heatmap Clustering of PIKK Candidates
For the intersected significant PIKK candidates, expression matrices are extracted from variance-stabilized data ($V$). Row-wise $Z$-scores are computed and clustered using:
- **McQuitty's Method (WPGMA - Weighted Pair Group Method with Arithmetic Mean)**:
  When clusters $\mathcal{C}_A$ and $\mathcal{C}_B$ merge into $\mathcal{C}_{AB}$, the distance to another cluster $\mathcal{C}_C$ is the unweighted arithmetic mean:
  $$D(\mathcal{C}_{AB}, \mathcal{C}_C) = \frac{D(\mathcal{C}_A, \mathcal{C}_C) + D(\mathcal{C}_B, \mathcal{C}_C)}{2}$$
  McQuitty's linkage prevents clusters of unequal size from exerting disproportionate influence, preserving intermediate branches in structured signaling networks.
- **Average Linkage (UPGMA)**: Evaluates the weighted pairwise average distance across all sample pairs, serving as a complementary, robust clustering standard.

---

## 4. Synthesis of Methodological Parameters and Comparative Framework

```
Table 2: Comparative Methodological Parameters and Alternative Algorithmic Frameworks
┌──────────────────────────────────────┬─────────────────────────────┬─────────────────────────────┐
│ Computational Parameter / Step       │ Implemented Standard        │ Conventional Alternative    │
├──────────────────────────────────────┼─────────────────────────────┼─────────────────────────────┤
│ Library Normalization Algorithm      │ TMM (calcNormFactors)       │ RPKM / FPKM / Total Counts  │
│ Reference Library Selection          │ Mean 75th percentile sample │ Arbitrary single library    │
│ M-Value Trimming Thresholds          │ Lower/Upper 30% (|M| = 0.3) │ Untrimmed mean log-ratio    │
│ A-Value Trimming Thresholds          │ Lower/Upper 5% (|A| = 0.05) │ Untrimmed absolute mean     │
│ Negative Binomial Model Family       │ Quasi-Likelihood GLM        │ Standard Likelihood Ratio   │
│ Hypothesis Testing Distribution      │ Exact F-distribution (QLF)  │ Asymptotic Chi-Square (LRT) │
│ Variance Shrinkage Mechanism         │ Robust Empirical Bayes      │ Standard Empirical Bayes    │
│ Exploratory Ordination Metrics       │ Leading logFC MDS & PCA     │ Raw Count PCA               │
│ Prior Count for log2(CPM) Ordination │ c = 2                       │ Naive log2(x + 1)           │
│ Sample Distance Metric               │ Pearson Dissimilarity (1-r) │ Raw Euclidean distance      │
│ Primary Significance Floor           │ FDR ≤ 0.05 & |log2FC| ≥ 1.0 │ p ≤ 0.05 unadjusted         │
│ Extended PIKK Intersection Floor     │ FDR ≤ 0.05 & |log2FC| ≥ 0.5 │ Arbitrary top N ranking     │
│ Cluster Linkage for DE Heatmaps      │ Ward's minimum variance     │ Single / Centroid linkage   │
│ Cluster Linkage for PIKK Heatmaps    │ McQuitty (WPGMA) & Average  │ Complete linkage            │
└──────────────────────────────────────┴─────────────────────────────┴─────────────────────────────┘
```

---

## 5. Ready-to-Publish Methods Section (Manuscript Reference Text)

> *Below is a formal, publication-ready draft formatted in standard third-person past tense, suitable for direct incorporation into the **Materials and Methods** section of an academic manuscript.*

### 2.X Differential Gene Expression Analysis via the edgeR Quasi-Likelihood Pipeline
Differential gene expression profiling between primary oral cavity squamous cell carcinoma (OSCC; $n = 224$) and non-malignant oral mucosal tissues ($n = 16$) was performed using the `edgeR` Bioconductor framework (`master_edgeR_OSCC.R`) [6, 7]. Matched raw count matrices encompassing 16,905 robustly expressed protein-coding genes (pre-filtered via `filterByExpr` relative to the minimum library size of the normal control cohort) were compiled into a `DGEList` object. Non-malignant oral tissue was designated as the explicit reference baseline (`ref = "Normal"`). 

To eliminate composition bias resulting from asymmetric transcriptional outputs between malignant and non-malignant tissues, library sizes were normalized using the Trimmed Mean of M-values (TMM) method (`calcNormFactors`) [4]. TMM normalization factors were derived by calculating the precision-weighted mean of log-expression ratios after asymmetrically trimming the top and bottom $30\%$ of log-fold changes ($M$) and the top and bottom $5\%$ of absolute expression intensities ($A$). Effective library sizes ($N^* = N \times f$) were incorporated into all subsequent statistical models.

### 2.X Exploratory Dimensionality Reduction and Quality Control
Unsupervised sample relationships and global data architectures were evaluated prior to model fitting. Multidimensional scaling (MDS) was conducted across the top 500 leading log-fold change dimensions. To decouple mean-variance dependence for linear ordination, count matrices were converted to $\log_2$-transformed counts per million ($\log_2\text{CPM}$) with an offset of $c = 2$. Principal component analysis (PCA) was executed on the centered and standardized $\log_2\text{CPM}$ space, and the proportion of variance explained by leading orthogonal components was quantified via scree eigen-decomposition. Pairwise sample dissimilarity was computed using Pearson correlation distance ($1 - r$), followed by agglomerative hierarchical clustering with complete linkage to confirm sample cluster cohesion.

### 2.X Negative Binomial Dispersion Estimation and Generalized Linear Modeling
Biological overdispersion across patient replicates was modeled under a negative binomial framework:
$$\text{Var}(Y_{gi}) = \mu_{gi} + \phi_g \mu_{gi}^2$$
Gene-wise dispersion parameters ($\phi_g$) were estimated using the empirical Bayes method (`estimateDisp`), incorporating a robust prior (`robust = TRUE`) to shrink individual tagwise dispersions toward a global mean-dispersion trend curve without distortion from high-variance hypervariable genes [3]. 

Generalized linear models (GLMs) were fitted using the quasi-likelihood (QL) framework (`glmQLFit`) [2]. The QL formulation introduces a gene-specific error variance parameter ($s_g^2$) reflecting dispersion uncertainty, which is squeezed toward a global trend using empirical Bayes moderation. Hypothesis testing for differential expression between Tumor and Normal phenotypes was executed using the Quasi-Likelihood $F$-test (`glmQLFTest`) [1]. Unlike asymptotic likelihood ratio tests, the QLF statistic follows an exact $F$-distribution with denominator degrees of freedom reflecting finite sample size, providing rigorous type-I error rate control in cohorts with unbalanced reference groups. Nominal $p$-values were adjusted for multiple testing across all 16,905 genes using the Benjamini–Hochberg False Discovery Rate (FDR) procedure [5]. Genes satisfying an $\text{FDR} \le 0.05$ and an absolute biological effect size of $|\log_2\text{FC}| \ge 1.0$ ($\ge 2.0$-fold change) were designated as significantly differentially expressed genes (DEGs).

### 2.X Targeted PIKK Pathway Universe Curation and Two-Tier Molecular Intersection
To delineate alterations within the phosphatidylinositol 3-kinase-related kinase (PIKK) signaling axis, an annotated reference universe comprising 204 functional gene records was curated (`intersect_PIKK.R`). The reference space encompassed the six human catalytic PIKK kinases—Ataxia Telangiectasia Mutated ($ATM$), ATM and Rad3-Related ($ATR$), DNA-Dependent Protein Kinase Catalytic Subunit ($PRKDC$), Mechanistic Target of Rapamycin ($MTOR$), Transformation/Transcription Domain-Associated Protein ($TRRAP$), and Suppressor with Morphogenetic Effect on Genitalia 1 ($SMG1$)—alongside their upstream sensors, downstream checkpoint effectors, and stoichiometric regulatory complex members.

Differential expression metrics were mapped onto the curated PIKK universe via set-theoretic inner joins using standardized uppercase HGNC symbols. Intersected candidates were evaluated across two formal statistical tiers:
1. **Canonical Strict Tier**: $\text{FDR} \le 0.05$ and $|\log_2\text{FC}| \ge 1.0$ ($\ge 2.0$-fold change), identifying high-magnitude transcriptional alterations.
2. **Extended Complex Tier**: $\text{FDR} \le 0.05$ and $|\log_2\text{FC}| \ge 0.5$ ($\ge 1.41$-fold change), capturing stoichiometric expression shifts within multi-protein regulatory complexes.

For unsupervised visualization of PIKK candidates, variance-stabilized expression values were converted to row-standardized $Z$-scores and subjected to hierarchical clustering using McQuitty's weighted pair-group method with arithmetic mean (WPGMA) and average linkage (`pheatmap`) [8].

---

## 6. Comprehensive Software, Dependency, and Environment Manifest

All differential modeling and pathway intersection procedures were executed in the R statistical computing environment (R version $\ge 4.2.0$). Below is the inventory of software dependencies:

```
Table 3: Software Dependencies and Computational Infrastructure
┌─────────────────┬──────────────┬────────────────────────────────────────────────────────┐
│ Software/Package│ Origin       │ Primary Methodological Role                            │
├─────────────────┼──────────────┼────────────────────────────────────────────────────────┤
│ edgeR           │ Bioconductor │ TMM normalization, empirical Bayes dispersion, QLF GLM │
│ limma           │ Bioconductor │ Linear modeling infrastructure, removeBatchEffect      │
│ ggplot2         │ CRAN         │ High-resolution vector graphic visualization           │
│ ggrepel         │ CRAN         │ Dynamic non-overlapping text labeling in volcano plots │
│ pheatmap        │ CRAN         │ Clustered heatmaps (Ward.D2, McQuitty, Average linkage)│
│ RColorBrewer    │ CRAN         │ Colorimetric palettes (RdBu, Set1) for heatmaps/PCA    │
│ matrixStats     │ CRAN         │ High-performance row-wise variance and dispersion math │
│ gridExtra       │ CRAN         │ Multi-panel figure compilation and grid alignment      │
│ dplyr / tidyr   │ CRAN         │ Tidy data frame manipulation, grouping, pivoting       │
│ readr / tibble  │ CRAN         │ Format-preserving I/O and column-name preservation     │
│ cairo_pdf       │ Base system  │ Cairo-based anti-aliased vector graphic PDF rendering  │
└─────────────────┴──────────────┴────────────────────────────────────────────────────────┘
```

---

## 7. Peer-Reviewed Scientific Bibliography

1. **Chen Y, Lun ATL, Smyth GK.** (2016). From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. *F1000Research*, 5:1438. DOI: [10.12688/f1000research.8987.2](https://doi.org/10.12688/f1000research.8987.2).
2. **Lund SP, Nettleton D, McCarthy DJ, Smyth GK.** (2012). Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates. *Statistical Applications in Genetics and Molecular Biology*, 11(5):Article 8. DOI: [10.1515/1544-6115.1826](https://doi.org/10.1515/1544-6115.1826).
3. **Phipson B, Lee S, Majewski IJ, Alexander WS, Smyth GK.** (2016). Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression. *Annals of Applied Statistics*, 10(2):946–963. DOI: [10.1214/16-AOAS920](https://doi.org/10.1214/16-AOAS920).
4. **Robinson MD, Oshlack A.** (2010). A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology*, 11(3):R25. DOI: [10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25).
5. **Benjamini Y, Hochberg Y.** (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B (Methodological)*, 57(1):289–300. DOI: [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).
6. **Robinson MD, McCarthy DJ, Smyth GK.** (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26(1):139–140. DOI: [10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616).
7. **McCarthy DJ, Chen Y, Smyth GK.** (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. *Nucleic Acids Research*, 40(10):4288–4297. DOI: [10.1093/nar/gks042](https://doi.org/10.1093/nar/gks042).
8. **Kolde R.** (2019). *pheatmap: Pretty Heatmaps*. R package version 1.0.12. CRAN: [https://CRAN.R-project.org/package=pheatmap](https://CRAN.R-project.org/package=pheatmap).
9. **Love MI, Huber W, Anders S.** (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12):550. DOI: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8).
10. **The Cancer Genome Atlas Network.** (2015). Comprehensive genomic characterization of head and neck squamous cell carcinomas. *Nature*, 517(7536):576–582. DOI: [10.1038/nature14129](https://doi.org/10.1038/nature14129).
11. **Blackford AN, Jackson SP.** (2017). ATM, ATR, and DNA-PK: The Trinity at the Heart of the DNA Damage Response. *Molecular Cell*, 66(6):801–817. DOI: [10.1016/j.molcel.2017.05.015](https://doi.org/10.1016/j.molcel.2017.05.015).
12. **Saxton RA, Sabatini DM.** (2017). mTOR Signaling in Growth, Metabolism, and Disease. *Cell*, 168(6):960–976. DOI: [10.1016/j.cell.2017.02.004](https://doi.org/10.1016/j.cell.2017.02.004).
