# Tumor Microenvironment Deconvolution, Immune Infiltration Modeling, and Therapeutic Druggability Profiling (Section 2.6) — Comprehensive Methodology Guide

**Project**: PIKK Pathway-Associated Biomarker Discovery and Oncogenic Characterization in Oral Squamous Cell Carcinoma (OSCC)  
**Target Manuscript Reference**: Materials and Methods / Bioinformatic Processing Pipeline  
**Primary Dataset**: The Cancer Genome Atlas Head and Neck Squamous Cell Carcinoma (TCGA-HNSC), Curated Oral Cavity Squamous Cell Carcinoma Cohort ($n = 240$ primary biospecimens; $n = 224$ Primary Solid Tumor, $n = 16$ Solid Tissue Normal)  
**Deconvolution Mixture Cohort**: $n = 224$ primary solid tumor transcriptomes converted to non-negative linear abundance space  
**Quality-Controlled Deconvolution Cohort**: $n = 137$ primary solid tumors passing empirical permutation significance ($P_{\text{perm}} < 0.05$) in CIBERSORTx  
**ESTIMATE Cellularity Cohort**: $n = 224$ primary solid tumors evaluated for ImmuneScore, StromalScore, and ESTIMATE TumorPurity ($n = 137$ overlapping with the QC-filtered deconvolution cohort)  
**Target Gene Panels**:
1. **Core Diamond Panel** ($n = 7$ genes): $PLK1$, $CDK2$, $TOPBP1$, $RAD51$, $FANCI$, $KAT2B$, $DEPTOR$ (representing the key prognostic and pathway-representative drivers across ATR, ATM, TRRAP, and mTOR axes).
2. **Comprehensive Hub Gene Space** ($n = 24$ genes): Full set of topologically prioritized network hubs ($PLK1$, $CDK2$, $TOPBP1$, $RAD51$, $FANCI$, $KAT2B$, $DEPTOR$, $EGFR$, $BRCA1$, $BRCA2$, $PARP1$, $H2AX$, $RUVBL1$, $CHEK1$, $CHEK2$, $AURKA$, $AURKB$, $E2F1$, $EXO1$, $CDC45$, $MRGBP$, $RPTOR$, $MLST8$, $RICTOR$).  
**Scripts Evaluated**:
1. `01_pikk_druggability_analysis.py` — *DGIdb v5.0 GraphQL Knowledgebase Mining, Regulatory Approval Tier Harmonization, and Multi-Signal Composite Druggability Scoring*
2. `02_pikk_druggability_plots.R` — *Weighted Sub-Signal Decomposition, Bipartite Target-Drug Interaction Cartography, and Translational Actionability Matrix Generation*
3. `03_prepare_cibersortx_input.R` — *Log-to-Linear Abundance Transformation, Tumor Filtering, Non-Negative Bounding, and Remote Deconvolution Packaging*
4. `04_immune_tme_analysis.R` — *CIBERSORTx LM22 High-Resolution Deconvolution Analysis, Permutation QC Filtering, ESTIMATE Tumor Purity Confounder Control, and Master Translational Synthesis*

---

## 1. Executive Summary & End-to-End Workflow

Identifying prognostic biomarkers through differential expression, network topology, and survival modeling reveals the molecular regulators of malignancy. However, translating these discoveries into clinical utility requires addressing two fundamental translational questions:
1. **Tumor Microenvironment (TME) Modulation & Immune Evasion**: Does the hyper-activation of PIKK-driven checkpoint and DNA repair machinery remodel the surrounding tumor immune microenvironment, promoting cytotoxic T-cell exclusion or immunosuppressive macrophage polarization?
2. **Therapeutic Tractability & Pharmacogenomic Actionability**: Are these prioritized drivers directly druggable by existing FDA-approved antineoplastics, clinical trial investigational agents, or selective chemical probes, enabling immediate clinical translation or repurposing?

This document formalizes the publication-grade computational methodology for **Tumor Microenvironment (TME) Deconvolution, Orthogonal Tumor Purity Modeling, and Multi-Signal Therapeutic Druggability Profiling (Section 2.6)**.

The computational architecture comprises four integrated analytical phases:
1. **Linear-Space Transformation & Remote CIBERSORTx Ingestion (`03_prepare_cibersortx_input.R`)**: Converts variance-stabilized $\log_2\text{CPM}$ expression into non-negative linear space, format-compliant with the Stanford CIBERSORTx server, deconvolving 22 human hematopoietic lineages (LM22 signature matrix).
2. **Permutation-Filtered Immune Deconvolution Analysis (`04_immune_tme_analysis.R`)**: Applies strict empirical permutation quality control ($P_{\text{perm}} < 0.05$) to exclude unreliable deconvolution mixtures ($n = 137$ retained), systematically evaluates non-parametric Spearman correlations between PIKK Diamond Panel targets and all 22 immune cell fractions (154 simultaneous hypothesis tests), and computes median-split Wilcoxon rank-sum tests for key immune lineages.
3. **Orthogonal ESTIMATE Purity Confounder Modeling (`04_immune_tme_analysis.R`)**: Runs the ESTIMATE algorithm on full-transcriptome bulk profiles ($n = 224$) to calculate ImmuneScore, StromalScore, and cosine-transformed TumorPurity, verifying that observed biomarker associations are not confounding artifacts of tumor cellularity.
4. **Pharmacogenomic Interrogation & Multi-Signal Druggability Scoring (`01_pikk_druggability_analysis.py` & `02_pikk_druggability_plots.R`)**: Automates GraphQL queries against the Drug-Gene Interaction Database (DGIdb v5.0), extracts clinical regulatory approval tiers, and computes a continuous composite druggability score integrating approval maturity, interaction confidence, source database diversity, and drug count.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                                 UPSTREAM INPUT ARTIFACTS                               │
│  • Primary Expression Matrix: logCPM_matrix.tsv (Tumor = 224, Normal = 16)             │
│  • Curated Clinical Metadata: matched_metadata_used.tsv (Deduplicated TCGA-HNSC OSCC)  │
│  • Prognostic Hub Manifest: hub_gene_characterisation_table.tsv                        │
│  • Core Diamond Target Panel: PLK1, CDK2, TOPBP1, RAD51, FANCI, KAT2B, DEPTOR        │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                      ┌─────────────────────┴─────────────────────┐
                      ▼                                           ▼
┌───────────────────────────────────────────┐ ┌───────────────────────────────────────────┐
│  TRACK A: PHARMACOGENOMIC DRUGGABILITY    │ │  TRACK B: TME DECONVOLUTION & CELLULARITY │
│  (01_pikk_druggability_analysis.py)       │ │  (03_prepare_cibersortx_input.R)          │
│  • Query DGIdb v5.0 GraphQL API           │ │  • Filter to primary oral tumors (n = 224)│
│    (genes: 24 Hubs + 7 Diamond Targets)   │ │  • Linearize: y = max(2^x - 1, 0)         │
│  • Two-pass querying:                     │ │  • Export: cibersortx_input_matrix.tsv    │
│    - Pass 1: Gene interaction records     │ └─────────────────────┬─────────────────────┘
│    - Pass 2: Curated drug approvals/MOA   │                       │
│  • Multi-Signal Composite Formulation:    │                       ▼
│    Comp = 0.35*App + 0.30*Int +           │ ┌───────────────────────────────────────────┐
│           0.20*Src + 0.15*Drugs           │ │  STANFORD CIBERSORTx REMOTE DECONVOLUTION │
│  • Categorical Clinical Maturity Tiers    │ │  • Signature: LM22 (22 immune phenotypes) │
│  • Export composite & raw TSVs            │ │  • Algorithm: ν-Support Vector Regression │
└─────────────────────┬─────────────────────┘ │  • Permutations: B = 500 null shuffles   │
                      │                       │  • Deliverable: cibersortx_results.csv    │
                      ▼                       └─────────────────────┬─────────────────────┘
┌───────────────────────────────────────────┐                       │
│  DRUGGABILITY CARTOGRAPHY                 │                       ▼
│  (02_pikk_druggability_plots.R)           │ ┌───────────────────────────────────────────┐
│  • Stacked sub-signal breakdown plots     │ │  QC FILTERING & IMMUNE PROFILING          │
│  • Bipartite target-drug interaction maps │ │  (04_immune_tme_analysis.R)               │
│  • Faceted lollipop plots of lead agents  │ │  • Permutation QC: P_perm < 0.05 (n = 137)│
│  • Translational Actionability Matrices   │ │  • Spearman correlations: 7 genes x 22    │
└─────────────────────┬─────────────────────┘ │    immune phenotypes (154 tests, BH FDR)  │
                      │                       │  • Median-split Wilcoxon rank-sum tests   │
                      │                       │    (CD8+ T-cells, M2, Tregs, M1)          │
                      │                       └─────────────────────┬─────────────────────┘
                      │                                             │
                      │                       ┌─────────────────────┘
                      │                       ▼
                      │       ┌───────────────────────────────────────────┐
                      │       │  ORTHOGONAL ESTIMATE PURITY CONTROLS      │
                      │       │  (04_immune_tme_analysis.R)               │
                      │       │  • Filter to 141 Immune & 141 Stromal     │
                      │       │    curated signature genes (Yoshihara)    │
                      │       │  • Compute ImmuneScore & StromalScore     │
                      │       │  • Cosine transformation: TumorPurity     │
                      │       │  • Linear regression & Spearman tests vs  │
                      │       │    PLK1/CDK2 (purity-confounder control)  │
                      │       └─────────────────────┬─────────────────────┘
                      │                             │
                      └──────────────────────┬──────┘
                                             │
                                             ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  MASTER TRANSLATIONAL SYNTHESIS & REPORTING (04_immune_tme_analysis.R)                 │
│  • Integration of Differential Expression, PPI Centrality, Survival Independence,     │
│    Immune Remodeling Correlations, Tumor Purity Independence, and Druggability Scores  │
│  • Master Deliverable: clinical_relevance_biomarker_table.tsv                          │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. Biological Context and Study Rationale

### 2.1 The Immune Landscape of OSCC: Replication Stress and TME Remodeling
Oral Squamous Cell Carcinoma (OSCC) arises in a complex, immunologically active mucosal environment. The genomic hallmark of HPV-negative OSCC—near-universal $TP53$ loss coupled with $CDK2NA$ inactivation—drives persistent replication stress and endogenous DNA damage. 

In classical immunology, damaged cells accumulate cytosolic DNA fragments (via micronuclei collapse or stalled replication fork breakdown), triggering the cyclic GMP-AMP synthase (cGAS)–stimulator of interferon genes (STING) pathway and eliciting type I interferon release. This pathway recruits CD8+ cytotoxic T lymphocytes and promotes anti-tumor immunity ("hot" immune phenotype). 

However, aggressive OSCC tumors frequently rewire DNA damage response (DDR) and cell cycle pathways to evade immune surveillance:
- **Mitotic Checkpoint Override ($PLK1$, $CDK2$)**: Rapid checkpoint transit and hyper-efficient homologous recombination repair prevent the accumulation of cytosolic double-stranded DNA, silencing cGAS-STING signaling.
- **Immunosuppressive Cell Recruitment**: Tumors displaying elevated PIKK signaling remodel the tumor microenvironment (TME) toward an immunosuppressive state ("cold" or immune-excluded phenotype), characterized by the recruitment of M2-polarized tumor-associated macrophages (TAMs), regulatory T cells ($T_{\text{regs}}$), and myeloid-derived suppressor cells (MDSCs).
- **Extracellular Stroma & Physical Barrier Formation**: Coordinated TRRAP-mediated chromatin remodeling and mTOR-driven metabolic shifts stimulate cancer-associated fibroblasts (CAFs), establishing dense desmoplastic stroma that physically excludes CD8+ T cells from the tumor core.

Characterizing the immune landscape associated with key PIKK drivers determines whether targeted inhibition of these kinases can overcome immune exclusion and convert "cold" OSCC tumors into immunologically responsive lesions susceptible to immune checkpoint blockade (e.g., anti-PD-1 / anti-PD-L1).

### 2.2 Methodological Rationale: Justification of Selected Analytical Paradigms

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                   METHODOLOGICAL COMPARISON: TME & DRUGGABILITY PARADIGMS              │
├───────────────────────────┬──────────────────────────────┬─────────────────────────────┤
│ Dimension                 │ Selected Pipeline Paradigm   │ Conventional Alternative    │
├───────────────────────────┼──────────────────────────────┼─────────────────────────────┤
│ Immune Deconvolution      │ CIBERSORTx LM22              │ Single-Sample GSEA (ssGSEA) │
│                           │ (ν-Support Vector Regression)│ or Mean Marker Expression   │
│ Rationale                 │ Resolves overlapping shared  │ Prone to extreme gene-gene  │
│                           │ lineage transcripts across   │ correlation bias; cannot    │
│                           │ 22 immune cell phenotypes;   │ infer relative fractional   │
│                           │ provides empirical p-values  │ proportions of mixtures     │
├───────────────────────────┼──────────────────────────────┼─────────────────────────────┤
│ Quality Control           │ Permutation Significance     │ Unfiltered Sample Ingestion │
│                           │ Filter (P_perm < 0.05)       │ (analyzing all samples)     │
│ Rationale                 │ Discards deconvolution fits  │ Retains computational noise;│
│                           │ that do not differ from      │ distorted mixtures produce  │
│                           │ random biological noise      │ spurious immune correlations│
├───────────────────────────┼──────────────────────────────┼─────────────────────────────┤
│ Cellularity Control       │ Orthogonal ESTIMATE Modeling │ No Purity Assessment        │
│                           │ (Immune/Stromal/TumorPurity) │ (uncontrolled bulk assays)  │
│ Rationale                 │ Verifies that biomarker      │ Vulnerable to cellularity   │
│                           │ expression is cell-intrinsic │ confounding (high purity    │
│                           │ rather than a purity artifact│ mimics biomarker induction) │
├───────────────────────────┼──────────────────────────────┼─────────────────────────────┤
│ Druggability Evaluation   │ Multi-Signal Continuous      │ Binary Literature Search    │
│                           │ Scoring (DGIdb v5.0 GraphQL) │ ("Druggable" Yes/No)        │
│ Rationale                 │ Synthesizes regulatory tier, │ Subjective; ignores clinical│
│                           │ interaction confidence, and  │ trial pipelines and database│
│                           │ database evidence breadth    │ provenance/evidence quality │
└───────────────────────────┴──────────────────────────────┴─────────────────────────────┘
```

#### 1. Why $\nu$-SVR Deconvolution (CIBERSORTx) Over Single-Gene Signatures?
Bulk RNA sequencing measures an aggregate transcriptional average of millions of malignant, stromal, and immune cells. Using naive marker genes (e.g., $CD8A$ for cytotoxic T cells or $FOXP3$ for $T_{\text{regs}}$) is confounded because many surface markers are shared across hematopoietic lineages (e.g., $CD4$ on monocytes and macrophages, $CD8$ on NK cell subsets). 

The pipeline uses **CIBERSORTx** (Newman et al., *Nat Methods* 2015 [1]; *Nat Biotechnol* 2019 [2]), which applies **$\nu$-support vector regression ($\nu$-SVR)**. This machine-learning technique minimizes linear loss functions while penalizing regression weights, making it robust against high collinearity among cell-type expression signatures.

#### 2. Why Strict Empirical Permutation QC is Mandatory
Deconvolution algorithms attempt to solve an underdetermined system of linear equations. If a tumor biospecimen has experienced severe RNA degradation, contains aberrant uncharacterized cell populations, or exhibits massive aneuploidy, the resulting cell fractions can represent computational noise rather than true biology. 

Enforcing an empirical permutation filter:
$$P_{\text{perm}} < 0.05 \quad (\text{filtering } 224 \to 137 \text{ tumors})$$
guarantees that every retained tumor possesses a signature deconvolution correlation significantly exceeding that of a randomly scrambled transcriptome ($B = 500$ permutations).

#### 3. Why Orthogonal ESTIMATE Modeling is Necessary
In bulk tumor profiling, an apparent negative correlation between a gene (e.g., $PLK1$) and immune infiltration (e.g., CD8+ T cells) could arise from a trivial technical confounder: tumor purity. If high-grade tumors have higher cellularity (95% cancer cells, 5% immune cells) compared to lower-grade tumors (50% cancer cells, 50% immune cells), any gene overexpressed in cancer cells will correlate negatively with immune signatures simply due to cell dilution. 

The pipeline deploys **ESTIMATE** (Yoshihara et al., *Nat Commun* 2013 [3]), an independent single-sample method that infers tumor purity from stroma- and immune-specific gene signatures. Interrogating the relationship between candidate biomarkers and ESTIMATE TumorPurity proves whether biomarker expression is truly independent of sample cellularity.

#### 4. Why Multi-Signal Druggability Scoring?
Traditional drug tractability classifications rely on binary designations (e.g., classifying a target as "undruggable" if no approved drug exists). In modern oncology, many targets lack FDA-approved indications but possess potent phase I/II clinical trial inhibitors or selective chemical probes. 

The pipeline implements a **Multi-Signal Composite Druggability Score**:
$$\text{Composite Score} = 0.35 \cdot S_{\text{app}} + 0.30 \cdot S_{\text{int}} + 0.20 \cdot S_{\text{src}} + 0.15 \cdot S_{\text{drugs}}$$
This continuous score ($0.0 \to 1.0$) integrates regulatory approval status, interaction confidence scores, source database diversity (e.g., ChEMBL, DrugBank, CIViC, PharmGKB), and candidate drug counts, systematically identifying repurposing opportunities.

---

## 3. Step-by-Step Computational Methodology

### 3.1 Step 1: Pharmacogenomic Ingestion & Multi-Signal Druggability Scoring
**Script**: `01_pikk_druggability_analysis.py`  
**Input Artifacts**:
- Core Diamond Panel gene list: `PLK1`, `CDK2`, `TOPBP1`, `RAD51`, `FANCI`, `KAT2B`, `DEPTOR`.
- Comprehensive Hub Gene manifest: 24 topological hub genes.
**Output Artifacts**:
- `pikk_druggability_composite_scores.tsv`: Comprehensive scores and clinical maturity tiers for all evaluated genes.
- `dgidb_raw_gene_interactions.tsv`: Granular edge-level table of all drug-target interactions, confidence scores, and literature sources.
- `dgidb_druggability_run.log`: Execution log tracking GraphQL network requests and batch pagination.

#### 3.1.1 Knowledgebase Interrogation via the DGIdb GraphQL API
Drug-gene interactions are retrieved programmatically from the **Drug-Gene Interaction Database (DGIdb v5.0)** (Freshour et al., *Nucleic Acids Res* 2021 [4]) via its GraphQL API endpoint:
$$\text{Endpoint: } \texttt{https://dgidb.org/api/graphql}$$

To handle rate limits and ensure robust retrieval, network communication is managed using an HTTP session adapter with an exponential backoff retry strategy ($3$ retries, backoff factor $1.5$, handling HTTP 429/500/502/503/504 errors).

A two-pass querying strategy is implemented:
- **Pass 1 (Gene Interaction Query)**: Queries the `genes(names: $names)` node to retrieve all physical and functional drug associations, interaction confidence scores, mechanisms of action (MOA), interaction types (e.g., inhibitor, antagonist, antibody), and reporting source databases.
- **Pass 2 (Curated Drug Regulatory Query)**: Extracts all unique drug chemical entities identified in Pass 1, removes salt suffixes (e.g., "HYDROCHLORIDE", "MALEATE"), and queries the `drugs(names: $names)` node in chunks of $50$ to extract canonical FDA approval flags and clinical disease indications.

#### 3.1.2 Mathematical Formulation of Druggability Sub-Signals
For each evaluated gene $g$, four orthogonal evidence signals are computed:

##### 1. Regulatory Approval Tier Subscore ($S_{\text{app}}$)
Maps the highest clinical development milestone achieved by any interacting drug into a bounded continuous score:
$$S_{\text{app}}(g) = \max_{d \in \mathcal{D}(g)} \text{Tier}(d)$$
Where $\text{Tier}(d)$ is assigned using an empirical hierarchy:
$$\text{Tier}(d) = \begin{cases} 
1.00 & \text{if FDA Approved / Companion Diagnostic} \\
0.90 & \text{if NCCN / Clinical Guideline Supported} \\
0.75 & \text{if Phase III Clinical Trial} \\
0.60 & \text{if Phase II / Phase Ib/II Clinical Trial} \\
0.45 & \text{if Phase I Clinical Trial} \\
0.30 & \text{if Observational Study / Case Series} \\
0.15 & \text{if Preclinical (Biochemical, Cell Line, PDX)} \\
0.00 & \text{if No Known Drug / Novel Target}
\end{cases}$$
If a drug possesses a verified positive approval flag in DGIdb Pass 2, $S_{\text{app}}(g)$ defaults to $1.00$.

##### 2. Interaction Evidence Confidence Subscore ($S_{\text{int}}$)
DGIdb computes interaction confidence scores based on publication citations and database evidence. To prevent outlier interactions from dominating, the maximum observed interaction score $I_{\max}(g)$ is normalized using log-scaling capped at $100$:
$$S_{\text{int}}(g) = \frac{\ln\left( 1 + \min(I_{\max}(g), 100) \right)}{\ln(1 + 100)}$$

##### 3. Candidate Drug Count Subscore ($S_{\text{drugs}}$)
Measures the depth of existing chemical matter targeting gene $g$. The total number of unique interacting drugs $N_{\text{drugs}}(g)$ is log-transformed and capped at $20$ molecules:
$$S_{\text{drugs}}(g) = \frac{\ln\left( 1 + \min(N_{\text{drugs}}(g), 20) \right)}{\ln(1 + 20)}$$

##### 4. Source Database Diversity Subscore ($S_{\text{src}}$)
Quantifies evidentiary consensus across independent pharmacological databases. Let $\mathcal{S}_{\text{known}}$ denote the set of $15$ curated oncology and pharmacology databases:
$$\mathcal{S}_{\text{known}} = \{\text{PharmGKB, CIViC, OncoKB, CGI, DoCM, FDA, CKB-CORE, Clearity, DrugBank, ChEMBL, TTD, TdgClinical, TALC, NCI, MyCancerGenome}\}$$
$$S_{\text{src}}(g) = \frac{|\mathcal{S}_{\text{observed}}(g) \cap \mathcal{S}_{\text{known}}|}{|\mathcal{S}_{\text{known}}|}$$

#### 3.1.3 Composite Druggability Score Formulation
The final **Composite Druggability Score** $S_{\text{comp}}(g) \in [0.0, 1.0]$ is computed as a weighted linear combination:
$$S_{\text{comp}}(g) = w_{\text{app}} S_{\text{app}}(g) + w_{\text{int}} S_{\text{int}}(g) + w_{\text{src}} S_{\text{src}}(g) + w_{\text{drugs}} S_{\text{drugs}}(g)$$
Where the empirical weighting vector $\mathbf{w}$ prioritizes clinical maturity and evidence quality:
$$w_{\text{app}} = 0.35, \quad w_{\text{int}} = 0.30, \quad w_{\text{src}} = 0.20, \quad w_{\text{drugs}} = 0.15 \quad \left( \sum w_i = 1.0 \right)$$

#### 3.1.4 Categorical Clinical Maturity Tier Assignment
Based on the computed subscores, every target is assigned to a discrete clinical maturity tier:
$$\text{Maturity Tier} = \begin{cases}
\text{"FDA Approved"} & \text{if } S_{\text{app}} \ge 1.00 \\
\text{"Phase III Clinical"} & \text{if } 0.70 \le S_{\text{app}} < 1.00 \\
\text{"Phase I/II Clinical"} & \text{if } 0.40 \le S_{\text{app}} < 0.70 \\
\text{"Preclinical / Experimental"} & \text{if } 0.00 < S_{\text{app}} < 0.40 \lor N_{\text{drugs}} > 0 \\
\text{"Undruggable / Novel"} & \text{if } S_{\text{app}} = 0.00 \land N_{\text{drugs}} = 0
\end{cases}$$

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        DRUGGABILITY PIPELINE PARAMETER MANIFEST                        │
├────────────────────────────────────────────────────────┬───────────────────────────────┤
│ Processing Parameter / Dimension                       │ Value / Configuration         │
├────────────────────────────────────────────────────────┼───────────────────────────────┤
│ Target Knowledgebase                                   │ DGIdb v5.0 (GraphQL API)      │
│ Target Evaluated Genes                                 │ G = 24 genes (incl. 7 Diamond)│
│ Core Diamond Panel Genes                               │ PLK1, CDK2, TOPBP1, RAD51,    │
│                                                        │ FANCI, KAT2B, DEPTOR          │
│ Known Curated Source Databases (|S_known|)             │ 15 databases                  │
│ Maximum Drug Count Normalization Cap                   │ N_max = 20 drugs              │
│ Maximum Interaction Score Normalization Cap            │ I_max = 100.0                 │
│ Composite Weight Distribution                          │ App: 0.35, Int: 0.30,         │
│                                                        │ Src: 0.20, Drugs: 0.15        │
│ Network HTTP Request Retry Configuration               │ 3 retries, backoff factor 1.5 │
│ Network Batch Chunk Size                               │ 50 drugs per batch            │
└────────────────────────────────────────────────────────┴───────────────────────────────┘
```

---

### 3.2 Step 2: Druggability Cartography & Bipartite Visualizations
**Script**: `02_pikk_druggability_plots.R`  
**Input Artifacts**: `pikk_druggability_composite_scores.tsv`, `dgidb_raw_gene_interactions.tsv`.  
**Output Artifacts**:
- `01_druggability_composite_score_breakdown.pdf` and `.png`: Stacked bar chart illustrating weighted sub-signal contributions.
- `02_drug_target_interaction_network_bipartite.pdf` and `.png`: Faceted lollipop plot mapping lead small molecules to targets.
- `03_therapeutic_actionability_heatmap.pdf` and `.png`: Categorical actionability matrix for the Diamond Panel.

#### 3.2.1 Sub-Signal Stacked Decomposition
To visualize the multi-signal composition of each target, the four weighted terms:
$$\{0.35 \cdot S_{\text{app}}, \; 0.30 \cdot S_{\text{int}}, \; 0.20 \cdot S_{\text{src}}, \; 0.15 \cdot S_{\text{drugs}}\}$$
are decomposed into stacked horizontal bar plots (`ggplot2`) for all Diamond targets and top-scoring hubs ($S_{\text{comp}} > 0.40$), with labels displaying composite scores and maturity tiers.

#### 3.2.2 Bipartite Target-Drug Interaction Curation
For lead target visualization:
1. Candidate drugs are filtered to remove raw chemical database identifiers (e.g., `CHEMBL:...`) and lengthy IUPAC chemical descriptions ($> 30$ characters).
2. Trailing chemical salt suffixes (e.g., "HYDROCHLORIDE", "MALEATE") are removed.
3. Up to six informative antineoplastics or trial inhibitors are selected per target.
4. For undruggable novel targets (e.g., $TOPBP1$, $FANCI$), placeholder indicators are inserted ("No Direct Small Molecule [Novel Target]") to maintain comprehensive panel reporting.
5. Visualized as a multi-panel faceted lollipop plot colored by clinical regulatory status: *FDA Approved* (red `#E64B35`), *Investigational / Trial* (blue `#3C5488`), and *Undruggable / Novel* (grey `#8491B4`).

---

### 3.3 Step 3: Linear-Space Matrix Transformation & CIBERSORTx Preprocessing
**Script**: `03_prepare_cibersortx_input.R`  
**Input Artifacts**: `matched_metadata_used.tsv`, `logCPM_matrix.tsv`.  
**Output Artifacts**: `cibersortx_input_matrix.tsv` (formatted linear expression matrix).

#### 3.3.1 Linear-Space Conversion Formulation
CIBERSORTx assumes that bulk RNA expression represents a **linear convolution** of constituent cell-type profiles:
$$\mathbf{m} = \sum_{k=1}^K f_k \mathbf{s}_k$$
Where $\mathbf{m}$ is the observed bulk tissue expression vector, $\mathbf{s}_k$ is the signature vector of cell type $k$, and $f_k \ge 0$ is the relative cell fraction ($\sum f_k = 1$). 

Because normalized input counts were previously transformed to $\log_2(\text{CPM} + 1)$ for differential expression and survival modeling, passing logarithmic data directly into CIBERSORTx violates the linear mixing assumption. 

Every entry $x_{gj}$ in the expression matrix is converted back to **non-negative linear space**:
$$y_{gj} = \max\left( 2^{x_{gj}} - 1, \; 0 \right)$$
Where $y_{gj}$ represents linear CPM abundance for gene $g$ in tumor sample $j$.

#### 3.3.2 Sample Filtering and File Serialization
1. Restrict columns to verified Primary Solid Tumors (`Condition == "Tumor"`, $n = 224$).
2. Ensure gene symbols are unique by deduplicating identifiers (`distinct(GeneSymbol)`).
3. Export the formatted linear matrix as a tab-delimited file (`cibersortx_input_matrix.tsv`) with genes in rows and sample identifiers in columns.

---

### 3.4 Step 4: CIBERSORTx Immune Deconvolution & Empirical Permutation QC
**Execution Platform**: Stanford CIBERSORTx Cloud Infrastructure (`https://cibersortx.stanford.edu`) [2]  
**Analysis Script**: `04_immune_tme_analysis.R`  
**Input Artifacts**: `cibersortx_input_matrix.tsv`, `cibersortx_results.csv`.  
**Output Artifacts**:
- `ciber_qc`: Filtered internal deconvolution matrix ($n = 137$ samples passing QC).
- `07_cibersortx_rmse_quality_control.pdf` and `.png`: Deconvolution Root Mean Squared Error (RMSE) diagnostic plots.

#### 3.4.1 $\nu$-Support Vector Regression ($\nu$-SVR) Formulation
High-resolution cellular deconvolution was conducted using the **LM22 leukocyte gene signature matrix** ($547$ signature genes resolving $22$ human immune cell phenotypes):
$$\mathcal{C}_{\text{LM22}} = \{\text{B cells naive}, \text{B cells memory}, \text{Plasma cells}, \text{T cells CD8}, \text{T cells CD4 naive}, \dots, \text{Neutrophils}\}$$

For each bulk tumor mixture vector $\mathbf{y}_j \in \mathbb{R}^{547}$, CIBERSORTx solves a constrained **$\nu$-support vector regression** optimization problem:
$$\min_{\mathbf{w}_j, b_j, \boldsymbol{\xi}_j, \boldsymbol{\xi}_j^*} \frac{1}{2} \|\mathbf{w}_j\|^2 + C \left( \nu \epsilon_j + \frac{1}{G} \sum_{g=1}^G (\xi_{jg} + \xi_{jg}^*) \right)$$
Subject to:
$$\left( \mathbf{w}_j^T \mathbf{S}_g + b_j \right) - y_{jg} \le \epsilon_j + \xi_{jg}$$
$$y_{jg} - \left( \mathbf{w}_j^T \mathbf{S}_g + b_j \right) \le \epsilon_j + \xi_{jg}^*$$
$$\xi_{jg} \ge 0, \quad \xi_{jg}^* \ge 0, \quad \epsilon_j \ge 0$$
Where $\mathbf{S}_g \in \mathbb{R}^{22}$ is the signature row vector for gene $g$ across the 22 immune cell types, $\mathbf{w}_j \in \mathbb{R}^{22}$ is the regression weight vector, and $\xi, \xi^*$ are slack variables.

Non-negative constraints are enforced, and estimated immune cell fractions $f_{jk}$ are calculated by normalizing the positive regression coefficients:
$$f_{jk} = \frac{\max(w_{jk}, 0)}{\sum_{m=1}^{22} \max(w_{jm}, 0)}, \quad \sum_{k=1}^{22} f_{jk} = 1$$

#### 3.4.2 Empirical Permutation Test and Quality Control Filtering
To evaluate whether the inferred cell fractions reflect biological signal rather than fitting noise, CIBERSORTx conducts an **empirical permutation test** ($B = 500$ permutations).

Let $r_{\text{obs}, j}$ denote the Pearson correlation coefficient between the observed bulk mixture vector $\mathbf{y}_j$ and the model-reconstructed expression profile $\hat{\mathbf{y}}_j = \mathbf{S} \mathbf{f}_j$:
$$r_{\text{obs}, j} = \text{Corr}(\mathbf{y}_j, \hat{\mathbf{y}}_j)$$

In each permutation $b \in \{1, \dots, B\}$, the gene labels of $\mathbf{y}_j$ are randomly shuffled, and $\nu$-SVR deconvolution is recomputed to obtain the null correlation $r_{b, j}^*$. The empirical permutation $p$-value is defined as:
$$P_{\text{perm}, j} = \frac{1}{B} \sum_{b=1}^B \mathbb{I}(r_{b, j}^* \ge r_{\text{obs}, j})$$

Samples failing empirical significance ($P_{\text{perm}} \ge 0.05$) exhibit deconvolution profiles indistinguishable from random noise and are filtered out:
$$\mathcal{S}_{\text{QC}} = \{ j \mid P_{\text{perm}, j} < 0.05 \}$$
Applying this filter retains **$n = 137$ out of $224$ primary oral tumors** (61.16% retention rate), ensuring high analytical fidelity.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        CIBERSORTx DECONVOLUTION QC MANIFEST                            │
├────────────────────────────────────────────────────────┬───────────────────────────────┤
│ Cohort / Processing Dimension                          │ Quantitative Metric           │
├────────────────────────────────────────────────────────┼───────────────────────────────┤
│ Total Deconvolved Primary Oral Tumors                  │ N = 224 tumors                │
│ Reference Signature Matrix                             │ LM22 (547 signature genes)    │
│ Immune Lineages Deconvolved                            │ K = 22 cell phenotypes        │
│ Permutation Shuffles for Significance Testing          │ B = 500 permutations          │
│ Permutation Quality Threshold                          │ P_perm < 0.05                 │
│ Samples Passing Permutation QC Filter                  │ n = 137 tumors (61.16%)       │
│ Samples Discarded Due to Poor Fit (P_perm ≥ 0.05)      │ n = 87 tumors (38.84%)        │
│ Final Matched Immune Infiltration Cohort               │ n = 137 primary oral tumors   │
└────────────────────────────────────────────────────────┴───────────────────────────────┘
```

---

### 3.5 Step 5: Non-Parametric Immune Correlation & Median-Split Stratification
**Script**: `04_immune_tme_analysis.R`  
**Input Artifacts**: Matched QC-filtered deconvolution table ($n = 137$), matched expression matrix.  
**Output Artifacts**:
- `pikk_immune_correlation_matrix.tsv`: Master correlation matrix reporting Spearman $\rho$, nominal $p$-values, and BH-adjusted $q$-values across all 154 pairs.
- `04_cibersortx_immune_correlation_heatmap.pdf` and `.png`: Comprehensive tile plot of all 22 immune cell types vs. Diamond Panel targets.
- `05_stratified_immune_boxplots_plk1_cdk2.pdf` and `.png`: Median-split boxplots with Wilcoxon rank-sum test annotations.

#### 3.5.1 Multi-Lineage Spearman Rank Correlation Analysis
To interrogate monotonic relationships without assuming linear scaling or normal bivariate distributions, the **Spearman Rank Correlation Coefficient** ($\rho$) is computed across all 7 Diamond Panel genes and 22 LM22 immune cell types ($M = 154$ hypothesis tests):
$$\rho(g, k) = 1 - \frac{6 \sum_{i=1}^N d_{i, gk}^2}{N(N^2 - 1)}$$
Where $d_{i, gk} = \text{rank}(x_{gi}) - \text{rank}(f_{ik})$ is the difference in ranks between gene expression and deconvolved cell fraction for patient $i$, and $N = 137$.

For sample sizes $N > 100$, statistical significance is evaluated using the asymptotic normal approximation with standard tie corrections:
$$z = \rho \sqrt{\frac{N - 2}{1 - \rho^2}} \sim \mathcal{N}(0, 1)$$
$$P_{\text{nominal}} = 2 \left[ 1 - \Phi(|z|) \right]$$

To control false discoveries across all 154 simultaneous tests, $p$-values are adjusted using the **Benjamini–Hochberg False Discovery Rate (FDR)** procedure:
$$q_{(i)} = \min_{m \ge i} \left( \frac{154 \cdot P_{(m)}}{m} \right)$$
In the correlation heatmap (`04_cibersortx_immune_correlation_heatmap.pdf`), cell types are arranged by descending mean absolute correlation ($\overline{|\rho|}$), with significance asterisks denoting nominal significance: `*` $p < 0.05$, `**` $p < 0.01$, `***` $p < 0.001$.

#### 3.5.2 Median-Split Stratification and Wilcoxon Rank-Sum Testing
To provide intuitive clinical comparisons, patients are dichotomized into **High vs. Low expression groups** based on cohort medians of key prognostic drivers ($PLK1$, $CDK2$):
$$\mathcal{G}_{\text{High}}(g) = \{i \mid x_{gi} \ge \text{median}(x_g)\}, \quad \mathcal{G}_{\text{Low}}(g) = \{i \mid x_{gi} < \text{median}(x_g)\}$$

Differences in immune cell infiltration fractions between High and Low strata are evaluated using the **Two-Sample Wilcoxon Rank-Sum Test (Mann–Whitney $U$ Test)**:
$$U = W_1 - \frac{n_1 (n_1 + 1)}{2}$$
Key immune subsets evaluated in four-panel composite boxplots:
1. **CD8+ Cytotoxic T Lymphocytes** (`T cells CD8`): Primary antitumor effector cells.
2. **M2 Polarized Macrophages** (`Macrophages M2`): Pro-tumoral, immunosuppressive TAMs.
3. **Regulatory T Cells** (`T cells regulatory (Tregs)`): Immunosuppressive T cells that inhibit cytotoxic activity.
4. **M1 Inflammatory Macrophages** (`Macrophages M1`): Pro-inflammatory, antitumor macrophages.

---

### 3.6 Step 6: ESTIMATE Tumor Purity, Stromal, and Immune Confounder Modeling
**Script**: `04_immune_tme_analysis.R`  
**Input Artifacts**: `estimate_scores.tsv` (pre-computed scores via R `estimate` library).  
**Output Artifacts**: `06_estimate_purity_immune_stromal_control.pdf` and `.png` (6-panel confounder control grid).

#### 3.6.1 Orthogonal ESTIMATE Architecture
To prove that observed immune correlations are not artifacts of tumor cellularity, the pipeline employs **ESTIMATE (Estimation of STromal and Immune cells in MAlignant Tumor tissues using Expression data)** (Yoshihara et al., *Nat Commun* 2013 [3]). 

ESTIMATE operates completely independently of CIBERSORTx, using single-sample Gene Set Enrichment Analysis (ssGSEA) across the entire transcriptome:
- **Immune Signature**: $141$ curated genes specific to hematopoietic cells.
- **Stromal Signature**: $141$ curated genes specific to stromal and extracellular matrix fibroblasts.

For each bulk tumor profile, empirical cumulative distribution functions generate normalized enrichment scores:
$$\text{ImmuneScore}_j = \sum_{g \in \mathcal{S}_{\text{immune}}} \Delta(g), \quad \text{StromalScore}_j = \sum_{g \in \mathcal{S}_{\text{stromal}}} \Delta(g)$$
$$\text{ESTIMATEScore}_j = \text{ImmuneScore}_j + \text{StromalScore}_j$$

#### 3.6.2 Non-Linear Cosine Transformation for Tumor Purity
Tumor purity represents the percentage of malignant epithelial cells in the bulk specimen. ESTIMATE computes purity through a non-linear cosine transformation fitted to DNA copy number-derived purity:
$$\text{TumorPurity}_j = \cos\left( 0.6049872010 + 0.0001469884 \cdot \text{ESTIMATEScore}_j \right)$$
The score is bounded within the biological interval $[0, 1]$:
$$\text{TumorPurity}_j = \max\left( 0, \; \min\left( 1, \; \text{TumorPurity}_j \right) \right)$$

#### 3.6.3 Confounder Regression Analysis
To evaluate purity confounding, bivariate linear regression models and Spearman rank correlations are fitted for leading prognostic drivers ($PLK1$, $CDK2$):
$$\text{Metric}_j = \beta_0 + \beta_1 \cdot x_{gj} + \varepsilon_j$$
Where $\text{Metric} \in \{\text{TumorPurity}, \text{ImmuneScore}, \text{StromalScore}\}$.
- A flat slope ($\beta_1 \approx 0, p > 0.05$) between gene expression and TumorPurity proves that elevated biomarker expression is an intrinsic biological property of malignant cells rather than an artifact of high tumor cell content.

---

### 3.7 Step 7: Master Translational Biomarker Synthesis
**Script**: `04_immune_tme_analysis.R`  
**Input Artifacts**: Upstream differential expression metrics, multivariate Cox survival statistics, CIBERSORTx correlation coefficients, ESTIMATE purity correlations, and DGIdb druggability scores.  
**Output Artifacts**: `clinical_relevance_biomarker_table.tsv` (comprehensive translational synthesis table).

#### 3.7.1 Multi-Dimensional Evidence Integration
The final synthesis compiles seven distinct analytical dimensions into a unified, publication-ready matrix (`clinical_relevance_biomarker_table.tsv`):
1. **Target Identification & PIKK Branch**: Gene symbol and upstream kinase branch ($ATM$, $ATR$, $TRRAP$, $mTOR$).
2. **Transcriptional Dysregulation**: EdgeR Quasi-Likelihood $\log_2\text{FC}$, FDR, and directionality in OSCC.
3. **Clinical Survival Relevance**: Kaplan–Meier log-rank $p$-value, multivariable adjusted hazard ratio ($\text{adj.HR}$), multivariable $p$-value, and independence status (*Independent*, *Confounded by Stage*, *Borderline*, *Silenced Suppressor*).
4. **Oncogenic Mechanism**: Explicit molecular role in replication stress tolerance, mitotic override, or chromatin remodeling.
5. **Therapeutic Druggability**: DGIdb composite druggability score, clinical maturity tier, total drug count, approved drug count, and approved/candidate drug manifests.
6. **Immune TME Impact**: Deconvolved Spearman correlations ($\rho$) and BH-adjusted $q$-values for cytotoxic CD8+ T cells, immunosuppressive M2 macrophages, and regulatory T cells ($T_{\text{regs}}$).
7. **Cellularity Independence**: Spearman correlation ($\rho$) and $p$-value against ESTIMATE TumorPurity.

---

## 4. Synthesis & Comparative Matrix

The following structured matrix compares the computational methods, sample sizes, parameter settings, and primary deliverables across all sub-components of the TME and druggability pipeline:

| Analytical Step | Computational Script | Primary Input Data / Sample Size | Core Method / Mathematical Model | Key Parameters / Filter Thresholds | Primary Output Deliverables |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Pharmacogenomic Ingestion** | `01_pikk_druggability_analysis.py` | 24 Hub genes; 7 Diamond Panel genes | DGIdb v5.0 GraphQL REST queries (Two-pass query) | Batch size = 50; retry factor = 1.5; timeout = 20s | `pikk_druggability_composite_scores.tsv`, `dgidb_raw_gene_interactions.tsv` |
| **Composite Druggability Scoring** | `01_pikk_druggability_analysis.py` | Raw interactions, approval status, source databases | Weighted Multi-Signal Model ($S_{\text{comp}} \in [0, 1]$) | Weights: App 0.35, Int 0.30, Src 0.20, Drugs 0.15 | Categorical Maturity Tiers (*FDA Approved*, *Phase I/II*, *Preclinical*, *Novel*) |
| **Druggability Cartography** | `02_pikk_druggability_plots.R` | Composite scores and interaction tables | Stacked bar decomposition; Bipartite mapping | Faceted lollipops; clean drug filtering ($< 30$ chars, no salt) | `01_druggability_composite_score_breakdown`, `02_drug_target_interaction_network_bipartite` |
| **Linear Deconvolution Prep** | `03_prepare_cibersortx_input.R` | $n = 224$ primary tumors; $\log_2\text{CPM}$ matrix | Non-negative linear transform: $y = \max(2^x - 1, 0)$ | Filter `Condition == "Tumor"`; deduplicate genes | `cibersortx_input_matrix.tsv` (224 tumors) |
| **Remote Deconvolution** | Stanford CIBERSORTx Cloud | $n = 224$ linear mixture matrix | $\nu$-Support Vector Regression ($\nu$-SVR) | LM22 signature ($547$ genes); $B = 500$ permutations | `cibersortx_results.csv` (22 cell fractions per tumor) |
| **Deconvolution QC Filtering** | `04_immune_tme_analysis.R` | $n = 224$ raw deconvolution outputs | Empirical permutation hypothesis test | Threshold: $P_{\text{perm}} < 0.05$; $n = 137$ tumors retained | `ciber_qc` internal matrix; `07_cibersortx_rmse_quality_control` |
| **Immune Correlation Modeling** | `04_immune_tme_analysis.R` | $n = 137$ QC tumors; 7 Diamond genes $\times$ 22 cells | Spearman rank correlation ($\rho$) | $M = 154$ simultaneous tests; Benjamini–Hochberg FDR | `pikk_immune_correlation_matrix.tsv`, `04_cibersortx_immune_correlation_heatmap` |
| **Stratified Immune Boxplots** | `04_immune_tme_analysis.R` | $n = 137$ tumors; CD8, M2, Tregs, M1 fractions | Two-sample Wilcoxon rank-sum test | Median-split stratification: High vs. Low expression | `05_stratified_immune_boxplots_plk1_cdk2.pdf/.png` |
| **ESTIMATE Purity Controls** | `04_immune_tme_analysis.R` | $n = 224$ tumors (intersected to $n = 137$) | ssGSEA (141 immune, 141 stromal genes); Cosine Purity | Cosine purity transform; linear regression with LOESS | `estimate_scores.tsv`, `06_estimate_purity_immune_stromal_control` |
| **Master Biomarker Synthesis** | `04_immune_tme_analysis.R` | Integrated outputs from all project phases | Multi-dimensional matrix compilation | Complete synthesis of DGEA, Survival, TME, Purity, Drugs | `clinical_relevance_biomarker_table.tsv` |

---

## 5. Ready-to-Publish Methods Section (Manuscript Reference Text)

### 2.6 Tumor Microenvironment Deconvolution, Immune Infiltration Modeling, and Therapeutic Druggability Profiling

#### 2.6.1 Bulk Transcriptomic Deconvolution via CIBERSORTx
To characterize the tumor immune microenvironment (TME) of oral cavity squamous cell carcinomas and evaluate its remodeling by PIKK pathway drivers, high-resolution cellular deconvolution was conducted using CIBERSORTx [1, 2]. Prior to deconvolution, normalized bulk RNA sequencing expression values ($\log_2(\text{CPM} + 1)$) for primary solid tumors ($n = 224$) were converted to non-negative linear space ($y = \max(2^x - 1, 0)$) to adhere strictly to linear mixture modeling assumptions. Deconvolution was executed via the Stanford CIBERSORTx high-performance cloud infrastructure utilizing the validated LM22 leukocyte gene signature matrix ($547$ signature genes resolving $22$ distinct human hematopoietic phenotypes) [1]. 

The relative abundance of each immune cell lineage was quantified via $\nu$-support vector regression ($\nu$-SVR) with disabled quantile normalization. Statistical significance of the deconvolution fit was determined through an empirical permutation test ($B = 500$ permutations). To eliminate computational noise from poorly fitting samples, a rigorous quality control threshold of $P_{\text{perm}} < 0.05$ was enforced, retaining $n = 137$ primary oral tumors for downstream correlation and stratification analyses. Deconvolution fit accuracy was verified by assessing root mean squared error (RMSE) distributions against biomarker expression.

#### 2.6.2 Correlation Analysis with Immune Subsets and Stratified Infiltration Modeling
Non-parametric Spearman rank correlation coefficients ($\rho$) were calculated to assess monotonic associations between the expression of the 7 core PIKK Diamond Panel targets ($PLK1$, $CDK2$, $TOPBP1$, $RAD51$, $FANCI$, $KAT2B$, $DEPTOR$) and the estimated fractions of all $22$ LM22 immune cell types ($154$ simultaneous hypothesis tests). Statistical significance was evaluated using the asymptotic normal approximation with standard tie corrections, and nominal $p$-values were adjusted for multiple testing across the hypothesis space using the Benjamini–Hochberg False Discovery Rate (FDR) procedure [5]. 

To evaluate discrete clinical shifts in key immune compartments, tumors were dichotomized into High and Low expression groups based on the cohort median for leading independent prognostic drivers ($PLK1$, $CDK2$). Inter-group differences in cytotoxic CD8+ T lymphocytes, M2-polarized macrophages, regulatory T cells ($T_{\text{regs}}$), and M1 macrophages were evaluated using the two-sample Wilcoxon rank-sum test (Mann–Whitney $U$ test).

#### 2.6.3 Orthogonal Estimation of Tumor Purity and Stromal Infiltration (ESTIMATE)
To confirm that observed immune correlations reflected true biological immunomodulation rather than technical confounding by tumor cellularity, an orthogonal analysis was conducted using the ESTIMATE algorithm in R [3]. Single-sample Gene Set Enrichment Analysis (ssGSEA) was performed across the full transcriptome using two curated signatures ($141$ immune genes and $141$ stromal genes) to derive normalized ImmuneScore and StromalScore values. Tumor purity was calculated using the validated non-linear cosine transformation:
$$\text{TumorPurity} = \cos(0.6049872010 + 0.0001469884 \cdot \text{ESTIMATEScore})$$
bounded within the interval $[0, 1]$. Linear regression models and Spearman rank correlation tests were performed between hub gene expression and ESTIMATE TumorPurity to ensure that biomarker associations were independent of cancer cell density.

#### 2.6.4 Pharmacogenomic Interrogation and Composite Druggability Scoring
To evaluate the therapeutic actionability of prioritized hub genes, comprehensive drug-gene interaction profiles were obtained programmatically from the Drug-Gene Interaction Database (DGIdb v5.0) via its GraphQL API [4]. A two-pass querying architecture was deployed: the first pass retrieved all physical and functional drug associations, interaction confidence scores, mechanisms of action, and database sources; the second pass fetched curated regulatory approval statuses (FDA approved) and clinical disease indications across $15$ curated oncology and pharmacology databases.

To prioritize targets based on translational maturity, a continuous Composite Druggability Score ($S_{\text{comp}} \in [0.0, 1.0]$) was formulated as a weighted linear combination of four normalized evidentiary signals:
$$S_{\text{comp}} = 0.35 \cdot S_{\text{app}} + 0.30 \cdot S_{\text{int}} + 0.20 \cdot S_{\text{src}} + 0.15 \cdot S_{\text{drugs}}$$
Where $S_{\text{app}}$ represents the regulatory approval tier (1.0 for FDA approved; 0.90 for clinical guidelines; 0.75 for Phase III; 0.60 for Phase II; 0.45 for Phase I; 0.15 for preclinical; 0.0 for novel), $S_{\text{int}}$ is the log-normalized maximum interaction confidence score (capped at $100$), $S_{\text{src}}$ is the proportion of matched source databases out of $15$ known repositories, and $S_{\text{drugs}}$ is the log-normalized count of unique interacting molecules (capped at $20$). Targets were assigned to categorical maturity tiers (*FDA Approved*, *Phase III Clinical*, *Phase I/II Clinical*, *Preclinical / Experimental*, or *Undruggable / Novel*). Bipartite target-drug interaction maps were generated to visualize repurposing candidates and approved antineoplastics.

---

## 6. Comprehensive Software, Dependency, and Environment Manifest

All TME deconvolution preparations, immune modeling, ESTIMATE purity controls, druggability scoring, and publication-quality figure generation were executed across Python (v3.10+) and R (v4.2.0+) environments:

| Software / Package | Language / Source | Version Tested | Primary Computational Purpose in Pipeline |
| :--- | :--- | :--- | :--- |
| **`requests`** | Python / PyPI | v2.31.0 | Automated HTTP communication, GraphQL query dispatching, and exponential backoff retry management for DGIdb v5.0 API queries. |
| **`pandas`** | Python / PyPI | v2.1.4 | Tabular data manipulation, GraphQL JSON parsing, drug salt cleaning, and composite score computation. |
| **`urllib3`** | Python / PyPI | v2.1.0 | Connection pooling, retry backoff configuration, and timeout enforcement for pharmacogenomic web scraping. |
| **`estimate`** | R / R-Forge | v1.0.13 | Orthogonal tumor microenvironment modeling: calculation of ImmuneScore, StromalScore, ESTIMATEScore, and cosine-transformed TumorPurity [3]. |
| **`ggplot2`** | R / CRAN | v3.4.4 | Publication graphics: stacked druggability bar charts, bipartite faceted lollipop plots, and multi-panel scatter grids. |
| **`patchwork`** | R / CRAN | v1.2.0 | Multi-panel figure composition, alignment, and hierarchical annotation for composite boxplot and scatter assemblies. |
| **`RColorBrewer`** | R / CRAN | v1.1-3 | Diverging and sequential color palettes for correlation heatmaps and regulatory status indicators. |
| **`scales`** | R / CRAN | v1.3.0 | Coordinate scaling, scientific axis formatting, and out-of-bounds (`oob = squish`) aesthetic handling. |
| **`dplyr`** | R / CRAN | v1.1.4 | Tidy manipulation, filtering of CIBERSORTx permutations ($P_{\text{perm}} < 0.05$), and median-split cohort assignment. |
| **`tidyr`** | R / CRAN | v1.3.1 | Data restructuring (`pivot_longer`, `pivot_wider`) for correlation matrices and stacked signal breakdowns. |
| **`readr`** | R / CRAN | v2.1.5 | High-speed, type-safe serialization of TSV and CSV input/output matrices. |
| **`tibble`** | R / CRAN | v3.2.1 | Column-to-row conversions and strict data frame handling. |

---

## 7. Peer-Reviewed Scientific Bibliography

1. **Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., Hoang, C. D., Diehn, M., & Alizadeh, A. A.** (2015). Robust enumeration of cell subsets from tissue expression profiles. *Nature Methods*, 12(5), 453–457. [DOI: 10.1038/nmeth.3337](https://doi.org/10.1038/nmeth.3337)
2. **Newman, A. M., Steen, C. B., Liu, C. L., Gentles, A. J., Chaudhuri, A. A., Scherer, F., Khodadoust, M. S., Esfahani, M. S., Luca, B. A., Steiner, D., Diehn, M., & Alizadeh, A. A.** (2019). Determining cell type abundance and expression from bulk tissues with digital cytometry. *Nature Biotechnology*, 37(7), 773–782. [DOI: 10.1038/s41587-019-0114-2](https://doi.org/10.1038/s41587-019-0114-2)
3. **Yoshihara, K., Shahmoradgoli, M., Martínez, E., Vegesna, R., Kim, H., Torres-Garcia, W., Treviño, V., Shen, H., Laird, P. W., Levine, D. A., Carter, S. L., Getz, G., Stemke-Hale, K., Mills, G. B., & Verhaak, R. G.** (2013). Inferring tumour purity and stromal and immune cell admixture from expression data. *Nature Communications*, 4, 2612. [DOI: 10.1038/ncomms3612](https://doi.org/10.1038/ncomms3612)
4. **Freshour, S. L., Kiwala, S., Cotto, K. C., Coffman, A. C., McMichael, J. F., Song, J. J., Griffith, M., & Griffith, O. L.** (2021). Integration of the Drug–Gene Interaction Database (DGIdb 4.0) with open crowdsource efforts. *Nucleic Acids Research*, 49(D1), D1144–D1151. [DOI: 10.1093/nar/gkaa1084](https://doi.org/10.1093/nar/gkaa1084)
5. **Benjamini, Y., & Hochberg, Y.** (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B (Methodological)*, 57(1), 289–300. [DOI: 10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
6. **Zar, J. H.** (2010). *Biostatistical Analysis* (5th ed.). Pearson Prentice Hall. ISBN: 978-0131008465.
7. **Chen, B., Khodadoust, M. S., Liu, C. L., Newman, A. M., & Alizadeh, A. A.** (2018). Profiling tumor infiltrating immune cells with CIBERSORT. In *Cancer Systems Biology: Methods and Protocols* (pp. 243–259). Springer. [DOI: 10.1007/978-1-4939-7493-1_12](https://doi.org/10.1007/978-1-4939-7493-1_12)
8. **Cancer Genome Atlas Network.** (2015). Comprehensive genomic characterization of head and neck squamous cell carcinomas. *Nature*, 517(7536), 576–582. [DOI: 10.1038/nature14129](https://doi.org/10.1038/nature14129)
9. **Kassambara, A.** (2020). *ggpubr: 'ggplot2' Based Publication Ready Plots*. R package version 0.4.0. [CRAN: ggpubr](https://CRAN.R-project.org/package=ggpubr)
