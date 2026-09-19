# Survival Association Analysis of PIKK-Associated Hub Genes (Section 2.5) — Comprehensive Methodology Guide

**Project**: PIKK Pathway-Associated Biomarker Discovery and Oncogenic Characterization in Oral Squamous Cell Carcinoma (OSCC)  
**Target Manuscript Reference**: Materials and Methods / Bioinformatic Processing Pipeline  
**Primary Dataset**: The Cancer Genome Atlas Head and Neck Squamous Cell Carcinoma (TCGA-HNSC), Curated Oral Cavity Squamous Cell Carcinoma Cohort ($n = 240$ primary biospecimens; $n = 224$ Primary Solid Tumor, $n = 16$ Solid Tissue Normal)  
**Clinical Survival Cohort**: $n = 222$ unique primary oral carcinoma patients with complete vital status and verified positive follow-up duration ($n = 98$ overall mortality events, $n = 124$ right-censored observations)  
**Multivariable Adjusted Cohort**: $n = 192$ patients ($n = 87$ mortality events, $n = 105$ right-censored observations) possessing complete multi-covariate clinical profiles (Patient Age, AJCC Pathologic Stage, Histologic Grade, and Gender)  
**Investigated Gene Space**: 24 prioritized PIKK-associated topological hub genes ($CDK2$, $PLK1$, $RAD51$, $TOPBP1$, $EGFR$, $BRCA1$, $BRCA2$, $PARP1$, $H2AX$, $RUVBL1$, $FANCI$, $EXO1$, $POLD1$, $KAT2A$, $MSH6$, $CHEK2$, $CDC45$, $PRKDC$, $TTI1$, $KAT2B$, $RPTOR$, $E2F1$, $RICTOR$, $CHEK1$) derived from multi-centrality network prioritization  
**Scripts Evaluated**:
1. `survival_analysis_hub_genes.R` — *Prognostic Association Pipeline: Kaplan–Meier Survival Stratification, Maximally Selected Rank Statistics Cutpoint Optimization, Univariate and Multivariable Cox Proportional Hazards Modeling, Schoenfeld Residual Assumption Verification, and Harrell's Concordance Metric Evaluation*

---

## 1. Executive Summary & End-to-End Workflow

In cancer systems biology, identifying topologically central regulatory nodes ("hub genes") within protein-protein interaction networks reveals the structural and functional backbones of malignant transformation. However, topological centrality does not automatically translate into clinical relevance. To establish translational and therapeutic utility, prioritized biomarkers must be interrogated against clinical longitudinal endpoints to determine whether their transcriptomic expression profiles reliably stratify patient survival and independently predict mortality risk beyond standard clinicopathologic staging.

This document establishes the publication-grade computational and statistical methodology for the **Survival Association and Clinical Prognostic Modeling pipeline (Section 2.5)** of the OSCC PIKK project.

The analytical pipeline evaluates candidate hub genes across complementary non-parametric and semi-parametric survival frameworks:
1. **Clinical Cohort Ingestion & Survival Variable Construction**: Curates overall survival endpoints ($T_{\text{OS}}$, event indicator $\delta$) from primary clinical metadata, deduplicates multi-sample patient records, and standardizes clinicopathologic covariates into ordered and binary risk tiers.
2. **Optimal Cutpoint Stratification & Kaplan–Meier Estimation**: Applies **maximally selected rank statistics** (`survminer::surv_cutpoint`) to identify biologically driven, data-adaptive expression cutpoints with minimum subgroup proportion constraints (`minprop = 0.25`), avoiding arbitrary median binning. Fits non-parametric product-limit survival curves and evaluates inter-group survival differences via two-sample log-rank tests.
3. **Continuous Univariate Cox Proportional Hazards Regression**: Quantifies the unadjusted hazard ratio (HR) per unit increase in normalized $\log_2\text{CPM}$ expression, preserving the continuous expression distribution to prevent information loss.
4. **Multivariable Cox Proportional Hazards Modeling**: Fits simultaneous semi-parametric models adjusting for established clinical prognostic determinants (Age at diagnosis, AJCC pathologic tumor stage, histologic grade, and patient gender), testing whether hub genes confer independent prognostic value.
5. **Proportional Hazards Assumption Verification**: Interrogates the time-invariance of regression coefficients using scaled Schoenfeld residuals and Grambsch–Therneau global and covariate-specific diagnostic tests (`cox.zph`).
6. **Model Discrimination & Concordance Quantification**: Quantifies predictive performance for all univariable and multivariable survival models using Harrell's Concordance Index ($C$-index) to measure improvements in risk discrimination.
7. **Multiple Testing Correction & Multi-Tiered Significance Synthesis**: Controls the family-wise False Discovery Rate (FDR) across all evaluated genes using the Benjamini–Hochberg procedure, synthesizing findings into unified prognostic classifications.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                                 UPSTREAM INPUT ARTIFACTS                               │
│  • Prioritized Hub Gene Manifest: hub_gene_characterisation_table.tsv (24 hub genes)   │
│  • Normalized Transcript Abundance: logCPM_matrix.tsv (Tumor = 224, Normal = 16)       │
│  • Matched Clinical Metadata: matched_metadata_used.tsv (Vital status, days, stage)   │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 1: CLINICAL COHORT HARMONIZATION & ENDPOINT ASSEMBLY (survival_analysis_hub_genes)│
│  • Filter to primary solid oral tumors (Sample_Type_Code == 1; n = 224)               │
│  • Deduplicate to unique patient cases (distinct `Case ID`)                            │
│  • Construct Overall Survival (OS) time in years:                                      │
│      - Dead: days_to_death / 365.25 (event = 1)                                        │
│      - Alive: days_to_last_follow_up / 365.25 (event = 0, right-censored)              │
│  • Exclude non-positive or missing follow-up durations (Survival Cohort: n = 222)     │
│  • Standardize clinical covariates:                                                    │
│      - Age at diagnosis (continuous, years)                                            │
│      - Pathologic Stage: Early (I–II) vs Advanced (III–IV)                             │
│      - Histologic Grade: Low Grade (G1–G2) vs High Grade (G3–G4)                       │
│      - Gender: Male vs Female                                                          │
│  • Multi-covariate complete case cohort assembly (Multivariate Cohort: n = 192)       │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 2: NON-PARAMETRIC STRATIFICATION & OPTIMAL CUTPOINT SELECTION                   │
│  • Maximally selected rank statistics (surv_cutpoint, minprop = 0.25)                  │
│  • Maximize standardized log-rank statistic across expression grid:                    │
│      M = max |T(μ)| subject to 0.25 ≤ n_High / N ≤ 0.75                                │
│  • Dichotomize patients into High vs Low expression groups per hub gene                │
│  • Compute Kaplan–Meier survival functions: S(t) = ∏ (1 - d_i / Y_i)                   │
│  • Hypothesis testing via Two-Sample Log-Rank Test (survdiff, df = 1)                  │
│  • Multiple testing correction via Benjamini–Hochberg FDR                              │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 3: CONTINUOUS UNIVARIATE COX PROPORTIONAL HAZARDS MODELING                      │
│  • Model hazard rate as continuous function: h(t | x) = h_0(t) exp(β * x_gene)        │
│  • Maximize Cox partial likelihood over failure times                                  │
│  • Estimate unadjusted Hazard Ratios (HR = exp(β)) and 95% Wald confidence intervals   │
│  • Compute Wald test statistics: z = β / SE(β), p-value from N(0, 1)                   │
│  • Adjust for multiple hypothesis testing across 24 hubs (BH method)                   │
│  • Visualize effect sizes via high-resolution univariate forest plots                  │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 4: MULTIVARIABLE COX REGRESSION & CLINICAL COVARIATE ADJUSTMENT                 │
│  • Model: h(t | Z) = h_0(t) exp(β_gene * x_gene + β_age * Age + β_stage * Stage +     │
│                                 β_grade * Grade + β_gender * Gender)                   │
│  • Evaluate independent prognostic significance (adj.HR = exp(β_gene))                 │
│  • Assess omnibus model fit via Likelihood Ratio Test (Score logtest p-value)          │
│  • Generate multi-parameter full forest plots for top-performing biomarkers            │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 5: MODEL DIAGNOSTICS & DISCRIMINATIVE PERFORMANCE EVALUATION                    │
│  • Test Proportional Hazards assumption via Scaled Schoenfeld Residuals:               │
│      E[r*_ik] ≈ β_k(t_i) - β_k (Grambsch–Therneau correlation test, cox.zph)           │
│  • Generate diagnostic plots of β(t) over time with LOESS smoothing curves             │
│  • Quantify discriminative ability via Harrell's Concordance Index (C-index ± SE):     │
│      C = P(Predicted Risk_i > Predicted Risk_j | T_i < T_j, event_i = 1)               │
│  • Compare discriminative gain: C_index(Multivariate) vs C_index(Univariate)           │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 6: MULTI-TIERED SURVIVAL SYNTHESIS & TRANSLATIONAL REPORTING                    │
│  • Assemble master prognostic integration table (survival_summary_table.tsv)           │
│  • Multi-tiered classification: Univariate vs Multivariate vs KM-only prognostic hubs  │
│  • Multi-test significance visualization: -log10(p) summary heatmap across KM,        │
│    Univariate Cox, and Multivariate Cox models                                         │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. Biological Context and Study Rationale

### 2.1 The Genomic Landscape of OSCC and the Necessity of Biomarker Stratification
Oral Squamous Cell Carcinoma (OSCC) represents the predominant malignancy of the oral cavity, displaying a notoriously poor 5-year overall survival rate (~50%) that has remained stagnant over recent decades despite advances in multimodal surgical and chemoradiotherapeutic regimens. Pathologically and genomically, HPV-negative OSCC is characterized by near-universal inactivation of $TP53$ (>80%), frequent homozygous deletion or promoter hypermethylation of $CDK2NA$ ($p16^{\text{INK4A}}$), and high copy number alterations. 

The loss of the G1/S checkpoint forces malignant oral epithelial cells to replicate with damaged templates, generating intense oncogene-induced **replication stress** and persistent DNA double-strand breaks (DSBs). To avert catastrophic mitotic failure, these cancer cells become addicted to the sensory and signaling networks governed by the **PIKK kinase superfamily** ($ATM$, $ATR$, $PRKDC$, $MTOR$, $TRRAP$, $SMG1$). 

Key downstream interactors and cell-cycle effectors—such as $PLK1$, $CDK2$, $TOPBP1$, $RAD51$, and $BRCA1$—serve as indispensable linchpins that facilitate replication fork stabilization, restart collapsed forks, and enforce the intra-S and G2/M checkpoints. If the transcriptional hyper-activation of these repair and checkpoint hubs enables tumors to withstand endogenous genotoxic stress and resist DNA-damaging therapeutics (e.g., cisplatin, ionizing radiation), high expression of these genes should directly correlate with aggressive clinical behavior, treatment failure, and accelerated patient mortality.

### 2.2 Methodological Rationale: Justification of Analytical Paradigms

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                    METHODOLOGICAL COMPARISON: SURVIVAL MODELING CHOICES                │
├───────────────────────────┬──────────────────────────────┬─────────────────────────────┤
│ Dimension                 │ Selected Pipeline Paradigm   │ Conventional Alternative    │
├───────────────────────────┼──────────────────────────────┼─────────────────────────────┤
│ Cutpoint Strategy         │ Maximally Selected Rank      │ Arbitrary Median Split      │
│                           │ Statistics (surv_cutpoint)   │ (50% / 50% fixed cut)       │
│ Rationale                 │ Data-adaptive, identifies    │ Arbitrary; masks non-linear │
│                           │ true biological threshold;   │ biological effects and can  │
│                           │ bounded by minprop = 0.25    │ split homogeneous clusters  │
├───────────────────────────┼──────────────────────────────┼─────────────────────────────┤
│ Variable Modeling         │ Dual-Engine: Continuous Cox  │ Categorical-only or         │
│                           │ + Dichotomous Kaplan–Meier   │ Continuous-only modeling    │
│ Rationale                 │ Preserves full variance and  │ Dichotomization loses power;│
│                           │ statistical power (Cox), while│ continuous modeling lacks   │
│                           │ providing intuitive clinical │ intuitive risk curves and   │
│                           │ risk stratifications (KM)    │ median survival times       │
├───────────────────────────┼──────────────────────────────┼─────────────────────────────┤
│ Clinical Confounding      │ Multivariable Adjustment for │ Univariate Analysis Only    │
│                           │ Stage, Grade, Age, Gender    │ (unadjusted biomarker calls)│
│ Rationale                 │ Demonstrates true independent│ Susceptible to confounding; │
│                           │ biomarker utility beyond     │ biomarker may simply reflect│
│                           │ standard anatomic staging    │ advanced clinical stage     │
├───────────────────────────┼──────────────────────────────┼─────────────────────────────┤
│ Model Verification        │ Scaled Schoenfeld Residual   │ Unchecked Proportional      │
│                           │ Diagnostics (cox.zph)        │ Hazards Assumption          │
│ Rationale                 │ Mathematically proves hazard │ Violated assumptions yield  │
│                           │ ratios remain constant over  │ misleading, time-biased, or │
│                           │ follow-up time               │ invalid effect estimates    │
├───────────────────────────┼──────────────────────────────┼─────────────────────────────┤
│ Discriminative Power      │ Harrell's Concordance Index  │ P-value Significance Alone  │
│                           │ (C-index ± asymptotic SE)    │                             │
│ Rationale                 │ Quantifies actual predictive │ Statistical significance    │
│                           │ discrimination and gain over │ does not equate to clinical │
│                           │ clinical staging baselines   │ predictive accuracy         │
└───────────────────────────┴──────────────────────────────┴─────────────────────────────┘
```

#### 1. Why Maximally Selected Rank Statistics Over Median Splitting?
Conventional survival analyses in oncology routinely dichotomize continuous gene expression by the cohort median. While superficially simple, median dichotomization suffers from major methodological flaws:
- It assumes that the true biological threshold separating high-risk from low-risk tumors always lies at the 50th percentile, an assumption with no physiological basis.
- In bimodally or right-skewed distributed oncogenes (e.g., $PLK1$, $EGFR$), a median split divides biologically identical patients into opposing groups, attenuating statistical power.
- The pipeline deploys **maximally selected rank statistics** (`survminer::surv_cutpoint`), which rigorously scans the candidate expression continuum to locate the cutpoint that maximizes the standardized two-sample log-rank statistic.
- To prevent statistical instability caused by extreme outliers (e.g., isolating 2 patients against 220), the algorithm enforces an explicit minimum proportion constraint:
  $$\min(n_{\text{Low}}, n_{\text{High}}) \ge 0.25 \cdot N$$
  This guarantees that both risk strata maintain at least 25% of the total cohort ($n \ge 55$ patients), ensuring robust survival curves and balanced degrees of freedom.

#### 2. Why Continuous Cox Modeling Alongside Kaplan–Meier?
Dichotomizing continuous biological variables inherently causes information loss, converting rich quantitative expression gradients into coarse binary states. Consequently, the pipeline operates a **dual-engine survival framework**:
- **Continuous Univariate Cox Regression**: Retains the complete, uncompressed expression gradient ($\log_2\text{CPM}$), evaluating the incremental increase in hazard per unit change in expression.
- **Categorical Kaplan–Meier Product-Limit Estimation**: Translates the continuous risk into actionable, discretized clinical risk cohorts (High vs. Low expression), yielding tangible clinical metrics such as median survival time and survival probability at discrete yearly milestones.

#### 3. Why Multivariable Adjustment is Essential
In clinical oncology, anatomic staging (AJCC Pathologic Stage) remains the primary determinant of prognosis. If a hub gene is upregulated purely as an epiphenomenon of massive tumor burden or advanced regional lymph node metastasis, it possesses no true autonomous diagnostic value. 

Fitting multivariable Cox models incorporating:
- **Age at Diagnosis** (accounting for age-related mortality and immunosenescence),
- **AJCC Pathologic Stage** (adjusting for primary tumor extent and nodal involvement),
- **Histologic Tumor Grade** (adjusting for cellular differentiation and architectural atypia), and
- **Patient Gender** (adjusting for baseline demographic disparities),
guarantees that any identified prognostic biomarker provides **independent additive information** that cannot be explained away by conventional clinicopathologic staging alone.

---

## 3. Step-by-Step Computational Methodology

### 3.1 Step 1: Clinical Cohort Harmonization & Survival Endpoint Construction
**Script**: `survival_analysis_hub_genes.R`  
**Input Artifacts**:
- `matched_metadata_used.tsv`: Curated TCGA-HNSC clinical metadata table.
- `logCPM_matrix.tsv`: Variance-stabilized, batch-corrected $\log_2\text{CPM}$ expression matrix.
- `hub_gene_characterisation_table.tsv`: Table of 24 topologically prioritized hub genes.
**Output Artifacts**:
- `surv_df`: Internal analysis data frame containing harmonized clinical endpoints and hub gene expression ($n = 222$).
- `surv_df_mv`: Complete-case subset for multivariable regression ($n = 192$).

#### 3.1.1 Clinical Inclusion Criteria and Cohort Deduplication
To ensure survival metrics reflect primary oncogenic biology without confounding from metastatic clones or normal tissue baselines, strict filtration rules are applied:
1. **Histological Filter**: Restrict exclusively to Primary Solid Tumor samples (`Condition == "Tumor"` and `Sample_Type_Code == "01"`).
2. **Case Deduplication**: For patients possessing multiple sequenced tumor vials, deduplication is enforced by preserving the first primary tumor vial per patient:
   $$\mathcal{P}_{\text{tumor}} = \text{distinct}(\mathcal{D}_{\text{metadata}}, \text{`Case ID`})$$
   Reducing the raw tumor biospecimens ($n = 224$) to unique individual patients.

#### 3.1.2 Formal Mathematical Definition of the Overall Survival (OS) Endpoint
The primary longitudinal outcome is defined as **Overall Survival (OS)**, parameterized as a right-censored time-to-event variable:
- **Survival Time ($T_i$)**: The elapsed time from initial pathologic diagnosis to either death or last confirmed clinical follow-up, expressed in decimal years:
  $$T_i = \begin{cases} \frac{\text{demographic.days\_to\_death}_i}{365.25} & \text{if } \text{vital\_status}_i = \text{"Dead"} \\[8pt] \frac{\text{diagnoses.days\_to\_last\_follow\_up}_i}{365.25} & \text{if } \text{vital\_status}_i = \text{"Alive"} \end{cases}$$
- **Event Indicator ($\delta_i$)**: A binary failure indicator denoting mortality:
  $$\delta_i = \begin{cases} 1 & \text{if } \text{vital\_status}_i = \text{"Dead"} \\ 0 & \text{if } \text{vital\_status}_i = \text{"Alive" (Right-Censored)} \end{cases}$$
- **Validity Exclusion**: Patients with missing vital status, non-positive follow-up durations ($T_i \le 0$), or missing survival tracking are excluded:
  $$\mathcal{P}_{\text{valid}} = \{ i \in \mathcal{P}_{\text{tumor}} \mid \delta_i \in \{0, 1\} \land T_i > 0 \}$$
  This filtering retains exactly $n = 222$ patients with complete longitudinal tracking ($n = 98$ death events, $n = 124$ right-censored).

#### 3.1.3 Clinical Covariate Harmonization & Classification
To support robust multivariable modeling, clinical covariates are mapped into standardized, clinically relevant factor levels:
1. **Patient Age at Diagnosis ($x_{\text{age}}$)**: Modeled as a continuous variable (years):
   $$x_{\text{age}} = \text{demographic.age\_at\_index}$$
2. **AJCC Pathologic Tumor Stage ($x_{\text{stage}}$)**: Binned into early versus advanced disease:
   $$x_{\text{stage}} = \begin{cases} \text{"Early (I–II)"} & \text{if } \text{stage} \in \{\text{Stage I, Stage IA, Stage IB, Stage II}\} \\ \text{"Advanced (III–IV)"} & \text{if } \text{stage} \in \{\text{Stage III, Stage IVA, Stage IVB, Stage IVC}\} \end{cases}$$
   Reference category: `Early (I–II)`.
3. **Histologic Tumor Grade ($x_{\text{grade}}$)**: Binned into low versus high histologic grade:
   $$x_{\text{grade}} = \begin{cases} \text{"Low Grade (G1–G2)"} & \text{if } \text{grade} \in \{\text{G1, G2}\} \\ \text{"High Grade (G3–G4)"} & \text{if } \text{grade} \in \{\text{G3, G4}\} \end{cases}$$
   Reference category: `Low Grade (G1–G2)`.
4. **Patient Gender ($x_{\text{gender}}$)**: Modeled as a binary demographic factor (`male`, `female`). Reference category: `male`.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        CLINICAL COHORT HARMONIZATION METRIC MANIFEST                   │
├────────────────────────────────────────────────────────┬───────────────────────────────┤
│ Cohort / Processing Dimension                          │ Quantitative Value            │
├────────────────────────────────────────────────────────┼───────────────────────────────┤
│ Total Curated TCGA-HNSC Oral Cavity Tumor Samples      │ N = 224                       │
│ Unique Primary Tumor Patients (Case ID Deduplicated)   │ N = 224                       │
│ Excluded: Zero / Negative Follow-up Duration           │ n = 2 patients                │
│ Verified Overall Survival (OS) Analytical Cohort       │ N = 222 patients              │
│   • Observed Failure Events (Patient Deaths, δ = 1)    │ n = 98 deaths (44.14%)        │
│   • Right-Censored Observations (Alive, δ = 0)         │ n = 124 censored (55.86%)     │
│ Complete-Case Multivariable Regression Cohort          │ N = 192 patients              │
│   • Observed Failure Events in Multivariable Cohort    │ n = 87 deaths (45.31%)        │
│   • Right-Censored Observations in Multivariable Cohort│ n = 105 censored (54.69%)     │
│ Excluded from Multivariable Models (Missing Stage/Grade)│ n = 30 patients               │
│ Candidate Hub Genes Evaluated                          │ G = 24 hub genes              │
└────────────────────────────────────────────────────────┴───────────────────────────────┘
```

---

### 3.2 Step 2: Non-Parametric Survival Function Estimation & Maximally Selected Rank Cutpoint Optimization
**Script**: `survival_analysis_hub_genes.R`  
**Input Artifacts**: `surv_df`, candidate hub gene expression vectors.  
**Output Artifacts**:
- `km_logrank_summary.tsv`: Master table of optimal cutpoints, subgroup sample sizes, log-rank $\chi^2$ statistics, nominal $p$-values, and FDR-adjusted $q$-values.
- Individual Kaplan–Meier plots: `km_plots/KM_{gene}.pdf` and `.png` (300 DPI) with risk tables and confidence intervals.
- `km_combined_grid.pdf` and `.png`: Multi-panel comparative Kaplan–Meier visualization grid.

#### 3.2.1 Maximally Selected Rank Statistics Formulation
Rather than bisecting expression at the sample median, data-driven cutpoint optimization is executed using **maximally selected rank statistics** (Hothorn & Lausen, 2003 [3]; implemented in `survminer::surv_cutpoint`).

Let $x_{gi}$ denote the normalized $\log_2\text{CPM}$ expression of gene $g$ in patient $i$. Let $\mu$ represent a candidate cutpoint along the expression range. Dichotomizing the cohort at threshold $\mu$ partitions patients into two groups:
$$\mathcal{G}_{\text{Low}}(\mu) = \{i \mid x_{gi} < \mu\}, \quad \mathcal{G}_{\text{High}}(\mu) = \{i \mid x_{gi} \ge \mu\}$$

For each candidate $\mu$, the standardized two-sample log-rank statistic $T(\mu)$ is computed:
$$T(\mu) = \frac{\sum_{j=1}^K \left( d_{1j}(\mu) - \frac{Y_{1j}(\mu) d_j}{Y_j} \right)}{\sqrt{\sum_{j=1}^K V_j(\mu)}}$$
Where:
- $K$ is the number of distinct death times observed in the cohort.
- $d_j$ is the total number of deaths occurring at time $t_j$.
- $Y_j$ is the total number of patients at risk immediately prior to time $t_j$.
- $d_{1j}(\mu)$ and $Y_{1j}(\mu)$ are the deaths and risk counts within $\mathcal{G}_{\text{High}}(\mu)$.
- $V_j(\mu)$ is the hypergeometric variance under the null hypothesis of no survival difference:
  $$V_j(\mu) = \frac{Y_{1j}(\mu) Y_{0j}(\mu) d_j (Y_j - d_j)}{Y_j^2 (Y_j - 1)}$$

The optimal cutpoint $\mu^*$ is selected as the value that maximizes the absolute standardized log-rank statistic:
$$\mu^* = \arg\max_{\mu} |T(\mu)|$$

#### 3.2.2 Minimum Proportion Search Constraint
To prevent overfitting to small boundary cohorts, the search space for $\mu$ is bounded by an explicit minimum proportion constraint (`minprop = 0.25`):
$$0.25 \le \frac{|\mathcal{G}_{\text{High}}(\mu)|}{N} \le 0.75$$
This mathematical constraint ensures that:
$$n_{\text{High}} \ge 55 \quad \text{and} \quad n_{\text{Low}} \ge 55 \quad (\text{for } N = 222)$$
If numerical convergence for optimal cutpoint estimation fails, the pipeline automatically falls back to the median expression value:
$$\mu_{\text{fallback}} = \text{median}(x_g)$$

#### 3.2.3 Kaplan–Meier Product-Limit Estimator
For both stratified groups ($\mathcal{G} \in \{\text{Low}, \text{High}\}$), the continuous survival probability function $S(t) = \mathbb{P}(T > t)$ is estimated non-parametrically using the **Kaplan–Meier product-limit estimator** (Kaplan & Meier, 1958 [1]):
$$\hat{S}(t) = \prod_{t_i \le t} \left( 1 - \frac{d_i}{Y_i} \right)$$
Where $t_i$ ranges over all distinct failure times up to time $t$.

The standard error of the estimated survival probability is derived using **Greenwood's formula**:
$$\widehat{\text{SE}}(\hat{S}(t)) = \hat{S}(t) \sqrt{ \sum_{t_i \le t} \frac{d_i}{Y_i (Y_i - d_i)} }$$
Log-log transformed 95% confidence intervals are computed as:
$$\text{CI}_{0.95}(\hat{S}(t)) = \hat{S}(t)^{\exp\left( \pm 1.96 \cdot \frac{\widehat{\text{SE}}(\hat{S}(t))}{\hat{S}(t) \ln(\hat{S}(t))} \right)}$$

#### 3.2.4 Two-Sample Log-Rank Hypothesis Testing
The global null hypothesis of equivalent survival curves between High and Low expression groups:
$$H_0: S_{\text{High}}(t) = S_{\text{Low}}(t) \quad \forall t \ge 0$$
is tested using the **Mantel–Cox Log-Rank Test** via `survival::survdiff`. The log-rank test statistic $\chi_{\text{LR}}^2$ is computed as:
$$\chi_{\text{LR}}^2 = \frac{\left[ \sum_{j=1}^K \left( d_{\text{High}, j} - e_{\text{High}, j} \right) \right]^2}{\sum_{j=1}^K V_j} \sim \chi^2(1)$$
Where $e_{\text{High}, j} = Y_{\text{High}, j} \frac{d_j}{Y_j}$ is the expected number of deaths in the High expression group under $H_0$. 

The exact two-sided nominal $p$-value is evaluated from the central chi-square distribution with $\text{df} = 1$:
$$p_{\text{logrank}} = 1 - F_{\chi^2(1)}(\chi_{\text{LR}}^2)$$

To control false discovery rates across all 24 evaluated hub genes, nominal $p$-values are corrected using the Benjamini–Hochberg procedure:
$$q_{\text{logrank}, (i)} = \min_{k \ge i} \left( \frac{24 \cdot p_{\text{logrank}, (k)}}{k} \right)$$

---

### 3.3 Step 3: Continuous Univariate Cox Proportional Hazards Regression
**Script**: `survival_analysis_hub_genes.R`  
**Input Artifacts**: `surv_df`, candidate hub gene expression vectors.  
**Output Artifacts**:
- `univariate_cox_results.tsv`: Complete univariate regression summary reporting coefficients ($\beta$), unadjusted hazard ratios (HR), 95% confidence intervals, standard errors, Wald $z$-scores, nominal $p$-values, BH-adjusted $q$-values, and concordance indices.
- `univariate_forest_plot.pdf` and `.png`: High-resolution forest plot illustrating univariate hazard ratios, 95% CIs, and significance annotations.

#### 3.3.1 Semi-Parametric Cox Proportional Hazards Formulation
To assess the intrinsic prognostic effect of hub gene expression without the arbitrary thresholding inherent to categorization, a continuous semi-parametric **Cox Proportional Hazards Model** (Cox, 1972 [2]) is fitted for each gene:
$$h(t \mid x_{gi}) = h_0(t) \exp(\beta_g \cdot x_{gi})$$
Where:
- $h(t \mid x_{gi})$ is the hazard rate of death at time $t$ for patient $i$ with normalized expression $x_{gi}$.
- $h_0(t)$ is the non-parametric baseline hazard function representing the underlying failure rate when $x_{gi} = 0$.
- $\beta_g$ is the regression coefficient representing the change in log-hazard per 1-unit increase in $\log_2\text{CPM}$.

#### 3.3.2 Partial Likelihood Estimation
Because the baseline hazard $h_0(t)$ is left completely unspecified, $\beta_g$ is estimated by maximizing the **Cox Partial Likelihood** over all $m$ observed distinct death events:
$$L(\beta_g) = \prod_{i=1}^m \frac{\exp(\beta_g \cdot x_{g(i)})}{\sum_{j \in \mathcal{R}(t_{(i)})} \exp(\beta_g \cdot x_{gj})}$$
Where:
- $t_{(1)} < t_{(2)} < \dots < t_{(m)}$ denote the ordered failure times.
- $x_{g(i)}$ is the expression level of the patient who experienced the event at time $t_{(i)}$.
- $\mathcal{R}(t_{(i)}) = \{j \mid T_j \ge t_{(i)}\}$ represents the set of all individuals at risk immediately prior to $t_{(i)}$.

Taking the natural logarithm yields the log-partial likelihood:
$$\ell(\beta_g) = \ln L(\beta_g) = \sum_{i=1}^m \left[ \beta_g \cdot x_{g(i)} - \ln\left( \sum_{j \in \mathcal{R}(t_{(i)})} \exp(\beta_g \cdot x_{gj}) \right) \right]$$
The maximum partial likelihood estimate $\hat{\beta}_g$ is computed iteratively via Newton–Raphson optimization.

#### 3.3.3 Hazard Ratio and Wald Confidence Intervals
The unadjusted **Hazard Ratio (HR)** reflects the relative increase in the instantaneous risk of death for every doubling of transcript abundance (per unit increase in $\log_2\text{CPM}$):
$$\text{HR}_g = \exp(\hat{\beta}_g)$$
The standard error $\text{SE}(\hat{\beta}_g)$ is obtained from the observed Fisher information matrix:
$$\text{SE}(\hat{\beta}_g) = \left( -\frac{\partial^2 \ell}{\partial \beta_g^2} \right)^{-1/2}$$

The 95% Wald confidence interval for the hazard ratio is computed as:
$$\text{CI}_{0.95}(\text{HR}_g) = \left[ \exp\left( \hat{\beta}_g - 1.96 \cdot \text{SE}(\hat{\beta}_g) \right), \; \exp\left( \hat{\beta}_g + 1.96 \cdot \text{SE}(\hat{\beta}_g) \right) \right]$$

The statistical significance of the association is evaluated using the **Wald test**:
$$z_g = \frac{\hat{\beta}_g}{\text{SE}(\hat{\beta}_g)} \sim \mathcal{N}(0, 1)$$
$$p_{\text{Wald}, g} = 2 \left[ 1 - \Phi(|z_g|) \right]$$
Where $\Phi(\cdot)$ denotes the standard normal cumulative distribution function. False discovery rates across the 24 hub genes are controlled via the Benjamini–Hochberg method ($q_{\text{uni}, g}$).

---

### 3.4 Step 4: Multivariable Cox Proportional Hazards Regression & Clinical Covariate Adjustment
**Script**: `survival_analysis_hub_genes.R`  
**Input Artifacts**: `surv_df_mv` ($n = 192$ patients with complete covariates).  
**Output Artifacts**:
- `multivariate_cox_results.tsv`: Master multivariable results reporting adjusted hazard ratios ($\text{adj.HR}$), 95% CIs, covariate-specific coefficients, and model omnibus test $p$-values.
- `multivariate_forest_plot.pdf` and `.png`: Comparative forest plot depicting adjusted biomarker hazard ratios alongside clinical covariate controls.
- `multivariate_full_forest_best_gene.pdf` and `.png`: Comprehensive covariate forest plot for the leading prognostic hub gene generated via `forestmodel`.

#### 3.4.1 Full Multivariable Proportional Hazards Formulation
To isolate the independent prognostic contribution of each hub gene from established clinical risk factors, a distinct multivariable Cox model is constructed for each candidate biomarker across the complete-case cohort ($n = 192$):
$$h(t \mid \mathbf{Z}_i) = h_0(t) \exp\left( \beta_{\text{gene}} x_{gi} + \beta_{\text{age}} z_{i, \text{age}} + \beta_{\text{stage}} z_{i, \text{stage}} + \beta_{\text{grade}} z_{i, \text{grade}} + \beta_{\text{gender}} z_{i, \text{gender}} \right)$$
Where the covariate vector $\mathbf{Z}_i = (x_{gi}, z_{i, \text{age}}, z_{i, \text{stage}}, z_{i, \text{grade}}, z_{i, \text{gender}})^T$ contains:
- $x_{gi}$: Continuous $\log_2\text{CPM}$ expression of hub gene $g$.
- $z_{i, \text{age}}$: Continuous patient age at diagnosis (in years).
- $z_{i, \text{stage}}$: Binary indicator for advanced anatomic stage ($1 = \text{Advanced [Stages III–IV]}, 0 = \text{Early [Stages I–II]}$).
- $z_{i, \text{grade}}$: Binary indicator for high histologic grade ($1 = \text{High Grade [G3–G4]}, 0 = \text{Low Grade [G1–G2]}$).
- $z_{i, \text{gender}}$: Binary indicator for patient sex ($1 = \text{Female}, 0 = \text{Male}$).

#### 3.4.2 Adjusted Hazard Ratio (adj.HR) and Covariate Effects
The **adjusted Hazard Ratio ($\text{adj.HR}$)** for the hub gene:
$$\text{adj.HR}_g = \exp(\hat{\beta}_{\text{gene}})$$
quantifies the independent multiplicative change in mortality hazard per 1-unit increase in gene expression, holding age, tumor stage, histologic grade, and gender strictly constant.

The standard error $\text{SE}(\hat{\beta}_{\text{gene}})$ is extracted from the diagonal element of the inverted Hessian matrix:
$$\text{Var}(\hat{\boldsymbol{\beta}}) = \left( -\nabla^2 \ell(\hat{\boldsymbol{\beta}}) \right)^{-1}$$
Adjusted 95% confidence intervals are computed:
$$\text{CI}_{0.95}(\text{adj.HR}_g) = \left[ \exp\left( \hat{\beta}_{\text{gene}} - 1.96 \cdot \text{SE}(\hat{\beta}_{\text{gene}}) \right), \; \exp\left( \hat{\beta}_{\text{gene}} + 1.96 \cdot \text{SE}(\hat{\beta}_{\text{gene}}) \right) \right]$$

#### 3.4.3 Omnibus Model Goodness-of-Fit (Likelihood Ratio Test)
The global explanatory power of each 5-covariate multivariable model is evaluated using the **Likelihood Ratio Test (Score / Log-Rank Test)**:
$$\Lambda = 2 \left[ \ell(\hat{\boldsymbol{\beta}}) - \ell(\mathbf{0}) \right]$$
Where $\ell(\mathbf{0})$ is the log-partial likelihood of the null model without any predictors. Under the null hypothesis that all regression coefficients are simultaneously zero ($H_0: \boldsymbol{\beta} = \mathbf{0}$):
$$\Lambda \sim \chi^2(p), \quad p = 5 \text{ degrees of freedom}$$
A model log-test $p$-value $< 0.05$ confirms that the combined biomarker-clinical model explains significant variance in patient survival.

---

### 3.5 Step 5: Proportional Hazards Assumption Verification via Schoenfeld Residuals
**Script**: `survival_analysis_hub_genes.R`  
**Input Artifacts**: Multivariable Cox model object (`fm_cox`) for leading prognostic candidate.  
**Output Artifacts**:
- `cox_ph_diagnostics.pdf` and `.png`: Diagnostic panel displaying scaled Schoenfeld residuals against follow-up time for all model covariates, overlaid with LOESS smoothing splines and 95% confidence bands.

#### 3.5.1 Mathematical Rationale of the Proportional Hazards Assumption
The fundamental mathematical assumption of the Cox model is that the hazard ratio between any two individuals is **time-invariant**:
$$\frac{h(t \mid \mathbf{Z}_A)}{h(t \mid \mathbf{Z}_B)} = \frac{h_0(t) \exp(\boldsymbol{\beta}^T \mathbf{Z}_A)}{h_0(t) \exp(\boldsymbol{\beta}^T \mathbf{Z}_B)} = \exp\left( \boldsymbol{\beta}^T (\mathbf{Z}_A - \mathbf{Z}_B) \right) = \text{constant}$$
If the biological effect of a biomarker wanes over time (e.g., strong prognostic impact in early years that attenuates due to secondary mutations or therapy), $\beta(t)$ becomes time-dependent, invalidating standard time-invariant hazard ratios.

#### 3.5.2 Scaled Schoenfeld Residual Formulation
To formally test for time-varying coefficients, the pipeline computes **scaled Schoenfeld residuals** (Schoenfeld, 1982 [4]; Grambsch & Therneau, 1994 [5]; via `survival::cox.zph`).

For each observed failure time $t_i$ ($\delta_i = 1$), the Schoenfeld residual vector for covariate $k$ is defined as the difference between the observed covariate value of the failing individual and the risk-weighted average of all individuals in the risk set:
$$r_{ik} = z_{ik} - \bar{z}_k(t_i)$$
Where:
$$\bar{z}_k(t_i) = \sum_{j \in \mathcal{R}(t_i)} w_j(t_i) z_{jk}, \quad w_j(t_i) = \frac{\exp(\hat{\boldsymbol{\beta}}^T \mathbf{Z}_j)}{\sum_{m \in \mathcal{R}(t_i)} \exp(\hat{\boldsymbol{\beta}}^T \mathbf{Z}_m)}$$

The **scaled Schoenfeld residual** vector $\mathbf{r}_i^*$ scales the raw residuals by the inverse of the estimated covariance matrix:
$$\mathbf{r}_i^* \approx m \cdot \left[ \mathcal{I}(\hat{\boldsymbol{\beta}}) \right]^{-1} \mathbf{r}_i + \hat{\boldsymbol{\beta}}$$
Where $m$ is the total number of events. Under the proportional hazards assumption, the expected value of the scaled Schoenfeld residual is asymptotically zero:
$$\mathbb{E}[\mathbf{r}_i^*] \approx \boldsymbol{\beta}(t_i) - \hat{\boldsymbol{\beta}}$$

#### 3.5.3 Grambsch–Therneau Correlation Test
The pipeline tests whether $\boldsymbol{\beta}(t)$ varies as a linear or monotonic function of time:
$$\beta_k(t) = \beta_k + \theta_k g(t)$$
Where $g(t)$ is a transformed time scale (e.g., Kaplan–Meier survival time transformation $g(t) = 1 - \hat{S}(t)$). 

The null hypothesis of proportional hazards ($H_0: \theta_k = 0$) is tested via the generalized correlation between $\mathbf{r}_{\cdot k}^*$ and $g(t)$:
$$\chi_{\text{zph}, k}^2 = \frac{\left( \sum_{i=1}^m (g(t_i) - \bar{g}) r_{ik}^* \right)^2}{\text{Var}\left( \sum_{i=1}^m (g(t_i) - \bar{g}) r_{ik}^* \right)} \sim \chi^2(1)$$
A non-significant $p$-value ($p > 0.05$) across individual covariates and the global omnibus test confirms that the proportional hazards assumption is mathematically valid.

---

### 3.6 Step 6: Discriminative Performance Assessment via Harrell's Concordance Index ($C$-index)
**Script**: `survival_analysis_hub_genes.R`  
**Input Artifacts**: Fitted univariable and multivariable Cox model objects.  
**Output Artifacts**:
- `concordance_summary.tsv`: Table comparing univariate and multivariate $C$-index values and standard errors across all 24 hub genes.
- `concordance_comparison.pdf` and `.png`: Horizontal bar plot illustrating $C$-index gains from unadjusted to multivariable models relative to random chance ($C = 0.5$).

#### 3.6.1 Mathematical Formulation of the Concordance Index
Statistical significance ($p < 0.05$) indicates an association with survival but does not quantify how accurately a model predicts the relative order of patient deaths. To measure discriminative capacity, the methodology calculates **Harrell's Concordance Index ($C$-index)** (Harrell et al., 1982 [6]; 1996 [7]).

For any randomly selected pair of patients $(i, j)$ where at least one patient experienced a mortality event, the pair is **informative** if the earlier observation is an event ($T_i < T_j$ and $\delta_i = 1$). 

Let $\hat{\eta}_i = \hat{\boldsymbol{\beta}}^T \mathbf{Z}_i$ denote the predicted linear risk score from the Cox model. The pair is **concordant** if the patient who died earlier was predicted to have a higher risk score ($\hat{\eta}_i > \hat{\eta}_j$). The $C$-index is formally defined as:
$$C = \frac{\sum_{i \ne j} \mathbb{I}(T_i < T_j) \cdot \mathbb{I}(\hat{\eta}_i > \hat{\eta}_j) \cdot \delta_i}{\sum_{i \ne j} \mathbb{I}(T_i < T_j) \cdot \delta_i}$$
Where $\mathbb{I}(\cdot)$ is the indicator function. If tied risk scores occur ($\hat{\eta}_i = \hat{\eta}_j$), a value of $0.5$ is added to the numerator.

#### 3.6.2 Standard Error and Confidence Interval Derivation
The standard error of the $C$-index ($\text{SE}(C)$) is estimated using the asymptotic variance formula derived from infinitesimal jackknife / U-statistic variance estimators:
$$\text{CI}_{0.95}(C) = \left[ C - 1.96 \cdot \text{SE}(C), \; C + 1.96 \cdot \text{SE}(C) \right]$$

#### 3.6.3 Interpretive Scale for Risk Discrimination
- **$C = 0.50$**: Random prediction (no discriminative ability; dashed baseline).
- **$0.51 \le C \le 0.60$**: Poor to marginal discrimination.
- **$0.61 \le C \le 0.70$**: Acceptable to good discriminative ability (typical of clinical oncology models).
- **$0.71 \le C \le 0.80$**: Strong discriminative ability.
- **$C > 0.80$**: Exceptional discrimination.
- **$C = 1.00$**: Perfect concordance (patient with higher risk score always dies first).

Comparing $C_{\text{multivariate}}$ versus $C_{\text{univariate}}$ directly quantifies the **incremental discriminative gain** provided by combining the molecular biomarker with standard staging parameters.

---

### 3.7 Step 7: Multi-Testing False Discovery Rate Control & Multi-Tiered Survival Synthesis
**Script**: `survival_analysis_hub_genes.R`  
**Input Artifacts**: Tabular results from Steps 2, 3, and 4.  
**Output Artifacts**:
- `survival_summary_table.tsv`: Master survival synthesis table integrating topological metrics, expression fold changes, cutpoints, log-rank tests, univariate Cox, and multivariate Cox statistics.
- `survival_significance_heatmap.pdf` and `.png`: Summary tile plot depicting $-\log_{10}(p)$ values across Log-Rank, Univariate Cox, and Multivariate Cox models with significance star annotations.

#### 3.7.1 Multi-Testing False Discovery Rate (FDR) Control
Across the pipeline, multiple simultaneous statistical tests are performed across $M = 24$ hub genes. To prevent the accumulation of Type I errors (false-positive biomarker discovery), nominal $p$-values are adjusted within each analytical family using the **Benjamini–Hochberg (BH) step-up procedure**:
1. Sort nominal $p$-values ascendingly: $p_{(1)} \le p_{(2)} \le \dots \le p_{(M)}$.
2. Compute adjusted false discovery rate $q$-values:
   $$q_{(i)} = \min_{k \ge i} \left( \frac{M \cdot p_{(k)}}{k} \right)$$
3. Enforce statistical significance at the standard discovery threshold:
   $$\text{FDR} \le 0.05 \quad (\text{or } q < 0.10 \text{ for exploratory trends})$$

#### 3.7.2 Multi-Tiered Hierarchical Survival Classification
To provide the client and manuscript readers with clear, transparent translational guidance, every hub gene is categorized into a rule-based prognostic tier in `survival_summary_table.tsv`:
$$\text{Prognostic Tier} = \begin{cases} \text{"Yes (univariate)"} & \text{if } q_{\text{uni}} < 0.05 \\[4pt] \text{"Yes (multivariate only)"} & \text{if } q_{\text{uni}} \ge 0.05 \land q_{\text{mv}} < 0.05 \\[4pt] \text{"Yes (KM only)"} & \text{if } q_{\text{logrank}} < 0.05 \land q_{\text{uni}} \ge 0.05 \land q_{\text{mv}} \ge 0.05 \\[4pt] \text{"No"} & \text{otherwise} \end{cases}$$

#### 3.7.3 Significance Synthesis Heatmap Architecture
The multi-method significance landscape is synthesized into a publication tile plot:
- **$x$-axis**: Evaluated survival test ($\text{Log-Rank [KM]}$, $\text{Univariate Cox}$, $\text{Multivariate Cox}$).
- **$y$-axis**: Hub genes ordered by ascending univariate Cox $p$-value.
- **Fill Color Gradient**: Transformed statistical significance:
  $$\text{Score} = -\log_{10}(p_{\text{value}})$$
  Centered at midpoint $1.301$ ($-\log_{10}(0.05)$).
- **Symbolic Overlays**: Standard scientific significance stars:
  $$\text{Annotation} = \begin{cases} \text{"***"} & \text{if } p < 0.001 \\ \text{"**"} & \text{if } p < 0.01 \\ \text{"*"} & \text{if } p < 0.05 \\ \text{""} & \text{if } p \ge 0.05 \end{cases}$$

---

## 4. Synthesis & Comparative Matrix

The following structured matrix summarizes the mathematical parameters, sample sizes, censoring proportions, and analytical outputs governing each step of the survival analysis pipeline:

| Analytical Step | Primary Computational Script | Input Data / Cohort Size | Mathematical Test / Estimator | Core Parameter / Threshold Settings | Primary Deliverable Artifacts |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Cohort & Endpoint Harmonization** | `survival_analysis_hub_genes.R` | $n = 224$ primary oral tumors; metadata table | Data curation; patient deduplication | Deduplicate `Case ID`; OS time $T > 0$; $T = \text{days} / 365.25$ | Harmonized survival cohort ($N = 222$, events = 98); Complete covariate cohort ($N = 192$, events = 87) |
| **Optimal Cutpoint Determination** | `survival_analysis_hub_genes.R` | $N = 222$ OS cohort; $\log_2\text{CPM}$ expression | Maximally Selected Rank Statistics (`surv_cutpoint`) | Search constraint: `minprop = 0.25` ($n \ge 55$ per subgroup); fallback = median | Optimal cutoff vector $\mu_g^*$; High vs. Low expression factors |
| **Stratified Survival Modeling** | `survival_analysis_hub_genes.R` | $N = 222$; High vs. Low factor per gene | Kaplan–Meier product-limit; Two-sample Log-Rank | Greenwood CI (95%); Mantel–Cox $\chi^2$ ($\text{df} = 1$); BH FDR | `km_logrank_summary.tsv`, `KM_{gene}.pdf/.png`, `km_combined_grid.pdf/.png` |
| **Continuous Univariate Cox Modeling** | `survival_analysis_hub_genes.R` | $N = 222$; continuous $\log_2\text{CPM}$ expression | Semi-parametric Cox proportional hazards | Partial likelihood maximization; Newton–Raphson; Wald test ($\text{df} = 1$) | `univariate_cox_results.tsv`, `univariate_forest_plot.pdf/.png` |
| **Multivariable Cox Adjustment** | `survival_analysis_hub_genes.R` | $N = 192$; complete clinical covariates | Multivariable Cox regression ($p = 5$ predictors) | Covariates: Age, Stage (Early/Adv), Grade (Low/High), Gender; Omnibus LRT | `multivariate_cox_results.tsv`, `multivariate_forest_plot.pdf/.png` |
| **Full Multivariable Forest Modeling** | `survival_analysis_hub_genes.R` | $N = 192$; leading prognostic hub gene | Multi-parameter hazard modeling (`forestmodel`) | Full parameter visualization; reference baselines: Male, Early Stage, Low Grade | `multivariate_full_forest_best_gene.pdf/.png` |
| **Proportional Hazards Diagnostics** | `survival_analysis_hub_genes.R` | Multivariable Cox model object (`fm_cox`) | Scaled Schoenfeld residuals (`cox.zph`) | Grambsch–Therneau correlation test; Kaplan–Meier time transform $g(t)$ | `cox_ph_diagnostics.pdf/.png` (residual plots with LOESS smoothing) |
| **Model Discrimination Evaluation** | `survival_analysis_hub_genes.R` | Univariate vs. Multivariate model pairs | Harrell's Concordance Index ($C$-statistic) | Infinitesimal jackknife $\text{SE}(C)$; null baseline $C = 0.50$ | `concordance_summary.tsv`, `concordance_comparison.pdf/.png` |
| **Prognostic Synthesis & Integration** | `survival_analysis_hub_genes.R` | Integrated results from all modeling steps | Multi-tier classification; $-\log_{10}(p)$ heatmap | Tier rules: Univariate, Multivariate, KM-only; BH FDR $q < 0.05$ | `survival_summary_table.tsv`, `survival_significance_heatmap.pdf/.png` |

---

## 5. Ready-to-Publish Methods Section (Manuscript Reference Text)

### 2.5 Survival Association and Prognostic Modeling of PIKK Hub Genes

#### 2.5.1 Patient Cohort and Clinical Endpoint Definition
To evaluate the clinical prognostic relevance of the 24 topologically prioritized PIKK-associated hub genes, clinical annotations and longitudinal tracking data were obtained from The Cancer Genome Atlas Head and Neck Squamous Cell Carcinoma (TCGA-HNSC) curated oral cavity cohort [8]. Analysis was strictly restricted to primary solid tumor specimens ($n = 224$); normal mucosal tissues and metastatic lesions were excluded. Multiple sequenced vials from individual patients were deduplicated to unique biological cases by retaining the primary diagnostic vial (`distinct(Case ID)`). 

The primary clinical endpoint was Overall Survival (OS), defined as the time interval in decimal years from initial histological diagnosis to death from any cause (failure event, $\delta = 1$) or to the date of last confirmed follow-up (right-censored, $\delta = 0$). Patients lacking valid longitudinal tracking or exhibiting zero/negative survival durations ($T \le 0$) were excluded, yielding a verified survival analytical cohort of $n = 222$ patients ($n = 98$ mortality events, $n = 124$ censored observations; event rate: 44.1%). 

For multivariable modeling, clinical covariates were harmonized into standardized biological categories: patient age at diagnosis (continuous, years), AJCC pathologic tumor stage (Early [Stages I–II] vs. Advanced [Stages III–IV]), histologic tumor grade (Low Grade [G1–G2] vs. High Grade [G3–G4]), and patient gender (Male vs. Female). Complete-case restriction across all four clinical covariates defined a multivariable modeling cohort of $n = 192$ patients ($n = 87$ mortality events, $n = 105$ censored observations).

#### 2.5.2 Stratified Survival Analysis and Data-Adaptive Cutpoint Determination
Rather than imposing arbitrary median dichotomization, optimal stratification thresholds for continuous normalized transcript abundance ($\log_2\text{CPM}$) were determined using maximally selected rank statistics (`survminer::surv_cutpoint` in R) [3, 9]. The algorithm evaluated standardized two-sample log-rank statistics across candidate expression values, selecting the cutpoint that maximized inter-group survival discrimination while enforcing a minimum subgroup proportion constraint (`minprop = 0.25`) to guarantee that both High and Low expression cohorts maintained at least 25% of the total cohort ($n \ge 55$ patients). In cases where maximally selected rank optimization did not converge, the cohort median was utilized as a default threshold.

Non-parametric cumulative survival functions were estimated using the Kaplan–Meier product-limit method [1], with variance and 95% confidence intervals calculated via Greenwood’s formula. Survival differences between High and Low expression groups were evaluated using the two-sample Mantel–Cox log-rank test (`survival::survdiff`, $\text{df} = 1$) [10]. Resulting $p$-values were adjusted for multiple testing across all 24 candidate genes using the Benjamini–Hochberg False Discovery Rate (FDR) procedure [11].

#### 2.5.3 Univariate and Multivariable Cox Proportional Hazards Regression
To assess prognostic associations across the continuous expression spectrum without the information loss inherent to binning, semi-parametric Cox proportional hazards regression models were fitted using `survival::coxph` [2, 10]. Univariate models evaluated the unadjusted hazard ratio (HR) per 1-unit increase in $\log_2\text{CPM}$ transcript abundance, with 95% Wald confidence intervals and two-sided $p$-values derived from the partial likelihood information matrix. Multiple testing across the candidate gene panel was controlled using Benjamini–Hochberg FDR adjustment.

To establish whether candidate hub genes provided independent prognostic utility beyond standard clinicopathologic risk factors, multivariable Cox proportional hazards models were constructed. Each model simultaneously incorporated the continuous biomarker expression alongside patient age at diagnosis (continuous), AJCC pathologic stage (Early [I–II] vs. Advanced [III–IV]), histologic grade (Low [G1–G2] vs. High [G3–G4]), and gender (Male vs. Female). Adjusted hazard ratios ($\text{adj.HR}$), 95% confidence intervals, and Wald test $p$-values were computed for each predictor. Global multivariable model significance was verified using the omnibus Likelihood Ratio Test against the null model ($\text{df} = 5$). Visualizations of hazard ratios and confidence intervals across models were generated using `ggplot2` and `forestmodel` [12].

#### 2.5.4 Proportional Hazards Assumption Verification and Discriminative Performance Evaluation
The validity of the proportional hazards assumption was formally assessed for multivariable models using scaled Schoenfeld residuals and Grambsch–Therneau correlation tests against transformed survival time (`survival::cox.zph`) [4, 5, 10]. Diagnostic plots of scaled residuals over time were generated with LOESS smoothing splines; absence of statistically significant non-zero slopes ($p > 0.05$) confirmed the temporal invariance of regression coefficients.

Model discriminative capacity was quantified using Harrell’s Concordance Index ($C$-index) [6, 7], which calculates the proportion of all informative patient pairs where predicted linear risk scores concord with the observed chronological sequence of death events. Asymptotic standard errors were estimated using infinitesimal jackknife variance estimators. Incremental discriminative gains achieved by integrating molecular biomarkers into clinical staging were determined by comparing multivariable model concordance against baseline univariate models. Multi-tiered survival significance was summarized in integrated cross-method matrices and $-\log_{10}(p)$ significance heatmaps.

---

## 6. Comprehensive Software, Dependency, and Environment Manifest

All survival association analyses, cutpoint optimizations, regression modeling, diagnostic testing, and publication visualizations were executed within the R statistical computing environment (R version $\ge 4.2.0$). The following specialized libraries were utilized:

| Package Name | Origin / Repository | Version Tested | Primary Computational Purpose in Pipeline |
| :--- | :--- | :--- | :--- |
| **`survival`** | CRAN | v3.5-7 | Core survival analysis engine: `Surv` object construction, `survfit` Kaplan–Meier estimation, `survdiff` log-rank testing, `coxph` continuous univariable and multivariable proportional hazards regression, and `cox.zph` Schoenfeld residual diagnostics [10]. |
| **`survminer`** | CRAN | v0.4.9 | Data-adaptive cutpoint optimization via `surv_cutpoint` (maximally selected rank statistics); generation of publication-grade `ggsurvplot` figures with integrated number-at-risk tables and confidence intervals [9]. |
| **`forestmodel`** | CRAN | v0.6.2 | Automated compilation of multi-parameter multivariable forest plots depicting adjusted hazard ratios, 95% CIs, and covariate reference levels (`forest_model`) [12]. |
| **`ggplot2`** | CRAN | v3.4.4 | High-resolution publication graphics generation, custom univariate and multivariate forest plots, and survival significance heatmaps. |
| **`patchwork`** | CRAN | v1.2.0 | Multi-panel plot composition, alignment, and hierarchical annotation for combined Kaplan–Meier survival grids. |
| **`gridExtra`** | CRAN | v2.3 | Low-level grid alignment and tabular figure layout rendering. |
| **`grid`** | Base R | v4.2.0 | Core graphic primitive manipulation, viewports, and multi-panel device layouts. |
| **`scales`** | CRAN | v1.3.0 | Coordinate scaling, scientific formatting, and aesthetic axis transformations. |
| **`RColorBrewer`** | CRAN | v1.1-3 | Colorblind-safe diverging and sequential color palettes for heatmaps and risk strata. |
| **`dplyr`** | CRAN | v1.1.4 | Tidy data manipulation, clinical filtering, patient case deduplication, and tabular transformation. |
| **`tidyr`** | CRAN | v1.3.1 | Data pivoting (`pivot_longer`, `pivot_wider`) for multi-test significance heatmaps. |
| **`readr`** | CRAN | v2.1.5 | Fast, type-safe serialization of TSV and CSV input/output artifacts. |
| **`stringr`** | CRAN | v1.5.1 | Regular expression processing for clinical stage harmonization and covariate cleaning. |

---

## 7. Peer-Reviewed Scientific Bibliography

1. **Kaplan, E. L., & Meier, P.** (1958). Nonparametric estimation from incomplete observations. *Journal of the American Statistical Association*, 53(282), 457–481. [DOI: 10.1080/01621459.1958.10501452](https://doi.org/10.1080/01621459.1958.10501452)
2. **Cox, D. R.** (1972). Regression models and life-tables (with discussion). *Journal of the Royal Statistical Society: Series B (Methodological)*, 34(2), 187–220. [DOI: 10.1111/j.2517-6161.1972.tb00899.x](https://doi.org/10.1111/j.2517-6161.1972.tb00899.x)
3. **Hothorn, T., & Lausen, B.** (2003). On the exact distribution of maximally selected rank statistics. *Computational Statistics & Data Analysis*, 43(2), 121–137. [DOI: 10.1016/S0167-9473(02)00225-6](https://doi.org/10.1016/S0167-9473(02)00225-6)
4. **Schoenfeld, D.** (1982). Partial residuals for the proportional hazards regression model. *Biometrika*, 69(1), 239–241. [DOI: 10.1093/biomet/69.1.239](https://doi.org/10.1093/biomet/69.1.239)
5. **Grambsch, P. M., & Therneau, T. M.** (1994). Proportional hazards tests and diagnostics based on weighted residuals. *Biometrika*, 81(3), 515–526. [DOI: 10.1093/biomet/81.3.515](https://doi.org/10.1093/biomet/81.3.515)
6. **Harrell, F. E., Califf, R. M., Pryor, D. B., Lee, K. L., & Rosati, R. A.** (1982). Evaluating the yield of medical tests. *JAMA*, 247(18), 2543–2546. [DOI: 10.1001/jama.1982.03320430047030](https://doi.org/10.1001/jama.1982.03320430047030)
7. **Harrell, F. E., Lee, K. L., & Mark, D. B.** (1996). Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. *Statistics in Medicine*, 15(4), 361–387. [DOI: 10.1002/(SICI)1097-0258(19960229)15:4<361::AID-SIM168>3.0.CO;2-4](https://doi.org/10.1002/(SICI)1097-0258(19960229)15:4<361::AID-SIM168>3.0.CO;2-4)
8. **Cancer Genome Atlas Network.** (2015). Comprehensive genomic characterization of head and neck squamous cell carcinomas. *Nature*, 517(7536), 576–582. [DOI: 10.1038/nature14129](https://doi.org/10.1038/nature14129)
9. **Kassambara, A., Kosinski, M., & Biecek, P.** (2021). *survminer: Drawing Survival Curves using 'ggplot2'*. R package version 0.4.9. [CRAN: survminer](https://CRAN.R-project.org/package=survminer)
10. **Therneau, T. M.** (2023). *A Package for Survival Analysis in R*. R package version 3.5-7. [CRAN: survival](https://CRAN.R-project.org/package=survival)
11. **Benjamini, Y., & Hochberg, Y.** (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B (Methodological)*, 57(1), 289–300. [DOI: 10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
12. **Kennedy, N.** (2022). *forestmodel: Forest Plots from Regression Models*. R package version 0.6.2. [CRAN: forestmodel](https://CRAN.R-project.org/package=forestmodel)
