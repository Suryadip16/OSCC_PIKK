# Functional Enrichment and Pathway Analysis (Section 2.3) — Comprehensive Methodology Guide

**Project**: PIKK Pathway-Associated Biomarker Discovery and Oncogenic Characterization in Oral Squamous Cell Carcinoma (OSCC)  
**Target Manuscript Reference**: Materials and Methods / Bioinformatic Processing Pipeline  
**Primary Dataset**: The Cancer Genome Atlas Head and Neck Squamous Cell Carcinoma (TCGA-HNSC), Curated Oral Cavity Cohort ($n = 240$ primary samples; $n = 224$ Primary Solid Tumor, $n = 16$ Solid Tissue Normal)  
**Background Universe**: 16,905 Tested Protein-Coding Genes (pre-filtered via `filterByExpr` and modeled in edgeR QLF) mapped to unique NCBI Entrez Gene Identifiers  
**Target Foreground**: Prioritized PIKK-Response Candidate Gene Set  
**Knowledge Bases Evaluated**: Gene Ontology (GO: Biological Process, Molecular Function, Cellular Component), Kyoto Encyclopedia of Genes and Genomes (KEGG), and Reactome Pathway Knowledgebase  
**Scripts Evaluated**:
1. `ora_edgeR.R` — *Over-Representation Analysis (ORA), Signed Quasi-Likelihood Metric Construction, and Unbiased Gene Set Enrichment Analysis (GSEA)*

---

## 1. Executive Summary & End-to-End Workflow

While differential expression analysis isolates discrete transcripts displaying statistically significant alterations between malignant and non-malignant tissues, cellular phenotypes are driven by higher-order functional modules, multi-protein complexes, and interconnected signaling cascades. Functional enrichment analysis translates unstructured gene lists into mechanistic insights, contextualizing oncogenic alterations within established biochemical and physiological pathways.

This document establishes the rigorous, publication-grade methodology governing the **functional annotation, pathway over-representation, and genome-wide gene set enrichment workflows** implemented for the OSCC PIKK biomarker discovery program.

The analytical pipeline implemented in `ora_edgeR.R` utilizes a dual, complementary bioinformatics strategy:
1. **Targeted Over-Representation Analysis (ORA)**: Tests whether the prioritized PIKK-response gene signature displays statistically significant over-representation within specific biological processes, molecular functions, cellular components, and pathway databases relative to an empirical background universe.
2. **Unbiased Gene Set Enrichment Analysis (GSEA)**: Evaluates global, coordinated pathway shifts across the entire transcriptome ranked by a signed quasi-likelihood test statistic, capturing biological processes driven by distributed, cumulative expression changes that circumvent rigid significance thresholds.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                                 UPSTREAM INPUT ARTIFACTS                               │
│  • Primary Foreground: Significant PIKK-Response Intersection Candidates               │
│  • Statistical Universe: edgeR All-Genes Differential Results (16,905 tested features)│
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 1: IDENTIFIER HARMONIZATION & BACKGROUND DEFINITION                             │
│  • Automated HGNC symbol standardization (case-normalization & deduplication)          │
│  • Cross-database Entrez ID mapping via clusterProfiler::bitr and org.Hs.eg.db         │
│  • Definition of empirical assay-specific background universe (N_universe)             │
│  • Elimination of false-positive inflation caused by whole-genome backgrounds          │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 2: OVER-REPRESENTATION ANALYSIS (ORA) VIA HYPERGEOMETRIC MODELING               │
│  • One-tailed Fisher's exact / Hypergeometric testing across 5 ontology spaces:        │
│      - GO Biological Process (BP), Molecular Function (MF), Cellular Component (CC)    │
│      - KEGG Pathway Database (organism = "hsa")                                        │
│      - Reactome Knowledgebase (organism = "human")                                     │
│  • Multiple testing correction via Benjamini-Hochberg (BH) FDR (p ≤ 0.05, q ≤ 0.20)     │
│  • Calculation of Rich Factor, Fold Enrichment, GeneRatio, and BgRatio metrics         │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 3: SIGNED QUASI-LIKELIHOOD METRIC DERIVATION FOR GSEA                           │
│  • Extraction of fold change (log2FC) and quasi-likelihood F-statistic (F)             │
│  • Construction of continuous signed ranking statistic: s_g = sign(log2FC) × sqrt(F)  │
│  • Ranking of entire background transcriptome from top-activated to top-repressed     │
│  • Disambiguation of multi-mapped Entrez IDs via maximum absolute statistic            │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 4: CONTINUOUS GENE SET ENRICHMENT ANALYSIS (GSEA)                               │
│  • Weighted Kolmogorov-Smirnov running sum calculation across GO-BP, KEGG, & Reactome │
│  • Gene set size filtering: 10 ≤ |S| ≤ 500 (mitigating broad or idiosyncratic bias)    │
│  • Normalized Enrichment Score (NES) calculation via permutation-based null models    │
│  • Phenotypic pathway bifurcation: Upregulated (NES > 0) vs Downregulated (NES < 0)   │
│  • Core enrichment (leading-edge) gene extraction                                     │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 5: MULTI-TIER PUBLICATION VISUALIZATION SUITE                                   │
│  • ORA Visualizations: Dot plots, Fold Enrichment bar charts, Magma bubble plots      │
│  • GSEA Visualizations: Multi-panel running score graphs (gseaplot2), Ridgeplots,     │
│    and direction-split dot plots                                                       │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. Biological Context and Study Rationale

### 2.1 Biological Imperative in Oral Carcinogenesis
Oral Squamous Cell Carcinoma (OSCC) is characterized by severe chromosomal instability, complex structural rearrangements, and widespread transcriptional reprogramming driven by tobacco carcinogens, alcohol abuse, and chronic mucosal irritation. In this context, functional enrichment analysis is critical to decipher how disruptions within the **PIKK (Phosphatidylinositol 3-kinase-related kinase) family signaling network** ($ATM$, $ATR$, $PRKDC$, $MTOR$, $TRRAP$, and $SMG1$) propagate across downstream cellular programs.

Transcriptional deregulation of PIKK nodes in OSCC directly impinges upon:
1. **Double-Strand Break (DSB) Repair and Checkpoint Control**: ATM- and PRKDC-mediated non-homologous end joining (NHEJ) and homologous recombination repair (HRR), preventing catastrophic mitotic death under persistent genotoxic stress.
2. **Replication Fork Protection and Intra-S Checkpoint Signaling**: ATR-mediated stabilization of stalled replication forks, preventing double-strand break conversion during oncogene-driven accelerated S-phase progression.
3. **Translational Competence, Autophagy, and Metabolic Rewiring**: mTORC1/mTORC2 signaling, coordinating nutrient uptake, ribosome biogenesis, lipid synthesis, and suppression of macroautophagy.
4. **Epigenetic Plasticity and Chromatin Remodeling**: TRRAP-containing histone acetyltransferase complexes (e.g., TIP60/NuA4, SAGA), modulating chromatin accessibility at DNA damage foci and oncogenic promoter elements.
5. **Transcript Quality Surveillance**: SMG1-dependent nonsense-mediated mRNA decay (NMD), degrading aberrant transcripts harboring premature termination codons and modulating stress-granule dynamics.

Evaluating these interconnected biological networks requires broad ontology profiling that goes beyond individual gene annotations to uncover higher-order oncogenic signatures.

### 2.2 Dual-Paradigm Rationale: Over-Representation vs. Gene Set Enrichment
To ensure analytical robustness, the study deploys two fundamentally distinct, complementary computational paradigms:

```
Table 1: Methodological Comparison of Enrichment Paradigms
┌──────────────────────────────────────┬─────────────────────────────┬─────────────────────────────┐
│ Methodological Dimension             │ Over-Representation (ORA)   │ Gene Set Enrichment (GSEA)  │
├──────────────────────────────────────┼─────────────────────────────┼─────────────────────────────┤
│ Input Data Architecture              │ Discrete, thresholded gene  │ Complete continuous ranked  │
│                                      │ list (foreground subset)    │ list (all tested genes)     │
│ Primary Statistical Model            │ Hypergeometric distribution │ Weighted Kolmogorov-Smirnov │
│                                      │ (Fisher's exact test)       │ running sum statistic       │
│ Hypothesis Tested                    │ Are pathway genes over-     │ Are pathway genes non-      │
│                                      │ represented in foreground?  │ randomly distributed in rank│
│ Sensitivity to Modest Shifts         │ Blind to sub-threshold      │ Highly sensitive to small,  │
│                                      │ genes (|log2FC| < cutoff)   │ coordinated pathway shifts  │
│ Directional Stratification           │ Evaluated on pre-filtered   │ Naturally bifurcates into   │
│                                      │ sets or pooled foregrounds  │ NES > 0 vs NES < 0 programs │
│ Primary Metric of Significance       │ Hypergeometric p-value / BH │ Permutation FDR /           │
│                                      │ adjusted p-value (p.adjust) │ Normalized Enrichment Score │
│ Primary Quantitative Magnitude       │ Fold Enrichment & Rich      │ Normalized Enrichment Score │
│                                      │ Factor                      │ (NES) & Leading Edge Size   │
└──────────────────────────────────────┴─────────────────────────────┴─────────────────────────────┘
```

1. **Why ORA?** ORA provides rigorous, hypothesis-driven verification that our prioritized PIKK-response gene signature is specifically concentrated in targeted DDR, checkpoint, and metabolic modules rather than representing stochastic background noise.
2. **Why GSEA?** Complex metabolic cascades, multi-subunit chromatin remodeling complexes, and broad immune programs often operate through subtle, concurrent shifts across dozens of pathway components. Imposing hard significance cutoffs (e.g., $|\log_2\text{FC}| \ge 1.0$) causes such pathways to fail ORA detection. GSEA circumvents thresholding bias by analyzing the complete ranked continuum of 16,905 genes.

---

## 3. Step-by-Step Computational Methodology

### 3.1 Step 1: Identifier Mapping, Harmonization, and Background Definition
**Script**: `ora_edgeR.R`  
**Input Artifacts**:
- `ext_sig_file`: Significant PIKK-response intersection table (`OSCC_PIKK_extended_edgeR_DEG_intersection_significant.tsv`).
- `edgeR_all_file`: Complete differential expression statistics from edgeR (`edgeR_all_results_tumor_vs_normal.tsv`).
**Output Artifacts**:
- `gene_mapping_summary_edgeR.tsv`: Record of mapped versus unmapped identifiers across foreground and background universes.

#### 3.1.1 Biological and Methodological Rationale
Biological knowledge bases (Gene Ontology, KEGG, and Reactome) index functional terms using NCBI Entrez Gene Identifiers, whereas RNA-seq alignment pipelines and primary data tables report HUGO Gene Nomenclature Committee (HGNC) gene symbols. Mismatches frequently occur due to historical symbol aliasing, deprecated nomenclature, or multi-mapping.

Furthermore, a critical source of false-positive error in ORA studies is the **improper definition of the statistical background universe**. Default software configurations frequently use the entire theoretical human genome ($\approx 20,000–40,000$ annotated genes) as the background. However, solid oral epithelial tissues express only a subset of the genome; thousands of genes (e.g., neuronal, olfactory, or cardiac-specific genes) have zero probability of being detected in an oral mucosal biopsy. Using a whole-genome background artificially inflates the background ratio, producing spuriously small $p$-values and high false-discovery rates (Huang, Sherman, & Lempicki, *Nucleic Acids Res* 2009 [1]). 

To maintain statistical integrity, the methodology explicitly restricts the background universe to the **16,905 protein-coding genes robustly expressed and statistically evaluated within the edgeR QLF pipeline**.

#### 3.1.2 Algorithmic Implementation
1. **Symbol Harmonization**:
   Gene symbol strings across foreground and background data tables are coerced to uppercase to eliminate case-sensitivity mismatches:
   $$\text{gene\_symbol} \leftarrow \text{toupper}(\text{gene\_name})$$
2. **Cross-Database Identifier Translation (`symbol_to_entrez`)**:
   Identifier translation is performed via `clusterProfiler::bitr` querying the organism-level annotation database `org.Hs.eg.db` (Yu et al., 2012 [2]):
   $$\mathcal{M} = \text{bitr}(\mathcal{G}_{\text{symbols}}, \text{fromType} = \text{"SYMBOL"}, \text{toType} = \text{"ENTREZID"}, \text{OrgDb} = \text{org.Hs.eg.db})$$
3. **Deduplication and Redundancy Filtering**:
   Where a gene symbol maps to multiple Entrez identifiers, the primary mapping is retained, and duplicate entries are discarded:
   $$\mathcal{M}_{\text{unique}} = \mathcal{M}[! \text{duplicated}(\mathcal{M}\$\text{SYMBOL}), ]$$
4. **Universe Set Formalization**:
   - Foreground Entrez Set: $\mathcal{S}_{\text{foreground}} = \{e \in \text{Entrez} \mid e \in \text{Mapped}(\mathcal{G}_{\text{PIKK\_sig}})\}$.
   - Background Universe Set: $\mathcal{U}_{\text{background}} = \{e \in \text{Entrez} \mid e \in \text{Mapped}(\mathcal{G}_{\text{edgeR\_all}})\}$.
   - Condition enforced: $\mathcal{S}_{\text{foreground}} \subseteq \mathcal{U}_{\text{background}}$.

---

### 3.2 Step 2: Over-Representation Analysis (ORA) via the Hypergeometric Model
**Script**: `ora_edgeR.R`  
**Input Artifacts**: Mapped foreground vector ($\mathcal{S}_{\text{foreground}}$) and background universe vector ($\mathcal{U}_{\text{background}}$).  
**Output Artifacts**:
- Tables: `ORA_GO_BP_edgeR_extended.tsv`, `ORA_GO_MF_edgeR_extended.tsv`, `ORA_GO_CC_edgeR_extended.tsv`, `ORA_KEGG_edgeR_extended.tsv`, `ORA_Reactome_edgeR_extended.tsv`.
- Publication Plots: Dot plots (`_dotplot.pdf/.png`), bar charts (`_barplot.pdf/.png`), and magma bubble plots (`_bubble.pdf/.png`).

#### 3.2.1 Mathematical Formulation of the Hypergeometric Model
Let the background universe contain $N$ total genes ($N = |\mathcal{U}_{\text{background}}|$), of which $M$ genes are annotated to a specific functional pathway or biological term $\mathcal{P}$. Suppose the experimental foreground contains $n$ prioritized genes ($n = |\mathcal{S}_{\text{foreground}}|$), of which $k$ genes are annotated to pathway $\mathcal{P}$.

```
Contingency Matrix for Pathway P:
┌──────────────────────────────┬─────────────────────────────┬─────────────────────────────┬─────────────┐
│ Category                     │ In Pathway P                │ Not in Pathway P            │ Total       │
├──────────────────────────────┼─────────────────────────────┼─────────────────────────────┼─────────────┤
│ In Foreground Set            │ k (Overlap)                 │ n - k                       │ n           │
│ In Background (Not in Set)   │ M - k                       │ (N - M) - (n - k)           │ N - n       │
├──────────────────────────────┼─────────────────────────────┼─────────────────────────────┼─────────────┤
│ Total Universe               │ M                           │ N - M                       │ N           │
└──────────────────────────────┴─────────────────────────────┴─────────────────────────────┴─────────────┘
```

Under the null hypothesis ($H_0$), the foreground set is assumed to be a random sample drawn without replacement from the background universe. The probability of observing exactly $k$ genes belonging to pathway $\mathcal{P}$ follows the Hypergeometric distribution:
$$\mathbb{P}(X = k) = \frac{\displaystyle\binom{M}{k} \binom{N - M}{n - k}}{\displaystyle\binom{N}{n}}$$
Where:
- $\binom{M}{k}$ is the number of ways to select $k$ pathway genes from the $M$ available pathway genes.
- $\binom{N - M}{n - k}$ is the number of ways to select the remaining $n - k$ non-pathway genes from the $N - M$ non-pathway background genes.
- $\binom{N}{n}$ is the total number of ways to draw $n$ genes from the background universe $N$.

The statistical significance (one-tailed cumulative $p$-value) representing the probability of observing $k$ or more pathway genes by chance is:
$$p = \mathbb{P}(X \ge k) = \sum_{i = k}^{\min(n, M)} \frac{\displaystyle\binom{M}{i} \binom{N - M}{n - i}}{\displaystyle\binom{N}{n}}$$

#### 3.2.2 Quantitative Enrichment Metrics
To assess biological magnitude alongside statistical significance, four standardized metrics are computed:
1. **Gene Ratio**: The fraction of foreground genes that belong to pathway $\mathcal{P}$:
   $$\text{GeneRatio} = \frac{k}{n}$$
2. **Background Ratio**: The fraction of background universe genes annotated to pathway $\mathcal{P}$:
   $$\text{BgRatio} = \frac{M}{N}$$
3. **Fold Enrichment**: The ratio of the observed pathway proportion in the foreground relative to the background:
   $$\text{FoldEnrichment} = \frac{\text{GeneRatio}}{\text{BgRatio}} = \frac{k / n}{M / N} = \frac{k \cdot N}{n \cdot M}$$
   A $\text{FoldEnrichment} > 1.0$ indicates positive over-representation.
4. **Rich Factor**: The proportion of genes in pathway $\mathcal{P}$ that are successfully captured within the experimental foreground:
   $$\text{RichFactor} = \frac{k}{M}$$

#### 3.2.3 Biological Knowledge Bases Interrogated
1. **Gene Ontology (GO)**:
   - **Biological Process (BP)**: Identifies multi-step biological objectives (e.g., double-strand break repair, replication fork processing).
   - **Molecular Function (MF)**: Captures elemental biochemical activities (e.g., ATP-dependent 3'-5' DNA helicase activity, protein kinase binding).
   - **Cellular Component (CC)**: Localizes active complexes to subcellular structures (e.g., site of double-strand break, replication fork, TORC1 complex).
2. **Kyoto Encyclopedia of Genes and Genomes (KEGG)**:
   - Evaluates high-level pathway wiring diagrams (`organism = "hsa"`), covering cell cycle checkpoints, Fanconi anemia pathways, and central carbon metabolism.
3. **Reactome Pathway Knowledgebase**:
   - Utilizes a reaction-centric hierarchical ontology (`organism = "human"`) capturing discrete biochemical transitions, signal transduction cascades, and post-translational regulation.

#### 3.2.4 False Discovery Rate Control
To correct for the thousands of simultaneous hypotheses tested across each ontology, nominal $p$-values are adjusted using the Benjamini–Hochberg (BH) step-up procedure:
$$P_{\text{adj}(i)} = \min_{j \ge i} \left( \frac{K \cdot p_{(j)}}{j} \right)$$
Where $K$ is the total number of tested pathway terms within the ontology. Terms are considered significantly enriched if:
$$P_{\text{adj}} \le 0.05 \quad \text{and} \quad q \le 0.20$$

---

### 3.3 Step 3: Derivation of the Signed Quasi-Likelihood Ranking Metric for GSEA
**Script**: `ora_edgeR.R`  
**Input Artifacts**: Complete edgeR differential expression results (`edgeR_all_results_tumor_vs_normal.tsv`).  
**Output Artifacts**: `GSEA_ranked_gene_list_edgeR.tsv`.

#### 3.3.1 Mathematical Formulation of the Ranking Statistic
GSEA requires a single, continuous, ordered numerical metric ($s_g$) for every gene $g$ in the transcriptome, ranking genes from the most strongly upregulated in tumor tissue to the most strongly downregulated.

In linear models, ranking purely by log-fold change ($\log_2\text{FC}$) is inadequate because it ignores statistical confidence, allowing noisy, low-count transcripts with erratic variance to dominate the top ranks. Conversely, ranking purely by $p$-value or $F$-statistic strips away the **direction of biological regulation**.

To combine statistical confidence and effect magnitude, we construct the **Signed Quasi-Likelihood Statistic**:
$$s_g = \text{sign}(\log_2\text{FC}_g) \times \sqrt{F_g}$$
Where:
- $\log_2\text{FC}_g$ is the empirical log-fold change from the edgeR GLM fit.
- $F_g$ is the Quasi-Likelihood $F$-statistic from `glmQLFTest`.
- $\text{sign}(\cdot)$ is the signum function:
  $$\text{sign}(x) = \begin{cases} +1 & \text{if } x > 0 \\ 0 & \text{if } x = 0 \\ -1 & \text{if } x < 0 \end{cases}$$

#### 3.3.2 Equivalence to the Signed Wald $z$-Score
Because the quasi-likelihood $F$-statistic tests a single numerator degree of freedom ($k = 1$ for `Tissue_TypeTumor`), $F_g$ is mathematically identical to the square of a quasi-$t$ statistic:
$$F_g = t_g^2 \implies \sqrt{F_g} = |t_g|$$
Multiplying by $\text{sign}(\log_2\text{FC}_g)$ yields:
$$s_g = \text{sign}(\log_2\text{FC}_g) \times |t_g| = t_g$$
This transformation yields a signed, variance-standardized test statistic analogous to the Wald $z$-statistic in DESeq2, providing an optimal input for weighted Kolmogorov–Smirnov scoring (Subramanian et al., *PNAS* 2005 [3]).

#### 3.3.3 List Formatting and Entrez Disambiguation
1. Gene symbols are mapped to unique Entrez identifiers via `org.Hs.eg.db`.
2. Where multiple gene symbols map to a single Entrez record, the entry with the highest absolute ranking score is preserved:
   $$\text{geneList}(e) = \arg\max_{g \in \text{Entrez}^{-1}(e)} |s_g|$$
3. The final ranked vector $L = \{g_1, g_2, \dots, g_N\}$ is sorted in descending order:
   $$s_{g_1} \ge s_{g_2} \ge \dots \ge s_{g_N}$$

---

### 3.4 Step 4: Continuous Gene Set Enrichment Analysis (GSEA) Random Walk Algorithm
**Script**: `ora_edgeR.R`  
**Input Artifacts**: Sorted ranked gene list ($L$)  
**Output Artifacts**:
- Tables: `GSEA_GO_BP_edgeR_allgenes.tsv`, `GSEA_KEGG_edgeR_allgenes.tsv`, `GSEA_Reactome_edgeR_allgenes.tsv`.
- Subsets: Direction-split tables (`_upregulated_terms.tsv`, `_downregulated_terms.tsv`).
- Visualizations: Dot plots, ridgeplots, and triple-panel running enrichment plots (`gseaplot2`).

#### 3.4.1 Mathematical Formulation of the Weighted Random Walk
Let $L = \{g_1, g_2, \dots, g_N\}$ be the genome-wide ranked gene list, and let $S$ denote an annotated biological gene set containing $N_H$ genes ($N_H = |S \cap L|$).

The GSEA algorithm executes a weighted random walk down the ordered list $L$, computing an accumulated running sum from rank $i = 1$ to $N$:
1. **Fraction of Genes in $S$ (Hits)**:
   The accumulated fraction of genes in gene set $S$, weighted by their ranking statistic $|s_j|^p$ (with weighting exponent $p = 1$):
   $$P_{\text{hit}}(S, i) = \sum_{\substack{g_j \in S \\ j \le i}} \frac{|s_j|^p}{N_R}, \quad \text{where } N_R = \sum_{g_j \in S} |s_j|^p$$
2. **Fraction of Genes Not in $S$ (Misses)**:
   The accumulated fraction of genes outside gene set $S$:
   $$P_{\text{miss}}(S, i) = \sum_{\substack{g_j \notin S \\ j \le i}} \frac{1}{N - N_H}$$
3. **Enrichment Score (ES)**:
   The Enrichment Score ($\text{ES}$) corresponds to the maximum deviation of the running difference from zero:
   $$\text{ES}(S) = \arg\max_{\delta} \left\{ P_{\text{hit}}(S, i) - P_{\text{miss}}(S, i) \right\}_{i=1}^N$$
   - **Positive ES ($\text{ES} > 0$)**: Indicates that genes of set $S$ cluster at the top of the ranked list (upregulated in primary oral carcinoma).
   - **Negative ES ($\text{ES} < 0$)**: Indicates that genes of set $S$ cluster at the bottom of the ranked list (repressed in oral carcinoma).

```
GSEA Running Score Architecture:
Running Sum
  ^
  │        /\  <- Peak = Maximum Deviation (Enrichment Score, ES)
  │       /  \
  │      /    \
  │_____/      \_______________________
0 ┼────────────────────────────────────> Ranked List Index (1 to N)
  │
Hits:   || | |  ||||  |  |   |    |    (Core Enrichment / Leading Edge)
Rank:  [Top Upregulated] ----------> [Top Downregulated]
```

#### 3.4.2 Normalized Enrichment Score (NES) and Permutation Testing
Because raw Enrichment Scores correlate with gene set size, scores cannot be compared directly across pathways.
1. **Empirical Permutation**:
   A null distribution of enrichment scores ($\text{ES}_{\text{null}}$) is generated by permuting gene labels across the ranked list ($B = 1,000$ permutations).
2. **Normalization**:
   The **Normalized Enrichment Score (NES)** scales the observed $\text{ES}$ by the mean of the positive or negative null permutations:
   $$\text{NES}(S) = \begin{cases} \frac{\text{ES}(S)}{\text{mean}(\text{ES}_{\text{null}} \mid \text{ES}_{\text{null}} > 0)} & \text{if } \text{ES}(S) \ge 0 \\ -\frac{\text{ES}(S)}{\text{mean}(\text{ES}_{\text{null}} \mid \text{ES}_{\text{null}} < 0)} & \text{if } \text{ES}(S) < 0 \end{cases}$$
3. **Leading-Edge (Core Enrichment) Subset**:
   For pathways with $\text{NES} > 0$, the leading-edge subset comprises all genes in $S$ encountered from rank 1 up to the peak index $i_{\text{peak}} = \arg\max_i (P_{\text{hit}} - P_{\text{miss}})$. These genes represent the primary biological drivers accounting for the enrichment signal.

#### 3.4.3 Filtering and Significance Thresholds
- **Gene Set Size Constraints**: Gene sets with fewer than 10 genes or more than 500 genes are excluded (`minGSSize = 10`, `maxGSSize = 500`) to prevent bias from overly specific or broadly non-specific categories.
- **Significance Floor**: Pathways are designated as significantly enriched at:
  $$\text{FDR } (q\text{-value}) \le 0.05$$

---

### 3.5 Step 5: Multi-Tier Publication Visualization Architectures
**Script**: `ora_edgeR.R`  
**Output Artifacts**: Comprehensive vector graphics (PDF) and publication raster graphics (PNG, 300 DPI).

#### 3.5.1 ORA Visualizations
1. **Enrichment Dot Plots (`dotplot`)**:
   Plots the top enriched terms along the $y$-axis, displaying **Gene Ratio** on the $x$-axis, point size scaled to **Gene Count** ($k$), and point fill mapped to statistical significance ($-\log_{10}(P_{\text{adj}})$).
2. **Fold Enrichment Bar Charts (`plot_ora_bar`)**:
   Visualizes the top 20 terms ordered by significance, plotting **Fold Enrichment** on the $x$-axis with fill mapped to $-\log_{10}(P_{\text{adj}})$ using the Viridis Plasma color scale.
3. **Magma Bubble Plots (`plot_ora_bubble`)**:
   Displays the top 20 terms plotting **Rich Factor** on the $x$-axis, circle diameter proportional to gene overlap count, and color mapped to $-\log_{10}(P_{\text{adj}})$ using the Viridis Magma scale.

#### 3.5.2 GSEA Visualizations
1. **Triple-Panel GSEA Running Score Graphs (`gseaplot2`)**:
   - **Upper Panel**: The running enrichment score curve ($P_{\text{hit}} - P_{\text{miss}}$), identifying the peak score ($\text{ES}$).
   - **Middle Panel**: Vertical tick marks indicating the positions of gene set members within the ranked list, demarcating the leading-edge core enrichment window.
   - **Lower Panel**: The rank metric profile ($s_g$) across the entire transcriptome.
2. **Ridgeplots (`ridgeplot`)**:
   Visualizes the frequency density distribution of the ranked statistic ($s_g$) for all constituent genes across the top enriched pathways, illustrating whether the term is driven by uniform moderate shifts or bimodal gene subsets.
3. **Phenotypic Split Dot Plots**:
   Separates pathways into distinct plots based on activation status: **Upregulated in Tumor ($\text{NES} > 0$)** versus **Downregulated in Tumor ($\text{NES} < 0$)**, preventing opposing biological processes from obscuring functional interpretation.

---

## 4. Synthesis of Methodological Parameters and Comparative Framework

```
Table 2: Operational Parameters Across Enrichment Analysis Protocols
┌──────────────────────────────────────┬─────────────────────────────┬─────────────────────────────┐
│ Operational Parameter / Setting      │ ORA Workflow (ora_edgeR.R)  │ GSEA Workflow (ora_edgeR.R) │
├──────────────────────────────────────┼─────────────────────────────┼─────────────────────────────┤
│ Input Gene Set                       │ Significant PIKK Candidates │ Complete Ranked Background  │
│ Ranking Statistic Formula            │ N/A (Discrete binary set)   │ s_g = sign(log2FC) × sqrt(F)│
│ Ranking Statistic Interpretation     │ N/A                         │ Signed Quasi-Likelihood z   │
│ Statistical Hypothesis Model         │ Hypergeometric Distribution │ Weighted Random Walk (KS)   │
│ Weighting Parameter (p)              │ N/A                         │ p = 1.0 (Linear weighting)  │
│ Minimum Gene Set Size (minGSSize)    │ N/A                         │ 10 genes                    │
│ Maximum Gene Set Size (maxGSSize)    │ N/A                         │ 500 genes                   │
│ Background Universe Constraint       │ 16,905 tested edgeR genes   │ 16,905 tested edgeR genes   │
│ Cross-Mapping Database               │ org.Hs.eg.db (NCBI Entrez)  │ org.Hs.eg.db (NCBI Entrez)  │
│ Ontology Databases Queried           │ GO (BP, MF, CC), KEGG,      │ GO (BP), KEGG, Reactome     │
│                                      │ Reactome                    │                             │
│ Multi-Testing Correction             │ Benjamini-Hochberg (BH) FDR │ Benjamini-Hochberg (BH) FDR │
│ Significance Cutoff                  │ p.adjust ≤ 0.05, q ≤ 0.20   │ p.adjust ≤ 0.05 (FDR)       │
│ Directional Stratification           │ Pre-stratified subsets      │ Explicit NES > 0 vs NES < 0 │
│ Primary Output Graphics              │ Dot, Bar, Magma Bubble      │ gseaplot2, Ridge, Split Dot │
└──────────────────────────────────────┴─────────────────────────────┴─────────────────────────────┘
```

---

## 5. Ready-to-Publish Methods Section (Manuscript Reference Text)

> *Below is a formal, publication-ready draft formatted in standard third-person past tense, suitable for direct incorporation into the **Materials and Methods** section of an academic manuscript.*

### 2.X Functional Annotation and Pathway Over-Representation Analysis (ORA)
To delineate the biological processes, molecular functions, and signaling pathways significantly enriched within the prioritized PIKK-response signature, Over-Representation Analysis (ORA) was performed using the `clusterProfiler` Bioconductor package (`ora_edgeR.R`) [2]. Gene symbols were converted to unique NCBI Entrez Gene Identifiers using the human genome annotation database `org.Hs.eg.db`. To prevent false-positive inflation associated with default whole-genome backgrounds, the statistical background universe was restricted strictly to the 16,905 protein-coding genes robustly expressed and statistically evaluated within our differential expression model [1].

Statistical over-representation was evaluated across five reference knowledge bases: the Gene Ontology (GO) Biological Process (BP), Molecular Function (MF), and Cellular Component (CC) ontologies, the Kyoto Encyclopedia of Genes and Genomes (KEGG; organism = "hsa") pathway repository [4], and the Reactome Knowledgebase (organism = "human") [5]. Overlap significance was calculated using the one-tailed hypergeometric test (Fisher's exact test). Multi-testing correction was applied across all terms using the Benjamini–Hochberg False Discovery Rate (FDR) method [6]. Ontological categories satisfying an adjusted $p$-value ($P_{\text{adj}}$) $\le 0.05$ and a $q$-value $\le 0.20$ were designated as significantly enriched. Enrichment magnitude was quantified using Fold Enrichment and Rich Factor metrics, visualised via high-resolution dot plots, fold-enrichment bar charts, and magma bubble plots.

### 2.X Unbiased Genome-Wide Gene Set Enrichment Analysis (GSEA)
To complement discrete over-representation testing and capture coordinated, genome-wide pathway shifts across the full dynamic range of the transcriptome, Gene Set Enrichment Analysis (GSEA) was performed (`gseGO`, `gseKEGG`, and `gsePathway`) [2, 3]. To integrate biological effect size with statistical precision, an empirical signed ranking metric ($s_g$) was derived from the edgeR quasi-likelihood output for each gene:
$$s_g = \text{sign}(\log_2\text{FC}_g) \times \sqrt{F_g}$$
Where $\log_2\text{FC}_g$ represents the model log-fold change and $F_g$ represents the Quasi-Likelihood $F$-statistic, yielding a signed test statistic mathematically equivalent to the signed Wald $z$-score. Multi-mapped Entrez identifiers were resolved by selecting the maximum absolute ranking statistic, producing a sorted, non-redundant ranked transcriptome vector.

Weighted Kolmogorov–Smirnov random walk statistics were computed with a weighting exponent of $p = 1.0$. Gene set boundaries were constrained to biological modules containing between 10 and 500 genes (`minGSSize = 10`, `maxGSSize = 500`), mitigating statistical bias from overly broad or idiosyncratic terms. Empirical $p$-values were derived from 1,000 permutations and adjusted using the Benjamini–Hochberg procedure. Significant pathways ($\text{FDR} \le 0.05$) were partitioned into phenotypically activated programs ($\text{NES} > 0$) and repressed programs ($\text{NES} < 0$). Leading-edge (core enrichment) genes accounting for the peak enrichment score were identified. Pathway architectures were visualised using three-panel GSEA running score graphs (`gseaplot2`), expression density ridgeplots (`enrichplot`), and direction-stratified dot plots.

---

## 6. Comprehensive Software, Dependency, and Environment Manifest

All functional enrichment procedures were executed in the R statistical computing environment (R version $\ge 4.2.0$). Below is the inventory of software dependencies:

```
Table 3: Software Dependencies and Computational Infrastructure
┌─────────────────┬──────────────┬────────────────────────────────────────────────────────┐
│ Software/Package│ Origin       │ Primary Methodological Role                            │
├─────────────────┼──────────────┼────────────────────────────────────────────────────────┤
│ clusterProfiler │ Bioconductor │ Core statistical engine for ORA, GSEA, and dotplots [2]│
│ org.Hs.eg.db    │ Bioconductor │ Genome-wide annotation mapping for Homo sapiens        │
│ ReactomePA      │ Bioconductor │ Reactome pathway over-representation and GSEA [5]      │
│ enrichplot      │ Bioconductor │ Visualization suite for gseaplot2, ridgeplots, dotplots│
│ ggplot2         │ CRAN         │ High-resolution vector graphic and bubble plot engine  │
│ viridis         │ CRAN         │ Colorblind-friendly scientific palettes (Magma, Plasma)│
│ stringr         │ CRAN         │ Automated text wrapping and sanitization for plot axes │
│ dplyr / tibble  │ CRAN         │ Data frame manipulation, sorting, and column filtering │
│ readr           │ CRAN         │ High-performance TSV table serialization               │
│ cairo_pdf       │ Base system  │ Cairo-based anti-aliased vector graphic PDF rendering  │
└─────────────────┴──────────────┴────────────────────────────────────────────────────────┘
```

---

## 7. Peer-Reviewed Scientific Bibliography

1. **Huang DW, Sherman BT, Lempicki RA.** (2009). Bioinformatics enrichment tools: paths toward the comprehensive functional analysis of large gene lists. *Nucleic Acids Research*, 37(1):1–13. DOI: [10.1093/nar/gkn923](https://doi.org/10.1093/nar/gkn923).
2. **Yu G, Wang LG, Han Y, He QY.** (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS: A Journal of Integrative Biology*, 16(5):284–287. DOI: [10.1089/omi.2011.0118](https://doi.org/10.1089/omi.2011.0118).
3. **Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, Mesirov JP.** (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proceedings of the National Academy of Sciences USA*, 102(43):15545–15550. DOI: [10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102).
4. **Kanehisa M, Goto S.** (2000). KEGG: kyoto encyclopedia of genes and genomes. *Nucleic Acids Research*, 28(1):27–30. DOI: [10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27).
5. **Yu G, He QY.** (2016). ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. *Molecular BioSystems*, 12(2):477–479. DOI: [10.1039/c5mb00663e](https://doi.org/10.1039/c5mb00663e).
6. **Benjamini Y, Hochberg Y.** (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B (Methodological)*, 57(1):289–300. DOI: [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).
7. **The Gene Ontology Consortium.** (2021). The Gene Ontology resource: enriching a GOld mine. *Nucleic Acids Research*, 49(D1):D325–D334. DOI: [10.1093/nar/gkaa1113](https://doi.org/10.1093/nar/gkaa1113).
8. **Gillespie M, Jassal B, Stephan R, Milacic M, Rothfels K, Senff-Ribeiro A, Shamovsky V, Sousa da Silva C, Fabregat A, Hermjakob H, D'Eustachio P.** (2022). The reactome pathway knowledgebase 2022. *Nucleic Acids Research*, 50(D1):D687–D692. DOI: [10.1093/nar/gkab1028](https://doi.org/10.1093/nar/gkab1028).
9. **Chen Y, Lun ATL, Smyth GK.** (2016). From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. *F1000Research*, 5:1438. DOI: [10.12688/f1000research.8987.2](https://doi.org/10.12688/f1000research.8987.2).
10. **Blackford AN, Jackson SP.** (2017). ATM, ATR, and DNA-PK: The Trinity at the Heart of the DNA Damage Response. *Molecular Cell*, 66(6):801–817. DOI: [10.1016/j.molcel.2017.05.015](https://doi.org/10.1016/j.molcel.2017.05.015).
