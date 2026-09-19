# Protein-Protein Interaction (PPI) Network Analysis and Hub Gene Identification (Section 2.4) — Comprehensive Methodology Guide

**Project**: PIKK Pathway-Associated Biomarker Discovery and Oncogenic Characterization in Oral Squamous Cell Carcinoma (OSCC)  
**Target Manuscript Reference**: Materials and Methods / Bioinformatic Processing Pipeline  
**Primary Dataset**: The Cancer Genome Atlas Head and Neck Squamous Cell Carcinoma (TCGA-HNSC), Curated Oral Cavity Cohort ($n = 240$ primary samples; $n = 224$ Primary Solid Tumor, $n = 16$ Solid Tissue Normal)  
**Input Gene Spaces**: Significant PIKK-Associated Differentially Expressed Genes (DEGs) from edgeR QLF, Curated PIKK Gene Universe ($n = 204$ genes), and Normalized Expression Space ($\log_2\text{CPM}$)  
**Interaction Data Repository**: Search Tool for the Retrieval of Interacting Genes/Proteins (STRING database, v12.0; species = 9606, *Homo sapiens*; High Confidence Interaction Score $\ge 0.700$)  
**Scripts Evaluated**:
1. `05_ppi_network_analysis.R` — *High-Confidence PPI Network Construction, Multi-Centrality Topological Characterization, Sub-Network Modular Decomposition, and Louvain Community Detection*
2. `06_hub_gene_identification.R` — *Tri-Metric Composite Hub Scoring, Non-Parametric Clinical Covariate Association Testing, and Biological Pathway Role Classification*

---

## 1. Executive Summary & End-to-End Workflow

Differential expression and functional enrichment analyses identify altered transcripts and biological terms, but cellular behaviors are governed by physical macromolecular protein assemblies and interconnected biochemical networks. In oncogenesis, cancer cells systematically rewire protein-protein interaction (PPI) networks to tolerate genomic instability, escape cell cycle checkpoints, and fuel autonomous growth. Identifying topologically central "hub" genes—nodes that command high connectivity or occupy critical communication bottlenecks—uncovers master regulators of the oncogenic state that represent primary prognostic biomarkers and synthetic lethal vulnerabilities.

This document formalizes the publication-grade methodology governing the **Protein-Protein Interaction (PPI) Network Topology and Multi-Dimensional Hub Gene Identification pipeline** for the OSCC PIKK characterization project.

The analytical architecture comprises two sequential, highly integrated computational phases:
1. **Network Topology & Modular Community Analysis (`05_ppi_network_analysis.R`)**: Retrieves high-confidence physical and functional interactome data from STRING (score $\ge 0.700$), builds unweighted and weighted undirected graph structures in `igraph`, quantifies five orthogonal mathematical centrality metrics, extracts nested functional sub-networks (first-order interactome and DDR-specific module), and detects densely interconnected communities via Louvain modularity optimization.
2. **Composite Hub Gene Discovery & Clinical Association Modeling (`06_hub_gene_identification.R`)**: Prioritizes key regulator genes using a standardized tri-metric composite score integrating degree centrality, betweenness centrality, and biological expression magnitude ($|\log_2\text{FC}|$), systematically interrogates correlations with clinical covariates (AJCC pathologic tumor stage, vital status, age, gender) using non-parametric statistics, and categorizes candidates into defined functional pathway roles.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                                 UPSTREAM INPUT ARTIFACTS                               │
│  • Significant PIKK DEGs: OSCC_PIKK_extended_edgeR_DEG_intersection_significant.tsv    │
│  • Curated Reference Space: PIKK_related_gene_universe.csv (204 genes)                 │
│  • Normalized Expression Matrix: logCPM_matrix.tsv (Tumor = 224, Normal = 16)          │
│  • Clinical Metadata: matched_metadata_used.tsv (AJCC Stage, Vital Status, Age, Gender)│
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 1: INTERACTION RETRIEVAL & ALIAS HARMONIZATION (05_ppi_network_analysis.R)      │
│  • Query STRING v12.0 REST API (species = 9606, Homo sapiens; caller = OSCC_PIKK_DDR) │
│  • High-confidence filtering: combined interaction score ≥ 0.700                       │
│  • Gene symbol alias harmonization (e.g., H2AX ↔ H2AFX, MRE11 ↔ MRE11A)                │
│  • Edge deduplication, undirected pair canonicalization, and self-loop removal         │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 2: GRAPH CONSTRUCTION & TOPOLOGICAL METRIC COMPUTATION (05_ppi_network_analysis)│
│  • Multi-scale graph instantiation in igraph:                                          │
│      - Core Network: DEGs only (both endpoints significant)                            │
│      - First-Order Sub-network: DEGs + PIKK-universe interactors                       │
│      - DDR-Focused Sub-network: ATR, ATM, and PRKDC signaling modules                  │
│  • Mathematical computation of 5 centrality metrics:                                   │
│      - Degree, Betweenness, Closeness, Eigenvector, and Local Clustering Coefficient   │
│  • Node attribute embedding: log2FC, FDR, regulation status, PIKK family group         │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 3: MODULAR COMMUNITY DETECTION (05_ppi_network_analysis.R)                      │
│  • Louvain modularity optimization algorithm (seed = 42)                               │
│  • Two-phase greedy modularity (Q) maximization                                        │
│  • Identification of dominant PIKK kinase family per topological module                │
│  • Serialization of Cytoscape-ready edge/node TSVs and GraphML containers              │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 4: TRI-METRIC COMPOSITE HUB GENE PRIORITIZATION (06_hub_gene_identification.R)  │
│  • Top-N rank selection: Top 15 by Degree ∪ Top 15 by Betweenness                      │
│  • Multi-dimensional standardization: Z-score transformation of Degree, Betweenness,   │
│    and absolute effect magnitude (|log2FC|)                                            │
│  • Composite Hub Score formulation: Composite = Z(deg) + Z(btwn) + Z(|log2FC|)         │
│  • Unbiased ranking of master transcriptional and topological hubs                     │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 5: CLINICAL COVARIATE ASSOCIATION & ROLE CLASSIFICATION (06_hub_gene_ident.)    │
│  • Extraction of hub expression across matched primary tumor cohort (n = 224)          │
│  • Non-parametric hypothesis testing:                                                  │
│      - AJCC Pathologic Stage (I–IV): Kruskal-Wallis rank-sum test                      │
│      - Overall Vital Status (Alive vs Dead): Wilcoxon rank-sum test                    │
│      - Patient Gender (Male vs Female): Wilcoxon rank-sum test                         │
│      - Patient Age at Diagnosis: Spearman rank-order correlation                       │
│  • Multi-testing correction via Benjamini-Hochberg FDR within each clinical variable   │
│  • Rule-based categorical classification into 8 distinct biological pathway roles      │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
                                            ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  PHASE 6: NETWORK VISUALIZATION & DOWNSTREAM BRIDGING                                  │
│  • High-resolution vector network layouts using ggraph and tidygraph                   │
│  • Export of Master Hub Characterization Table for clinical survival modeling          │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. Biological Context and Study Rationale

### 2.1 The Centrality-Lethality Rule and Biological Network Resilience in OSCC
In structural network biology, the **Centrality-Lethality Rule** (Jeong et al., *Nature* 2001 [1]) posits that the phenotypic importance of a protein correlates with its topological position within the interactome. Highly connected proteins ("hubs") and critical path intermediaries ("bottlenecks") are far more likely to be biologically essential than peripheral nodes. 

In Oral Squamous Cell Carcinoma (OSCC), carcinogen-induced damage generates profound replication stress and double-strand breaks. Malignant mucosal epithelial cells must continuously coordinate DNA repair, checkpoint surveillance, and metabolic resource allocation to survive. This adaptation is orchestrated by the **PIKK kinase superfamily** ($ATM$, $ATR$, $PRKDC$, $MTOR$, $TRRAP$, and $SMG1$).

Rather than acting as isolated enzymes, PIKK family members function as scaffolds and central coordinators within multi-protein macromolecular complexes:
- **ATM & PRKDC**: Organize high-order repair hubs at DNA double-strand breaks, interacting physically with the MRN complex ($MRE11$-$RAD50$-$NBN$) and Ku heterodimer ($XRCC5$/$XRCC6$) to dictate pathway choice between homologous recombination (HR) and non-homologous end joining (NHEJ).
- **ATR**: Coordinates with $ATRIP$, the 9-1-1 clamp ($RAD9A$-$RAD1$-$HUS1$), and $TOPBP1$ to stabilize stalled replication forks and prevent fork collapse.
- **MTOR**: Operates within multi-subunit TORC1 ($RPTOR$, $MLST8$) and TORC2 ($RICTOR$, $MAPKAP1$) signaling complexes, serving as an information integration hub linking DNA damage to metabolic capacity.
- **TRRAP**: Acts as a massive pseudokinase scaffold within the TIP60/NuA4 and SAGA complexes, recruiting histone acetyltransferases ($KAT5$, $KAT2A$) to open chromatin structure at DNA damage sites and active transcriptional promoters.
- **SMG1**: Forms the SURF and DECID surveillance complexes with $UPF1$, $UPF2$, $SMG7$, and $SMG9$ to execute nonsense-mediated mRNA decay (NMD) and safeguard transcriptome fidelity.

Identifying which PIKK components and functional interactors serve as network hubs provides direct insight into the structural dependencies sustaining the malignant state in OSCC.

### 2.2 Methodological Rationale: Why Tri-Metric Composite Scoring?
Conventional bioinformatic workflows frequently identify hub genes based on a single parameter, such as degree centrality alone ($k > 10$) or differential expression fold change alone. These simplistic approaches suffer from well-documented systematic biases:
1. **Degree Centrality Bias**: Prioritizes "sticky" or promiscuous structural proteins (e.g., chaperones, ribosomal subunits) that interact non-specifically with hundreds of partners, regardless of whether their expression is altered in the disease state.
2. **Betweenness Centrality Bias**: Identifies inter-modular bridges that may connect disparate cellular processes but lack dense local connectivity or robust transcriptional regulation.
3. **Expression Fold Change Bias**: Identifies highly responsive downstream transcripts (e.g., secreted cytokines, structural keratins) that have no regulatory authority over upstream oncogenic signaling.

To overcome these limitations, the methodology implements a **Tri-Metric Composite Hub Scoring Framework**:
$$\text{Composite Score} = Z(\text{Degree}) + Z(\text{Betweenness}) + Z(|\log_2\text{FC}|)$$
This mathematical formulation guarantees that prioritized hub genes fulfill three biological criteria:
- **High Local Connectivity ($Z_{\text{deg}}$)**: Regulates or physically complexes with numerous direct interaction partners.
- **High Global Communication Authority ($Z_{\text{btwn}}$)**: Controls the flow of molecular signals between distinct functional modules as a key network bottleneck.
- **Robust Biological Modulation ($Z_{\text{FC}}$)**: Undergoes significant, non-random transcriptional dysregulation between malignant and healthy tissues.

---

## 3. Step-by-Step Computational Methodology

### 3.1 Step 1: Protein-Protein Interaction Ingestion & High-Confidence Filtering
**Script**: `05_ppi_network_analysis.R`  
**Input Artifacts**:
- `pikk_deg_file`: Significant PIKK-associated DEGs (`OSCC_PIKK_extended_edgeR_DEG_intersection_significant.tsv`).
- `pikk_univ_file`: Curated PIKK gene universe (`PIKK_related_gene_universe.csv`).
**Output Artifacts**:
- `STRING_interactions_all_PIKK_universe.tsv`: Master table of high-confidence interactome edges.

#### 3.1.1 Interactome Knowledgebase Interrogation
Protein-protein interactions were retrieved from the **STRING database (Search Tool for the Retrieval of Interacting Genes/Proteins, version 12.0)** (Szklarczyk et al., *Nucleic Acids Res* 2023 [2]) via the STRING REST API using species identifier 9606 (*Homo sapiens*):
$$\text{Query URL: } \texttt{https://string-db.org/api/tsv/network?identifiers=\{genes\}\&species=9606\&required\_score=700}$$

#### 3.1.2 High-Confidence Score Thresholding
The STRING database calculates a probabilistic combined interaction score ($S_{\text{combined}} \in [0, 1000]$) by integrating seven independent evidentiary channels:
$$S_{\text{combined}} = 1 - \prod_{c \in \mathcal{C}} (1 - S_c)$$
Where $\mathcal{C} = \{\text{neighborhood}, \text{fusion}, \text{co-occurrence}, \text{co-expression}, \text{experimental}, \text{database}, \text{text-mining}\}$.

To eliminate false-positive associations and computationally noisy text-mining artifacts, an empirical **high-confidence threshold** is enforced:
$$S_{\text{combined}} \ge 700 \quad (\ge 0.700 \text{ on a unit scale})$$
Interactions satisfying this threshold possess strong biochemical support from physical binding assays, crystallographic complexes, or curated pathway databases.

#### 3.1.3 Alias Harmonization & Edge Sanitization
1. **Identifier Harmonization**: Gene symbols diverging between modern HGNC designations and STRING canonical identifiers are remapped using an automated lookup dictionary:
   $$H2AX \longrightarrow H2AFX, \quad MRE11 \longrightarrow MRE11A$$
2. **Self-Loop Removal**: Autocrine or homodimeric self-interactions ($v_i = v_j$) are excluded:
   $$\mathcal{E}_{\text{clean}} = \{ (u, v) \in \mathcal{E} \mid u \ne v \}$$
3. **Canonical Undirected Edge Deduplication**:
   Because STRING exports bi-directional directed entries ($A \to B$ and $B \to A$), edges are sorted lexicographically:
   $$\text{edge\_key} = \min(u, v) \mathbin{\Vert} \text{"::"} \mathbin{\Vert} \max(u, v)$$
   Unique undirected edges are preserved:
   $$\mathcal{E}_{\text{undirected}} = \text{distinct}(\mathcal{E}_{\text{clean}}, \text{edge\_key})$$

---

### 3.2 Step 2: Graph Representation & Topological Centrality Quantitation
**Script**: `05_ppi_network_analysis.R`  
**Input Artifacts**: Cleaned STRING edge list and node annotation tables.  
**Output Artifacts**:
- `network_topology_metrics.tsv`: Comprehensive topology table reporting all 5 centrality metrics for every node.
- Graph objects: `g_core` (DEGs only), `g_sub` (First-order interactome), and `g_ddr` (DDR-specific modules).

#### 3.2.1 Graph Formulation
Let $\mathcal{G} = (\mathcal{V}, \mathcal{E})$ denote an unweighted, undirected graph where $\mathcal{V}$ is the set of $N = |\mathcal{V}|$ protein nodes, and $\mathcal{E}$ is the set of $M = |\mathcal{E}|$ physical interaction edges. The network structure is encoded by the symmetric adjacency matrix $\mathbf{A} \in \{0, 1\}^{N \times N}$:
$$A_{uv} = \begin{cases} 1 & \text{if } (u, v) \in \mathcal{E} \\ 0 & \text{otherwise} \end{cases}$$

#### 3.2.2 Mathematical Formulations of Topological Centrality Metrics

##### 1. Degree Centrality ($k_v$)
Quantifies the direct local connectivity of node $v$, reflecting its immediate interaction density:
$$k_v = \sum_{u \in \mathcal{V} \setminus \{v\}} A_{vu} = |\mathcal{N}(v)|$$
Where $\mathcal{N}(v) = \{u \in \mathcal{V} \mid (v, u) \in \mathcal{E}\}$ represents the set of direct neighbors of node $v$.

##### 2. Betweenness Centrality ($C_B(v)$)
Measures the frequency with which node $v$ falls on the shortest path (geodesic) between all other node pairs, quantifying its role as a network communication bottleneck:
$$C_B(v) = \sum_{\substack{s \ne v \ne t \\ s, t \in \mathcal{V}}} \frac{\sigma_{st}(v)}{\sigma_{st}}$$
Where:
- $\sigma_{st}$ is the total number of shortest paths connecting node $s$ to node $t$.
- $\sigma_{st}(v)$ is the number of those shortest paths that pass through node $v$.
To allow direct comparison across networks of differing sizes, normalized betweenness ($C_B'(v)$) scales $C_B(v)$ by the maximum theoretical paths in an undirected graph:
$$C_B'(v) = \frac{2 \cdot C_B(v)}{(N - 1)(N - 2)}$$

##### 3. Closeness Centrality ($C_C(v)$)
Quantifies the reciprocal of the sum of shortest-path distances from node $v$ to all other nodes in the network, measuring how rapidly information propagates from node $v$ to the rest of the interactome:
$$C_C(v) = \frac{N - 1}{\sum_{u \in \mathcal{V} \setminus \{v\}} d(v, u)}$$
Where $d(v, u)$ is the shortest topological distance (number of edges) between node $v$ and node $u$.

##### 4. Eigenvector Centrality ($x_v$)
Measures the influence of a node by accounting not only for its immediate degree, but also for the centrality of its neighbors. A node connected to a few highly central hubs achieves a higher score than one connected to many peripheral leaves:
$$x_v = \frac{1}{\lambda} \sum_{u \in \mathcal{N}(v)} x_u = \frac{1}{\lambda} \sum_{u \in \mathcal{V}} A_{vu} x_u$$
In matrix notation:
$$\mathbf{A} \mathbf{x} = \lambda \mathbf{x}$$
By the Perron–Frobenius theorem, eigenvector centrality corresponds to the principal eigenvector $\mathbf{x}$ associated with the largest positive eigenvalue $\lambda_{\max}$ of $\mathbf{A}$.

##### 5. Local Clustering Coefficient / Transitivity ($C_v$)
Measures the cliquishness of a node's local neighborhood—the probability that two neighbors of node $v$ are also connected to each other (Watts & Strogatz, *Nature* 1998 [3]):
$$C_v = \frac{2 e_v}{k_v (k_v - 1)} = \frac{\sum_{j, k \in \mathcal{V}} A_{vj} A_{jk} A_{kv}}{k_v (k_v - 1)}$$
Where $e_v$ is the actual number of edges existing among the neighbors of $v$, and $k_v (k_v - 1) / 2$ is the maximum possible edges. For isolated nodes or nodes with $k_v \le 1$, $C_v$ is formally defined as 0.

---

### 3.3 Step 3: Modular Community Detection via the Louvain Modularity Algorithm
**Script**: `05_ppi_network_analysis.R`  
**Input Artifacts**: Core network graph `g_core`.  
**Output Artifacts**: `community_assignments.tsv`, `PPI_core_network.graphml`, Cytoscape edge/node tables.

#### 3.3.1 Mathematical Formulation of Modularity ($Q$)
Biological networks are partitioned into modular functional sub-units (communities) where intra-module edge density significantly exceeds inter-module density. Community detection was executed using the **Louvain modularity maximization algorithm** (Blondel et al., *J Stat Mech* 2008 [4]) with a fixed random seed (`set.seed(42)`):
$$Q = \frac{1}{2m} \sum_{i, j \in \mathcal{V}} \left[ A_{ij} - \frac{k_i k_j}{2m} \right] \delta(c_i, c_j)$$
Where:
- $m = |\mathcal{E}|$ is the total number of edges in the network ($2m = \sum_i k_i$).
- $A_{ij}$ is the observed adjacency between nodes $i$ and $j$.
- $\frac{k_i k_j}{2m}$ represents the expected number of edges between nodes $i$ and $j$ in a random null network preserving the empirical degree distribution (Newman-Girvan null model).
- $c_i$ denotes the community assignment of node $i$.
- $\delta(c_i, c_j)$ is the Kronecker delta:
  $$\delta(c_i, c_j) = \begin{cases} 1 & \text{if } c_i = c_j \text{ (same community)} \\ 0 & \text{otherwise} \end{cases}$$

#### 3.3.2 Two-Phase Algorithmic Iteration
1. **Phase 1 (Local Greedy Optimization)**: Each node begins in its own distinct community. For each node $i$, the algorithm computes the gain in modularity ($\Delta Q$) obtained by moving $i$ into the community of each neighbor $j$:
   $$\Delta Q = \left[ \frac{\Sigma_{\text{in}} + 2 k_{i,\text{in}}}{2m} - \left( \frac{\Sigma_{\text{tot}} + k_i}{2m} \right)^2 \right] - \left[ \frac{\Sigma_{\text{in}}}{2m} - \left( \frac{\Sigma_{\text{tot}}}{2m} \right)^2 - \left( \frac{k_i}{2m} \right)^2 \right]$$
   Where $\Sigma_{\text{in}}$ is the sum of edge weights within community $C$, $\Sigma_{\text{tot}}$ is the total edge weight incident to nodes in $C$, $k_i$ is the degree of node $i$, and $k_{i,\text{in}}$ is the edge weight from $i$ to community $C$. Node $i$ is relocated to the community providing maximal $\Delta Q > 0$. This process iterates until local convergence is reached.
2. **Phase 2 (Community Aggregation)**: A new coarse-grained network is constructed whose vertices represent the communities discovered in Phase 1. Edges between communities represent the sum of weights connecting constituent nodes, and internal edges become self-loops.
3. Phases 1 and 2 repeat hierarchically until no further modularity increase is mathematically achievable.

#### 3.3.3 Modular Characterization
For each discovered community, the **dominant PIKK kinase group** is determined by majority representation, and the **top local hub** is identified by maximal within-community degree:
$$\text{Dominant Group}(C) = \arg\max_{G \in \{\text{ATR, ATM, PRKDC, MTOR, TRRAP, SMG1}\}} \sum_{v \in C} \mathbb{I}(\text{PIKK\_Group}_v = G)$$

---

### 3.4 Step 4: Tri-Metric Composite Hub Scoring and Candidate Prioritization
**Script**: `06_hub_gene_identification.R`  
**Input Artifacts**: `network_topology_metrics.tsv`  
**Output Artifacts**:
- `hub_genes_top_degree.tsv`: Top 15 ranked genes by Degree centrality.
- `hub_genes_top_betweenness.tsv`: Top 15 ranked genes by Betweenness centrality.
- `hub_genes_union_composite.tsv`: Final prioritized union set ranked by composite hub score.

#### 3.4.1 Union Hub Set Formulation
To prevent bias toward connectivity alone or bottleneck routing alone, candidates are drawn from the set-theoretic union of the upper extremes of both distributions ($N_{\text{degree}} = 15$, $N_{\text{betweenness}} = 15$):
$$\mathcal{H}_{\text{deg}} = \arg\max_{\mathcal{S} \subset \mathcal{V}, |\mathcal{S}|=15} \sum_{v \in \mathcal{S}} k_v$$
$$\mathcal{H}_{\text{btwn}} = \arg\max_{\mathcal{S} \subset \mathcal{V}, |\mathcal{S}|=15} \sum_{v \in \mathcal{S}} C_B'(v)$$
$$\mathcal{H}_{\text{union}} = \mathcal{H}_{\text{deg}} \cup \mathcal{H}_{\text{btwn}}$$

#### 3.4.2 Standardization and Composite Score Formulation
Because Degree, Betweenness, and Fold Change operate on completely divergent numerical scales and variance distributions, raw metrics cannot be directly summed. Each metric is standardized to a standard normal distribution ($Z$-score) across all candidate nodes in $\mathcal{H}_{\text{union}}$:
$$Z_{\text{deg}}(v) = \frac{k_v - \bar{k}}{s_k}$$
$$Z_{\text{btwn}}(v) = \frac{C_B'(v) - \bar{C}_B'}{s_{C_B'}}$$
$$Z_{\text{FC}}(v) = \frac{|\log_2\text{FC}_v| - \overline{|\log_2\text{FC}|}}{s_{|\log_2\text{FC}|}}$$
Where $\bar{x}$ and $s_x$ represent the sample mean and standard deviation of metric $x$ across $\mathcal{H}_{\text{union}}$.

The **Composite Hub Score** ($S_{\text{composite}}(v)$) is defined as the unweighted linear sum:
$$S_{\text{composite}}(v) = Z_{\text{deg}}(v) + Z_{\text{btwn}}(v) + Z_{\text{FC}}(v)$$
This score awards top ranks to nodes that simultaneously maximize topological degree, bottleneck control, and transcriptional effect magnitude.

---

### 3.5 Step 5: Non-Parametric Clinical Covariate Association Testing
**Script**: `06_hub_gene_identification.R`  
**Input Artifacts**: Matched tumor expression matrix (`logCPM_matrix.tsv`, restricted to $n = 224$ primary oral tumors) and clinical annotations (`matched_metadata_used.tsv`).  
**Output Artifacts**: `hub_gene_clinical_correlation_tests.tsv`.

#### 3.5.1 Biological and Methodological Rationale
Clinical attributes in observational oncology cohorts (e.g., AJCC pathologic tumor stage, overall vital status, patient age, gender) frequently display non-normal distributions, ordinal staging boundaries, and sample size imbalances. Applying standard parametric models (e.g., Pearson correlation or Student's $t$-test) violates distributional assumptions. 

Consequently, the methodology deploys **non-parametric statistical hypothesis tests** that operate on rank transformations of variance-stabilized expression ($\log_2\text{CPM}$).

#### 3.5.2 Mathematical Formulation of Clinical Hypothesis Tests

##### 1. AJCC Pathologic Tumor Stage (Ordinal / Multi-Group Factor)
To evaluate whether hub gene expression correlates with clinical tumor progression, pathologic stages are harmonized into four ordered biological categories: Stage I, Stage II, Stage III, and Stage IV (pooling sub-stages IVA, IVB, and IVC).

Expression differences across the $k = 4$ independent clinical stage groups are evaluated using the **Kruskal–Wallis Rank-Sum Test** (Kruskal & Wallis, 1952 [5]):
$$H = \frac{12}{N(N + 1)} \sum_{j=1}^k \frac{R_j^2}{n_j} - 3(N + 1)$$
Where:
- $N$ is the total number of tumor samples with documented clinical stage.
- $k$ is the number of stage levels ($k = 4$).
- $n_j$ is the number of samples belonging to stage $j$.
- $R_j$ is the sum of ranks for stage $j$ after pooling all $N$ observations and ranking them in ascending order.
Under the null hypothesis that the expression medians across all stages are identical, $H$ is asymptotically distributed as a chi-square random variable with $k - 1$ degrees of freedom:
$$H \sim \chi^2(k - 1) = \chi^2(3)$$

##### 2. Overall Vital Status and Patient Gender (Binary Factors)
For binary clinical comparisons (Vital Status: *Alive* vs *Dead*; Gender: *Male* vs *Female*), expression differences are evaluated using the **Two-Sample Wilcoxon Rank-Sum Test (Mann–Whitney $U$ Test)** (Mann & Whitney, 1947 [6]):
$$U = W_1 - \frac{n_1(n_1 + 1)}{2}$$
Where $W_1$ is the sum of ranks assigned to group 1 (sample size $n_1$), and $n_2$ is the sample size of group 2. For larger sample sizes, the standardized $z$-statistic is computed:
$$z = \frac{U - \frac{n_1 n_2}{2}}{\sqrt{\frac{n_1 n_2 (n_1 + n_2 + 1)}{12}}}$$
Two-sided nominal $p$-values are computed from the standard normal cumulative distribution function $\Phi(z)$.

##### 3. Patient Age at Diagnosis (Continuous Variable)
To test for monotonic expression trends associated with patient aging without assuming linear scaling, we compute the non-parametric **Spearman Rank-Order Correlation Coefficient** ($\rho$):
$$\rho = 1 - \frac{6 \sum_{i=1}^n d_i^2}{n(n^2 - 1)}$$
Where $d_i = \text{rank}(y_i) - \text{rank}(\text{age}_i)$ is the difference between the expression rank and chronological age rank for patient $i$. Significance is evaluated using the exact $t$-transformation:
$$t = \rho \sqrt{\frac{n - 2}{1 - \rho^2}} \sim t(n - 2)$$

#### 3.5.3 Multiple Testing Correction and Dual Significance Reporting
To prevent false-positive reporting across the multiplicity of hub gene tests, nominal $p$-values are adjusted within each clinical covariate family using the Benjamini–Hochberg False Discovery Rate (FDR) procedure:
$$P_{\text{adjusted}(i)} = \min_{j \ge i} \left( \frac{M \cdot p_{(j)}}{j} \right)$$
Where $M$ is the number of hub genes tested against that specific clinical variable.

To provide transparency for translational manuscript drafting, results are annotated across two formal statistical thresholds:
1. **Nominal Significance**: $P_{\text{value}} < 0.05$ (exploratory biological threshold).
2. **FDR Significance**: $P_{\text{adjusted}} < 0.05$ (conservative, high-stringency threshold).

---

### 3.6 Step 6: Biological Pathway Role Classification
**Script**: `06_hub_gene_identification.R`  
**Input Artifacts**: Prioritized hub gene table.  
**Output Artifacts**: Final integrated hub gene characterization table (`hub_gene_characterisation_table.tsv`).

#### 3.6.1 Curated Biological Role Taxonomy
To contextualize the topological hub genes within cellular cancer hallmarks, every candidate gene is classified into an explicit biological pathway category using a rule-based hierarchy grounded in DDR and PIKK signaling:

```
Table 1: Functional Classification Scheme for PIKK Pathway Hub Genes
┌───────────────────────────────┬────────────────────────────────────────────────────────────────────────┐
│ Biological Pathway Role       │ Primary Biochemical & Cellular Hallmarks                               │
├───────────────────────────────┼────────────────────────────────────────────────────────────────────────┤
│ Core PIKK Kinase              │ The 6 primary catalytic enzymes (ATM, ATR, PRKDC, MTOR, TRRAP, SMG1)   │
├───────────────────────────────┼────────────────────────────────────────────────────────────────────────┤
│ Checkpoint Regulator          │ Cell cycle checkpoint transducers, clamp loaders, and DNA sensors      │
│                               │ (CHEK1, CHEK2, TP53, TP53BP1, MDC1, CLSPN, TOPBP1, RAD17, 9-1-1 clamp)│
├───────────────────────────────┼────────────────────────────────────────────────────────────────────────┤
│ Cell Cycle Regulator          │ Direct mitotic drivers, cyclin-dependent kinases, and replication      │
│                               │ initiators (PLK1, AURKA, AURKB, CCNE1, CDK2, CDC45, CDC7, CDKN2A)      │
├───────────────────────────────┼────────────────────────────────────────────────────────────────────────┤
│ DNA Repair Effector           │ Homologous recombination (HR), non-homologous end joining (NHEJ),      │
│                               │ Fanconi anemia, and excision repair (BRCA1/2, RAD51, PARP1, XRCC1/4/5) │
├───────────────────────────────┼────────────────────────────────────────────────────────────────────────┤
│ Signalling Effector           │ Downstream PI3K/AKT/mTOR cascades, nutrient-sensing G-proteins, and    │
│                               │ metabolic adaptors (RPTOR, RICTOR, DEPTOR, RHEB, EIF4EBP1, RPS6KB1)   │
├───────────────────────────────┼────────────────────────────────────────────────────────────────────────┤
│ Chromatin/Transcription       │ Histone acetyltransferases, SAGA/NuA4 subunits, and chromatin          │
│                               │ remodelers (KAT2A, KAT5, EP400, INO80, RUVBL1/2, TAF complex members)  │
├───────────────────────────────┼────────────────────────────────────────────────────────────────────────┤
│ mRNA Surveillance             │ Nonsense-mediated decay (NMD) machinery and exon junction complex      │
│                               │ components (UPF1, UPF2, UPF3A, SMG5, SMG6, SMG7, SMG8, SMG9, ETF1)    │
├───────────────────────────────┼────────────────────────────────────────────────────────────────────────┤
│ Other                         │ Peripheral or unassigned metabolic/signaling interactors               │
└───────────────────────────────┴────────────────────────────────────────────────────────────────────────┘
```

---

## 4. Synthesis of Methodological Parameters and Comparative Framework

```
Table 2: Comparative Methodological Parameters in Network Analysis and Hub Discovery
┌──────────────────────────────────────┬─────────────────────────────┬─────────────────────────────┐
│ Computational Parameter / Metric     │ Implemented Methodology     │ Conventional Alternative    │
├──────────────────────────────────────┼─────────────────────────────┼─────────────────────────────┤
│ PPI Database Source                  │ STRING v12.0 (REST API)     │ BioGRID / IntAct (raw)      │
│ Confidence Score Threshold           │ Combined Score ≥ 0.700      │ Medium Confidence (≥ 0.400) │
│ Interaction Types Captured           │ Physical & Functional       │ Binary physical binding only│
│ Hub Selection Strategy               │ Tri-Metric Composite Score  │ Degree Centrality cutoff    │
│ Centrality Normalization             │ Z-score standardization     │ Unscaled raw rank summing   │
│ Multi-Metric Component Weights       │ Equal weights (1:1:1)       │ Ad-hoc arbitrary weighting  │
│ Community Detection Algorithm        │ Louvain Modularity (Q)      │ Walktrap / Fast Greedy      │
│ Random Seed for Modularity           │ Fixed (seed = 42)           │ Unseeded (non-reproducible) │
│ Clinical Tumor Staging Model         │ Kruskal-Wallis (Ordinal)    │ Linear regression (ANOVA)   │
│ Binary Clinical Association Tests    │ Wilcoxon Rank-Sum (U-test)  │ Student's t-test            │
│ Continuous Age Correlation           │ Spearman Rank Correlation   │ Pearson linear correlation  │
│ Multiple Testing Correction          │ Benjamini-Hochberg FDR      │ Bonferroni (over-stringent) │
└──────────────────────────────────────┴─────────────────────────────┴─────────────────────────────┘
```

---

## 5. Ready-to-Publish Methods Section (Manuscript Reference Text)

> *Below is a formal, publication-ready draft formatted in standard third-person past tense, suitable for direct incorporation into the **Materials and Methods** section of an academic manuscript.*

### 2.X Protein-Protein Interaction (PPI) Network Construction and Centrality Analysis
To investigate the macromolecular organization and topological architecture of the prioritized PIKK signaling axis in Oral Squamous Cell Carcinoma (OSCC), protein-protein interaction (PPI) networks were constructed (`05_ppi_network_analysis.R`) using the `igraph` and `tidygraph` frameworks [7, 8]. Physical and functional interactions among significant PIKK-associated differentially expressed genes (DEGs) and the broader PIKK reference universe were retrieved from the STRING database (version 12.0; species = 9606, *Homo sapiens*) via the STRING REST API [2]. Interactions were filtered using an empirical high-confidence threshold of combined score $\ge 0.700$. Gene symbol discrepancies between HGNC designations and STRING canonical identifiers were harmonized using an automated alias mapping dictionary. Redundant undirected edges and homodimeric self-loops were excluded.

Network structures were represented as undirected graphs $\mathcal{G} = (\mathcal{V}, \mathcal{E})$. Five orthogonal topological centrality parameters were computed for every vertex in the interactome:
1. **Degree Centrality ($k_v$)**: Quantifying direct local connection density.
2. **Normalized Betweenness Centrality ($C_B'(v)$)**: Quantifying the proportion of all shortest paths traversing node $v$.
3. **Normalized Closeness Centrality ($C_C(v)$)**: Measuring the reciprocal of geodesic path distances.
4. **Eigenvector Centrality ($x_v$)**: Quantifying connection to neighboring topological hubs.
5. **Local Clustering Coefficient ($C_v$)**: Quantifying local neighbourhood cliquishness and transitivity [3].

Modular functional sub-networks, including the first-order extended interactome and the core DNA damage response (DDR) module ($ATM$, $ATR$, and $PRKDC$), were systematically extracted and mapped.

### 2.X Modular Community Detection via Louvain Modularity Maximization
Densely interconnected functional modules within the core interactome were identified using the Louvain modularity maximization algorithm (`cluster_louvain`) with a fixed initialization seed (`seed = 42`) [4]. The algorithm optimized the modularity objective function ($Q$), which compares observed intra-module edge density against the Newman-Girvan null model preserving empirical degree sequences. Dominant functional assignments for each community were determined by majority representation of specific PIKK kinase groups. Network topologies and community partitions were formatted for interactive exploration in Cytoscape [9] and exported as GraphML files.

### 2.X Tri-Metric Composite Hub Gene Prioritization
To prevent bias associated with evaluating single centrality parameters alone, master regulatory hub genes were identified using a standardized **Tri-Metric Composite Scoring Framework** (`06_hub_gene_identification.R`). The top 15 genes ranked by degree centrality and the top 15 genes ranked by normalized betweenness centrality were merged to form a union hub candidate set ($\mathcal{H}_{\text{union}}$). For each candidate gene, degree, betweenness, and absolute biological fold change magnitude ($|\log_2\text{FC}|$) were standardized to standard normal $Z$-scores across the candidate space:
$$S_{\text{composite}}(v) = Z(\text{Degree}_v) + Z(\text{Betweenness}_v) + Z(|\log_2\text{FC}_v|)$$
Candidates were ranked in descending order based on their composite hub score, integrating local physical connectivity, communication bottleneck control, and robust transcriptional dysregulation.

### 2.X Non-Parametric Clinical Covariate Association and Pathway Role Classification
To assess the clinical relevance of prioritized hub genes, expression levels ($\log_2\text{CPM}$) were interrogated across matched primary tumor specimens ($n = 224$) in relation to clinicopathologic variables. Given the non-normal distribution and ordinal boundaries of clinical parameters, non-parametric statistical hypothesis testing was systematically executed:
- **AJCC Pathologic Tumor Stage (Stages I–IV)**: Tested using the Kruskal–Wallis rank-sum test [5].
- **Overall Vital Status (Alive vs Dead) and Gender (Male vs Female)**: Tested using the two-sample Wilcoxon rank-sum test (Mann–Whitney $U$ test) [6].
- **Patient Age at Diagnosis**: Evaluated using the Spearman rank-order correlation coefficient ($\rho$).

Multi-testing adjustment was applied within each clinical covariate family using the Benjamini–Hochberg False Discovery Rate (FDR) procedure [10]. Both nominal significance ($P_{\text{value}} < 0.05$) and FDR-adjusted significance ($P_{\text{adjusted}} < 0.05$) were reported. All hub genes were subsequently classified into eight curated biological pathway roles (Core PIKK Kinase, Checkpoint Regulator, Cell Cycle Regulator, DNA Repair Effector, Signalling Effector, Chromatin/Transcription, mRNA Surveillance, or Other) based on established DDR molecular literature [11, 12].

---

## 6. Comprehensive Software, Dependency, and Environment Manifest

All network analyses, centrality calculations, and clinical association tests were executed in the R statistical computing environment (R version $\ge 4.2.0$). Below is the inventory of software dependencies:

```
Table 3: Software Dependencies and Computational Infrastructure
┌─────────────────┬──────────────┬────────────────────────────────────────────────────────┐
│ Software/Package│ Origin       │ Primary Methodological Role                            │
├─────────────────┼──────────────┼────────────────────────────────────────────────────────┤
│ igraph          │ CRAN         │ Graph data structures, centrality math, Louvain [7]    │
│ tidygraph       │ CRAN         │ Tidy graph manipulation and relational data management │
│ ggraph          │ CRAN         │ Publication-grade force-directed vector network plots  │
│ ggplot2         │ CRAN         │ Statistical scatterplots, boxplots, and bar charts     │
│ ggrepel         │ CRAN         │ Dynamic collision-free text placement on network plots │
│ pheatmap        │ CRAN         │ Hierarchical clustering heatmaps of centrality metrics │
│ RColorBrewer    │ CRAN         │ Qualitative and diverging color palettes               │
│ scales          │ CRAN         │ Graphical scaling, transformations, and axis formatting│
│ cowplot         │ CRAN         │ Publication grid alignment and multi-panel compilation │
│ reshape2        │ CRAN         │ Melt and cast operations for wide-to-long conversions  │
│ dplyr / tidyr   │ CRAN         │ Relational table joins, groupings, and data filtering  │
│ readr / tibble  │ CRAN         │ Format-preserving I/O and column-name preservation     │
│ cairo_pdf       │ Base system  │ Cairo-based anti-aliased vector graphic PDF rendering  │
└─────────────────┴──────────────┴────────────────────────────────────────────────────────┘
```

---

## 7. Peer-Reviewed Scientific Bibliography

1. **Jeong H, Mason SP, Barabási AL, Oltvai ZN.** (2001). Lethality and centrality in protein networks. *Nature*, 411(6833):41–42. DOI: [10.1038/35075138](https://doi.org/10.1038/35075138).
2. **Szklarczyk D, Kirsch R, Koutrouli M, Nastou K, Mehryary F, Hachilif R, Gable AL, Fang T, Doncheva NT, Pyysalo S, Bork P, Jensen LJ, von Mering C.** (2023). The STRING database in 2023: protein-protein association networks with increased coverage, integration of open access data and new computational tools. *Nucleic Acids Research*, 51(D1):D638–D646. DOI: [10.1093/nar/gkac1000](https://doi.org/10.1093/nar/gkac1000).
3. **Watts DJ, Strogatz SH.** (1998). Collective dynamics of 'small-world' networks. *Nature*, 393(6684):440–442. DOI: [10.1038/30918](https://doi.org/10.1038/30918).
4. **Blondel VD, Guillaume JL, Lambiotte R, Lefebvre E.** (2008). Fast unfolding of communities in large networks. *Journal of Statistical Mechanics: Theory and Experiment*, 2008(10):P10008. DOI: [10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008).
5. **Kruskal WH, Wallis WA.** (1952). Use of ranks in one-criterion variance analysis. *Journal of the American Statistical Association*, 47(260):583–621. DOI: [10.1080/01621459.1952.10483441](https://doi.org/10.1080/01621459.1952.10483441).
6. **Mann HB, Whitney DR.** (1947). On a test of whether one of two random variables is stochastically larger than the other. *Annals of Mathematical Statistics*, 18(1):50–60. DOI: [10.1214/aoms/1177730491](https://doi.org/10.1214/aoms/1177730491).
7. **Csardi G, Nepusz T.** (2006). The igraph software package for complex network research. *InterJournal, Complex Systems*, 1695(5):1–9. URI: [https://igraph.org](https://igraph.org).
8. **Pedersen TL.** (2022). *tidygraph: A Tidy API for Graph Manipulation*. R package version 1.2.2. CRAN: [https://CRAN.R-project.org/package=tidygraph](https://CRAN.R-project.org/package=tidygraph).
9. **Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, Amin N, Schwikowski B, Ideker T.** (2003). Cytoscape: a software environment for integrated models of biomolecular interaction networks. *Genome Research*, 13(11):2498–2504. DOI: [10.1101/gr.1239303](https://doi.org/10.1101/gr.1239303).
10. **Benjamini Y, Hochberg Y.** (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B (Methodological)*, 57(1):289–300. DOI: [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).
11. **Blackford AN, Jackson SP.** (2017). ATM, ATR, and DNA-PK: The Trinity at the Heart of the DNA Damage Response. *Molecular Cell*, 66(6):801–817. DOI: [10.1016/j.molcel.2017.05.015](https://doi.org/10.1016/j.molcel.2017.05.015).
12. **Saxton RA, Sabatini DM.** (2017). mTOR Signaling in Growth, Metabolism, and Disease. *Cell*, 168(6):960–976. DOI: [10.1016/j.cell.2017.02.004](https://doi.org/10.1016/j.cell.2017.02.004).
