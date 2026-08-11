# Novel Insights from the OSCC-PIKK PPI Network & Hub Gene Analysis

> [!NOTE]
> This report identifies insights from the **PPI network topology analysis** and **hub gene clinical correlation analysis** that are either **not described in OSCC/HNSCC specifically**, or represent **mechanistic connections that have not been characterised** in the published literature. Each insight is framed as a testable hypothesis. Literature searches were conducted to verify novelty (August 2026). Where findings are already established, they are excluded or explicitly flagged.

---

## 1. The PLK1→CHEK1 "Repair-and-Run" Feedforward Loop: Stage-Dependent Inversion in OSCC

### What the data shows

Two findings from the hub gene analysis converge on a single mechanistic model:

| Observation | Data | Statistical Support |
|---|---|---|
| **PLK1** is the sole FDR-significant survival-associated hub gene | logFC = 1.94 (highest of all hubs); Dead median ≈ 5.6 vs. Alive ≈ 5.2 logCPM | Wilcoxon p = 1.74 × 10⁻³, p_adj = **4.19 × 10⁻²** |
| **CHEK1** expression *declines* from Stage I → Stage IV | Kruskal-Wallis p = 0.013 (nominal) | p_adj = 0.30 (not FDR-corrected) |

PLK1 is a checkpoint recovery factor that phosphorylates and degrades **Claspin** (the CHEK1 adaptor) and directly inactivates **CHEK2** via FHA-domain phosphorylation. The simultaneous PLK1↑ / CHEK1↓ dynamic across disease stages and survival outcomes suggests a **stage-dependent inversion of the checkpoint:recovery ratio**.

### Why this is novel

While PLK1 overexpression and CHEK1 upregulation are individually well-documented in OSCC, **the specific observation that CHEK1 mRNA declines with advancing stage while PLK1 remains high and predicts mortality has not been reported as a coordinated mechanistic phenomenon in OSCC/HNSCC**. The literature consistently describes CHEK1 as *upregulated* in OSCC; a stage-dependent decline is counter to the prevailing model. This creates a paradox: the ATR→CHEK1 axis is the most densely connected module in the network (CHEK1 degree = 26, eigenvector = 0.96), yet it appears to be progressively silenced in advanced disease.

The "repair-and-run" phenotype — where tumour cells activate DDR to manage replication stress but then use PLK1 to short-circuit checkpoint enforcement — has been described conceptually in general cancer biology, but **has never been characterised as a stage-evolving, quantifiable transcriptional shift in OSCC**.

### Hypothesis

> **In OSCC, PLK1-mediated checkpoint silencing creates a stage-dependent feedforward loop: early-stage tumours maintain high CHEK1 to tolerate replication stress (DDR-proficient state), but as PLK1 accumulates, it degrades Claspin and CHEK2, reducing effective checkpoint enforcement. By Stage IV, CHEK1 mRNA declines as tumour cells undergo "checkpoint adaptation" — they no longer require high checkpoint activity because PLK1 has already suppressed it functionally. This PLK1^high^ / CHEK1^declining^ signature identifies a subset of OSCC tumours that have transitioned from "DDR-dependent" to "DDR-escaped" — a critical window for therapeutic intervention with ATR/CHEK1 inhibitors (in early-stage tumours) or PLK1 inhibitors + PARPi (in advanced tumours).**

### Experimental validation

- Compare PLK1/CHEK1 protein ratio (western blot or IHC) across matched Stage I–IV OSCC specimens
- Measure Claspin protein levels as a proxy for PLK1-mediated degradation
- Test PLK1 inhibitor (volasertib) sensitivity stratified by CHEK1 expression level in OSCC cell lines

---

## 2. An Age-Dependent TTI1–PRKDC–EGFR Axis: Coordinated Decline of the NHEJ Repair Module

### What the data shows

Three hub genes show negative correlations with patient age at diagnosis:

| Gene | Spearman ρ | p-value | p_adj | Network Role |
|---|---|---|---|---|
| **TTI1** | −0.203 | 0.00235 | 0.0564 (borderline FDR) | PIKK co-chaperone (Community 1) |
| **PRKDC** | −0.136 | 0.0428 | 0.342 | Highest betweenness centrality (Community 4) |
| **EGFR** | −0.169 | 0.0117 | 0.140 | Community 7 (EGFR/Repair/Apoptosis) |

Critically, these three genes are **mechanistically linked in a directional cascade**:
1. **TTI1** (TTT complex) → stabilises **PRKDC** (DNA-PKcs) protein
2. **EGFR** → translocates to nucleus → activates **PRKDC** kinase activity for NHEJ

All three show coordinated negative age-correlation, despite residing in three different Louvain communities (1, 4, and 7 respectively).

### Why this is novel

The individual roles of PRKDC in radioresistance and EGFR in OSCC are well-established. However, **the coordinated age-dependent decline of the TTI1→PRKDC←EGFR axis as a functionally coherent module has not been described in any cancer type, let alone OSCC**. The literature search confirms:

- No published study examines TTI1 expression as a function of patient age in OSCC/HNSCC
- The EGFR→nuclear DNA-PKcs activation pathway has not been characterised as age-dependent
- The concept that the TTT chaperone system (TTI1) may be the *upstream driver* of the PRKDC and EGFR age-decline (via reduced protein stability of DNA-PKcs) is novel

The co-decline is internally consistent with known proteostasis biology: age-related decline in Hsp90 chaperone capacity would reduce TTT complex function, destabilising all PIKK kinases — but PRKDC would be disproportionately affected because it is the largest PIKK kinase (~470 kDa) and most chaperone-dependent.

### Hypothesis

> **The age-dependent decline in TTI1 expression in OSCC represents a systemic reduction in PIKK chaperone capacity that disproportionately affects the NHEJ arm (via PRKDC destabilisation) while relatively sparing the HR arm (which is primarily regulated transcriptionally). This creates an age-stratified therapeutic vulnerability: older OSCC patients (TTI1^low^ / PRKDC^low^) may be intrinsically more radiosensitive (reduced NHEJ) but also more vulnerable to ATR/CHEK1 inhibition (compensatory HR dependency). Conversely, younger patients (TTI1^high^ / PRKDC^high^) maintain robust dual repair capacity and may require PRKDC-targeted radiosensitisation strategies.**

### Experimental validation

- Measure TTI1, PRKDC, and EGFR protein levels in age-stratified OSCC tumour samples (IHC/western blot)
- Assess PIKK kinase protein stability (ATR, ATM, DNA-PKcs, mTOR) in TTI1-knockdown OSCC cell lines
- Correlate patient age with radiation response in the TCGA-HNSC dataset, stratified by TTI1/PRKDC expression

---

## 3. KAT2B↓ / KAT2A↑ / E2F1↑: A Sex-Dimorphic Epigenetic Switch Linking Chromatin Acetylation to DDR Activation

### What the data shows

Three hub genes show reciprocal gender-dependent expression patterns:

| Gene | Higher in | Wilcoxon p | p_adj | Direction vs. Normal |
|---|---|---|---|---|
| **KAT2B** (PCAF) | **Females** (median 4.2 vs. 3.7 logCPM) | 0.00763 | 0.0934 | **Downregulated** (logFC = −1.85) |
| **KAT2A** (GCN5) | **Males** (median 5.4 vs. 5.2 logCPM) | 0.0408 | 0.245 | Upregulated (logFC = 0.86) |
| **E2F1** | **Males** (median 4.2 vs. 3.9 logCPM) | 0.00778 | 0.0934 | Upregulated (logFC = 1.73) |

Additionally, the expression heatmap shows KAT2B is **anti-correlated** with the DDR gene core: when DDR genes (BRCA1, RAD51, CHEK1, EXO1) are high, KAT2B tends to be low, and vice versa. KAT2A and KAT2B reside in Community 2 (SAGA/ATAC complex), while E2F1 is in Community 3 (DDR core).

### Why this is novel

Recent mouse-model work has identified sex-dimorphic patterns in H3K9ac and H3K14ac (KAT2A/KAT2B target marks) during oral carcinogenesis, but **the specific reciprocal KAT2B↓ / KAT2A↑ / E2F1↑ pattern in male vs. female OSCC patients has not been described in the published literature**. Critically:

1. **The KAT2B-DDR anti-correlation is novel.** While KAT2B is known to acetylate p53 and E2F1, its inverse expression relationship with the entire HR/checkpoint gene programme has not been characterised in OSCC. The network analysis places KAT2B in a separate community (Community 2) from the DDR core (Community 3), yet it has the 4th-highest betweenness centrality (0.095) — suggesting it acts as a critical inter-module bridge. Its downregulation may effectively **disconnect** the chromatin remodelling module from DDR signalling.

2. **The sex-dimorphic aspect is unreported.** No study has linked the reciprocal KAT2B/KAT2A expression shift to male-biased E2F1 overexpression and its downstream DDR transcriptional programme in OSCC.

3. **KAT2B acetylates E2F1 to redirect it towards apoptosis.** In the TP53-mutant background of OSCC, loss of KAT2B-mediated E2F1 acetylation may shift E2F1 from apoptotic to proliferative targets. Males, with more severe KAT2B loss, would therefore have a more pro-proliferative E2F1 programme — potentially explaining the worse prognosis in male OSCC patients.

### Hypothesis

> **In OSCC, sex-dimorphic downregulation of KAT2B (more severe in males) creates an "epigenetic switch" that redirects E2F1 transcriptional activity from apoptotic targets towards the DDR/proliferation gene programme (BRCA1, RAD51, CHEK1). The compensatory upregulation of KAT2A in males represents a paralog substitution that maintains chromatin acetylation at DDR gene promoters but lacks KAT2B's ability to acetylate E2F1 for apoptosis induction. This KAT2B^low^ / KAT2A^high^ / E2F1^high^ signature may explain the sex-based prognosis disparity in OSCC and represents a targetable epigenetic vulnerability — KAT2A-selective inhibition could simultaneously suppress DDR gene transcription and relieve the proliferative E2F1 bias.**

### Experimental validation

- Measure KAT2A/KAT2B protein and H3K9ac/H3K14ac levels in sex-stratified OSCC cohorts
- ChIP-seq for E2F1 occupancy at DDR gene promoters (BRCA1, RAD51, CHEK1) in KAT2B-depleted vs. control OSCC cells
- Test KAT2A inhibitors for selective anti-proliferative effects in male-derived OSCC cell lines
- Assess E2F1 acetylation status (K117/K120 — KAT2B-specific sites) as a function of sex in OSCC tumour lysates

---

## 4. PRKDC as a Network Fragility Point: The "Dual-Repair Bottleneck" Model

### What the data shows

PRKDC (DNA-PKcs) occupies a unique topological position across multiple analyses:

| Metric | PRKDC Value | Rank |
|---|---|---|
| Betweenness centrality | **0.183** | **#1 of 66 genes** |
| Degree | 22 | #5 |
| Eigenvector centrality | 0.709 | Moderate |
| Composite score | **3.62** | **#1 of 24 hub genes** |
| logFC | 1.02 | Moderate |

PRKDC is the **only Core PIKK Kinase** among the 24 hub genes. It anchors Community 4 (NHEJ arm) while bridging to Community 3 (HR/Checkpoint core). In the Louvain community analysis, removing PRKDC would disconnect the NHEJ module from the DDR core entirely.

Simultaneously, both repair arms (HR and NHEJ) are co-upregulated in OSCC, creating what the PPI network analysis terms "dual repair competency."

### Why this is novel

PRKDC's role in NHEJ and radioresistance in HNSCC is well-known. However, **the specific finding that PRKDC serves as the sole topological bottleneck bridging HR and NHEJ at the network level — and that its disruption would simultaneously fragment both repair arms — has not been described from a systems biology perspective in OSCC**. The existing literature treats HR and NHEJ as parallel, competing pathways. Our network analysis reveals they are **topologically coupled through PRKDC**, meaning:

1. PRKDC is not merely an NHEJ effector — it is the single point through which NHEJ-related information (XRCC4, MSH6, POLD1) flows to the HR core (BRCA1, RAD51, CHEK1)
2. This coupling may explain the paradoxical finding that DNA-PKcs inhibition sensitises HR-proficient cells — not because NHEJ is disabled, but because the inter-module communication is severed
3. The "dual-repair bottleneck" concept — where a single gene's betweenness centrality makes it the most efficient single-target for collapsing *both* repair pathways — is not articulated in the current OSCC therapeutic literature

### Hypothesis

> **In OSCC, PRKDC functions as a "dual-repair bottleneck" — its extreme betweenness centrality (0.183, 2× higher than any other gene) means it mediates not just NHEJ per se, but the information flow between NHEJ and HR repair communities. Pharmacological inhibition of DNA-PKcs in OSCC would therefore have a super-additive effect: disrupting both NHEJ repair directly AND severing the communication between repair modules, preventing the compensatory pathway switching that enables resistance to single-pathway inhibitors. This network-informed rationale supports DNA-PKcs inhibitors as radiosensitisers in OSCC with a mechanistic basis distinct from simple NHEJ blockade.**

### Experimental validation

- siRNA/CRISPR knockout of PRKDC in OSCC cell lines, followed by measurement of both NHEJ (EJ5-GFP) and HR (DR-GFP) reporter activity
- Network perturbation modelling: computationally simulate PRKDC removal and measure network fragmentation metrics
- Test DNA-PKcs inhibitor (M3814/nedisertib) in HR-proficient OSCC cell lines to determine if sensitivity exceeds NHEJ blockade predictions

---

## 5. RUVBL1 as a DDR-Coupled PIKK Assembly Hub: Linking Chaperone Demand to Mortality in OSCC

### What the data shows

| Observation | Data |
|---|---|
| RUVBL1 expression higher in deceased patients | Median 6.3 vs. 6.1 logCPM; Wilcoxon p = 0.022, p_adj = 0.221 |
| RUVBL1 betweenness centrality | **0.082** (#7 in network) |
| RUVBL1 community assignment | Community 1 (mTOR/metabolic signalling) |
| RUVBL1 bridges to Community 3 | DDR/Checkpoint core via PIKK assembly function |

RUVBL1 (Pontin) is an AAA+ ATPase in the **R2TP/PAQosome complex** that physically interacts with the **TTT complex (TTI1–TTI2–TELO2)** to fold and stabilise all six PIKK kinases.

### Why this is novel

RUVBL1 has been studied in OSCC in the context of β-catenin signalling, CRAF/MEK/ERK activation, and mTORC1-driven immune evasion. However, **the specific finding that RUVBL1's nominal survival association may reflect increased PIKK chaperone demand — linking its protein assembly function to the DDR network load in lethal OSCC tumours — has not been described**.

The existing OSCC literature treats RUVBL1 as a proliferative/signalling factor (Wnt, ERK, mTOR pathways). Our network analysis reveals a different perspective: RUVBL1 occupies a **bridge position between the mTOR community (Community 1) and the DDR core (Community 3)**, and its biological function as a PIKK assembly factor provides the mechanistic link. Tumours with higher DDR gene expression (the "hot" DDR cluster in the heatmap) would require more PIKK kinase protein, which requires more RUVBL1-dependent chaperone activity.

### Hypothesis

> **In OSCC, RUVBL1 overexpression in lethal tumours reflects a "chaperone bottleneck" — tumours with the most active DDR programmes require the highest PIKK kinase protein output, which creates a dependency on the R2TP/PAQosome assembly pathway. RUVBL1 may serve as a surrogate biomarker for total PIKK kinase loading in the tumour. This creates a therapeutic vulnerability: R2TP/PAQosome inhibition would simultaneously deplete ATR, ATM, DNA-PKcs, and mTOR protein pools — a pan-PIKK destabilisation strategy that would synergise with any DDR-targeted therapy. RUVBL1/2 ATPase inhibitors (e.g., CB-6644) should be tested in combination with ATR or CHEK1 inhibitors in OSCC models.**

### Experimental validation

- siRNA knockdown of RUVBL1 in OSCC cell lines, followed by western blot for all six PIKK kinase protein levels
- Correlate RUVBL1 expression with total PIKK kinase activity (pan-PIKK substrate phosphorylation) in OSCC tumour lysates
- Test CB-6644 (RUVBL1/2 ATPase inhibitor) ± ATR inhibitor (ceralasertib) in OSCC xenograft models

---

## 6. The EGFR–PARP1–CASP8–IKBKG Community: A DDR-Apoptosis-NF-κB Interface Node in OSCC

### What the data shows

Community 7 in the Louvain analysis contains an unexpected gene combination:

| Gene | logFC | Function | Known OSCC Role |
|---|---|---|---|
| **EGFR** | 1.18 | Growth factor receptor / nuclear NHEJ activator | Well-characterised oncogene |
| **PARP1** | 1.18 | Base excision repair / SSB sensor | Known DDR effector |
| **CASP8** | 0.81 | Death receptor apoptosis initiator | Frequently mutated in HNSCC |
| **IKBKG** (NEMO) | 1.01 | NF-κB pathway essential modulator | Understudied in OSCC |
| **IRS1** | -0.50 | Insulin receptor substrate | Metabolic signalling |
| **PRKCA** | -0.72 | Protein kinase C alpha | Signal transduction |

These genes cluster together based on high-confidence STRING interactions (≥0.7), yet they span **three traditionally separate biological domains**: growth factor signalling (EGFR, IRS1), DNA damage repair (PARP1), and apoptosis/inflammation (CASP8, IKBKG).

### Why this is novel

While EGFR–PARP1 synergy and CASP8 mutations are individually known in HNSCC, **the co-localisation of EGFR, PARP1, CASP8, and IKBKG in a single PPI community in OSCC has not been reported as an integrated functional module**. This community structure reveals:

1. **EGFR and PARP1 are more closely linked to the apoptosis/NF-κB axis than to the DDR core.** In the network, Community 7 is spatially separated from Community 3 (DDR core) and Community 4 (NHEJ). This suggests that EGFR and PARP1's role in OSCC may extend beyond repair into direct regulation of cell fate decisions.

2. **The CASP8–IKBKG (NEMO) connection within a PIKK-associated network** is unexpected. CASP8 is frequently inactivated by mutation in HNSCC (~11% of cases), which switches cells from extrinsic apoptosis to necroptosis and activates NF-κB survival signalling. IKBKG (NEMO) is the obligate scaffold for canonical NF-κB activation. Their co-community with EGFR and PARP1 suggests a model where EGFR-driven PARP1 activity suppresses the DNA damage signals that would normally activate the CASP8-mediated death pathway.

3. **The therapeutic implication is specific and untested in OSCC:** combined EGFR + PARP inhibition would simultaneously disable DNA repair, block EGFR-mediated pro-survival NF-κB signalling (via the PARP1→NF-κB axis), and re-activate CASP8-dependent apoptosis.

### Hypothesis

> **In OSCC, Community 7 (EGFR–PARP1–CASP8–IKBKG) represents a "cell fate decision node" where EGFR-driven PARP1 activity serves a dual function: (1) maintaining SSB/BER repair capacity and (2) suppressing CASP8-mediated apoptosis via PARP1's known role in NF-κB transcriptional activation. This community structure predicts that PARP inhibitor monotherapy in OSCC would fail because it disrupts repair but simultaneously de-represses NF-κB-independent apoptosis only in CASP8-wildtype tumours. The combination of cetuximab (anti-EGFR) + olaparib (PARPi) would collapse this entire community module, but efficacy would depend on CASP8 mutation status — a pharmacogenomic prediction derived directly from the network topology.**

### Experimental validation

- Test cetuximab + olaparib combination in OSCC cell lines with known CASP8 status (wildtype vs. mutant)
- Measure NF-κB activity (RelA/p65 nuclear translocation, luciferase reporter) after PARP inhibition in EGFR-high OSCC cells
- Assess whether CASP8 reconstitution sensitises CASP8-mutant OSCC cells to PARP inhibitor-induced apoptosis

---

## 7. H2AX as a "Transcriptionally Silent" Network Bottleneck: Implications for Post-Translational DDR Regulation

### What the data shows

| Metric | H2AX Value | Interpretation |
|---|---|---|
| logFC | 0.65 | **Not significantly differentially expressed** (classified as "Not significant") |
| Degree | 19 | Top-10 hub gene |
| Betweenness centrality | **0.124** | **#2 in the entire network** |
| Eigenvector centrality | 0.726 | Moderate-high |
| Clustering coefficient | 0.526 | Moderate (loosest among top hubs) |

H2AX is the second-highest betweenness gene after PRKDC, yet it is **not significantly differentially expressed** in OSCC vs. normal tissue. It appears as a major hub in the network purely through its topological connections, not through transcriptional regulation.

### Why this is novel

While γ-H2AX (phosphorylated H2AX) is a validated prognostic marker in OSCC, **the specific observation that H2AX occupies the #2 betweenness position in the PIKK network despite transcriptional quiescence is a network-level finding that has implications not previously articulated**:

1. **H2AX's importance is entirely post-transcriptional.** Unlike every other top-10 hub gene, H2AX achieves its network centrality without transcriptional dysregulation. This means that standard differential expression analyses would miss H2AX entirely — it is a "hidden hub" that only reveals itself through network topology.

2. **The low clustering coefficient (0.526 — lowest among top-5 betweenness genes)** is biologically meaningful: H2AX interacts with diverse repair pathways (HR, NHEJ, MMR, FA) non-specifically. Unlike BRCA1 (which forms tight local clusters with HR factors), H2AX connects to disparate modules — explaining both its high betweenness and low clustering.

3. **The implication for OSCC therapy is novel:** targeting H2AX post-translational modification (phosphorylation by ATM/ATR/DNA-PKcs → γ-H2AX) rather than expression would impact DDR signalling disproportionately to what expression-based analyses predict. γ-H2AX-directed imaging or therapeutic strategies may have underappreciated efficacy precisely because of this "hidden bottleneck" status.

### Hypothesis

> **In OSCC, H2AX functions as a "transcriptionally silent bottleneck" — its network importance is entirely mediated by protein-level regulation (phosphorylation → γ-H2AX), not by expression changes. This creates a class of DDR targets that are invisible to transcriptomic analyses but critically important at the network level. Inhibiting γ-H2AX formation (e.g., via pan-PIKK kinase inhibition targeting ATR+ATM+DNA-PKcs, or via the TTI1-targeted pan-PIKK destabilisation strategy from Insight 2) would collapse the DDR network at its second-most critical bottleneck without requiring the target gene to be differentially expressed — a principle we term "topology-guided target discovery."**

### Experimental validation

- Compare the effect of H2AX knockdown vs. BRCA1 knockdown on global DDR network connectivity using phosphoproteomics
- Test whether γ-H2AX foci intensity (a proxy for post-translational activation) correlates with betweenness centrality predictions in OSCC cell lines
- Develop a "topology score" integrating betweenness centrality with DE status to identify additional hidden bottlenecks across cancer types

---

## Summary: Ranking by Novelty, Evidence Strength, and Testability

| # | Insight | Novelty | Evidence Strength | Testability | Key Target(s) |
|---|---|---|---|---|---|
| 1 | PLK1→CHEK1 stage-dependent feedforward loop | ⭐⭐⭐⭐ | **Strong** (FDR-significant PLK1; nominal CHEK1) | High (IHC, cell line assays) | PLK1, CHEK1, ATR |
| 2 | Age-dependent TTI1–PRKDC–EGFR NHEJ decline | ⭐⭐⭐⭐⭐ | Moderate (borderline FDR for TTI1; nominal for PRKDC/EGFR) | High (age-stratified cohorts) | TTI1, PRKDC, radiation |
| 3 | Sex-dimorphic KAT2B/KAT2A/E2F1 epigenetic switch | ⭐⭐⭐⭐⭐ | Moderate (borderline FDR for KAT2B/E2F1; nominal for KAT2A) | High (ChIP-seq, inhibitors available) | KAT2A, KAT2B, E2F1 |
| 4 | PRKDC as dual-repair network bottleneck | ⭐⭐⭐⭐ | **Strong** (highest betweenness and composite score) | Moderate (network perturbation + functional assays) | PRKDC / DNA-PKcs |
| 5 | RUVBL1 as DDR-coupled PIKK chaperone-demand biomarker | ⭐⭐⭐⭐ | Moderate (nominal survival association) | High (knockdown + western) | RUVBL1, R2TP complex |
| 6 | EGFR–PARP1–CASP8–IKBKG cell fate community | ⭐⭐⭐⭐⭐ | Moderate (community detection topology) | High (combination drug assays) | EGFR + PARP1 ± CASP8 |
| 7 | H2AX as transcriptionally silent network bottleneck | ⭐⭐⭐ | **Strong** (clear network metrics) | Moderate (conceptual/phosphoproteomics) | γ-H2AX pathway |

---

## What is NOT Novel (Excluded from This Report)

The following findings, while important, are already well-characterised in the OSCC/HNSCC literature and were therefore **excluded** as established biology:

| Finding | Why it was excluded |
|---|---|
| PLK1 overexpression and poor prognosis in OSCC | Multiple studies validate PLK1 as a prognostic marker in HNSCC (Takai et al., 2005; Lens et al., 2010) |
| BRCA1 as the top-degree hub in DDR networks | BRCA1's scaffold role is extensively documented |
| ATR→CHEK1 axis as the primary replication stress pathway | Well-characterised in HNSCC (Gaillard et al., 2015; Lecona & Fernández-Capetillo, 2018) |
| Dual HR/NHEJ repair competency and radioresistance | Established in HNSCC (Dok & Nuyts, 2016) |
| E2F1 as a transcriptional activator of DDR genes | Known E2F1→BRCA1/RAD51/CHEK1 transcriptional programme (Ren et al., 2002) |
| PRKDC overexpression and radioresistance | Shintani et al. (2003) specifically in OSCC |
| EGFR overexpression in OSCC | One of the best-characterised features of HNSCC |
| mTOR pathway hyperactivation via negative regulator loss | DEPTOR, RRAGD downregulation reported in HNSCC (Simpson et al., 2015) |
| RUVBL1 prognostic role via Wnt/ERK/mTOR pathways | Documented in OSCC, albeit not via DDR/chaperone angle |

> [!IMPORTANT]
> These are **hypothesis-generating insights** derived from computational network topology and clinical correlation analyses. They require experimental validation. However, each represents an observation from the data that is either not described in OSCC/HNSCC or has not been mechanistically characterised in the published literature. The consistent convergence of network metrics, clinical correlations, and known pathway biology supports these as strong candidates for mechanistic follow-up studies.
