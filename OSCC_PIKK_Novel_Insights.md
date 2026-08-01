# Novel Insights from the OSCC-PIKK ORA Results

> [!NOTE]
> Each insight below represents an observation from your data that is either **not described in OSCC/HNSCC** specifically, or represents a **mechanistic connection that has not been characterised** in the published literature. These are framed as testable hypotheses for further investigation or discussion.

---

## 1. The TRRAP/SAGA–Centrosome–DDR Axis: An Uncharacterised Triangle in OSCC

### What the data shows

Several centrosome-related terms are significantly enriched, but the gene membership is unexpected:

| Term | Genes | FE | p.adj |
|---|---|---|---|
| Centrosome duplication | **KAT2B, CDK2, BRCA2, BRCA1, KAT2A** | 16.14 | 3.08 × 10⁻⁴ |
| Negative regulation of centrosome duplication | **KAT2B, BRCA1, KAT2A** | 58.12 | 3.47 × 10⁻⁴ |
| Regulation of centrosome cycle | **KAT2B, AURKA, BRCA1, KAT2A** | 17.55 | 1.31 × 10⁻³ |
| Centriole replication | **KAT2B, CDK2, BRCA1, KAT2A** | 23.25 | 5.17 × 10⁻⁴ |
| Centriole assembly | **KAT2B, CDK2, BRCA1, KAT2A** | 20.66 | 7.59 × 10⁻⁴ |

### Why this is novel

The recurrence of **KAT2A (GCN5)** and **KAT2B (PCAF)** — the catalytic acetyltransferases of the SAGA and ATAC complexes — in centrosome biology terms alongside **BRCA1** is striking. While BRCA1's role in centrosome amplification control is established (Starita et al., 2004, *Journal of Cell Biology*), the idea that **TRRAP-associated acetyltransferases KAT2A/KAT2B coordinately regulate centrosome duplication in the context of DDR in OSCC is unreported**. The literature search confirms that the SAGA complex's role in centrosome biology is described as "an area of ongoing investigation" with no published OSCC-specific studies.

### Hypothesis

> **In OSCC, dysregulation of TRRAP-complex acetyltransferases (KAT2A/KAT2B) may decouple centrosome number control from the DNA damage checkpoint, permitting centrosome amplification and chromosomal instability via a BRCA1-dependent licensing mechanism that is distinct from their canonical transcriptional coactivator role.**

This would represent a non-canonical, non-transcriptional function of the SAGA/ATAC complex in OSCC — connecting chromatin acetylation, centrosome integrity, and the DDR in a single axis.

---

## 2. The AURKB–PRKDC–PARP1 Triad: Coordinated Suppression of cGAS-STING Innate Immune Sensing

### What the data shows

Two cGAS-STING related terms are enriched with the same three genes:

| Term | Genes | FE | p.adj |
|---|---|---|---|
| **Negative regulation of cGAS/STING signaling pathway** | AURKB, PRKDC, PARP1 | 53.65 | 4.34 × 10⁻⁴ |
| **cGAS/STING signaling pathway** | AURKB, PRKDC, PARP1 | 34.87 | 1.36 × 10⁻³ |

### Why this is novel

While PARP1 and DNA-PKcs are individually known to suppress cGAS-STING activation by efficiently repairing DNA before it leaks into the cytoplasm, the co-enrichment of **AURKB** in this module is unexpected. AURKB is primarily a mitotic kinase; its role in cGAS-STING regulation is "less characterised" per the literature. 

The critical insight is that these three genes represent **three different PIKK-associated functional modules** (DNA-PKcs = NHEJ, PARP1 = BER/SSB repair, AURKB = mitotic checkpoint) that are **converging on innate immune suppression**. This suggests:

### Hypothesis

> **OSCC tumours deploy a multi-layered, PIKK-dependent immunoevasion strategy where DNA-PKcs suppresses cGAS-STING at DSBs, PARP1 prevents cytosolic DNA accumulation from SSBs, and AURKB suppresses micronuclei-derived cGAS activation during aberrant mitosis. The coordinated dysregulation of all three in OSCC may represent a previously unrecognised "immune cold" phenotype driven by DDR hyperactivity rather than classical immune checkpoint upregulation.**

This has direct therapeutic implications: combined inhibition of PARP + DNA-PKcs (or PARP + Aurora B) may simultaneously unleash immune sensing and impair repair — a dual-hit strategy that has not been explored in OSCC.

---

## 3. TTI1 as a Pan-PIKK Vulnerability Hub — Unexplored in OSCC

### What the data shows

**TTI1** appears in 5 enriched GO terms:
- Cell cycle checkpoint signaling
- Regulation of cell cycle checkpoint
- Regulation of DNA damage checkpoint
- TOR signaling
- Regulation of TOR signaling

It is present in the PIKK gene universe as an **mTOR-associated** gene, yet the ORA maps it to checkpoint and DDR processes as well.

### Why this is novel

TTI1 is a component of the TTT complex (TELO2-TTI1-TTI2) that acts as a co-chaperone required for the stability and assembly of **all six PIKK kinases** (ATM, ATR, DNA-PKcs, mTOR, SMG1, TRRAP). The literature confirms TTI1 dysregulation in glioblastoma and NSCLC, but **no study has examined TTI1 expression or function in OSCC/HNSCC**. Its role in OSCC is explicitly described as "uncharacterised" in the search results.

### Hypothesis

> **TTI1 dysregulation in OSCC acts as a single-point "gain-of-function" amplifier that simultaneously boosts the protein stability of ATM, ATR, DNA-PKcs, and mTOR, providing the tumour with enhanced DDR capacity AND metabolic fitness. TTI1 may represent a uniquely efficient therapeutic target in OSCC because its inhibition would simultaneously destabilise all PIKK-dependent survival mechanisms — a pan-PIKK synthetic lethality.**

This is a highly targetable concept: existing research shows that ivermectin targets TELO2 (a TTT complex partner), suggesting pharmacological accessibility of this axis.

---

## 4. NuA4/TIP60 Chromatin Remodeling Subunits as Positive Regulators of HR in OSCC

### What the data shows

The "positive regulation of DSB repair via HR" term is enriched with a specific gene set:

| Term | Genes | FE | p.adj |
|---|---|---|---|
| Positive regulation of DSB repair via HR | **MRGBP, ACTL6A, RUVBL1, PARP1** | 22.14 | 6.09 × 10⁻⁴ |
| Positive regulation of DSB repair | MRGBP, ACTL6A, PRKDC, RUVBL1, PARP1 | 12.63 | 8.41 × 10⁻⁴ |
| Positive regulation of DNA repair | MRGBP, ACTL6A, PRKDC, EGFR, BRCA1, RUVBL1, PARP1, H2AX | 13.78 | 3.78 × 10⁻⁶ |

These are predominantly **chromatin remodelling subunits**: MRGBP, ACTL6A, and RUVBL1 are all components of the **NuA4/TIP60 histone acetyltransferase complex** (which also includes TRRAP).

### Why this is novel

While the NuA4/TIP60 complex is known to participate in DDR via H4K16 acetylation at DSB sites, **its specific role as a positive regulator of HR in the context of OSCC is not established**. The literature describes these proteins as "frequently upregulated in various cancers" but does not characterise their HR-promoting function in oral/head and neck malignancies specifically.

The novel angle is that the **TRRAP family genes from your curated PIKK universe are not merely transcriptional regulators — they are actively promoting HR-dependent repair in OSCC**, and this creates a specific vulnerability:

### Hypothesis

> **OSCC tumours exploit TRRAP-associated NuA4/TIP60 chromatin remodelling (via MRGBP, ACTL6A, RUVBL1) to create a permissive chromatin state at DSBs that favours HR over NHEJ. This "chromatin-priming" for HR may explain why some OSCC tumours are resistant to radiation (which induces DSBs repaired by NHEJ) but could be sensitised by NuA4/TIP60 complex disruption — a strategy that would simultaneously impair HR and transcriptional programmes.**

---

## 5. KDM1A/LSD1 as a Multi-Pathway PIKK-Network Epigenetic Integrator in OSCC

### What the data shows

**KDM1A (LSD1)** appears across a remarkably diverse set of enriched pathways:

| Pathway category | Terms where KDM1A appears |
|---|---|
| DSB repair | Double-strand break repair, DSB repair via HR, recombinational repair |
| Checkpoint signalling | Signal transduction in response to DNA damage, regulation of signal transduction by p53 class mediator |
| Radiation response | Response to ionizing radiation, cellular response to ionizing radiation, cellular response to gamma radiation |
| UV response | Response to UV, cellular response to UV |
| Apoptosis | Intrinsic apoptotic signalling in response to DNA damage |

KDM1A is classified in the PIKK universe under **PRKDC (DNA-PKcs)-associated** genes.

### Why this is novel

While KDM1A/LSD1 is known to be recruited to DSB sites and facilitate 53BP1/BRCA1 recruitment via H3K4me2 demethylation, **its simultaneous involvement in DDR, checkpoint signalling, p53 regulation, radiation response, and apoptotic signalling within a single PIKK-centric gene network in OSCC has not been described as an integrated phenomenon**.

Most studies examine LSD1 in one of these contexts (usually transcription or HR individually). The data here suggests it acts as a **master epigenetic coordinator** across the entire PIKK-dependent DDR programme.

### Hypothesis

> **In OSCC, KDM1A/LSD1 functions not merely as a chromatin modifier at DSB sites, but as a pan-DDR epigenetic integrator that coordinates DSB repair pathway choice (HR vs. NHEJ), checkpoint enforcement (via p53 regulation), and damage-induced apoptosis within the PIKK signalling network. Its inhibition in OSCC may therefore have pleiotropic anti-tumour effects that extend far beyond transcriptional suppression — disrupting repair, checkpoint, and survival simultaneously. This "epigenetic integrator" concept for LSD1 in the PIKK network of OSCC is uncharacterised.**

---

## 6. Meiotic Gene Programme Reactivation in Somatic OSCC — A DDR-Driven "Cancer/Testis"-like Phenotype

### What the data shows

Multiple meiosis-specific terms are enriched with high significance:

| Term | Genes | FE | p.adj |
|---|---|---|---|
| Meiotic cell cycle | AURKA, PLK1, EXO1, RAD51, CCNE1, CDK2, FANCA, BRCA2, MSH6, FANCM, H2AX, RBBP8 | 11.25 | 3.07 × 10⁻⁸ |
| Meiosis I cell cycle process | AURKA, PLK1, RAD51, CCNE1, BRCA2, MSH6, FANCM | 14.79 | 1.36 × 10⁻⁵ |
| Meiotic nuclear division | AURKA, PLK1, RAD51, CCNE1, FANCA, BRCA2, FANCM | 9.86 | 1.51 × 10⁻⁴ |
| Meiosis I | AURKA, PLK1, RAD51, CCNE1, BRCA2, FANCM | 13.41 | 1.28 × 10⁻⁴ |

Additionally, **Reactome: Meiotic recombination** (RAD51, CDK2, BRCA2, BRCA1, H2AX, RBBP8) is enriched at FE = 15.99, p.adj = 2.75 × 10⁻⁴.

### Why this is novel

The enrichment of meiotic terms in a somatic tumour analysis is not merely an artefact of shared genes between meiosis and mitotic DDR. While genes like RAD51 and BRCA2 participate in both contexts, the specific combination of **FANCM + MSH6 + EXO1** — which are critical for meiotic crossover resolution and mismatch-dependent heteroduplex rejection — suggests that OSCC may be engaging molecular machinery that is normally restricted to the germline.

This phenomenon — the aberrant expression of germline/meiosis-specific gene programmes in somatic tumours — is related to the concept of "cancer/testis" (CT) gene reactivation. However, **the specific activation of meiotic recombination machinery (rather than CT antigens) through PIKK-network dysregulation in OSCC is not described in the literature**.

### Hypothesis

> **OSCC tumours exhibit a "meiotic-like" recombination programme driven by the aberrant co-activation of PIKK-network DDR genes (FANCM, MSH6, EXO1, RAD51, BRCA2) that normally coordinate crossover resolution during meiosis I. This programme may provide OSCC cells with enhanced inter-homolog recombination capacity, generating genetic diversity (via loss of heterozygosity or mitotic crossovers) that accelerates tumour evolution and clonal adaptation to therapy. Unlike classical cancer-testis antigen expression, this represents a functional (not merely antigenic) germline programme reactivation that is mechanistically linked to the PIKK DDR network.**

---

## Summary: Ranking by Novelty and Testability

| # | Insight | Novelty Level | Testability | Therapeutic Angle |
|---|---|---|---|---|
| 1 | SAGA–Centrosome–DDR axis | ⭐⭐⭐⭐⭐ | High (immunofluorescence + siRNA) | SAGA inhibitors + anti-mitotics |
| 2 | AURKB–PRKDC–PARP1 cGAS-STING suppression | ⭐⭐⭐⭐⭐ | High (cGAS-STING reporter assays) | PARPi + AurkBi immune activation |
| 3 | TTI1 pan-PIKK hub | ⭐⭐⭐⭐ | Moderate (protein stability assays) | Ivermectin repurposing / TTT inhibitors |
| 4 | NuA4/TIP60 HR promotion | ⭐⭐⭐⭐ | High (DR-GFP HR reporter) | TIP60 inhibitors + radiotherapy |
| 5 | KDM1A pan-DDR integrator | ⭐⭐⭐ | High (LSD1 inhibitors are available) | LSD1i + PARPi / CHK1i |
| 6 | Meiotic programme reactivation | ⭐⭐⭐⭐⭐ | Moderate (expression of meiotic markers) | Targeting meiotic recombination |

> [!IMPORTANT]
> These are hypothesis-generating insights derived from computational enrichment analysis. They require experimental validation. However, the consistent gene overlaps, high fold-enrichments, and statistical significance support them as strong candidates for mechanistic follow-up studies in OSCC. Each represents a potential angle for an original research contribution.
