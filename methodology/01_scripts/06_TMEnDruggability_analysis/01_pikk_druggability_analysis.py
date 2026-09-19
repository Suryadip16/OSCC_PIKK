#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PIKK Pathway Druggability & Therapeutic Actionability Analysis in OSCC      ║
║  Script: 01_pikk_druggability_analysis.py                                    ║
║  Source: DGIdb v5.0 GraphQL API (https://dgidb.org/api/graphql)              ║
║  Purpose: Query drug interactions, clinical approval tiers, and compute      ║
║           multi-signal composite druggability scores for PIKK targets.       ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import os
import sys
import math
import time
import logging
from pathlib import Path
from typing import Iterable, NamedTuple
import pandas as pd

# ── Paths & Output Directories ────────────────────────────────────────────────
BASE_DIR = Path("C:/Users/dipak/Desktop/Riku/NyberMan/OSCC_PIKK_Project")
OUT_DIR = BASE_DIR / "05_results" / "localisation_outputs" / "tables"
LOG_DIR = BASE_DIR / "05_results" / "localisation_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

# Configure dual logging (Console + Persistent Log File)
log_file_path = LOG_DIR / "dgidb_druggability_run.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_file_path, mode="w", encoding="utf-8")
    ]
)
logger = logging.getLogger("DGIdb_Druggability")

# ── Target Gene Lists ─────────────────────────────────────────────────────────
# Core 7-gene Diamond Panel
DIAMOND_PANEL = ["PLK1", "CDK2", "TOPBP1", "RAD51", "FANCI", "KAT2B", "DEPTOR"]

# All 24 Hub Genes for broader comparative landscape
ALL_HUB_GENES = [
    "PLK1", "CDK2", "TOPBP1", "RAD51", "FANCI", "KAT2B", "DEPTOR", "EGFR",
    "BRCA1", "BRCA2", "PARP1", "H2AX", "RUVBL1", "CHEK1", "CHEK2", "AURKA",
    "AURKB", "E2F1", "EXO1", "CDC45", "MRGBP", "RPTOR", "MLST8", "RICTOR"
]

# ── DGIdb Configuration & GraphQL Queries ─────────────────────────────────────
DGIDB_URL = "https://dgidb.org/api/graphql"
DGIDB_REQUEST_TIMEOUT = 20.0
DGIDB_BATCH_SIZE = 50

DGIDB_GRAPHQL_QUERY = """
query($names: [String!]!) {
  genes(names: $names) {
    nodes {
      name
      interactions {
        drug { name conceptId }
        interactionScore
        interactionTypes { type directionality }
        interactionAttributes { name value }
        publications { pmid }
        sources { sourceDbName }
      }
    }
  }
}
"""

DGIDB_DRUGS_QUERY = """
query($names: [String!]!) {
  drugs(names: $names) {
    nodes {
      name
      approved
      drugAttributes {
        name
        value
      }
    }
  }
}
"""

# Regulatory Approval Mapping Tiers
APPROVAL_TIER_MAP = {
    "fda approved": 1.0,
    "fda-approved": 1.0,
    "fda approved - on companion diagnostic": 1.0,
    "fda approved - has companion diagnostic": 1.0,
    "guideline": 0.90,
    "phase iii": 0.75,
    "phase ii": 0.60,
    "phase ib/ii": 0.60,
    "phase i": 0.45,
    "clinical study": 0.30,
    "case reports/case series": 0.30,
    "preclinical": 0.15,
    "preclinical - biochemical": 0.15,
    "preclinical - cell culture": 0.15,
    "preclinical - cell line xenograft": 0.15,
    "preclinical - pdx": 0.15,
}

KNOWN_SOURCES = frozenset({
    "PharmGKB", "CIViC", "OncoKB", "CGI", "DoCM", "FDA", "CKB-CORE",
    "ClearityFoundationBiomarkers", "DrugBank", "ChemblInteractions",
    "TTD", "TdgClinicalTrial", "TALC", "NCI", "MyCancerGenome"
})

DRUG_SALT_SUFFIXES = (
    "DIHYDROCHLORIDE", "HYDROCHLORIDE", "DISODIUM", "SODIUM", "POTASSIUM",
    "CALCIUM", "MALEATE", "MESYLATE", "TOSYLATE", "BESYLATE", "SULFATE",
    "SULPHATE", "PHOSPHATE", "ACETATE", "CITRATE", "TARTRATE", "FUMARATE",
    "SUCCINATE", "FREE BASE", "FREE ACID", "MONOHYDRATE", "TRIHYDRATE"
)

# ── Helper Functions ──────────────────────────────────────────────────────────
class DrugFacts(NamedTuple):
    approved: bool
    indications: tuple[str, ...]

def tier_of(val: str) -> float | None:
    text = val.strip().lower()
    if not text:
        return None
    if text in APPROVAL_TIER_MAP:
        return APPROVAL_TIER_MAP[text]
    for k in sorted(APPROVAL_TIER_MAP, key=len, reverse=True):
        if text.startswith(k):
            return APPROVAL_TIER_MAP[k]
    return None

def strip_salt(name: str) -> str:
    upper = name.strip().upper()
    for s in DRUG_SALT_SUFFIXES:
        if upper.endswith(" " + s):
            return upper[: -(len(s) + 1)].strip()
    return upper

# ── Query DGIdb ───────────────────────────────────────────────────────────────
def run_dgidb_query(gene_list: list[str]) -> dict:
    try:
        import requests
        from requests.adapters import HTTPAdapter
        from urllib3.util.retry import Retry
    except ImportError:
        logger.error("The 'requests' package is required. Run: pip install requests")
        sys.exit(1)

    session = requests.Session()
    retry_cfg = Retry(total=3, backoff_factor=1.5, status_forcelist=[429, 500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retry_cfg))

    logger.info(f"Querying DGIdb for {len(gene_list)} genes...")
    out = {}

    # Pass 1: Gene Interactions
    try:
        resp = session.post(
            DGIDB_URL,
            json={"query": DGIDB_GRAPHQL_QUERY, "variables": {"names": gene_list}},
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            timeout=DGIDB_REQUEST_TIMEOUT
        )
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        logger.error(f"GraphQL request failed: {e}")
        return {}

    nodes = ((data.get("data") or {}).get("genes") or {}).get("nodes") or []
    all_drugs_to_fetch = set()

    for node in nodes:
        gene = (node.get("name") or "").strip().upper()
        if not gene:
            continue

        interactions = node.get("interactions") or []
        drugs_list = []
        interaction_scores = []
        approval_values = []
        source_names = set()
        moa_list = []
        interaction_types = []

        for inter in interactions:
            drug_name = ((inter.get("drug") or {}).get("name") or "").strip()
            if drug_name:
                drugs_list.append(drug_name)
                all_drugs_to_fetch.add(drug_name.upper())

            score = inter.get("interactionScore")
            if score is not None:
                try:
                    interaction_scores.append(float(score))
                except (ValueError, TypeError):
                    pass

            for attr in inter.get("interactionAttributes") or []:
                aname = (attr.get("name") or "").strip().lower()
                aval = (attr.get("value") or "").strip()
                if aname == "approval status" and aval:
                    approval_values.append(aval)
                if aname == "mechanism of action" and aval:
                    moa_list.append(aval)

            for src in inter.get("sources") or []:
                sname = (src.get("sourceDbName") or "").strip()
                if sname:
                    source_names.add(sname)

            for itype in inter.get("interactionTypes") or []:
                tval = (itype.get("type") or "").strip()
                if tval:
                    interaction_types.append(tval)

        drugs_unique = sorted(set(drugs_list))
        max_score = max(interaction_scores) if interaction_scores else 0.0
        matched_sources = source_names & KNOWN_SOURCES
        source_div = len(matched_sources) / len(KNOWN_SOURCES) if KNOWN_SOURCES else 0.0

        out[gene] = {
            "gene_symbol": gene,
            "raw_interactions": interactions,
            "drugs": drugs_unique,
            "n_interactions": len(interactions),
            "interaction_score": max_score,
            "approval_status": "; ".join(sorted(set(approval_values))),
            "source_databases": "; ".join(sorted(source_names)),
            "source_diversity": source_div,
            "mechanisms_of_action": "; ".join(sorted(set(moa_list))),
            "interaction_types": "; ".join(sorted(set(interaction_types))),
            "approved_drugs": [],
            "n_approved_drugs": 0,
            "indications": ""
        }

    # Pass 2: Curated Drug Approvals
    if all_drugs_to_fetch:
        logger.info(f"Fetching curated regulatory data for {len(all_drugs_to_fetch)} unique drugs...")
        drug_names = list(all_drugs_to_fetch)
        drug_facts = {}

        for i in range(0, len(drug_names), DGIDB_BATCH_SIZE):
            chunk = drug_names[i : i + DGIDB_BATCH_SIZE]
            try:
                d_resp = session.post(
                    DGIDB_URL,
                    json={"query": DGIDB_DRUGS_QUERY, "variables": {"names": chunk}},
                    headers={"Content-Type": "application/json", "Accept": "application/json"},
                    timeout=DGIDB_REQUEST_TIMEOUT
                )
                if d_resp.ok:
                    d_body = d_resp.json()
                    d_nodes = ((d_body.get("data") or {}).get("drugs") or {}).get("nodes") or []
                    for dn in d_nodes:
                        name = str(dn.get("name") or "").strip().upper()
                        if name:
                            inds = []
                            for attr in dn.get("drugAttributes") or []:
                                if (attr.get("name") or "").strip().lower() == "indication":
                                    inds.extend((attr.get("value") or "").split(";"))
                            drug_facts[name] = DrugFacts(
                                approved=bool(dn.get("approved")),
                                indications=tuple(sorted(set(x.strip() for x in inds if x.strip())))
                            )
            except Exception as e:
                logger.warning(f"Drug batch query error: {e}")

        # Update per-gene records
        for gene, rec in out.items():
            approved = []
            gene_indications = set()
            for d in rec["drugs"]:
                df = drug_facts.get(d.upper()) or drug_facts.get(strip_salt(d))
                if df:
                    if df.approved:
                        approved.append(d)
                    gene_indications.update(df.indications)

            rec["approved_drugs"] = sorted(set(approved))
            rec["n_approved_drugs"] = len(rec["approved_drugs"])
            rec["indications"] = "; ".join(sorted(gene_indications))

    session.close()
    return out

# ── Compute Composite Druggability Score ──────────────────────────────────────
def compute_druggability_scores(records: dict) -> pd.DataFrame:
    rows = []
    
    for gene, rec in records.items():
        n_drugs = len(rec["drugs"])
        inter_score = rec["interaction_score"]
        source_div = rec["source_diversity"]

        # Signal 1: Normalized Interaction Score (log-scaled, cap 100)
        sig_int = math.log1p(min(inter_score, 100.0)) / math.log1p(100.0)

        # Signal 2: Approval Status Tier
        best_tier = 0.0
        for part in rec["approval_status"].split(";"):
            t = tier_of(part)
            if t is not None:
                best_tier = max(best_tier, t)
        if rec["n_approved_drugs"] > 0:
            best_tier = 1.0  # Curated approved drug overrides text tier
        sig_app = best_tier

        # Signal 3: Drug Count (log-scaled, cap 20)
        sig_drugs = math.log1p(min(n_drugs, 20)) / math.log1p(20)

        # Signal 4: Source Diversity
        sig_src = min(max(source_div, 0.0), 1.0)

        # Composite Score (Weights: App 0.35, Int 0.30, Src 0.20, Drugs 0.15)
        w_app, w_int, w_src, w_drugs = 0.35, 0.30, 0.20, 0.15
        composite = (
            (w_app * sig_app) +
            (w_int * sig_int) +
            (w_src * sig_src) +
            (w_drugs * sig_drugs)
        )

        # Categorical Tier Label
        if sig_app >= 1.0:
            maturity = "FDA Approved"
        elif sig_app >= 0.70:
            maturity = "Phase III Clinical"
        elif sig_app >= 0.40:
            maturity = "Phase I/II Clinical"
        elif sig_app > 0.0 or n_drugs > 0:
            maturity = "Preclinical / Experimental"
        else:
            maturity = "Undruggable / Novel"

        is_diamond = "Yes" if gene in DIAMOND_PANEL else "No"

        rows.append({
            "gene_symbol": gene,
            "is_diamond_panel": is_diamond,
            "composite_druggability_score": round(composite, 4),
            "clinical_maturity": maturity,
            "approval_tier_score": round(sig_app, 4),
            "interaction_subscore": round(sig_int, 4),
            "source_diversity_subscore": round(sig_src, 4),
            "drug_count_subscore": round(sig_drugs, 4),
            "n_drugs": n_drugs,
            "n_approved_drugs": rec["n_approved_drugs"],
            "approved_drugs": "; ".join(rec["approved_drugs"][:10]),
            "top_candidate_drugs": "; ".join(rec["drugs"][:8]),
            "interaction_types": rec["interaction_types"],
            "mechanisms_of_action": rec["mechanisms_of_action"],
            "indications": rec["indications"][:150]
        })

    df = pd.DataFrame(rows).sort_values("composite_druggability_score", ascending=False)
    return df

# ── Main Execution ────────────────────────────────────────────────────────────
def main():
    print("\n" + "═"*70)
    print("  PIKK Pathway Druggability Analysis (DGIdb GraphQL Engine)")
    print("═"*70 + "\n")

    # 1. Query DGIdb
    all_genes = sorted(set(ALL_HUB_GENES + DIAMOND_PANEL))
    records = run_dgidb_query(all_genes)

    if not records:
        logger.error("No records returned. Check network connection.")
        sys.exit(1)

    # 2. Compute Scores
    df_scores = compute_druggability_scores(records)

    # 3. Export Tables
    score_file = OUT_DIR / "pikk_druggability_composite_scores.tsv"
    df_scores.to_csv(score_file, sep="\t", index=False)
    logger.info(f"Saved composite druggability scores: {score_file}")

    # Export raw interactions detail
    raw_rows = []
    for gene, rec in records.items():
        for inter in rec["raw_interactions"]:
            drug_name = ((inter.get("drug") or {}).get("name") or "").strip()
            score = inter.get("interactionScore")
            srcs = [str(s.get("sourceDbName") or "") for s in (inter.get("sources") or [])]
            raw_rows.append({
                "gene_symbol": gene,
                "drug_name": drug_name,
                "interaction_score": score,
                "sources": "; ".join(filter(None, srcs))
            })
    if raw_rows:
        raw_file = OUT_DIR / "dgidb_raw_gene_interactions.tsv"
        pd.DataFrame(raw_rows).to_csv(raw_file, sep="\t", index=False)
        logger.info(f"Saved raw interaction records: {raw_file}")

    # 4. Print Summary Table
    print("\n" + "─"*70)
    print("  PIKK DIAMOND PANEL DRUGGABILITY SUMMARY")
    print("─"*70)
    diamond_summary = df_scores[df_scores["is_diamond_panel"] == "Yes"][
        ["gene_symbol", "composite_druggability_score", "clinical_maturity", "n_drugs", "approved_drugs"]
    ]
    print(diamond_summary.to_string(index=False))
    print("─"*70 + "\n")
    logger.info("Druggability analysis complete!")

if __name__ == "__main__":
    main()
