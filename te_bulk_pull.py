"""
te_bulk_pull.py
Pulls hundreds of real TE-encoded proteins in bulk from UniProt (by Pfam
domain, across taxonomy), scores them with PLAAC, and writes the merged
per-protein results to an Excel file for easy review.

Statistical comparison and plotting live in te_bulk_analysis.py.
"""

import argparse # reads command-line flags
import json # caches UniProt results to/from disk as JSON
import logging # logs progress/info messages instead of using print()
import os # reads RSCRIPT_PATH env var override
import random # seeded random subsampling of UniProt pools
import sys # sys.exit() for fatal error paths
from pathlib import Path # cross-platform filesystem path handling

import pandas as pd # tabular data handling (DataFrames) and Excel output
import requests # makes the HTTP calls to the UniProt REST API

from plaac_utils import run_plaac, _parse_plaac_tsv, setup_logging # runs PLAAC, parses its TSV output, configures logging

log = logging.getLogger(__name__)

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"

# TE-associated Pfam domains, verified against the UniProt API on 2026-06-17
# (each returned real TE proteins on a spot check, not false-positive gene families):
#   PF03732 Retrotrans_gag  - Gag (LTR retrotransposon/endogenous retrovirus capsid)
#   PF00078 RVT_1           - Reverse transcriptase (LINE/LTR retrotransposons)
#   PF00665 rve             - Integrase core domain (LTR retrotransposons)
#   PF03004 DDE_Tnp_1       - Mutator/En-Spm-like DNA transposase
#   PF01498 Transposase_8   - IS4/hAT-like DNA transposase
#   PF03108 MULE            - MuDR/Foldback-like DNA transposase
# Note: no reliable Pfam domain was found specifically for Helitron Rep/helicase
# (PF14214 and similar guesses returned 0 UniProt hits) — Helitrons are
# under-represented here as a result.
TE_PFAM_DOMAINS = {
    "PF03732": "Gag",
    "PF00078": "RT",
    "PF00665": "Integrase",
    "PF03004": "Transposase_DDE",
    "PF01498": "Transposase_IS4",
    "PF03108": "Transposase_MULE",
}

# All paths are anchored to this file's location (not the working directory)
# so the script behaves the same regardless of the caller's cwd.
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"

CONFIG = {
    "cache_dir": str(DATA_DIR / "bulk_cache"),
    "fasta_dir": str(DATA_DIR / "bulk_fasta"),
    "plaac_output_dir": str(DATA_DIR / "bulk_plaac_output"),
    "plaac_jar": str(BASE_DIR / "tools" / "plaac" / "plaac.jar"),
    "plaac_plot_r": str(BASE_DIR / "tools" / "plaac" / "R" / "plaac_plot.r"),
    "results_dir": str(DATA_DIR / "results"),
    "prld_plots_dir": str(DATA_DIR / "results" / "plots_prld_positive"),
    "alpha": 0.5,
    "core_length": 60,
    "min_len": 50,
    "max_len": 3000,
    # Rscript path — defaults to PATH; set RSCRIPT_PATH if R isn't on PATH
    "rscript": os.environ.get("RSCRIPT_PATH", "Rscript"),
}


def _uniprot_fetch(query: str, size: int, cache_path: Path) -> list[dict]:
    """Fetch up to `size` UniProt entries matching query; cache to disk as JSON."""
    if cache_path.exists():
        with open(cache_path) as f:
            return json.load(f)

    r = requests.get(UNIPROT_API, params={
        "query": query,
        "fields": "accession,protein_name,organism_name,length,sequence",
        "format": "tsv",
        "size": size,
    }, timeout=60)
    r.raise_for_status()

    records = []
    lines = r.text.strip().split("\n")
    for line in lines[1:]:  # skip header
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        acc, name, organism, length, seq = parts[0], parts[1], parts[2], parts[3], parts[4]
        records.append({
            "accession": acc, "protein_name": name, "organism": organism,
            "length": int(length), "sequence": seq,
        })

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(records, f, indent=2)
    return records


def fetch_te_set(per_domain: int, fetch_pool: int, cache_dir: str, seed: int) -> pd.DataFrame:
    """Fetch up to `per_domain` proteins per TE Pfam domain (randomly subsampled
    from a larger fetched pool so results aren't just UniProt's default sort order)."""
    rng = random.Random(seed)
    rows = []
    for pfam_id, label in TE_PFAM_DOMAINS.items():
        cache_path = Path(cache_dir) / f"te_{pfam_id}.json"
        query = f"xref:pfam-{pfam_id} AND length:[{CONFIG['min_len']} TO {CONFIG['max_len']}]"
        pool = _uniprot_fetch(query, fetch_pool, cache_path)
        log.info("TE domain %s (%s): fetched pool of %d", pfam_id, label, len(pool))
        sample = rng.sample(pool, min(per_domain, len(pool)))
        for rec in sample:
            rows.append({**rec, "group": "TE", "domain": label, "pfam_id": pfam_id})
    df = pd.DataFrame(rows).drop_duplicates(subset="accession")
    log.info("TE set: %d unique proteins across %d domains", len(df), len(TE_PFAM_DOMAINS))
    return df


def write_fasta(df: pd.DataFrame, path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for i, row in enumerate(df.itertuples(), start=1):
            f.write(f">{row.group}_{i}_{row.accession}\n{row.sequence}\n")
    log.info("Wrote %d sequences to %s", len(df), path)


def score_with_plaac(df: pd.DataFrame, fasta_path: str, alpha: float, dry_run: bool) -> pd.DataFrame:
    tsv_path = run_plaac(fasta_path, CONFIG["plaac_jar"], alpha, CONFIG["core_length"],
                          CONFIG["plaac_output_dir"], dry_run=dry_run)
    if dry_run or tsv_path is None:
        return df

    plaac_df = _parse_plaac_tsv(tsv_path)
    if plaac_df is None:
        log.error("Could not parse PLAAC output; aborting analysis")
        sys.exit(1)

    # SEQid format is "{group}_{i}_{accession} ..."; recover accession by index.
    plaac_df["stable_id"] = plaac_df["SEQid"].str.split().str[0]
    plaac_df["accession"] = plaac_df["stable_id"].str.split("_").str[-1]
    plaac_df["prd_called"] = (plaac_df["PRDlen"].fillna(0).astype(float) > 0).astype(int)

    # Carry the full PLAAC output through (COREscore, COREstart/end, VITmaxrun,
    # FoldIndex/PAPA stats, HMM scores, etc.) rather than a hand-picked subset —
    # everything PLAAC computes is potentially useful for follow-up analysis.
    plaac_cols = [c for c in plaac_df.columns if c not in ("SEQid",)]
    merged = df.merge(plaac_df[plaac_cols], on="accession", how="left")
    return merged


def write_proteins_excel(df: pd.DataFrame, results_dir: str) -> None:
    """Write the same merged table to an .xlsx for quick spreadsheet review."""
    Path(results_dir).mkdir(parents=True, exist_ok=True)
    out_xlsx = str(Path(results_dir) / "bulk_proteins.xlsx")
    df.to_excel(out_xlsx, sheet_name="proteins", index=False)
    log.info("Wrote %d protein rows to %s", len(df), out_xlsx)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bulk-pull TE-encoded proteins from UniProt and score them "
                    "with PLAAC, writing the results to Excel for review."
    )
    parser.add_argument("--per-domain", type=int, default=100,
                        help="Max proteins to sample per TE Pfam domain")
    parser.add_argument("--fetch-pool", type=int, default=500,
                        help="Pool size fetched per query before local random subsampling")
    parser.add_argument("--alpha", type=float, default=CONFIG["alpha"])
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip-fetch", action="store_true",
                        help="Use only cached UniProt data; make no network calls")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(args.verbose)
    log.info("=== TE bulk-pull starting ===")

    if args.skip_fetch:
        cache_dir = Path(CONFIG["cache_dir"])
        missing = [p for p in TE_PFAM_DOMAINS if not (cache_dir / f"te_{p}.json").exists()]
        if missing:
            log.error("--skip-fetch set but cache is incomplete (missing: %s); run once without it first", missing)
            sys.exit(1)

    te_df = fetch_te_set(args.per_domain, args.fetch_pool, CONFIG["cache_dir"], args.seed)

    fasta_path = str(Path(CONFIG["fasta_dir"]) / "bulk_combined.fasta")
    write_fasta(te_df, fasta_path)

    scored = score_with_plaac(te_df, fasta_path, args.alpha, args.dry_run)
    if args.dry_run:
        log.info("[dry-run] Stopping before Excel write")
        return

    proteins = scored.drop(columns=["sequence"])
    write_proteins_excel(proteins, CONFIG["results_dir"])


if __name__ == "__main__":
    main()
