"""
te_clan_pull.py
Full-population UniProt pull for a Pfam clan's member families (as exported
from InterPro's "entry matching" TSV), skipping DUF-labeled (domain-of-
unknown-function) entries. Unlike te_bulk_pull.py, this pulls every matching
UniProt entry per Pfam ID (via the /stream endpoint, no pagination cap and no
random subsampling) rather than a fixed-size random sample -- appropriate
when the goal is an unbiased taxonomic census, not a bounded PLAAC input set.

No PLAAC scoring here by design (see conversation this was built from) --
this script only produces the protein pull + taxonomic annotation. Scoring
is a deliberately separate, later step once the mega-domain/clustering
questions this pull exposes are resolved.
"""

import argparse
import json
import logging
import re
import time
from pathlib import Path

import pandas as pd
import requests

from te_bulk_pull import CONFIG, BASE_DIR
from te_taxon_analysis import annotate_taxon
from plaac_utils import setup_logging

log = logging.getLogger(__name__)

UNIPROT_STREAM_API = "https://rest.uniprot.org/uniprotkb/stream"

# Manually curated from Pfam Clan CL0523 (entry-matching-CL0523.tsv), after
# excluding the 5 DUF-labeled entries in that export (PF03564, PF13961,
# PF23055, PF23309, PF30931). Split into two biologically distinct groups:
#
#   te_capsid            -- active retrotransposon Gag/capsid domains, still
#                            found as repeated TE copies in genomes.
#   domesticated_host_gene -- host genes domesticated from an ancestral Gag
#                            capsid millions of years ago (Arc, PEG10, PNMA,
#                            RTL1); fixed single-copy eutherian mammal genes,
#                            no longer transposons. Expect near-total mammal
#                            taxonomic skew here -- that's confirmation the
#                            classification is right, not a TE finding.
CLAN_CL0523_DOMAINS: dict[str, dict] = {
    "PF01021": {"label": "TYA_Ty1_capsid", "group": "te_capsid"},
    "PF03732": {"label": "PEG10_N-capsid", "group": "te_capsid"},
    "PF14223": {"label": "Retrotran_gag_2", "group": "te_capsid"},
    "PF14244": {"label": "Retrotran_gag_3", "group": "te_capsid"},
    "PF17241": {"label": "Ty5_gag", "group": "te_capsid"},
    "PF19259": {"label": "Ty3_capsid", "group": "te_capsid"},
    "PF14893": {"label": "PNMA_C", "group": "domesticated_host_gene"},
    "PF16297": {"label": "RTL1-8_LDOC", "group": "domesticated_host_gene"},
    "PF18162": {"label": "Arc_C-lobe", "group": "domesticated_host_gene"},
    "PF21395": {"label": "Arc_N-lobe", "group": "domesticated_host_gene"},
    "PF30901": {"label": "PEG10-RTL1_C-capsid", "group": "domesticated_host_gene"},  # 0 UniProt hits as of 2026-07-10
}

DUF_PATTERN = re.compile(r"\bDUF\d+\b|domain of unknown function|protein of unknown function", re.IGNORECASE)

DATA_DIR = BASE_DIR / "data"
CLAN_CACHE_DIR = DATA_DIR / "clan_cache"
CLAN_FASTA_DIR = DATA_DIR / "clan_fasta"
CLAN_RESULTS_DIR = DATA_DIR / "results"

MIN_LEN = 60   # >= PLAAC core_length (50 in te_bulk_pull.py's CONFIG is inconsistent with
               # core_length=60 -- a protein shorter than the core window can never be
               # called PrLD-positive, so it shouldn't be in scope here either)
MAX_LEN = 3000


def parse_entry_matching_tsv(tsv_path: str) -> pd.DataFrame:
    """Parse an InterPro 'entry matching' export and flag DUF-labeled rows."""
    df = pd.read_csv(tsv_path, sep="\t")
    df = df[df["Source Database"].str.lower() == "pfam"].copy()
    df["is_duf"] = df["Name"].str.contains(DUF_PATTERN, na=False)
    log.info("Parsed %d Pfam entries from %s (%d DUF, %d non-DUF)",
             len(df), tsv_path, df["is_duf"].sum(), (~df["is_duf"]).sum())
    return df


def _uniprot_stream_fetch(pfam_id: str, cache_path: Path) -> list[dict]:
    """Fetch the FULL set of UniProt entries for a Pfam ID via /stream (no
    pagination cap, no subsampling). Cached to disk as JSON."""
    if cache_path.exists():
        with open(cache_path) as f:
            return json.load(f)

    query = f"xref:pfam-{pfam_id} AND length:[{MIN_LEN} TO {MAX_LEN}]"
    r = requests.get(UNIPROT_STREAM_API, params={
        "query": query,
        "fields": "accession,protein_name,organism_name,length,sequence",
        "format": "tsv",
    }, timeout=600)
    r.raise_for_status()

    records = []
    lines = r.text.strip().split("\n")
    for line in lines[1:]:
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


def fetch_clan_set(domains: dict[str, dict], cache_dir: str, skip_fetch: bool) -> pd.DataFrame:
    """Full-population pull (no subsampling) for every domain in `domains`."""
    rows = []
    for pfam_id, meta in domains.items():
        cache_path = Path(cache_dir) / f"clan_{pfam_id}.json"
        if skip_fetch and not cache_path.exists():
            log.warning("MISS: %s not cached and --skip-fetch set; skipping", pfam_id)
            continue
        if not skip_fetch:
            log.info("Fetching full population for %s (%s) ...", pfam_id, meta["label"])
        records = _uniprot_stream_fetch(pfam_id, cache_path)
        log.info("Domain %s (%s, %s): %d proteins", pfam_id, meta["label"], meta["group"], len(records))
        for rec in records:
            rows.append({**rec, "domain": meta["label"], "pfam_id": pfam_id, "clan_group": meta["group"]})
        if not skip_fetch:
            time.sleep(0.3)

    df = pd.DataFrame(rows).drop_duplicates(subset=["accession", "domain"])
    log.info("Clan set: %d protein-domain rows (%d unique accessions) across %d domains",
              len(df), df["accession"].nunique(), len(domains))
    return df


def write_fasta(df: pd.DataFrame, path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    stable_ids = []
    with open(path, "w") as f:
        for i, row in enumerate(df.itertuples(), start=1):
            stable_id = f"CLAN_{i}_{row.accession}"
            stable_ids.append(stable_id)
            f.write(f">{stable_id}\n{row.sequence}\n")
    df["stable_id"] = stable_ids
    log.info("Wrote %d sequences to %s", len(df), path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Full-population UniProt pull for a Pfam clan's non-DUF member "
                    "families, tagged by domain/pfam_id/clan_group and annotated with "
                    "taxonomic lineage. No PLAAC scoring in this script."
    )
    parser.add_argument("--tsv", default=str(Path.home() / "Downloads" / "entry-matching-CL0523.tsv"),
                        help="InterPro 'entry matching' TSV export for the clan")
    parser.add_argument("--out-tag", default="clan_cl0523")
    parser.add_argument("--cache-dir", default=str(CLAN_CACHE_DIR))
    parser.add_argument("--fasta-dir", default=str(CLAN_FASTA_DIR))
    parser.add_argument("--results-dir", default=str(CLAN_RESULTS_DIR))
    parser.add_argument("--skip-fetch", action="store_true",
                        help="Use only cached UniProt data; make no network calls")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(args.verbose)
    log.info("=== TE clan pull starting ===")

    # Cross-check the hardcoded CLAN_CL0523_DOMAINS classification against the TSV's
    # own DUF flags, so a future re-export that changes DUF status doesn't go unnoticed.
    tsv_df = parse_entry_matching_tsv(args.tsv)
    tsv_non_duf = set(tsv_df.loc[~tsv_df["is_duf"], "Accession"])
    hardcoded = set(CLAN_CL0523_DOMAINS)
    if tsv_non_duf != hardcoded:
        log.warning("TSV non-DUF set differs from hardcoded CLAN_CL0523_DOMAINS: "
                     "in TSV only=%s, in hardcoded only=%s",
                     tsv_non_duf - hardcoded, hardcoded - tsv_non_duf)

    df = fetch_clan_set(CLAN_CL0523_DOMAINS, args.cache_dir, args.skip_fetch)
    if df.empty:
        log.error("No proteins fetched; aborting")
        return

    fasta_path = str(Path(args.fasta_dir) / f"{args.out_tag}_combined.fasta")
    write_fasta(df, fasta_path)

    df = annotate_taxon(df, args.cache_dir, args.out_tag)

    Path(args.results_dir).mkdir(parents=True, exist_ok=True)
    out_csv = Path(args.results_dir) / f"{args.out_tag}_proteins.csv"
    df.drop(columns=["sequence"]).to_csv(out_csv, index=False)
    log.info("Wrote %d protein-domain rows to %s (sequence excluded; see %s for sequences)",
              len(df), out_csv, fasta_path)


if __name__ == "__main__":
    main()
