"""
te_taxon_analysis.py
Reads the Excel file written by te_bulk_pull.py, looks up each accession's
NCBI/UniProt taxonomic lineage (cached to disk), classifies every protein
into a broad taxonomic group (Plant/Animal/Fungi/Virus/Bacteria/...) and,
for animals, a finer subgroup (Fish/Mammal/Insect/...), and writes:

  - the merged per-protein table with taxon columns
  - taxonomic composition (protein count per group)
  - domain x taxon composition crosstab (the domain/taxon confound check)
  - PrLD-positive rate by taxon
  - mean protein length by domain (the length-confound check)

This does not replace te_bulk_analysis.py's per-domain breakdown -- it's a
second, taxon-oriented breakdown of the same bulk-pull output.
"""

import argparse
import json
import logging
import re
import time
from pathlib import Path

import pandas as pd
import requests

from te_bulk_pull import CONFIG
from plaac_utils import setup_logging

log = logging.getLogger(__name__)

UNIPROT_ACCESSIONS_API = "https://rest.uniprot.org/uniprotkb/accessions"
LINEAGE_BATCH_SIZE = 100

ANIMAL_SUBGROUPS = [
    ("Actinopterygii", "Fish"), ("Mammalia", "Mammal"), ("Aves", "Bird"),
    ("Insecta", "Insect"), ("Arachnida", "Arachnid"), ("Amphibia", "Amphibian"),
    ("Reptilia", "Reptile"), ("Testudines", "Reptile"), ("Squamata", "Reptile"),
    ("Mollusca", "Mollusc"), ("Crustacea", "Crustacean"),
]

KINGDOM_LABELS = {"Viridiplantae": "Plant", "Metazoa": "Animal", "Fungi": "Fungi"}


def load_proteins(xlsx_path: str) -> pd.DataFrame:
    df = pd.read_excel(xlsx_path, sheet_name="proteins")
    log.info("Loaded %d proteins from %s", len(df), xlsx_path)
    return df


def fetch_lineages(accessions: list[str], cache_path: Path) -> pd.DataFrame:
    """Fetch NCBI/UniProt taxonomic lineage per accession, cached to disk as JSON."""
    cache: dict[str, str] = {}
    if cache_path.exists():
        with open(cache_path) as f:
            cache = json.load(f)

    def _flush():
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = cache_path.with_suffix(".json.tmp")
        with open(tmp, "w") as f:
            json.dump(cache, f, indent=2)
        tmp.replace(cache_path)  # atomic on both POSIX and Windows

    missing = [a for a in accessions if a not in cache]
    if missing:
        log.info("Fetching lineage for %d accessions (%d already cached)",
                  len(missing), len(accessions) - len(missing))
        for i in range(0, len(missing), LINEAGE_BATCH_SIZE):
            batch = missing[i:i + LINEAGE_BATCH_SIZE]

            for attempt in range(5):
                try:
                    r = requests.get(UNIPROT_ACCESSIONS_API, params={
                        "accessions": ",".join(batch),
                        "fields": "accession,lineage",
                        "format": "tsv", "size": 500,
                    }, timeout=60)
                    r.raise_for_status()
                    break
                except (requests.exceptions.RequestException,) as e:
                    if attempt == 4:
                        # Flush everything fetched so far before giving up -- a
                        # crash here must not lose already-completed lookups.
                        _flush()
                        log.error("Lineage fetch failed after retries at batch %d/%d "
                                  "(flushed %d entries fetched so far): %s",
                                  i, len(missing), len(cache), e)
                        raise
                    backoff = 2 ** attempt
                    log.warning("Lineage fetch batch %d/%d failed (attempt %d/5): %s -- retrying in %ds",
                                i, len(missing), attempt + 1, e, backoff)
                    time.sleep(backoff)

            lines = r.text.strip().split("\n")
            for line in lines[1:]:
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                cache[parts[0]] = parts[1]
            log.info("Lineage fetch: %d/%d accessions done", min(i + LINEAGE_BATCH_SIZE, len(missing)), len(missing))

            # Flush periodically so a later failure never loses more than one
            # checkpoint's worth of already-completed lookups.
            if (i // LINEAGE_BATCH_SIZE) % 20 == 0:
                _flush()
            time.sleep(0.3)

        _flush()

    rows = [{"accession": a, "lineage": cache.get(a)} for a in accessions]
    return pd.DataFrame(rows)


def classify_broad_group(lineage: str | float) -> str:
    if not isinstance(lineage, str):
        return "Unknown"
    if lineage.startswith("Viruses"):
        return "Virus"
    m = re.search(r"([A-Za-z]+) \(kingdom\)", lineage)
    if m:
        return KINGDOM_LABELS.get(m.group(1), m.group(1))
    if "Bacteria (domain)" in lineage:
        return "Bacteria"
    if "Archaea (domain)" in lineage:
        return "Archaea"
    if "Eukaryota (domain)" in lineage:
        return "Other Eukaryote"
    return "Unknown"


def classify_animal_subgroup(lineage: str | float) -> str | None:
    if not isinstance(lineage, str) or "Metazoa (kingdom)" not in lineage:
        return None
    for tag, label in ANIMAL_SUBGROUPS:
        if tag in lineage:
            return label
    return "Other Animal"


def annotate_taxon(df: pd.DataFrame, cache_dir: str, out_tag: str) -> pd.DataFrame:
    """Merge broad_group/animal_subgroup/taxon columns onto df (must have 'accession')."""
    accessions = df["accession"].dropna().unique().tolist()
    cache_path = Path(cache_dir) / f"{out_tag}_lineage.json"
    lineage_df = fetch_lineages(accessions, cache_path)

    df = df.merge(lineage_df, on="accession", how="left")
    df["broad_group"] = df["lineage"].apply(classify_broad_group)
    df["animal_subgroup"] = df["lineage"].apply(classify_animal_subgroup)
    df["taxon"] = df.apply(
        lambda r: r["animal_subgroup"] if r["broad_group"] == "Animal" and r["animal_subgroup"] else r["broad_group"],
        axis=1,
    )
    return df


def taxonomic_composition(df: pd.DataFrame) -> pd.DataFrame:
    comp = df.groupby("broad_group").size().rename("n").reset_index()
    comp["pct"] = (100 * comp["n"] / comp["n"].sum()).round(1)
    comp = comp.sort_values("n", ascending=False)
    log.info("\nTaxonomic composition:\n%s", comp.to_string(index=False))
    return comp


def domain_taxon_crosstab(df: pd.DataFrame) -> pd.DataFrame:
    ct = pd.crosstab(df["domain"], df["broad_group"])
    log.info("\nDomain x taxon composition (counts):\n%s", ct.to_string())
    return ct


def prld_rate_by_taxon(df: pd.DataFrame) -> pd.DataFrame:
    rate = (df.groupby("broad_group")["prd_called"].agg(["sum", "count"])
            .rename(columns={"sum": "positive", "count": "n"}))
    rate["rate_pct"] = (100 * rate["positive"] / rate["n"]).round(1)
    rate = rate.sort_values("rate_pct", ascending=False).reset_index()
    log.info("\nPrLD-positive rate by taxon:\n%s", rate.to_string(index=False))
    return rate


def length_by_domain(df: pd.DataFrame) -> pd.DataFrame:
    length = df.groupby("domain")["PROTlen"].agg(["mean", "median", "min", "max"]).round(0)
    length = length.sort_values("mean", ascending=False).reset_index()
    log.info("\nProtein length by domain:\n%s", length.to_string(index=False))
    return length


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Taxonomic breakdown of a te_bulk_pull.py run: composition, "
                    "domain x taxon crosstab (confound check), PrLD rate by taxon, "
                    "and length by domain (length-confound check)."
    )
    parser.add_argument("--out-tag", default="bulk",
                        help="Filename prefix matching the te_bulk_pull.py run to analyze")
    parser.add_argument("--xlsx", default=None,
                        help="Defaults to {results_dir}/{out-tag}_proteins.xlsx")
    parser.add_argument("--results-dir", default=CONFIG["results_dir"])
    parser.add_argument("--cache-dir", default=CONFIG["cache_dir"],
                        help="Where to cache per-accession lineage lookups")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    xlsx_path = args.xlsx or str(Path(CONFIG["results_dir"]) / f"{args.out_tag}_proteins.xlsx")

    setup_logging(args.verbose)
    log.info("=== TE taxon analysis starting ===")

    df = load_proteins(xlsx_path)
    df = annotate_taxon(df, args.cache_dir, args.out_tag)

    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    df.drop(columns=["sequence"], errors="ignore").to_csv(
        results_dir / f"{args.out_tag}_proteins_taxon.csv", index=False)
    taxonomic_composition(df).to_csv(
        results_dir / f"{args.out_tag}_taxon_composition.csv", index=False)
    domain_taxon_crosstab(df).to_csv(
        results_dir / f"{args.out_tag}_domain_taxon_crosstab.csv")
    prld_rate_by_taxon(df).to_csv(
        results_dir / f"{args.out_tag}_prld_rate_by_taxon.csv", index=False)
    length_by_domain(df).to_csv(
        results_dir / f"{args.out_tag}_length_by_domain.csv", index=False)

    log.info("Wrote taxon breakdown tables to %s", results_dir)


if __name__ == "__main__":
    main()
