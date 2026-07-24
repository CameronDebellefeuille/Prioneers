"""
te_clan_plaac_rescore.py
Rescores the clan_cl0523 protein pull (te_clan_pull.py) with PLAAC at a given
alpha, reproducing the recipe used for the existing alpha=0.5 run
(clan_cl0523_proteins_plaac_clustered.csv):

  1. Dedup the ~157k protein-domain rows to ~155k unique accessions (each
     accession maps to exactly one sequence in this pull -- verified: 0
     accessions have >1 distinct length across their domain rows).
  2. Run PLAAC once per unique accession at the given --alpha (core_length
     fixed at 60, matching every other PLAAC run in this repo).
  3. Merge the PLAAC score back onto every protein-domain row by accession,
     so a protein sampled under more than one clan domain carries the same
     score into each of its rows.
  4. Reattach the existing per-domain MMseqs2 cluster_id/cluster_size (95%
     identity, 80% coverage, clustered independently within each Pfam
     domain's own set) unchanged -- clustering is sequence-identity-based,
     not alpha-dependent, so it is not recomputed here.

--mode benchmark scores a seeded random subset of accessions, for measuring
PLAAC's real throughput at this input size and validating the reproduction
(the alpha=0.5 benchmark subset should reproduce the existing scored file's
values for the same accessions) before committing to the full run.
"""

import argparse
import logging
import random
import time
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from te_bulk_pull import CONFIG, BASE_DIR
from plaac_utils import run_plaac, _parse_plaac_tsv, setup_logging

log = logging.getLogger(__name__)

DATA_DIR = BASE_DIR / "data"
CLAN_FASTA = DATA_DIR / "clan_fasta" / "clan_cl0523_combined.fasta"
CLAN_CLUSTERED_CSV = DATA_DIR / "results" / "clan_cl0523_proteins_plaac_clustered.csv"
RESCORE_FASTA_DIR = DATA_DIR / "clan_fasta" / "rescore"
RESCORE_PLAAC_DIR = DATA_DIR / "clan_plaac_output"
RESULTS_DIR = DATA_DIR / "results"

CORE_LENGTH = 60  # fixed, matches every other PLAAC run in this repo
PLAAC_COLS = ["accession", "prd_called", "PRDlen", "COREscore", "PROTlen"]


def load_base_table() -> pd.DataFrame:
    """All non-PLAAC columns from the existing clustered file: taxon/domain/
    cluster metadata to reattach to freshly-scored rows."""
    df = pd.read_csv(CLAN_CLUSTERED_CSV)
    keep = [c for c in df.columns if c not in ("prd_called", "PRDlen", "COREscore", "PROTlen")]
    log.info("Loaded %d protein-domain rows (%d unique accessions) from %s",
              len(df), df["accession"].nunique(), CLAN_CLUSTERED_CSV)
    return df[keep]


def load_sequences() -> dict[str, str]:
    """stable_id -> sequence, from the clan combined FASTA."""
    seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(CLAN_FASTA), "fasta")}
    log.info("Loaded %d sequences from %s", len(seqs), CLAN_FASTA)
    return seqs


def unique_accession_table(base: pd.DataFrame) -> pd.DataFrame:
    """One row per accession (first occurrence), carrying its stable_id so the
    matching sequence can be looked up."""
    return base.drop_duplicates(subset="accession")[["accession", "stable_id"]]


def write_rescore_fasta(acc_table: pd.DataFrame, seqs: dict[str, str], path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(path, "w") as f:
        for row in acc_table.itertuples():
            seq = seqs.get(row.stable_id)
            if seq is None:
                log.warning("No sequence found for stable_id %s (accession %s); skipping",
                            row.stable_id, row.accession)
                continue
            f.write(f">{row.accession}\n{seq}\n")
            n += 1
    log.info("Wrote %d unique-accession sequences to %s", n, path)
    return n


def score(fasta_path: Path, alpha: float, output_dir: Path, timeout: int, heap: str) -> pd.DataFrame:
    tsv_path = run_plaac(str(fasta_path), CONFIG["plaac_jar"], alpha, CORE_LENGTH,
                          str(output_dir), timeout=timeout, heap=heap)
    if tsv_path is None:
        raise SystemExit("PLAAC run failed -- see the .err file next to the output TSV")
    plaac_df = _parse_plaac_tsv(tsv_path)
    if plaac_df is None:
        raise SystemExit("Could not parse PLAAC output")
    plaac_df["accession"] = plaac_df["SEQid"].str.split().str[0]
    plaac_df["prd_called"] = (plaac_df["PRDlen"].fillna(0).astype(float) > 0).astype(int)
    return plaac_df[PLAAC_COLS]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rescore the clan_cl0523 protein pull with PLAAC at a given "
                    "alpha, merging results back onto every protein-domain row "
                    "and reattaching the existing MMseqs2 clustering unchanged."
    )
    parser.add_argument("--alpha", type=float, required=True)
    parser.add_argument("--mode", choices=["benchmark", "full"], default="full")
    parser.add_argument("--subset-n", type=int, default=5000,
                        help="Benchmark mode only: number of unique accessions to sample")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--timeout", type=int, default=300,
                        help="PLAAC subprocess timeout in seconds -- raise this a lot for --mode full")
    parser.add_argument("--heap", default="6g", help="JVM -Xmx for the PLAAC process")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(args.verbose)
    alpha_tag = str(args.alpha).replace(".", "p")
    log.info("=== Clan PLAAC rescore starting (alpha=%s, mode=%s) ===", args.alpha, args.mode)

    base = load_base_table()
    seqs = load_sequences()
    acc_table = unique_accession_table(base)

    if args.mode == "benchmark":
        rng = random.Random(args.seed)
        sample_accessions = set(rng.sample(list(acc_table["accession"]), min(args.subset_n, len(acc_table))))
        acc_table = acc_table[acc_table["accession"].isin(sample_accessions)]
        base = base[base["accession"].isin(sample_accessions)]
        tag = f"benchmark_a{alpha_tag}"
    else:
        tag = f"a{alpha_tag}"

    fasta_path = RESCORE_FASTA_DIR / f"{tag}_unique.fasta"
    write_rescore_fasta(acc_table, seqs, fasta_path)

    output_dir = RESCORE_PLAAC_DIR / tag
    start = time.monotonic()
    scored = score(fasta_path, args.alpha, output_dir, args.timeout, args.heap)
    elapsed = time.monotonic() - start
    log.info("PLAAC scored %d unique accessions in %.1fs (%.4fs/accession)",
              len(scored), elapsed, elapsed / max(len(scored), 1))

    merged = base.merge(scored, on="accession", how="left")
    missing = merged["prd_called"].isna().sum()
    if missing:
        log.warning("%d/%d rows have no PLAAC score after merge (sequence lookup or PLAAC parse gap)",
                    missing, len(merged))

    out_csv = RESULTS_DIR / f"clan_cl0523_proteins_plaac_{tag}.csv"
    merged.to_csv(out_csv, index=False)
    log.info("Wrote %d rows to %s", len(merged), out_csv)
    log.info("=== Done: alpha=%s, %d accessions scored in %.1fs ===", args.alpha, len(scored), elapsed)


if __name__ == "__main__":
    main()
