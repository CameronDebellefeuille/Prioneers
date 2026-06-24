"""
te_bulk_analysis.py
Reads the Excel file written by te_bulk_pull.py and produces a TE-only
per-domain PrLD-positive breakdown and PLAAC per-residue plots.

No control group is involved here yet — that comparison will be added back
once the control-set pull is reimplemented.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd

from te_bulk_pull import CONFIG
from plaac_utils import run_plaac_plots, setup_logging

log = logging.getLogger(__name__)


def load_proteins(xlsx_path: str) -> pd.DataFrame:
    df = pd.read_excel(xlsx_path, sheet_name="proteins")
    log.info("Loaded %d proteins from %s", len(df), xlsx_path)
    return df


def per_domain_breakdown(df: pd.DataFrame) -> pd.DataFrame:
    """TE-only PrLD-positive count/rate per Pfam domain (no control rate)."""
    per_domain = (df.groupby("domain")["prd_called"]
                  .agg(["sum", "count"]).rename(columns={"sum": "positive", "count": "n"}))
    per_domain["rate"] = per_domain["positive"] / per_domain["n"]
    per_domain = per_domain.reset_index()
    log.info("\nPer-domain breakdown:\n%s", per_domain.to_string())
    return per_domain


def plot_prd_positive(df: pd.DataFrame, fasta_path: str, alpha: float,
                       core_length: int, plots_dir: str, rscript: str) -> None:
    """Generate a PLAAC per-residue PNG for every protein PLAAC actually called
    a PrLD in (prd_called==1), organized into its own folder."""
    positive_ids = df.loc[df["prd_called"] == 1, "stable_id"].dropna().tolist()
    if not positive_ids:
        log.info("No PrLD-positive proteins to plot.")
        return
    log.info("Plotting %d PrLD-positive proteins to %s", len(positive_ids), plots_dir)
    run_plaac_plots(fasta_path, CONFIG["plaac_jar"], CONFIG["plaac_plot_r"],
                    positive_ids, plots_dir, alpha=alpha, core_length=core_length,
                    rscript=rscript, plot_all=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze the bulk-pull Excel output: TE-only per-domain "
                    "PrLD breakdown and PLAAC per-residue plots."
    )
    parser.add_argument("--xlsx", default=str(Path(CONFIG["results_dir"]) / "bulk_proteins.xlsx"))
    parser.add_argument("--fasta", default=str(Path(CONFIG["fasta_dir"]) / "bulk_combined.fasta"))
    parser.add_argument("--alpha", type=float, default=CONFIG["alpha"])
    parser.add_argument("--results-dir", default=CONFIG["results_dir"])
    parser.add_argument("--prld-plots-dir", default=CONFIG["prld_plots_dir"])
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(args.verbose)
    log.info("=== TE bulk-pull analysis starting ===")

    df = load_proteins(args.xlsx)

    per_domain = per_domain_breakdown(df)
    Path(args.results_dir).mkdir(parents=True, exist_ok=True)
    per_domain.to_csv(Path(args.results_dir) / "per_domain_breakdown.csv", index=False)

    plot_prd_positive(df, args.fasta, args.alpha, CONFIG["core_length"],
                      args.prld_plots_dir, CONFIG["rscript"])


if __name__ == "__main__":
    main()
