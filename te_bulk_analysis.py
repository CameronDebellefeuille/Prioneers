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
    parser.add_argument("--out-tag", default="bulk",
                        help="Filename prefix matching the te_bulk_pull.py run to analyze")
    parser.add_argument("--xlsx", default=None,
                        help="Defaults to {results_dir}/{out-tag}_proteins.xlsx")
    parser.add_argument("--fasta", default=None,
                        help="Defaults to {fasta_dir}/{out-tag}_combined.fasta")
    parser.add_argument("--alpha", type=float, default=CONFIG["alpha"])
    parser.add_argument("--results-dir", default=CONFIG["results_dir"])
    parser.add_argument("--prld-plots-dir", default=None,
                        help="Defaults to {prld_plots_dir}_{out-tag}")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    xlsx_path = args.xlsx or str(Path(CONFIG["results_dir"]) / f"{args.out_tag}_proteins.xlsx")
    fasta_path = args.fasta or str(Path(CONFIG["fasta_dir"]) / f"{args.out_tag}_combined.fasta")
    prld_plots_dir = args.prld_plots_dir or f"{CONFIG['prld_plots_dir']}_{args.out_tag}"

    setup_logging(args.verbose)
    log.info("=== TE bulk-pull analysis starting ===")

    df = load_proteins(xlsx_path)

    per_domain = per_domain_breakdown(df)
    Path(args.results_dir).mkdir(parents=True, exist_ok=True)
    per_domain.to_csv(Path(args.results_dir) / f"{args.out_tag}_per_domain_breakdown.csv", index=False)

    plot_prd_positive(df, fasta_path, args.alpha, CONFIG["core_length"],
                      prld_plots_dir, CONFIG["rscript"])


if __name__ == "__main__":
    main()
