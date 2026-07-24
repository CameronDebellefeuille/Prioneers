"""
te_clan_taxon_heatmap.py
Builds a multi-panel taxonomic PrLD-positive-rate heatmap from the three
clan_cl0523 PLAAC rescores (alpha=0, 0.5, 1; core_length=60 fixed), produced
by te_clan_plaac_rescore.py.

Design decisions (agreed with Cameron before building this):
  - Unit of analysis is cluster-corrected: one representative row per
    (domain, cluster_id) MMseqs2 cluster (95% identity, 80% coverage), not
    raw protein-domain rows -- raw rows pseudoreplicate near-identical
    paralogs/strains (e.g. TYA_Ty1_capsid's raw rate is inflated by ~40
    near-identical S. cerevisiae Ty1 copies).
  - domesticated_host_gene (Arc/PNMA/RTL1 -- domesticated non-TE host genes)
    is plotted as a separate, clearly labeled reference panel, not blended
    into the main TE PrLD-by-taxon heatmap.
  - Metric is PrLD-positive rate (%), not raw counts -- raw counts would
    mostly replot this pull's taxonomic sampling skew (76% Plant), not
    biology.
  - Cells with cluster-corrected n < MIN_N are masked (grey, hatched) rather
    than shown as a spurious 0%/100% -- several taxa here have n in the
    single digits (Virus n=1, Bacillati/Pseudomonadati n<25).
"""

import argparse
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless: this script only writes a PNG, never shows a window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from te_bulk_pull import BASE_DIR
from plaac_utils import setup_logging

log = logging.getLogger(__name__)

DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = DATA_DIR / "results"

ALPHAS = [("a0p0", 0.0), ("a0p5", 0.5), ("a1p0", 1.0)]
CORE_LENGTH = 60
MIN_N = 20  # cluster-corrected n below this is masked as unreliable

GROUPS = [
    ("te_capsid", "TE capsid domains (active transposons)"),
    ("domesticated_host_gene", "Domesticated host genes (Arc / PNMA / RTL1)"),
]

# The "fine" taxon column (te_taxon_analysis.py's annotate_taxon()) splits Animal
# into Fish/Mammal/Insect/Arachnid/Mollusc/etc. via animal_subgroup. "broad"
# collapses back to broad_group, matching the 5-group convention used in
# clan_fig1-3 (Plant/Animal/Fungi/Other Eukaryote/Bacteria-or-other).
BROAD_GROUPS_KEPT = ["Plant", "Animal", "Fungi", "Other Eukaryote"]
BROAD_OTHER_LABEL = "Bacteria/other"


def load_alpha(tag: str, taxon_level: str) -> pd.DataFrame:
    path = RESULTS_DIR / f"clan_cl0523_proteins_plaac_{tag}.csv"
    df = pd.read_csv(path)
    if taxon_level == "broad":
        df["taxon"] = df["broad_group"].where(
            df["broad_group"].isin(BROAD_GROUPS_KEPT), BROAD_OTHER_LABEL)
    log.info("Loaded %s: %d rows", path.name, len(df))
    return df


def cluster_corrected_rate_table(df: pd.DataFrame, clan_group: str, min_n: int) -> pd.DataFrame:
    """One row per (domain, cluster_id); returns long table of taxon x domain
    PrLD-positive rate/n/reliability."""
    g = df[df["clan_group"] == clan_group]
    cc = g.drop_duplicates(subset=["domain", "cluster_id"])
    tab = (cc.groupby(["taxon", "domain"])["prd_called"]
           .agg(positive="sum", n="count").reset_index())
    tab["rate_pct"] = 100 * tab["positive"] / tab["n"]
    tab["reliable"] = tab["n"] >= min_n
    return tab


def order_by_total_n(tab: pd.DataFrame, col: str) -> list[str]:
    totals = tab.groupby(col)["n"].sum().sort_values(ascending=False)
    return totals.index.tolist()


def pivot_grid(tab: pd.DataFrame, taxon_order: list[str], domain_order: list[str]):
    rate = tab.pivot(index="taxon", columns="domain", values="rate_pct").reindex(
        index=taxon_order, columns=domain_order)
    n = tab.pivot(index="taxon", columns="domain", values="n").reindex(
        index=taxon_order, columns=domain_order).fillna(0)
    return rate, n


def draw_panel(ax, rate: pd.DataFrame, n: pd.DataFrame, vmax: float, title: str, cmap):
    reliable = n.values >= MIN_N
    masked_rate = np.ma.masked_where(~reliable, rate.values)

    im = ax.imshow(masked_rate, cmap=cmap, vmin=0, vmax=vmax, aspect="auto")

    n_rows, n_cols = rate.shape
    for i in range(n_rows):
        for j in range(n_cols):
            n_ij = n.values[i, j]
            rate_ij = rate.values[i, j]
            if n_ij == 0:
                ax.text(j, i, "–", ha="center", va="center", fontsize=7, color="0.6")
                continue
            if n_ij < MIN_N:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="0.88",
                                        hatch="///", edgecolor="0.6", linewidth=0.5))
                ax.text(j, i, f"n={int(n_ij)}", ha="center", va="center",
                        fontsize=6, color="0.35")
                continue
            # Text color flips to white on dark (high-rate) cells for legibility.
            text_color = "white" if rate_ij / vmax > 0.6 else "black"
            ax.text(j, i, f"{rate_ij:.1f}%\n(n={int(n_ij)})", ha="center", va="center",
                    fontsize=6.5, color=text_color)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(rate.columns, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(rate.index, fontsize=8)
    ax.set_title(title, fontsize=10)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return im


def build_figure(alpha_tables: dict, taxon_order: list[str], domain_orders: dict, out_path: Path):
    fig, axes = plt.subplots(2, 3, figsize=(19, 11.5))
    # Fixed margins, reserving a right-hand strip for colorbars -- deliberately
    # NOT using tight_layout() afterward, since it fights the colorbar axes
    # placed via fig.add_axes() below and was clipping the rightmost (alpha=1)
    # column's labels/values.
    fig.subplots_adjust(left=0.11, right=0.90, top=0.86, bottom=0.09, wspace=0.5, hspace=0.65)
    cmaps = {"te_capsid": plt.get_cmap("Blues"), "domesticated_host_gene": plt.get_cmap("Oranges")}

    for row_idx, (clan_group, row_label) in enumerate(GROUPS):
        domain_order = domain_orders[clan_group]
        vmax = max(
            np.nanmax(np.where(alpha_tables[tag][clan_group][1].values >= MIN_N,
                                alpha_tables[tag][clan_group][0].values, np.nan))
            for tag, _ in ALPHAS
        )
        vmax = max(vmax, 5.0)  # floor so a near-zero-rate row doesn't get an unreadably tight scale
        im = None
        for col_idx, (tag, alpha) in enumerate(ALPHAS):
            rate, n = alpha_tables[tag][clan_group]
            rate = rate.reindex(columns=domain_order)
            n = n.reindex(columns=domain_order).fillna(0)
            ax = axes[row_idx, col_idx]
            title = f"α = {alpha}, core_length = {CORE_LENGTH}"
            im = draw_panel(ax, rate, n, vmax, title, cmaps[clan_group])

        # Dedicated colorbar strip to the right of this row's 3 panels, sized
        # to that row's own vertical extent so it can't overlap the row below.
        pos_first = axes[row_idx, 0].get_position()
        pos_last = axes[row_idx, -1].get_position()
        cbar_ax = fig.add_axes([0.92, pos_last.y0, 0.014, pos_first.y1 - pos_last.y0])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label("PrLD-positive rate (%)\n(cluster-corrected)", fontsize=8)
        axes[row_idx, 0].set_ylabel(row_label, fontsize=9, fontweight="bold", labelpad=10)

    fig.suptitle(
        "PLAAC PrLD-positive rate by taxon and domain — Pfam Clan CL0523\n"
        "Unit of analysis: cluster-corrected (one representative per within-domain MMseqs2 cluster, 95% identity / 80% coverage)",
        fontsize=12, y=0.965,
    )
    fig.text(
        0.5, 0.015,
        "Hatched/grey cells: cluster-corrected n < %d (unreliable rate, shown but not colored). "
        "Bottom row (domesticated host genes) is a non-TE classification sanity check, not part of the TE PrLD claim." % MIN_N,
        ha="center", fontsize=8, color="0.3",
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    log.info("Wrote figure to %s", out_path)


def main() -> None:
    global MIN_N
    parser = argparse.ArgumentParser(
        description="Multi-panel taxonomic PrLD-positive-rate heatmap across "
                    "PLAAC alpha=0/0.5/1 for the clan_cl0523 pull."
    )
    parser.add_argument("--min-n", type=int, default=MIN_N,
                        help="Cluster-corrected n below which a cell is masked as unreliable")
    parser.add_argument("--taxon-level", choices=["fine", "broad"], default="fine",
                        help="'fine' = taxon column (Animal split into Fish/Mammal/Insect/...); "
                             "'broad' = broad_group collapsed to Plant/Animal/Fungi/Other Eukaryote/Bacteria-or-other, "
                             "matching the clan_fig1-3 convention")
    parser.add_argument("--out", default=None,
                        help="Defaults to clan_cl0523_taxon_prld_heatmap.png (fine) or "
                             "clan_cl0523_taxon_prld_heatmap_broad.png (broad)")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    MIN_N = args.min_n
    suffix = "" if args.taxon_level == "fine" else "_broad"
    out_path = Path(args.out) if args.out else RESULTS_DIR / f"clan_cl0523_taxon_prld_heatmap{suffix}.png"

    setup_logging(args.verbose)
    log.info("=== Clan taxon PrLD heatmap starting (min_n=%d, taxon_level=%s) ===", MIN_N, args.taxon_level)

    tables: dict[str, dict[str, pd.DataFrame]] = {}
    long_rows = []
    taxon_order_source = None
    domain_orders = {}

    for tag, alpha in ALPHAS:
        df = load_alpha(tag, args.taxon_level)
        tables[tag] = {}
        for clan_group, _ in GROUPS:
            tab = cluster_corrected_rate_table(df, clan_group, MIN_N)
            tables[tag][clan_group] = tab

            tagged = tab.copy()
            tagged.insert(0, "alpha", alpha)
            tagged.insert(1, "clan_group", clan_group)
            long_rows.append(tagged)

            if clan_group not in domain_orders:
                domain_orders[clan_group] = order_by_total_n(tab, "domain")
            if tag == "a0p5" and clan_group == "te_capsid":
                taxon_order_source = order_by_total_n(tab, "taxon")

    # Union taxon order across both groups, ranked by te_capsid's ordering (the
    # main analysis), with any taxon absent from te_capsid appended at the end.
    all_taxa = set()
    for chunk in long_rows:
        all_taxa |= set(chunk["taxon"].unique())
    taxon_order = [t for t in taxon_order_source if t in all_taxa]
    taxon_order += sorted(all_taxa - set(taxon_order))

    alpha_tables: dict[str, dict[str, tuple]] = {}
    for tag, alpha in ALPHAS:
        alpha_tables[tag] = {}
        for clan_group, _ in GROUPS:
            alpha_tables[tag][clan_group] = pivot_grid(
                tables[tag][clan_group], taxon_order, domain_orders[clan_group])

    out_csv = RESULTS_DIR / f"clan_cl0523_taxon_prld_rate_by_alpha{suffix}.csv"
    long_df = pd.concat(long_rows, ignore_index=True)
    long_df.to_csv(out_csv, index=False)
    log.info("Wrote backing table (%d rows) to %s", len(long_df), out_csv)

    build_figure(alpha_tables, taxon_order, domain_orders, out_path)


if __name__ == "__main__":
    main()
