"""
plant_prion_plaac_benchmark.py
Standalone benchmark, not part of the TE bulk-pull or clan pipelines: scores
three literature-validated prion/PrLD proteins (Arabidopsis ELF3 and LD, plus
yeast Sup35 as a positive-control anchor) with PLAAC at alpha=0/0.5/1, and
assembles a publication figure combining PLAAC's own per-residue score tracks
with a summary metrics panel.

Answers "does PLAAC recover known/suspected plant PrDs, and how alpha-
sensitive is the call?" -- not a TE screen.

EVD (ATCOPIA93) was scoped for this benchmark too but deliberately excluded:
it maps to TAIR/Araport locus AT5G17125, which has no UniProt entry and no
annotated CDS/translation in GenBank (flagged /pseudo, gene+mRNA feature only,
no protein sequence) -- getting a correct Gag ORF would mean hand-calling the
ORF and Gag/Pol boundary from raw genomic DNA, which isn't defensible for a
publication figure without a citable source.

Alpha direction (verified against plaac.jar's own -h text and the upstream
PLAAC CLI docs, github.com/whitehead/plaac/cli/README.md -- this repo's
README previously had it backwards): alpha=0 uses a background AA-frequency
model built from the INPUT FASTA itself; alpha=1 uses a fixed, external
S. cerevisiae proteome-wide background. With only 3 proteins in this
benchmark's FASTA, the alpha=0 "background" is a small, self-referential
estimate built from the very proteins being scored (including their own
PrLD-like stretches) -- weaker discrimination, more conservative calls, not
a "compositionally neutral" baseline. alpha=1's fixed external background is
the more stable estimate at this sample size. See the figure's own footer
note.
"""

import argparse
import json
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless: this script only writes a PNG, never shows a window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests

from te_bulk_pull import CONFIG, BASE_DIR
from plaac_utils import run_plaac, run_plaac_plots, _parse_plaac_tsv, setup_logging

log = logging.getLogger(__name__)

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"

# Literature-validated set, confirmed against UniProt 2026-07-22. Order here
# is fixed and deliberate: it drives the FASTA record order, which in turn
# drives the R plotter's indexed PNG order (plaac_plot_00001.png = record 1,
# etc.) -- so this list is the single source of truth for row identity
# throughout the script.
PROTEINS = [
    {"accession": "P05453", "label": "Sup35", "organism": "S. cerevisiae", "length": 685,
     "note": "Canonical yeast [PSI+] prion — PLAAC HMM training positive"},
    {"accession": "Q38796", "label": "LD", "organism": "A. thaliana", "length": 953,
     "note": "Yeast-assay-confirmed heritable prion behavior (Chakrabortee et al. 2016, PNAS)"},
    {"accession": "O82804", "label": "ELF3", "organism": "A. thaliana", "length": 695,
     "note": "Thermosensor PrD (Jung et al. 2020, Nature)"},
]

ALPHAS = [0.0, 0.5, 1.0]
CORE_LENGTH = 60  # matches every other PLAAC run in this repo

DATA_DIR = BASE_DIR / "data"
CACHE_DIR = DATA_DIR / "prion_cache"
FASTA_DIR = DATA_DIR / "prion_fasta"
PLAAC_OUTPUT_DIR = DATA_DIR / "prion_plaac_output"
RESULTS_DIR = DATA_DIR / "results" / "plant_prion_benchmark"

# dataviz skill's validated categorical palette, slots 1-3 (blue/aqua/yellow),
# fixed assignment to alpha values -- never cycled or reassigned.
ALPHA_COLORS = {0.0: "#2a78d6", 0.5: "#1baf7a", 1.0: "#eda100"}
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"
SURFACE = "#fcfcfb"


def fetch_sequences(cache_path: Path) -> dict[str, str]:
    """accession -> sequence, fetched from UniProt's FASTA endpoint and cached
    to disk. Hard-fails if a fetched length doesn't match the length recorded
    in PROTEINS at scoping time -- a cheap guard against silently scoring the
    wrong protein if an accession ever gets edited without checking."""
    cache: dict[str, str] = {}
    if cache_path.exists():
        with open(cache_path) as f:
            cache = json.load(f)

    for p in PROTEINS:
        acc = p["accession"]
        if acc in cache:
            continue
        r = requests.get(f"{UNIPROT_API}/{acc}.fasta", timeout=30)
        r.raise_for_status()
        lines = r.text.strip().split("\n")
        seq = "".join(lines[1:])
        cache[acc] = seq
        log.info("Fetched %s (%s): %d aa", acc, p["label"], len(seq))

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(cache, f, indent=2)

    for p in PROTEINS:
        seq = cache[p["accession"]]
        if len(seq) != p["length"]:
            raise SystemExit(
                f"{p['accession']} ({p['label']}): expected {p['length']} aa, got "
                f"{len(seq)} aa -- sequence identity mismatch, aborting rather than "
                f"scoring the wrong protein")
    return cache


def write_fasta(sequences: dict[str, str], path: Path) -> list[str]:
    """One record per protein, in PROTEINS order. Returns the stable_id list in
    that same order (see PROTEINS docstring note on why order matters)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    stable_ids = []
    with open(path, "w") as f:
        for p in PROTEINS:
            stable_id = f"{p['label']}_{p['accession']}"
            stable_ids.append(stable_id)
            f.write(f">{stable_id}\n{sequences[p['accession']]}\n")
    log.info("Wrote %d sequences to %s", len(PROTEINS), path)
    return stable_ids


def score_alpha(fasta_path: Path, alpha: float) -> pd.DataFrame:
    """Run PLAAC's summary scoring at one alpha; return a long-format frame
    with a prd_called column and an alpha column, matching the convention
    used by score_with_plaac() in te_bulk_pull.py."""
    alpha_tag = str(alpha).replace(".", "p")
    output_dir = PLAAC_OUTPUT_DIR / f"a{alpha_tag}"
    tsv_path = run_plaac(str(fasta_path), CONFIG["plaac_jar"], alpha, CORE_LENGTH, str(output_dir))
    if tsv_path is None:
        raise SystemExit(f"PLAAC run failed at alpha={alpha} -- see the .err file next to the output TSV")
    plaac_df = _parse_plaac_tsv(tsv_path)
    if plaac_df is None:
        raise SystemExit(f"Could not parse PLAAC output at alpha={alpha}")
    plaac_df["stable_id"] = plaac_df["SEQid"].str.split().str[0]
    plaac_df["prd_called"] = (plaac_df["PRDlen"].fillna(0).astype(float) > 0).astype(int)
    plaac_df["alpha"] = alpha
    return plaac_df


def build_summary(fasta_path: Path) -> pd.DataFrame:
    """Score all three alphas and concatenate into one long-format table."""
    label_map = {f"{p['label']}_{p['accession']}": p["label"] for p in PROTEINS}
    frames = [score_alpha(fasta_path, alpha) for alpha in ALPHAS]
    long_df = pd.concat(frames, ignore_index=True)
    long_df["label"] = long_df["stable_id"].map(label_map)
    return long_df


def make_plots(fasta_path: Path) -> dict[float, Path]:
    """Per-residue PLAAC + R plotting at each alpha, reusing run_plaac_plots()
    exactly as te_bulk_analysis.py does -- plot_all=True since we want every
    protein plotted at every alpha, not just PrLD-positive calls (that's the
    whole point of the alpha=0/1 stability comparison)."""
    plot_dirs = {}
    for alpha in ALPHAS:
        alpha_tag = str(alpha).replace(".", "p")
        plots_dir = RESULTS_DIR / f"plots_a{alpha_tag}"
        run_plaac_plots(str(fasta_path), CONFIG["plaac_jar"], CONFIG["plaac_plot_r"],
                        prd_ids=[], plots_dir=str(plots_dir), alpha=alpha,
                        core_length=CORE_LENGTH, rscript=CONFIG["rscript"], plot_all=True)
        plot_dirs[alpha] = plots_dir
    return plot_dirs


def draw_summary_panel(ax, summary_df: pd.DataFrame) -> None:
    """Grouped bar chart: PLAAC COREscore by protein x alpha, with the PrLD
    call (PRDlen > 0) marked above each bar rather than double-encoded into
    color -- alpha already owns the color channel here (dataviz skill: color
    follows one job at a time)."""
    labels = [p["label"] for p in PROTEINS]
    x = np.arange(len(labels))
    bar_width = 0.24
    max_score = summary_df["COREscore"].fillna(0).max()
    headroom = max(max_score * 0.08, 5)

    for i, alpha in enumerate(ALPHAS):
        sub = summary_df[summary_df["alpha"] == alpha].set_index("label").reindex(labels)
        offset = (i - 1) * (bar_width + 0.02)
        scores = sub["COREscore"].fillna(0).values
        called = sub["prd_called"].fillna(0).values
        ax.bar(x + offset, scores, width=bar_width, color=ALPHA_COLORS[alpha],
               label=f"α = {alpha}", zorder=3)
        for xi, score, is_called in zip(x + offset, scores, called):
            marker = "✓" if is_called else "–"
            marker_color = INK_PRIMARY if is_called else INK_MUTED
            ax.text(xi, score + headroom * 0.15, marker, ha="center", va="bottom",
                    fontsize=9, color=marker_color)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10.5, color=INK_PRIMARY)
    ax.set_ylabel("PLAAC COREscore", fontsize=10, color=INK_SECONDARY)
    ax.set_ylim(top=max_score + headroom)
    ax.tick_params(colors=INK_MUTED, length=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(BASELINE)
    ax.grid(axis="y", color=GRIDLINE, linewidth=1, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY,
              loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.set_title("PrLD call (✓ = PLAAC PRDlen > 0) and COREscore by protein × α",
                 fontsize=10, color=INK_SECONDARY, loc="left")


def assemble_figure(plot_dirs: dict[float, Path], summary_df: pd.DataFrame, out_path: Path) -> None:
    """Combines the per-alpha per-residue PNGs (already self-labeled with each
    protein's SEQid in the title by plaac_plot_util.r) into a protein x alpha
    grid, plus the summary panel below."""
    n_proteins = len(PROTEINS)
    n_alphas = len(ALPHAS)

    fig = plt.figure(figsize=(6.2 * n_alphas, 3.3 * n_proteins + 3.4), facecolor=SURFACE)
    gs = fig.add_gridspec(n_proteins + 1, n_alphas,
                          height_ratios=[1] * n_proteins + [1.3], hspace=0.25, wspace=0.03)

    for row, p in enumerate(PROTEINS):
        for col, alpha in enumerate(ALPHAS):
            ax = fig.add_subplot(gs[row, col])
            img_path = plot_dirs[alpha] / f"plaac_plot_{row + 1:05d}.png"
            ax.imshow(plt.imread(img_path))
            ax.axis("off")
            if row == 0:
                ax.set_title(f"α = {alpha}", fontsize=13, color=INK_PRIMARY, pad=10)
            if col == 0:
                ax.text(-0.03, 0.5, f"{p['label']}\n{p['note']}", transform=ax.transAxes,
                        ha="right", va="center", fontsize=8.5, color=INK_SECONDARY,
                        wrap=True, linespacing=1.4)

    ax_summary = fig.add_subplot(gs[n_proteins, :])
    draw_summary_panel(ax_summary, summary_df)

    fig.suptitle(
        "PLAAC PrLD scoring across α = 0 / 0.5 / 1 (core_length = 60)\n"
        "Plant PrLD candidates vs. a canonical yeast prion",
        fontsize=13.5, color=INK_PRIMARY, y=0.995)

    fig.text(
        0.5, 0.006,
        "α=0 background = this input FASTA's own composition (here, only 3 proteins — a small, "
        "self-referential estimate); α=1 background = a fixed external S. cerevisiae proteome-wide "
        "composition, independent of the query set. α=0 is therefore the less stable estimate at "
        "this sample size, not a compositionally neutral baseline.",
        ha="center", fontsize=7.8, color=INK_MUTED, wrap=True,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    log.info("Wrote figure to %s", out_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Score ELF3, LD, and Sup35 with PLAAC at alpha=0/0.5/1 and "
                    "assemble a combined per-residue + summary publication figure."
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(args.verbose)
    log.info("=== Plant prion PLAAC benchmark starting ===")

    sequences = fetch_sequences(CACHE_DIR / "sequences.json")
    fasta_path = FASTA_DIR / "plant_prion_benchmark.fasta"
    write_fasta(sequences, fasta_path)

    if args.dry_run:
        for alpha in ALPHAS:
            alpha_tag = str(alpha).replace(".", "p")
            run_plaac(str(fasta_path), CONFIG["plaac_jar"], alpha, CORE_LENGTH,
                      str(PLAAC_OUTPUT_DIR / f"a{alpha_tag}"), dry_run=True)
        log.info("[dry-run] Stopping before scoring/plotting/figure assembly")
        return

    summary_df = build_summary(fasta_path)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    summary_csv = RESULTS_DIR / "plant_prion_plaac_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    log.info("Wrote summary table to %s", summary_csv)

    plot_dirs = make_plots(fasta_path)

    out_fig = RESULTS_DIR / "plant_prion_plaac_figure.png"
    assemble_figure(plot_dirs, summary_df, out_fig)

    log.info("=== Done ===")


if __name__ == "__main__":
    main()
