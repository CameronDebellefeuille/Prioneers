"""
plaac_utils.py
Shared helpers for running PLAAC (https://github.com/whitehead/plaac) and its
R plotting script, and for parsing PLAAC's TSV output. Used by both
te_bulk_pull.py and te_bulk_analysis.py.
"""

import logging
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

log = logging.getLogger(__name__)


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
        level=level,
        handlers=[logging.StreamHandler(sys.stdout)],
    )


def run_plaac(fasta_path: str, plaac_jar: str, alpha: float,
              core_length: int, output_dir: str, dry_run: bool = False) -> str | None:
    """
    Run PLAAC on a protein FASTA. Returns path to output TSV, or None on failure.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    alpha_tag = str(alpha).replace(".", "p")
    out_tsv = str(Path(output_dir) / f"summary_a{alpha_tag}.tsv")

    jar = Path(plaac_jar)
    if not jar.exists():
        log.error("plaac.jar not found at %s — see setup_instructions.md", plaac_jar)
        return None

    cmd = [
        "java", "-jar", str(jar),
        "-i", fasta_path,
        "-a", str(alpha),
        "-c", str(core_length),
    ]
    log.info("PLAAC command: %s", " ".join(cmd))

    if dry_run:
        log.info("[dry-run] Would write output to %s", out_tsv)
        return out_tsv

    try:
        with open(out_tsv, "w") as f_out, open(out_tsv.replace(".tsv", ".err"), "w") as f_err:
            result = subprocess.run(cmd, stdout=f_out, stderr=f_err, timeout=300)
        if result.returncode != 0:
            log.error("PLAAC exited with code %d — check %s", result.returncode,
                      out_tsv.replace(".tsv", ".err"))
            return None
    except FileNotFoundError:
        log.error("'java' not found. Install Java JRE — see setup_instructions.md")
        return None
    except subprocess.TimeoutExpired:
        log.error("PLAAC timed out after 5 minutes")
        return None

    log.info("PLAAC output written to %s", out_tsv)
    return out_tsv


def run_plaac_plots(fasta_path: str, plaac_jar: str, plaac_plot_r: str,
                    prd_ids: list[str], plots_dir: str, alpha: float = 0.5,
                    core_length: int = 60, dry_run: bool = False,
                    rscript: str = "Rscript", plot_all: bool = False) -> None:
    """
    Generate per-residue PLAAC output (-p) then call Rscript to produce PNG plots.
    By default only plots proteins with a called PrLD; pass plot_all=True to
    plot every scored protein regardless of PrLD call.
    """
    if not plot_all and not prd_ids:
        log.info("No PrLD-positive proteins; skipping plots.")
        return

    plots_path = Path(plots_dir).resolve()
    plots_path.mkdir(parents=True, exist_ok=True)
    r_script = Path(plaac_plot_r).resolve()
    jar = Path(plaac_jar).resolve()

    if not r_script.exists():
        log.warning("plaac_plot.r not found at %s; skipping plots.", plaac_plot_r)
        return
    if not jar.exists():
        log.warning("plaac.jar not found; skipping plots.")
        return

    # Write a filtered FASTA: PrLD-positive proteins, or all scored proteins if plot_all
    filtered_fasta = plots_path / ("all_scored.fasta" if plot_all else "prd_positive.fasta")
    all_records = list(SeqIO.parse(fasta_path, "fasta"))
    records = all_records if plot_all else [r for r in all_records if r.id in set(prd_ids)]
    if not records:
        log.warning("Could not match PrLD IDs to FASTA records; skipping plots.")
        return

    with open(filtered_fasta, "w") as f:
        SeqIO.write(records, f, "fasta")

    if dry_run:
        log.info("[dry-run] Would run per-residue PLAAC and Rscript for %d proteins", len(records))
        return

    # Per-residue mode: PLAAC with -p all writes a per-residue table to stdout
    perres_tsv = plots_path / "perresidues.tsv"
    cmd_perres = ["java", "-jar", str(jar), "-i", str(filtered_fasta),
                  "-a", str(alpha), "-c", str(core_length), "-p", "all"]
    try:
        with open(perres_tsv, "w") as out:
            subprocess.run(cmd_perres, check=True, timeout=120, stdout=out, stderr=subprocess.PIPE, text=True)
    except Exception as e:
        log.warning("Per-residue PLAAC run failed: %s", e)
        return

    # plaac_plot.r does source("plaac_plot_util.r") with a relative path, so it
    # must run with its own directory as the working directory. Figname must
    # end in .png/.pdf; for png, one indexed file is written per protein.
    figname = plots_path / "plaac_plot.png"
    cmd_r = [rscript, str(r_script), str(perres_tsv), str(figname)]
    try:
        result = subprocess.run(cmd_r, capture_output=True, text=True, timeout=180, cwd=str(r_script.parent))
        if result.returncode != 0:
            log.warning("Rscript plot failed:\n%s", (result.stderr or result.stdout)[:1000])
        else:
            log.info("Plots written to %s/", plots_path)
    except FileNotFoundError:
        log.warning("'Rscript' not found — install R to generate plots.")


def _parse_plaac_tsv(tsv_path: str) -> pd.DataFrame | None:
    """Read PLAAC summary TSV; return DataFrame or None if unreadable."""
    try:
        df = pd.read_csv(tsv_path, sep="\t", comment="#")
        log.info("PLAAC TSV: %d rows, columns: %s", len(df), list(df.columns))
        return df
    except Exception as e:
        log.error("Failed to parse PLAAC TSV %s: %s", tsv_path, e)
        return None
