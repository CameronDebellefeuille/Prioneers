"""
te_plaac_pipeline.py
Screens TE-encoded proteins for prion-like domains using PLAAC.
"""

import argparse
import json
import logging
import os
import re
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import requests

# ---------------------------------------------------------------------------
# CONFIG — edit paths here or override via CLI flags
# ---------------------------------------------------------------------------
CONFIG = {
    "plan_doc": "TE_PLAAC_pipeline_plan.md",
    "db": "results/prioneers.db",
    "fasta_dir": "data/fasta",
    "cache_dir": "data/cache",
    "plaac_output_dir": "data/plaac_output",
    "plaac_jar": "tools/plaac/plaac.jar",
    "plaac_plot_r": "tools/plaac/R/plaac_plot.r",
    "plots_dir": "results/plots",
    "conda_env": "prioneers",
    "entrez_email": "21cd51@queensu.ca",
    "entrez_sleep": 0.4,        # seconds between NCBI calls
    "orf_min_nt": 300,          # minimum ORF length passed to getorf (nucleotides)
    "orf_min_aa": 100,          # secondary AA filter after translation
    "alpha": 0.5,
    "core_length": 60,
    # Rscript path — defaults to PATH; set RSCRIPT_PATH if R isn't on PATH
    "rscript": os.environ.get("RSCRIPT_PATH", "Rscript"),
}

# Keywords used to label/prefer ORFs that correspond to known TE domains
TE_DOMAIN_KEYWORDS = [
    "gag", "pol", "rt", "reverse transcriptase", "integrase", "rve",
    "transposase", "rep", "helicase", "orf1", "orf2", "orf1p", "orf2p",
    "endonuclease", "ribonuclease h", "rnase h",
]

# ---------------------------------------------------------------------------
# Manual retrieval overrides
# ---------------------------------------------------------------------------
# The pipeline plan table uses descriptive group labels ("TART / HeT-A",
# "CR1 / CR1-C", ...) that are NOT literal Dfam family names or NCBI search
# terms, so the generic name-search/Entrez-AND-query retrieval in
# fetch_sequences() returns nothing for them (Dfam's name filter requires an
# exact canonical name and rejects "/"; NCBI's AND-of-free-text query only
# matches if that exact compound string appears in a record).
#
# These entries were confirmed by hand against the Dfam API and NCBI
# (accession, organism, and CDS/product checked) and are fetched directly by
# accession, bypassing the fragile generic search entirely.
#
# Families with no entry here (Crwydryn, Gmr9) have no real deposited
# sequence under that name — Crwydryn isn't in Dfam/NCBI at all, and "Gmr9"
# collides with an unrelated bacterial gene abbreviation ("Gmr" cyclic-di-GMP
# phosphodiesterase) — both still require pulling sequences from their
# source papers' supplementary data.
FAMILY_RETRIEVAL_OVERRIDES: dict[str, list[dict]] = {
    "TART / HeT-A": [
        {"kind": "dfam", "accession": "DF000001711"},   # TART-A consensus (nt)
        {"kind": "dfam", "accession": "DF000001642"},   # HeT-A consensus (nt)
        {"kind": "ncbi_protein", "id": "542574", "protein_name": "Gag"},  # D. melanogaster HeT-A gag-like polyprotein (fragment)
    ],
    "CR1 / CR1-C": [
        {"kind": "dfam", "accession": "DF000000110"},   # CR1_Mam consensus (nt)
        {"kind": "ncbi_protein", "id": "347227", "protein_name": "ORF2-RT"},  # Gallus gallus CR1 reverse transcriptase (partial)
    ],
    "Atenspm2 (CACTA/EnSpm)": [
        {"kind": "ncbi_protein", "id": "BAE98974.1", "protein_name": "Transposase"},  # A. thaliana En/Spm-like transposon protein
    ],
    "Sore1 (potato LTR retro)": [
        # NB: the only deposited "SORE1" sequence is from Glycine max (soybean),
        # not potato — the plan's host-organism label appears to be wrong.
        {"kind": "ncbi_protein", "id": "BAB40834.1", "protein_name": "RT"},  # Glycine max gypsy-type retrotransposon SORE1, reverse transcriptase (partial)
    ],
    "Ty3/gypsy ogre elements": [
        {"kind": "ncbi_protein", "id": "ACX37123.1", "protein_name": "Gag-Pro"},  # Lathyrus oleraceus Ogre retrotransposon gag-pro
        {"kind": "ncbi_protein", "id": "AAQ82035.1", "protein_name": "Gag-Pol"},  # Pisum sativum Ogre retrotransposon gag/pol polyprotein (partial)
    ],
    "CRM1 / CRM4 (maize centromeric)": [
        {"kind": "ncbi_protein", "id": "ACY00046.1", "protein_name": "Polyprotein"},  # Zea mays CRM1 subgroup A/R2/R3/R4 polyprotein (partial)
    ],
}

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
        level=level,
        handlers=[logging.StreamHandler(sys.stdout)],
    )


log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Step 1 — Parse the source table from the plan document
# ---------------------------------------------------------------------------

# Disposition text → scope code
_SCOPE_MAP = {
    "in scope": "in_scope",
    "flag": "broad_group",
    "exclude": "non_coding",
    "broad group": "broad_group",
    "non-coding": "non_coding",
    "non_coding": "non_coding",
}

def _map_scope(disposition: str) -> tuple[str, str]:
    """Return (scope_code, scope_note) from a raw disposition cell."""
    lower = disposition.lower()
    if "exclude" in lower or "non-coding" in lower or "none" in lower:
        return "non_coding", disposition
    if "flag" in lower or "broad group" in lower or "uncertain" in lower or "non-autonomous" in lower:
        return "broad_group", disposition
    if "in scope" in lower:
        return "in_scope", disposition
    return "broad_group", disposition


def _infer_te_class(family_name: str) -> str:
    """Heuristically assign a TE class from the family name."""
    fn = family_name.lower()
    if any(x in fn for x in ["line-1", "l1", "cr1", "tart", "het-a"]):
        return "LINE"
    if any(x in fn for x in ["ty3", "gypsy", "ogre", "crm", "ltr", "crwydryn", "sore1", "gmr9", "malr"]):
        return "LTR-Ty3gypsy"
    if any(x in fn for x in ["cacta", "enspm", "atenspm"]):
        return "CACTA"
    if "helitron" in fn:
        return "Helitron"
    if any(x in fn for x in ["sine", "alu", "b1"]):
        return "SINE"
    if "mite" in fn:
        return "MITE"
    return "unknown"


def parse_table(plan_doc: str) -> list[dict]:
    """
    Parse the markdown element-classification table in the plan doc.
    Returns a list of dicts (one per unique family_name).
    """
    path = Path(plan_doc)
    if not path.exists():
        raise FileNotFoundError(f"Plan document not found: {plan_doc}")

    text = path.read_text(encoding="utf-8")

    # Find the markdown table between the header row and the next blank line after rows
    table_pattern = re.compile(
        r"\| Family.*?Disposition.*?\n"   # header
        r"\|[-| ]+\|\n"                   # separator
        r"((?:\|.*\|\n?)+)",              # rows
        re.IGNORECASE,
    )
    match = table_pattern.search(text)
    if not match:
        raise ValueError("Could not locate the element-classification table in the plan document.")

    rows_text = match.group(1)
    elements = []
    seen = set()

    for line in rows_text.strip().splitlines():
        cells = [c.strip() for c in line.split("|") if c.strip()]
        if len(cells) < 3:
            continue
        family_raw, proteins_raw, disposition_raw = cells[0], cells[1], cells[2]

        # Strip markdown bold markers
        family_name = re.sub(r"\*+", "", family_raw).strip().strip('"')
        encoded_proteins = re.sub(r"\*+", "", proteins_raw).strip()
        disposition = re.sub(r"\*+", "", disposition_raw).strip()

        if not family_name or family_name.startswith("-"):
            continue

        scope, scope_note = _map_scope(disposition)
        te_class = _infer_te_class(family_name)

        # De-duplicate by family_name
        key = family_name.lower()
        if key in seen:
            continue
        seen.add(key)

        elements.append({
            "family_name": family_name,
            "encoded_proteins": encoded_proteins,   # stored on the dict, used for ORF labelling
            "te_class": te_class,
            "host_organism": None,
            "derived_seq": None,
            "reference_num": None,
            "scope": scope,
            "scope_note": scope_note,
        })

    log.info("Parsed %d unique families from table (%d in_scope, %d broad_group, %d non_coding)",
             len(elements),
             sum(1 for e in elements if e["scope"] == "in_scope"),
             sum(1 for e in elements if e["scope"] == "broad_group"),
             sum(1 for e in elements if e["scope"] == "non_coding"))
    return elements


# ---------------------------------------------------------------------------
# Step 2 — Retrieve nucleotide sequences (Dfam + NCBI)
# ---------------------------------------------------------------------------

DFAM_API = "https://dfam.org/api"


def _dfam_search_family(name: str) -> str | None:
    """Return Dfam accession for a family name, or None."""
    try:
        r = requests.get(f"{DFAM_API}/families",
                         params={"name": name, "limit": 5},
                         timeout=20)
        r.raise_for_status()
        data = r.json()
        results = data.get("results", [])
        if results:
            acc = results[0].get("accession")
            log.debug("Dfam: '%s' → %s", name, acc)
            return acc
    except Exception as e:
        log.debug("Dfam search failed for '%s': %s", name, e)
    return None


def _dfam_get_sequence(accession: str) -> str | None:
    """Fetch consensus nucleotide sequence for a Dfam accession."""
    try:
        r = requests.get(f"{DFAM_API}/family/{accession}/sequence", timeout=30)
        r.raise_for_status()
        data = r.json()
        seq = data.get("sequence") or data.get("cons_seq") or data.get("consensus")
        if seq:
            return seq.upper().replace(" ", "").replace("\n", "")
    except Exception as e:
        log.debug("Dfam sequence fetch failed for %s: %s", accession, e)
    return None


def _ncbi_fetch_protein_by_id(ncbi_id: str, protein_name: str) -> dict | None:
    """Fetch one protein record directly by accession/GI (used for confirmed overrides)."""
    Entrez.email = CONFIG["entrez_email"]
    try:
        handle = Entrez.efetch(db="protein", id=ncbi_id, rettype="fasta", retmode="text")
        fasta_text = handle.read()
        handle.close()
        time.sleep(CONFIG["entrez_sleep"])
        rec = next(SeqIO.parse(__import__("io").StringIO(fasta_text), "fasta"))
        return {
            "accession": rec.id,
            "protein_name": protein_name,
            "sequence": str(rec.seq),
            "source": "NCBI",
            "seq_type": "aa",
        }
    except Exception as e:
        log.warning("NCBI direct protein fetch failed for id=%s: %s", ncbi_id, e)
        return None


def _ncbi_fetch_protein(family_name: str, encoded_proteins: str) -> list[dict]:
    """
    Search NCBI protein database for known proteins of this TE family.
    Returns list of {accession, protein_name, sequence (AA), source}.
    """
    Entrez.email = CONFIG["entrez_email"]
    results = []

    # Build search terms combining family name + expected protein types
    protein_hints = [p.strip() for p in re.split(r"[,/]", encoded_proteins) if p.strip() and p.strip() != "—"]
    if not protein_hints:
        protein_hints = ["transposase"]

    for hint in protein_hints:
        term = f'"{family_name}"[All Fields] AND "{hint}"[All Fields]'
        try:
            handle = Entrez.esearch(db="protein", term=term, retmax=3)
            record = Entrez.read(handle)
            handle.close()
            ids = record.get("IdList", [])
            time.sleep(CONFIG["entrez_sleep"])

            if not ids:
                log.debug("NCBI protein: no hits for '%s' + '%s'", family_name, hint)
                continue

            fetch_handle = Entrez.efetch(db="protein", id=ids[0], rettype="fasta", retmode="text")
            fasta_text = fetch_handle.read()
            fetch_handle.close()
            time.sleep(CONFIG["entrez_sleep"])

            for rec in SeqIO.parse(__import__("io").StringIO(fasta_text), "fasta"):
                results.append({
                    "accession": rec.id,
                    "protein_name": hint,
                    "sequence": str(rec.seq),
                    "source": "NCBI",
                    "seq_type": "aa",
                })
                log.debug("NCBI: fetched %s for '%s'/%s", rec.id, family_name, hint)
                break  # one protein per hint is enough
        except Exception as e:
            log.warning("NCBI fetch error for '%s'/'%s': %s", family_name, hint, e)
            time.sleep(CONFIG["entrez_sleep"])

    return results


def fetch_sequences(elements: list[dict], cache_dir: str, skip_fetch: bool) -> list[dict]:
    """
    For each in_scope element, attempt Dfam then NCBI retrieval.
    Attaches results to each element dict under key 'retrieved'.
    Returns the same list (mutated).
    Logs misses to data/retrieval_misses.txt.
    """
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    misses = []

    for elem in elements:
        if elem["scope"] != "in_scope":
            elem["retrieved"] = []
            continue

        fname = elem["family_name"]
        cache_key = re.sub(r"[^a-zA-Z0-9_-]", "_", fname)
        cache_file = cache / f"{cache_key}.json"

        if cache_file.exists():
            with open(cache_file) as f:
                elem["retrieved"] = json.load(f)
            log.info("Cache hit: %s (%d records)", fname, len(elem["retrieved"]))
            continue

        if skip_fetch:
            log.warning("MISS: '%s' not cached and --skip-fetch set; no network call made", fname)
            elem["retrieved"] = []
            misses.append(fname)
            continue

        log.info("Fetching sequences for: %s", fname)
        retrieved = []

        overrides = FAMILY_RETRIEVAL_OVERRIDES.get(fname)
        if overrides:
            # Confirmed-by-hand accessions; bypass the fragile generic search entirely.
            for ov in overrides:
                if ov["kind"] == "dfam":
                    nt_seq = _dfam_get_sequence(ov["accession"])
                    if nt_seq:
                        retrieved.append({
                            "accession": ov["accession"],
                            "protein_name": None,
                            "sequence": nt_seq,
                            "source": "Dfam",
                            "seq_type": "nt",
                        })
                        log.info("Dfam: got nt consensus %s for %s (%d bp)", ov["accession"], fname, len(nt_seq))
                elif ov["kind"] == "ncbi_protein":
                    rec = _ncbi_fetch_protein_by_id(ov["id"], ov["protein_name"])
                    if rec:
                        retrieved.append(rec)
                        log.info("NCBI: got %s (%s) for %s", rec["accession"], ov["protein_name"], fname)
        else:
            # --- Try NCBI protein first (annotated AA is preferred) ---
            ncbi_proteins = _ncbi_fetch_protein(fname, elem["encoded_proteins"])
            retrieved.extend(ncbi_proteins)

            # --- Try Dfam for nucleotide consensus ---
            dfam_acc = _dfam_search_family(fname)
            if dfam_acc:
                nt_seq = _dfam_get_sequence(dfam_acc)
                if nt_seq:
                    retrieved.append({
                        "accession": dfam_acc,
                        "protein_name": None,
                        "sequence": nt_seq,
                        "source": "Dfam",
                        "seq_type": "nt",
                    })
                    log.info("Dfam: got nt consensus for %s (%d bp)", fname, len(nt_seq))

        if not retrieved:
            log.warning("MISS: no sequence retrieved for '%s'", fname)
            misses.append(fname)

        elem["retrieved"] = retrieved

        # Cache the result (even empty, to avoid re-hitting network)
        with open(cache_file, "w") as f:
            json.dump(retrieved, f, indent=2)

    # Write miss log
    miss_path = Path(cache_dir).parent / "retrieval_misses.txt"
    with open(miss_path, "w") as f:
        f.write("# Families with no automatic sequence retrieval\n")
        f.write("# Add sequences manually to data/cache/<name>.json and re-run with --skip-fetch\n\n")
        for m in misses:
            f.write(m + "\n")

    if misses:
        log.warning("Retrieval misses (%d): %s — see data/retrieval_misses.txt", len(misses), ", ".join(misses))

    return elements


# ---------------------------------------------------------------------------
# Step 3 — Extract proteins (annotated AA or getorf on nt)
# ---------------------------------------------------------------------------

def _rank_orfs(orfs: list[dict]) -> list[dict]:
    """Prefer domain-keyword ORFs; within each tier, prefer longest."""
    keyword_hits = [o for o in orfs if o["protein_name"] != "ORF"]
    if keyword_hits:
        return keyword_hits
    if orfs:
        return [max(orfs, key=lambda o: len(o["sequence"]))]
    return []


def _biopython_find_orfs(nt_seq: str, family_name: str, min_nt: int, min_aa: int) -> list[dict]:
    """Pure-Biopython ORF finder — used when EMBOSS getorf is unavailable (e.g. Windows)."""
    from Bio.Seq import Seq as BioSeq

    seq = BioSeq(nt_seq.upper())
    orfs = []
    for strand, nuc in ((+1, seq), (-1, seq.reverse_complement())):
        for frame in range(3):
            trans = nuc[frame:].translate()
            aa_str = str(trans)
            start = 0
            while True:
                stop = aa_str.find("*", start)
                fragment = aa_str[start:stop] if stop != -1 else aa_str[start:]
                if len(fragment) * 3 >= min_nt and len(fragment) >= min_aa:
                    desc_lower = family_name.lower()
                    protein_name = next((k for k in TE_DOMAIN_KEYWORDS if k in desc_lower), "ORF")
                    orfs.append({
                        "accession": f"{family_name}_s{strand}_f{frame}",
                        "protein_name": protein_name,
                        "sequence": fragment,
                        "source": "biopython_orf",
                        "seq_type": "aa",
                    })
                if stop == -1:
                    break
                start = stop + 1
    log.debug("biopython_orf: %d raw ORFs for %s", len(orfs), family_name)
    return _rank_orfs(orfs)


def _run_getorf(nt_seq: str, family_name: str, conda_env: str,
                cache_dir: str) -> list[dict]:
    """
    Call EMBOSS getorf via conda run. Falls back to Biopython ORF finder if
    conda/getorf is unavailable (bioconda does not support Windows natively).
    Returns list of {protein_name, sequence (AA), source}.
    """
    safe_name = re.sub(r"[^a-zA-Z0-9_-]", "_", family_name)
    tmp_in = Path(cache_dir) / f"{safe_name}_nt.fna"
    tmp_out = Path(cache_dir) / f"{safe_name}_orfs.faa"

    tmp_in.write_text(f">{safe_name}\n{nt_seq}\n", encoding="utf-8")

    cmd = [
        "conda", "run", "-n", conda_env,
        "getorf",
        "-sequence", str(tmp_in),
        "-outseq", str(tmp_out),
        "-minsize", str(CONFIG["orf_min_nt"]),
        "-find", "1",
        "-auto",
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            log.warning("getorf failed for %s: %s", family_name, result.stderr[:300])
            log.info("Falling back to Biopython ORF finder for %s", family_name)
            return _biopython_find_orfs(nt_seq, family_name, CONFIG["orf_min_nt"], CONFIG["orf_min_aa"])
    except FileNotFoundError:
        log.info("conda/getorf not found — using Biopython ORF finder for %s", family_name)
        return _biopython_find_orfs(nt_seq, family_name, CONFIG["orf_min_nt"], CONFIG["orf_min_aa"])
    except subprocess.TimeoutExpired:
        log.warning("getorf timed out for %s", family_name)
        return []

    if not tmp_out.exists():
        return []

    orfs = []
    for rec in SeqIO.parse(str(tmp_out), "fasta"):
        aa_seq = str(rec.seq).replace("*", "")
        if len(aa_seq) < CONFIG["orf_min_aa"]:
            continue
        desc_lower = rec.description.lower()
        protein_name = next((k for k in TE_DOMAIN_KEYWORDS if k in desc_lower), None)
        if protein_name is None:
            protein_name = "ORF"
        orfs.append({
            "accession": rec.id,
            "protein_name": protein_name,
            "sequence": aa_seq,
            "source": "getorf",
            "seq_type": "aa",
        })

    return _rank_orfs(orfs)


def extract_proteins(elements: list[dict], conda_env: str, cache_dir: str) -> list[dict]:
    """
    For each in_scope element, produce a list of AA sequence dicts under 'proteins'.
    Priority: annotated NCBI AA > getorf on Dfam nt.
    """
    for elem in elements:
        if elem["scope"] != "in_scope":
            elem["proteins"] = []
            continue

        retrieved = elem.get("retrieved", [])
        aa_seqs = [r for r in retrieved if r["seq_type"] == "aa"]
        nt_seqs = [r for r in retrieved if r["seq_type"] == "nt"]

        proteins = list(aa_seqs)  # annotated proteins go straight in

        if not proteins and nt_seqs:
            log.info("Running getorf on Dfam nt for: %s", elem["family_name"])
            for nt in nt_seqs:
                orfs = _run_getorf(nt["sequence"], elem["family_name"], conda_env, cache_dir)
                proteins.extend(orfs)

        if not proteins:
            log.warning("No proteins extracted for '%s'", elem["family_name"])

        elem["proteins"] = proteins
        log.info("  %s -> %d protein(s)", elem["family_name"], len(proteins))

    return elements


# ---------------------------------------------------------------------------
# Step 4 — Write combined protein FASTA
# ---------------------------------------------------------------------------

def write_protein_fasta(elements: list[dict], fasta_dir: str) -> str:
    """Write all proteins to a combined FASTA; returns the output path."""
    Path(fasta_dir).mkdir(parents=True, exist_ok=True)
    out_path = str(Path(fasta_dir) / "combined.fasta")

    records = []
    for elem in elements:
        eid = elem.get("element_id")
        for i, prot in enumerate(elem.get("proteins", []), start=1):
            raw_name = prot.get("protein_name") or f"protein{i}"
            safe_name = re.sub(r"[^a-zA-Z0-9_-]", "_", raw_name)
            seq_id = f"{eid}_{safe_name}_{i}"
            rec = SeqRecord(
                Seq(prot["sequence"]),
                id=seq_id,
                description=f"{elem['family_name']} | {raw_name} | {prot.get('accession', '')}",
            )
            records.append(rec)

    with open(out_path, "w") as f:
        SeqIO.write(records, f, "fasta")

    log.info("Wrote %d protein sequences to %s", len(records), out_path)
    return out_path


# ---------------------------------------------------------------------------
# Step 5 — Run PLAAC
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Step 6 — PLAAC plots (optional)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Step 7 — SQLite schema + loading
# ---------------------------------------------------------------------------

DDL = """
CREATE TABLE IF NOT EXISTS elements (
    element_id    INTEGER PRIMARY KEY AUTOINCREMENT,
    family_name   TEXT NOT NULL,
    te_class      TEXT,
    host_organism TEXT,
    derived_seq   TEXT,
    reference_num TEXT,
    scope         TEXT NOT NULL,
    scope_note    TEXT
);

CREATE TABLE IF NOT EXISTS sequences (
    seq_id        INTEGER PRIMARY KEY AUTOINCREMENT,
    element_id    INTEGER NOT NULL REFERENCES elements(element_id),
    seq_type      TEXT NOT NULL,   -- 'nt' or 'aa'
    source        TEXT,
    accession     TEXT,
    protein_name  TEXT,
    length        INTEGER,
    sequence      TEXT
);

CREATE TABLE IF NOT EXISTS plaac_staging (
    id           INTEGER PRIMARY KEY AUTOINCREMENT,
    alpha        REAL,
    raw_line     TEXT
);

CREATE TABLE IF NOT EXISTS plaac_results (
    result_id   INTEGER PRIMARY KEY AUTOINCREMENT,
    seq_id      INTEGER REFERENCES sequences(seq_id),
    alpha       REAL,
    core_length INTEGER,
    COREscore   REAL,
    LLR         REAL,
    NLLR        REAL,
    PRDstart    INTEGER,
    PRDstop     INTEGER,
    PRDscore    REAL,
    PAPAprop    REAL,
    PAPAfi      REAL,
    prd_called  INTEGER,
    run_date    TEXT,
    raw_tsv     TEXT
);
"""


def init_db(db_path: str) -> sqlite3.Connection:
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.executescript(DDL)
    conn.commit()
    return conn


def load_elements(conn: sqlite3.Connection, elements: list[dict]) -> list[dict]:
    """Insert elements into DB; attach element_id back to each dict."""
    cur = conn.cursor()
    for elem in elements:
        cur.execute("""
            INSERT INTO elements (family_name, te_class, host_organism, derived_seq,
                                  reference_num, scope, scope_note)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """, (elem["family_name"], elem["te_class"], elem["host_organism"],
              elem["derived_seq"], elem["reference_num"], elem["scope"], elem["scope_note"]))
        elem["element_id"] = cur.lastrowid
    conn.commit()
    return elements


def load_sequences(conn: sqlite3.Connection, elements: list[dict]) -> dict[str, int]:
    """
    Insert sequences and return a mapping of FASTA stable_id → seq_id
    so PLAAC results can be linked back.
    """
    cur = conn.cursor()
    fasta_id_to_seq_id: dict[str, int] = {}

    for elem in elements:
        eid = elem["element_id"]
        # nt sequences
        for rec in elem.get("retrieved", []):
            if rec["seq_type"] != "nt":
                continue
            cur.execute("""
                INSERT INTO sequences (element_id, seq_type, source, accession,
                                       protein_name, length, sequence)
                VALUES (?, 'nt', ?, ?, NULL, ?, ?)
            """, (eid, rec.get("source"), rec.get("accession"),
                  len(rec["sequence"]), rec["sequence"]))

        # aa sequences
        for i, prot in enumerate(elem.get("proteins", []), start=1):
            raw_name = prot.get("protein_name") or f"protein{i}"
            safe_name = re.sub(r"[^a-zA-Z0-9_-]", "_", raw_name)
            stable_id = f"{eid}_{safe_name}_{i}"

            cur.execute("""
                INSERT INTO sequences (element_id, seq_type, source, accession,
                                       protein_name, length, sequence)
                VALUES (?, 'aa', ?, ?, ?, ?, ?)
            """, (eid, prot.get("source"), prot.get("accession"),
                  raw_name, len(prot["sequence"]), prot["sequence"]))
            fasta_id_to_seq_id[stable_id] = cur.lastrowid

    conn.commit()
    return fasta_id_to_seq_id


def _parse_plaac_tsv(tsv_path: str) -> pd.DataFrame | None:
    """Read PLAAC summary TSV; return DataFrame or None if unreadable."""
    try:
        df = pd.read_csv(tsv_path, sep="\t", comment="#")
        log.info("PLAAC TSV: %d rows, columns: %s", len(df), list(df.columns))
        return df
    except Exception as e:
        log.error("Failed to parse PLAAC TSV %s: %s", tsv_path, e)
        return None


_PLAAC_COL_MAP = {
    # Maps common PLAAC output column names → our DB column names
    "SEQname": "seq_name",
    "COREscore": "COREscore",
    "LLR": "LLR",
    "NLLR": "NLLR",
    "PRDstart": "PRDstart",
    "PRDstop": "PRDstop",
    "PRDscore": "PRDscore",
    "PAPAprop": "PAPAprop",
    "PAPAfi": "PAPAfi",
}


def load_plaac_results(conn: sqlite3.Connection,
                       tsv_path: str,
                       alpha: float,
                       core_length: int,
                       fasta_id_to_seq_id: dict[str, int]) -> list[str]:
    """
    Parse PLAAC TSV, insert into plaac_results, return list of seq names with prd_called=1.
    """
    from datetime import date
    df = _parse_plaac_tsv(tsv_path)
    if df is None:
        return []

    cur = conn.cursor()
    prd_positive_ids = []
    run_date = date.today().isoformat()

    for _, row in df.iterrows():
        # PLAAC column names vary slightly across versions; be flexible.
        # PLAAC's SEQname keeps the whole FASTA header, but fasta_id_to_seq_id
        # is keyed on just the leading id token, so look up by that token.
        # PLAAC's actual column is "SEQid", not "SEQname"
        seq_name = str(row.get("SEQid", row.get("SEQname", row.get("NAME", row.get("name", "")))))
        seq_id = fasta_id_to_seq_id.get(seq_name.split()[0] if seq_name else seq_name)

        core_score = float(row.get("COREscore", row.get("CORE", 0)) or 0)
        llr = float(row.get("LLR", 0) or 0)
        nllr = float(row.get("NLLR", 0) or 0)
        prd_start = row.get("PRDstart", row.get("PRDSTART"))
        # PLAAC's column is "PRDend", not "PRDstop"
        prd_stop = row.get("PRDend", row.get("PRDstop", row.get("PRDSTOP")))
        prd_score = float(row.get("PRDscore", row.get("PRDSCORE", 0)) or 0)
        papa_prop = float(row.get("PAPAprop", row.get("PAPAPROP", 0)) or 0)
        papa_fi = float(row.get("PAPAfi", row.get("PAPAFI", 0)) or 0)

        # PLAAC reports PRDstart=0/PRDend=-1 placeholders when no domain is found,
        # so presence of values isn't a valid positive test — PRDlen>0 is.
        prd_len = row.get("PRDlen", row.get("PRDLEN", 0))
        prd_called = 1 if (prd_len is not None and str(prd_len).strip() not in ("", "NA", "nan")
                           and float(prd_len) > 0) else 0

        raw_tsv = "\t".join(str(v) for v in row.values)

        cur.execute("""
            INSERT INTO plaac_results
              (seq_id, alpha, core_length, COREscore, LLR, NLLR,
               PRDstart, PRDstop, PRDscore, PAPAprop, PAPAfi,
               prd_called, run_date, raw_tsv)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (seq_id, alpha, core_length, core_score, llr, nllr,
              prd_start, prd_stop, prd_score, papa_prop, papa_fi,
              prd_called, run_date, raw_tsv))

        if prd_called:
            prd_positive_ids.append(seq_name)

    conn.commit()
    log.info("Loaded %d PLAAC results (alpha=%.1f); %d PrLD-positive",
             len(df), alpha, len(prd_positive_ids))
    return prd_positive_ids


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Screen TE-encoded proteins for prion-like domains using PLAAC."
    )
    parser.add_argument("--table", default=CONFIG["plan_doc"],
                        help="Path to plan document containing the element table")
    parser.add_argument("--db", default=CONFIG["db"])
    parser.add_argument("--alpha", type=float, default=CONFIG["alpha"])
    parser.add_argument("--conda-env", default=CONFIG["conda_env"],
                        dest="conda_env",
                        help="Name of Conda environment containing EMBOSS getorf")
    parser.add_argument("--skip-plots", action="store_true")
    parser.add_argument("--plot-all", action="store_true",
                        help="Plot every scored protein, not just PrLD-positive ones")
    parser.add_argument("--skip-fetch", action="store_true",
                        help="Use only cached data; make no network calls")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print PLAAC commands without executing them")
    parser.add_argument("--stability-check", action="store_true",
                        help="Also run PLAAC at alpha=0 and alpha=1")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(args.verbose)

    # Override CONFIG from args
    CONFIG["conda_env"] = args.conda_env
    CONFIG["alpha"] = args.alpha

    log.info("=== TE-PLAAC pipeline starting ===")

    # Step 1: Parse table
    elements = parse_table(args.table)

    # Step 2: Retrieve sequences
    elements = fetch_sequences(elements, CONFIG["cache_dir"], args.skip_fetch)

    # Step 3: Extract proteins
    elements = extract_proteins(elements, args.conda_env, CONFIG["cache_dir"])

    # Step 4: Init DB + load elements
    conn = init_db(args.db)
    elements = load_elements(conn, elements)

    # Build FASTA ID → seq_id mapping (must happen before write_protein_fasta
    # because stable IDs are constructed from element_id which comes from DB)
    fasta_id_to_seq_id = load_sequences(conn, elements)

    # Step 5: Write combined FASTA
    fasta_path = write_protein_fasta(elements, CONFIG["fasta_dir"])

    # Step 6: Run PLAAC at primary alpha
    alphas = [args.alpha]
    if args.stability_check:
        alphas += [0.0, 1.0]

    prd_positive_ids: list[str] = []
    for a in alphas:
        tsv = run_plaac(fasta_path, CONFIG["plaac_jar"], a,
                        CONFIG["core_length"], CONFIG["plaac_output_dir"],
                        dry_run=args.dry_run)
        if tsv and not args.dry_run:
            pos = load_plaac_results(conn, tsv, a, CONFIG["core_length"], fasta_id_to_seq_id)
            if a == args.alpha:
                prd_positive_ids = pos

    # Step 7: Plots
    if not args.skip_plots:
        run_plaac_plots(fasta_path, CONFIG["plaac_jar"], CONFIG["plaac_plot_r"],
                        prd_positive_ids, CONFIG["plots_dir"], alpha=args.alpha,
                        core_length=CONFIG["core_length"], dry_run=args.dry_run,
                        rscript=CONFIG.get("rscript", "Rscript"), plot_all=args.plot_all)

    conn.close()
    log.info("=== Pipeline complete. DB: %s ===", args.db)

    # Summary
    total_proteins = sum(len(e.get("proteins", [])) for e in elements if e["scope"] == "in_scope")
    log.info("  In-scope families scored: %d",
             sum(1 for e in elements if e["scope"] == "in_scope" and e.get("proteins")))
    log.info("  Total protein sequences: %d", total_proteins)
    log.info("  PrLD-positive at alpha=%.1f: %d", args.alpha, len(prd_positive_ids))


if __name__ == "__main__":
    main()
