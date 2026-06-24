# Prioneers

A bioinformatics pipeline to screen transposable element (TE)-encoded proteins
for prion-like domains (PrLDs) using [PLAAC](https://github.com/whitehead/plaac).

It pulls TE-encoded proteins from UniProt by TE-associated Pfam domain (Gag,
RT, Integrase, and several transposase families) across taxonomy, scores
them with PLAAC, and writes the results to Excel plus per-residue PNG plots
for every PrLD-positive hit.

## Layout

| Path | Purpose |
|---|---|
| `te_bulk_pull.py` | Fetches/caches proteins from UniProt per Pfam domain, scores them with PLAAC, writes `data/results/bulk_proteins.xlsx` |
| `te_bulk_analysis.py` | Reads that Excel file, writes a per-domain PrLD-positive rate breakdown, and plots every PrLD-positive protein |
| `plaac_utils.py` | Shared helpers (run PLAAC, parse its TSV output, generate per-residue plots) used by both scripts above |
| `tools/plaac/` | Bundled PLAAC jar and its R plotting scripts |

## Setup

Full step-by-step instructions (Java, R, conda) are in
[`setup_instructions.md`](setup_instructions.md). Short version:

```bash
conda create -n prioneers python=3.11 -y
conda activate prioneers
pip install -r requirements.txt
```

You'll also need:
- **Java JRE** (to run `tools/plaac/plaac.jar`)
- **R + Rscript** (to generate plots; only base R is required, no extra packages)

`tools/plaac/` is already included in this repo — no separate PLAAC download needed.

### Rscript not on PATH?

Both scripts call `Rscript` assuming it's on your `PATH`. If it isn't (common
on Windows installs), set the `RSCRIPT_PATH` environment variable to the full
path of your `Rscript` executable instead of editing the scripts:

```bash
# Windows (PowerShell)
$env:RSCRIPT_PATH = "C:\Program Files\R\R-4.6.0\bin\Rscript.exe"

# macOS/Linux
export RSCRIPT_PATH=/usr/local/bin/Rscript
```

## Usage

```bash
conda activate prioneers
python te_bulk_pull.py --skip-fetch --per-domain 100 --verbose
python te_bulk_analysis.py --verbose
```

- `te_bulk_pull.py` fetches/caches proteins per TE-associated Pfam domain from
  UniProt, scores them with PLAAC, and writes
  `data/results/bulk_proteins.xlsx`. Drop `--skip-fetch` on first run to
  populate the UniProt cache.
- `te_bulk_analysis.py` reads that Excel file, writes a per-domain
  PrLD-positive rate breakdown (`data/results/per_domain_breakdown.csv`), and
  plots every PrLD-positive protein into
  `data/results/plots_prld_positive/`.

Both scripts also accept `--alpha`, `--seed`, `--dry-run`, and `--verbose`;
run with `--help` for the full list.

## Methodology

This is the part of the pipeline most worth understanding before trusting the
numbers it produces.

### PLAAC parameters

PLAAC scores a protein by comparing its amino acid composition against a
prion-like HMM trained on known yeast prion domains. Two parameters control
that scoring (both set in `CONFIG` in `te_bulk_pull.py`, overridable via
`--alpha` on the CLI):

- **`alpha = 0.5`** — blends the background amino acid distribution PLAAC
  scores against. `alpha=1` uses the S. cerevisiae amino acid background frequencies;
  `alpha=0` uses the frequencies of the organisms being scored. Some organisms/datasets
  have unusual compositional bias (e.g. low-complexity-rich proteomes) that
  would inflate false positives at `alpha=1`, so I used
  `0.5` as a balanced default rather than either extreme. PLAAC's own
  documentation suggests re-checking PRD calls at `alpha=0.0` and `alpha=1.0`
  for stability if you need confidence in a specific result — pass
  `--alpha 0.0` / `--alpha 1.0` as separate runs to do that.
- **`core_length = 60`** — the minimum contiguous "core" window (in
  residues) PLAAC's Viterbi parse must find before it will call a region a
  prion-like domain (PRD), passed as PLAAC's `-c` flag.

A protein counts as PrLD-positive (`prd_called = 1`) when PLAAC's own output
reports `PRDlen > 0` — i.e. PLAAC's HMM actually found a qualifying domain,
not just a nonzero score.

### Pfam domains searched

`te_bulk_pull.py` queries UniProt for proteins annotated with one of 6
TE-associated Pfam domains (`TE_PFAM_DOMAINS` in the script):

| Pfam ID | Label | Domain |
|---|---|---|
| PF03732 | Gag | Retrotrans_gag — Gag (LTR retrotransposon/endogenous retrovirus capsid) |
| PF00078 | RT | RVT_1 — reverse transcriptase (LINE/LTR retrotransposons) |
| PF00665 | Integrase | rve — integrase core domain (LTR retrotransposons) |
| PF03004 | Transposase_DDE | DDE_Tnp_1 — Mutator/En-Spm-like DNA transposase |
| PF01498 | Transposase_IS4 | Transposase_8 — IS4/hAT-like DNA transposase |
| PF03108 | Transposase_MULE | MULE — MuDR/Foldback-like DNA transposase |

Each was spot-checked against the UniProt API to confirm it returns real TE
proteins rather than an unrelated gene family sharing the same domain name.
No reliable Pfam domain exists specifically for Helitron Rep/helicase
(candidates like PF14214 returned 0 UniProt hits), so Helitrons are
under-represented in this search.

UniProt queries are additionally constrained to `length:[50 TO 3000]`
(`min_len`/`max_len` in `CONFIG`) to exclude short fragments and unusually
long multi-domain outliers before sequences ever reach PLAAC.

### Pool, fetch, and random selection

The protein set per domain isn't just "the first N UniProt hits" — that
would bias toward whatever UniProt's default sort order favors (e.g.
well-annotated model organisms). Instead:

1. For each Pfam domain, fetch a **pool** of up to `--fetch-pool` (default
   500) matching UniProt entries and cache it to
   `data/bulk_cache/te_<PFAM_ID>.json`.
2. Randomly subsample `--per-domain` (default 100) proteins from that pool
   using a seeded RNG (`--seed`, default 42), so results are reproducible
   across runs as long as the cached pool and seed don't change.
3. After combining all domains' samples, drop any duplicate UniProt
   accessions (a protein can carry more than one TE-associated domain, e.g.
   a Gag-Pol polyprotein, and get pulled by two separate domain queries).

Re-running with `--skip-fetch` reuses the cached pool with no network calls;
the sample drawn from it is identical given the same seed. Re-running
without `--skip-fetch` re-fetches each domain's pool fresh from UniProt, so
the sample can change if UniProt's matching set has changed since the cache
was written.

## Notes

- Pipeline outputs (Excel files, PLAAC TSVs, PNG plots, UniProt JSON caches
  under `data/`) are gitignored — they're regenerated by running the scripts
  above, not checked into version control.
- The pipeline currently reports TE-only PrLD rates per domain; a non-TE
  control-set comparison is planned but not yet implemented.
