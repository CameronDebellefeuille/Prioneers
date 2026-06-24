# Plan: TE-Encoded Proteins → PLAAC → SQLite

*Preliminary pipeline for screening transposable-element proteins for prion-like amino acid composition (PrLDs).*

---

## 0. The one constraint that shapes everything

**PLAAC scores proteins, not DNA.** It detects prion-like *amino-acid* composition, so its input is protein FASTA. Two consequences:

- The right-hand column of the source table (satDNA, micro-/minisatellites) is **non-coding** and is metadata only — never sent to PLAAC.
- In the left-hand column, only **autonomous, protein-coding TEs** are PLAAC-able (they encode Gag, ORF1p/ORF2p, transposase, RT/integrase, Rep/helicase, etc.). Non-autonomous elements (SINEs/Alu, MITEs) encode no protein and are excluded.

**Background frequency setting:** run PLAAC at **α = 0.5** (the authors' recommended setting for cross-species scanning). α is exposed as a script parameter so the set can be re-scored at α = 0 and α = 1 to check call stability. Core length stays at the default 60 aa.

---

## 1. Scope (confirmed)

- Target = **TE-encoded proteins from autonomous families only**.
- Broad/generic group labels are **flagged in the database but not run** through PLAAC.
- Non-autonomous / non-coding entries are **recorded as metadata, excluded** from PLAAC.
- This is preliminary work; no attempt to force-retrieve every obscure element yet.

### Element classification (families de-duplicated; LINE-1 appears 4×, Alu 6×, MITE 4× in the table)

| Family (left column) | Encoded protein(s) | Disposition |
|---|---|---|
| LINE-1 | ORF1p, ORF2p (RT/EN) | **In scope** |
| TART / HeT-A | Gag (+ HeT-A ORF) | **In scope** |
| CR1 / CR1-C | ORF2-like (RT/EN) | **In scope** |
| Atenspm2 (CACTA/EnSpm) | Transposase | **In scope** |
| Crwydryn (rye LTR retro) | Gag-Pol | **In scope** (retrieval may need cited paper) |
| Sore1 (potato LTR retro) | Gag-Pol | **In scope** (may need cited paper) |
| Ty3/gypsy **ogre** elements | Gag-Pol (Tat/ogre clade) | **In scope** |
| Gmr9 / Gm ogre (soybean) | Gag-Pol | **In scope** (may need cited paper) |
| CRM1 / CRM4 (maize centromeric) | Gag-Pol | **In scope** |
| Helitrons | Rep / helicase | **In scope** |
| "LTR retrotransposons" (generic) | — | **Flag: broad group, do not run** |
| "Ty3/gypsy-retroelement" (generic) | — | **Flag: broad group, do not run** |
| pDv mobile element | uncertain | **Flag: verify coding status before running** |
| SGM-IS transposons | uncertain | **Flag: verify coding status before running** |
| MaLR retrotransposon | largely non-autonomous | **Flag: likely non-coding, default exclude** |
| SINE / SINE-like / Alu / SINE B1 | none | **Exclude (non-coding)** |
| MITE / MITE-like | none | **Exclude (non-coding)** |

> Expect good coverage for LINE-1, CR1, ogre, CRM, Helitron, TART/HeT-A via Dfam/NCBI, and manual digging in the cited papers for the named obscure elements (Crwydryn, Sore1, Gmr9, SGM-IS).

---

## 2. Sources (confirmed)

- **Dfam** (dfam.org) — primary for consensus sequences / HMMs / family records (open, has an API).
- **NCBI** nucleotide + protein via Entrez (Biopython) — for named accessions and annotated CDS/protein products.
- **Cited papers' supplementary data** — fallback for obscure named elements.
- Repbase: not available; not used.

Preference order for getting the actual **protein**: (1) an annotated protein/CDS product if one exists → (2) translate the annotated CDS → (3) call ORFs on the nucleotide consensus and translate.

---

## 3. Pipeline

1. **Parse the table** into structured metadata (one row per table entry).
2. **Classify** each entry: `in_scope` / `broad_group` / `non_coding`, with a short note. De-duplicate families for retrieval but keep every original row as metadata.
3. **Retrieve nucleotide sequences** for in-scope families (Dfam first, NCBI for named accessions). Store nt sequence + source + accession.
4. **Obtain proteins:** prefer annotated protein/CDS; otherwise call ORFs (EMBOSS `getorf` or NCBI ORFfinder) and translate. Filter spurious ORFs by length and, where feasible, by expected TE domains (e.g. Gag, RVT/RT, integrase `rve`, transposase, Rep/helicase).
5. **Write a combined protein FASTA** of all in-scope proteins, with stable IDs that encode `element_id` + protein name.
6. **Run PLAAC** (`java -jar plaac.jar`) on the protein FASTA at α = 0.5 → per-protein summary TSV.
7. **Generate plots:** run the jar a second time in per-residue mode for the proteins of interest (driven by a print-list), then call `Rscript plaac_plot.r` → PNG/PDF per protein.
8. **Parse PLAAC output** and **load everything into SQLite**.

---

## 4. Database (SQLite)

Single-file DB; clean to read from both R (`DBI` + `RSQLite`) and Python (`sqlite3` / pandas).

```
elements
  element_id      INTEGER PK
  family_name     TEXT       -- e.g. "LINE-1", "ogre"
  te_class        TEXT       -- LINE / LTR-Ty3gypsy / CACTA / Helitron / ...
  host_organism   TEXT
  derived_seq     TEXT       -- right-column satDNA/micro/minisat (metadata)
  reference_num   TEXT       -- citation number from the table
  scope           TEXT       -- in_scope | broad_group | non_coding
  scope_note      TEXT

sequences
  seq_id          INTEGER PK
  element_id      INTEGER FK -> elements.element_id
  seq_type        TEXT       -- nt | aa
  source          TEXT       -- Dfam | NCBI | paper
  accession       TEXT
  protein_name    TEXT       -- ORF1p / Gag / transposase / ...
  length          INTEGER
  sequence        TEXT

plaac_results
  result_id       INTEGER PK
  seq_id          INTEGER FK -> sequences.seq_id   (aa rows only)
  alpha           REAL
  core_length     INTEGER
  COREscore       REAL
  LLR             REAL
  NLLR            REAL
  PRDstart        INTEGER
  PRDstop         INTEGER
  PRDscore        REAL
  PAPAprop        REAL
  PAPAfi          REAL
  prd_called      INTEGER    -- boolean
  run_date        TEXT
  raw_tsv         TEXT       -- full original PLAAC line, for safety
```

> Suggestion: load the **entire** PLAAC TSV into a staging table with `pandas.to_sql` (so a future PLAAC version adding columns doesn't break the loader), then populate `plaac_results` from it. Keep `raw_tsv` so nothing is ever lost.

---

## 5. Dependencies

Hard requirements (independent of wrapper language):

- **Java JRE** — to run `plaac.jar`.
- **`plaac.jar`** + **`R/plaac_plot.r`** + **`R/plaac_plot_util.r`** from `github.com/whitehead/plaac`.
- **R + `Rscript`** — for the canonical plots.

Python (the orchestrator):

- `biopython` (Entrez/SeqIO), `requests` (Dfam API), `pandas`, `sqlite3` (stdlib), `argparse`, `logging`.
- **EMBOSS `getorf`** (or NCBI ORFfinder) for ORF calling, called via `subprocess`.

R (for downstream querying, since you're stronger in R): `DBI`, `RSQLite`.

---

## 6. Single Python script — structure

One `te_plaac_pipeline.py` with a config block up top and modular functions:

```
CONFIG: paths (plaac.jar, R scripts, db file, FASTA dir), alpha=0.5, core_len=60, ORF min length

parse_table(path)            -> list[dict]        # step 1
classify_elements(rows)      -> annotated rows    # step 2
fetch_sequences(rows)        -> nt seqs           # step 3 (Dfam + NCBI)
extract_proteins(nt_seqs)    -> aa seqs           # step 4 (annotated or getorf)
write_protein_fasta(aa_seqs) -> combined.fasta    # step 5
run_plaac(fasta, alpha)      -> summary.tsv       # step 6 (subprocess java)
run_plaac_plots(...)         -> png/pdf           # step 7 (subprocess Rscript)
parse_plaac_output(tsv)      -> dataframe         # step 8
load_to_sqlite(...)          -> db                # step 8

main() ties them together with argparse flags (e.g. --alpha, --skip-plots, --db)
```

Beginner-friendly notes for the build phase: keep network calls (Dfam/NCBI) cached to disk so re-runs don't re-download; rate-limit Entrez (set email + small sleep); log every retrieval miss so the manual-digging list is explicit.

---

## 7. Open items to settle while coding

- Confirm the Dfam family IDs / NCBI accessions for each in-scope family before bulk retrieval (worth a quick manual check so the script isn't guessing).
- Decide the ORF acceptance rule (min length + optional domain check) — affects which translated proteins reach PLAAC.
- For Helitrons specifically, retrieval of a clean autonomous Rep/helicase protein can be fiddly; flag if it stalls.
- Verify `pDv` and `SGM-IS` coding status before deciding whether they move from "flag" to "in scope."
