# Setup Instructions

## 1. Java JRE (required for PLAAC)

Install via winget (run in PowerShell as Administrator):
```powershell
winget install Microsoft.OpenJDK.21
```
Or download manually from https://adoptium.net/temurin/releases/

Verify: `java -version`

---

## 2. PLAAC jar + R scripts

Download the latest release from https://github.com/whitehead/plaac/releases

You need:
- `plaac.jar`
- `R/plaac_plot.r`
- `R/plaac_plot_util.r`

Place them in this layout inside the repo:
```
tools/
└── plaac/
    ├── plaac.jar
    └── R/
        ├── plaac_plot.r
        └── plaac_plot_util.r
```

---

## 3. R + Rscript (required for plots)

Download R from https://cran.r-project.org/bin/windows/base/

After installing, open R or RScript and run:
```r
install.packages(c("DBI", "RSQLite"))
```

Verify: `Rscript --version`

---

## 4. Conda environment with Python dependencies

### 4a. Install Miniconda (if not already installed)

In PowerShell (per-user install, no admin needed):
```powershell
winget install --id Anaconda.Miniconda3 -e --source winget
```
Then close and reopen your terminal so `conda` is on PATH. If `conda` still isn't found, run:
```powershell
& "$env:USERPROFILE\miniconda3\Scripts\conda.exe" init powershell
```
and restart the shell.

### 4b. Create the `prioneers` environment
EMBOSS is not available via bioconda on native Windows, so the environment only
needs Python + Biopython. The pipeline detects this automatically and falls back
to its built-in Biopython ORF finder (`_run_getorf` in `te_plaac_pipeline.py`)
instead of calling `getorf`.

```bash
conda create -n prioneers python=3.11 -y
conda activate prioneers
pip install -r requirements.txt
```

### 4c. Run the pipeline inside the environment
```bash
conda activate prioneers
python te_plaac_pipeline.py
```

Or use the `--conda-env prioneers` flag if calling from outside the environment.
On Windows, the `--conda-env` flag mainly just needs the env to exist with
Biopython installed — `getorf` is optional and will be skipped if missing.

---

## 5. Quick-check all dependencies

```bash
java -version
Rscript --version
conda activate prioneers && python -c "import Bio, requests, pandas; print('Python deps OK')"
```

All three should print version strings without errors.
