# Setup Instructions

## 1. Java JRE (required for PLAAC)

Install via winget (run in PowerShell as Administrator):
```powershell
winget install Microsoft.OpenJDK.21
```
Or download manually from https://adoptium.net/temurin/releases/

Verify: `java -version`

---

## 2. R + Rscript (required for plots)

Download R from https://cran.r-project.org/bin/windows/base/ (or your OS's
package manager). The bundled plotting scripts (`tools/plaac/R/`) use only
base R — no extra packages to install.

Verify: `Rscript --version`

If `Rscript` isn't on your `PATH`, set the `RSCRIPT_PATH` environment variable
to its full path instead of relying on `PATH` — see the README for details.

---

## 3. Python dependencies (conda)

```bash
conda create -n prioneers python=3.11 -y
conda activate prioneers
pip install -r requirements.txt
```

---

## 4. Quick-check all dependencies

```bash
java -version
Rscript --version
conda activate prioneers && python -c "import Bio, requests, pandas; print('Python deps OK')"
```

All three should print version strings without errors.

PLAAC itself (`tools/plaac/plaac.jar` and its R scripts) is already bundled in
this repo — no separate download needed.
