# README — `BRC02_Fig23_SpectraPCA_SNV_SG1st_SG2nd_v1.m`

## What this script does
This script generates **two double‑panel, publication‑ready figures** for the manuscript **Section 3.2**:

- **Fig. 2** — *Mean spectra by plant part (Part)*  
  **(A)** SNV  
  **(B)** SNV + Savitzky–Golay **1st derivative**

- **Fig. 3** — *PCA score space (PC1 vs PC2), coloured by plant part (Part)*  
  **(A)** SNV  
  **(B)** SNV + Savitzky–Golay **2nd derivative**

Key points:
- The **spectral mean curves** are computed **per Part level** (one curve per level).
- PCA is run on **autoscaled** spectra (**column‑wise mean‑centre + unit variance**), consistent with the OBJ2 compare‑preprocessing logic.
- Savitzky–Golay derivatives use **mirror padding** to reduce endpoint artefacts.

## Inputs
- Excel file: `Matriz_Brocoli_SUM_1nm_ASCII.xlsx`
- Sheet: `Matriz` (configurable)
- Spectral variables must be named: `nm_<integer>` (e.g., `nm_250 … nm_1800`)

Required metadata column:
- `Part` (or `Parte`) — used to group curves and colour PCA scores.

Optional metadata columns (only required if `USE_CORE_FILTER = true`):
- `Cultivar` / `Variety` / `Variedad`
- `Extraction` / `Extraccion` / `Extracción`

## Outputs
All outputs are written to:
- `FIG_3p2_OUT/` (created if missing)

Files generated:
- `Fig2_Spectra_SNV_vs_SG1st.png` and `.fig`
- `Fig3_PCA_SNV_vs_SG2nd.png` and `.fig`

Figures are exported at **400 dpi**.

## How to run
1. Place the script in your MATLAB working directory (or add it to the path).
2. Ensure the input Excel file is accessible (same folder, or edit `INPUT_XLSX`).
3. Edit the **USER PARAMETERS** block if needed (file name, sheet, filtering).
4. Run:

```matlab
BRC02_Fig23_SpectraPCA_SNV_SG1st_SG2nd_v1
```

## Key parameters (USER PARAMETERS)
- `INPUT_XLSX`, `SHEET_NAME`: input dataset location.
- `OUTDIR`: output folder.
- `USE_CORE_FILTER`:
  - `true`: subsets rows by `CULTIVAR_CORE` and `EXTRACTION_CORE`.
  - `false`: uses all rows.
- Savitzky–Golay:
  - `SG_POLY_ORDER` (default 2)
  - `SG_FRAME_LEN` (default 11, **must be odd**)
  - `DELTA_NM` (default 1; wavelength step in nm)
- Plot style:
  - `FONT_NAME`, `FS_AX`, `FS_LAB`

## Interpretation notes
- **Fig. 2** is intended as a **visual reference** of how preprocessing changes separability/shape across Parts; it is not a statistical test.
- **Fig. 3** shows how the **Part separation** appears in PCA space under two preprocessing choices. Percent variance is reported on axes.

## Reproducibility
- The PCA implementation uses deterministic SVD (`Algorithm='svd'`) and no random sampling.
- If the input file is unchanged and the same MATLAB version is used, figures should be reproducible up to minor rendering differences.

## Dependencies
- MATLAB toolboxes/functions used:
  - `detectImportOptions`, `readtable`
  - `sgolay` (Signal Processing Toolbox)
  - `pca` (Statistics and Machine Learning Toolbox)
  - `exportgraphics` (fallback to `print` if unavailable)

## Suggested repository placement
Recommended location within `BRC02 / Objective 2` materials:
- `02_objective2/figures/` (script)
- `02_objective2/figures/FIG_3p2_OUT/` (generated outputs; usually **not** committed)

If you later add more scripts for Objective 2, keep **one README per script** (this file), and optionally a **folder‑level** `README.md` summarising how scripts relate to each other.

