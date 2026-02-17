# Objective 3 — Tier3S PLSR screening (SNV+SG 1st vs 2nd derivative)

This script implements **Objective 3 (Tier3S)** for the broccoli UV–VIS–NIR chemometrics workflow. It performs an **audit-ready, repeated nested cross-validation (CV)** PLSR screening for **all Tier 3 endpoints** (i.e., all remaining **numeric, non-spectral, non-factor** variables after excluding Tier 1 and Tier 2) within the **Pathernon–UAE (Ultrasonido) core** subset.

Unlike the main Objective 3 script (Tier 1 + Tier 2; 4 preprocessing strategies), this Tier3S script compares **only two preprocessing pipelines** (both include SNV):

- **A. SNV + Savitzky–Golay 1st derivative** (poly = 2, window = 11, Δ = 1 nm)  
- **B. SNV + Savitzky–Golay 2nd derivative** (poly = 2, window = 11, Δ = 1 nm)

The goal is **fast, consistent screening** of many endpoints using the same nested-CV PLSR framework.

---

## Script

**Filename:** `BRC06c_OBJ3_TIER3S_PLSR_NestedCV_SNV_SG1st_vs_SG2nd_v1.m`

**Recommended location in this repository:**
- `matlab/obj3_supervised_models/drivers/`

(Older or superseded versions can go to `matlab/obj3_supervised_models/legacy/`.)

---

## Inputs

### 1) Data matrix (Excel)
- `INPUT_XLSX` (default): `Matriz_Brocoli_SUM_1nm_ASCII.xlsx`
- `SHEET_NAME` (default): `Matriz`

### 2) Required columns in the matrix
The script expects **factor columns** and **spectral columns**:

**Factor columns (names are detected robustly):**
- Cultivar/Variety: `Variedad` / `Variety` / `Cultivar`
- Plant part: `Parte` / `Part`
- Maturity: `Maduracion` / `Maduración` / `Maturity`
- Nitrogen treatment: `Aplicacion_N2` / `Aplicación_N2` / `N2` / `Nitrogen`
- Extraction: `Extraccion` / `Extracción` / `Extraction`
- Optional ID/code: `Codigo` / `Código` / `Code`

**Spectral columns:**
- Must start with `nm_` (e.g., `nm_250`, `nm_251`, …).  
- The script assumes a **1 nm grid** (`DELTA_NM = 1`). If your grid is different, update `DELTA_NM` accordingly.

### 3) Tier definitions (must match exact headers)
- **Tier 1 list:** `TIER1_LIST`
- **Tier 2 list:** `TIER2_LIST`

Tier3S is then defined automatically as **everything else numeric** (after excluding factors + spectra + Tier 1 + Tier 2), subject to QC.

---

## Tier3S endpoint definition and QC

Tier3S endpoints are built automatically from the core-filtered table (`Tcore`) using these rules:

1. Exclude all **spectral** variables (`nm_*`).
2. Exclude all **factor** columns (Variety, Part, Maturity, N2, Extraction, optional Code).
3. Exclude variables listed in `TIER1_LIST` and `TIER2_LIST`.
4. Keep only **numeric** variables.
5. Apply QC filters:
   - at least `MIN_N_NONNAN` usable rows (default: 20)
   - NaN fraction ≤ `MAX_NAN_FRAC` (default: 0.30)
   - non-zero variance

Optional exclusions can be applied via `EXCLUDE_REGEX` (e.g., to skip certain endpoint families by name prefix).

The script exports:
- `OBJ3S_Tier3S_Included.xlsx`
- `OBJ3S_Tier3S_Skipped.xlsx` (with reasons)

---

## Modelling protocol (audit-ready)

### Core subset
The analysis is restricted to:
- `CULTIVAR_CORE = "Pathernon"`
- `EXTRACTION_CORE = "Ultrasonido"`

### Cross-validation design
- **Outer CV:** K = 5 folds, **repeated R = 50** times
- **Stratification:** by `Maturity × N2`
- **Guardrail (recommended):** each outer **TRAIN** split must include **all Part levels** (`ENFORCE_PART_COVERAGE = true`)

### Inner CV (LV selection)
- **Inner CV:** K = 4
- LV selection: **min-RMSE** over `LV = 1..LVmax`
- `LVmax = min(LV_MAX_CAP, n_train − 2, p)` with `LV_MAX_CAP = 15`

### Preprocessing
For each outer fold (train/test):
- Apply SNV (row-wise, sample-wise)  
- Apply Savitzky–Golay derivative (1st or 2nd)  
- Optionally standardise predictors with **TRAIN statistics only** (`SCALE_X = true`)

### Metrics reported
For each endpoint:
- **R²**, **RMSE**, **MAE**, **Bias** (mean(y_pred − y_true)), **RPD** (= sd(y) / RMSE)
- Reported as:
  - **pooled** over all outer predictions (`*_pooled`)
  - **repeat-level** mean and SD across repeats (`*_mean`, `*_sd`)
  - `LV_median` (median latent variables used across outer fits)

---

## Outputs

By default, outputs are written to:
- `OUTDIR = ./OBJ3_TIER3S_OUT/` (relative to the current working directory)

For each preprocessing option (`SNV_SG1st`, `SNV_SG2nd`), the script writes:

### Core endpoint lists
- `OBJ3S_Tier3S_Included.xlsx`
- `OBJ3S_Tier3S_Skipped.xlsx`

### Metrics
- `OBJ3S_Summary_Tier3S_<Config>.xlsx`
- `OBJ3S_RepeatMetrics_Tier3S_<Config>.xlsx`

### Optional artefacts (controlled by switches)
- **Predictions per endpoint** (`WRITE_PREDICTIONS = true`):
  - `OBJ3S_Predictions_Tier3S_<Config>.xlsx` (one sheet per Y)
- **Raw VIP per outer fit** (`SAVE_VIPRAW = true`):
  - `OBJ3S_VIPraw_Tier3S_<Config>.mat`  
  (intended for downstream VIP-stability mapping)

---

## How to run

1. Open MATLAB and set the working directory to the repository root (or to `matlab/obj3_supervised_models/drivers/`).
2. Ensure `INPUT_XLSX` points to the correct file location.
3. Run the script:
   - Either press **Run** in the MATLAB editor, or execute:
     ```matlab
     BRC06c_OBJ3_TIER3S_PLSR_NestedCV_SNV_SG1st_vs_SG2nd_v1
     ```

> Tip: Do **not** commit `OBJ3_TIER3S_OUT/` outputs to Git. Add the output directory to `.gitignore` if you run inside the repo.

---

## Requirements

- MATLAB (tested conceptually with modern releases; adjust if needed)
- **Statistics and Machine Learning Toolbox** (`plsregress`, `cvpartition`)
- **Signal Processing Toolbox** (`sgolay`)

---

## Notes and cautions

- The script assumes a **uniform wavelength grid** consistent with `DELTA_NM` for derivative scaling.
- Tier3S can be computationally heavy if many endpoints pass QC. Use:
  - `WRITE_PREDICTIONS = false` (default) to reduce I/O
  - `SAVE_VIPRAW = true/false` depending on whether you need VIP stability downstream

---

## Citation / provenance (recommended)

If this repository is linked to a manuscript or Zenodo release, add the corresponding citation here (DOI, version tag, and license).
