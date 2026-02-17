# README — Objective 3 (BRC06) PLSR preprocessing benchmark (Tier 1 + Tier 2)

**Script:** `BRC06_OBJ3_PLSR_NestedCV_Tier1_Tier2_4PreprocCompare_v1.m`  
**Repo path (recommended):** `matlab/obj3_supervised_models/drivers/`  
**Objective:** Compare four spectral preprocessing strategies for **PLSR** prediction using **repeated nested cross-validation** for **Tier 1** endpoints and **Tier 2** class aggregates.

This README is **script-specific** (do not replace the general Objective 3 README).

---

## 1) What this script does

For each response variable (**Y**) in **Tier 1** and **Tier 2**, the script:

1. Loads an Excel matrix containing:
   - spectral predictors: columns named `nm_<integer>` (e.g., `nm_250`, …),
   - experimental factors (Part, Maturity, N2, Cultivar, Extraction),
   - response variables (Tier 1 endpoints and Tier 2 aggregates).

2. Filters a **core subset** (default in script):
   - `CULTIVAR_CORE = "Pathernon"`
   - `EXTRACTION_CORE = "Ultrasonido"`

3. Benchmarks **4 preprocessing configurations**:
   1. **RAW** — no preprocessing
   2. **SNV** — row-wise Standard Normal Variate
   3. **SNV + SG 1st derivative** — Savitzky–Golay (poly=2, window=11)
   4. **SNV + SG 2nd derivative** — Savitzky–Golay (poly=2, window=11)

4. Uses **repeated nested CV** (audit-friendly):
   - Outer CV: `K_OUTER = 5` folds, repeated `R_REPEATS = 50` times  
   - Outer stratification: `Maturity × N2`
   - Optional guardrail: each outer TRAIN split must include all **Part** levels
   - Inner CV (LV selection): `K_INNER = 4`, choose LV by **minimum RMSE**
   - LV upper bound: `LVmax = min(LV_MAX_CAP, nTrain-2, p)` (default `LV_MAX_CAP = 15`)

5. Exports performance metrics, out-of-fold predictions, and raw VIP vectors.

---

## 2) Requirements

### MATLAB toolboxes
- **Statistics and Machine Learning Toolbox** (`plsregress`, `cvpartition`)
- **Signal Processing Toolbox** (`sgolay`) for Savitzky–Golay derivatives

---

## 3) Input data expectations (Excel)

Default input in the script:
- `INPUT_XLSX = 'Matriz_Brocoli_SUM_1nm_ASCII.xlsx'`
- `SHEET_NAME = 'Matriz'`

Minimum required content:

### A) Spectral columns
- Must be named `nm_<integer>` (e.g., `nm_250`, `nm_251`, …).
- The script parses the numeric wavelengths and sorts columns ascending.

### B) Factor columns
The script finds these by robust matching (Spanish/English variants):
- Cultivar: `Variedad` / `Variety` / `Cultivar`
- Part: `Parte` / `Part`
- Maturity: `Maduracion` / `Maduración` / `Maturity`
- N2: `Aplicacion_N2` / `Aplicación_N2` / `N2` / `Nitrogen` (and close variants)
- Extraction: `Extraccion` / `Extracción` / `Extraction`

### C) Response columns
- Tier 1: entries in `TIER1_LIST`
- Tier 2: entries in `TIER2_LIST` (`SUM_*` aggregates)

**Y QC rules (script defaults):**
- Skip Y if `nanFrac > 0.30`, or `nNonNaN < 20`, or `std(Y)=0`.

---

## 4) How to run

From the repo root in MATLAB:

```matlab
run('matlab/obj3_supervised_models/drivers/BRC06_OBJ3_PLSR_NestedCV_Tier1_Tier2_4PreprocCompare_v1.m')
```

Outputs are written to `OUTDIR` (default: `./OBJ3_OUT_4CFG` under the current working directory).

---

## 5) Outputs created

Outputs are produced **per tier** and **per preprocessing config**.

### A) Per tier + per config
For each `Tier ∈ {Tier1, Tier2}` and `Config ∈ {RAW, SNV, SNV_SG1st, SNV_SG2nd}`:

1. `OBJ3_Summary_<Tier>_<Config>.xlsx`  
   One row per Y with pooled metrics and repeat-level mean ± SD:
   - `R2`, `RMSE`, `MAE`, `Bias`, `RPD`, and `LV_median`.

2. `OBJ3_RepeatMetrics_<Tier>_<Config>.xlsx`  
   Repeat-level metrics per Y (one row per repeat).

3. `OBJ3_Predictions_<Tier>_<Config>.xlsx`  
   One sheet per Y (long format): `row`, `repeat`, `fold`, `y_true`, `y_pred`, `LV`.

4. `OBJ3_VIPraw_<Tier>_<Config>.mat`  
   Raw VIP vectors per outer fit (for later VIP stability mapping).  
   Also stores wavelength vector `wl`, selected LVs, and the fold plan for audit.

### B) Per tier (all configs combined)
- `OBJ3_Compare_4CFG_<Tier>.xlsx`
  - `Summary_AllCFG` (long format)
  - `RepeatMetrics_AllCFG` (long format)
  - optional wide sheets if `unstack` succeeds:
    - `R2pooled_wide`, `RMSEpooled_wide`

---

## 6) Reproducibility notes

- For each Y, the script builds **one** outer fold plan and reuses it across the four configs (fair comparison).
- Folds are seeded (`BASE_SEED + repeat index`) for repeatability.
- The Part-coverage guardrail can fail on small/imbalanced designs; if so:
  1) check balance across `Maturity × N2` and `Part`,  
  2) reduce `K_OUTER`, or  
  3) set `ENFORCE_PART_COVERAGE = false`.

---

## 7) Runtime notes

This benchmark is compute-intensive (outer repeats × folds × inner LV search).  
For quick tests, temporarily reduce `R_REPEATS`, cap `LV_MAX_CAP`, or run only one tier.

---

## 8) Changelog

- **v1:** Initial script-specific README for the PLSR preprocessing benchmark (Tier 1 + Tier 2).
