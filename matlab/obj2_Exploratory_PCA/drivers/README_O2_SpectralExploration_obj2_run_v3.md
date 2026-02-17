# O2_SpectralExploration_obj2_run_v3 — Objective 2 (PCA Spectral Exploration)

This script generates the **Objective 2** outputs for the manuscript: **unsupervised spectral exploration via PCA** on **FX10 (VIS–NIR)** and **FX17 (SWIR)** sample-level mean spectra, comparing **RAW reflectance** vs **SNV** preprocessing. It produces **publication-ready figures** and an **Excel summary** (variance explained, diagnostic ANOVA on score PCs, and top loadings).

---

## What this script delivers (audit-ready)

### Figures (PNG + FIG)
Saved in: `Objetivo_2/`

**Main manuscript panels (Figure 3):**
- `Fig3A_FX10_raw_byDate` — FX10 RAW, coloured by **SamplingDate**
- `Fig3B_FX10_SNV_byIrrigation` — FX10 SNV, coloured by **IrrigationRegime**
- `Fig3C_FX17_raw_byDate` — FX17 RAW, coloured by **SamplingDate**
- `Fig3D_FX17_SNV_byIrrigation` — FX17 SNV, coloured by **IrrigationRegime**

**Additional recommended outputs (Figure 4 and Figure 5):**
- `Loadings_<Sensor>_<Preproc>` — PC1/PC2 loadings (per sensor × preprocessing)
- `MeanDiffSpectra_<Sensor>_<Preproc>_byIrrigation` — mean spectra by irrigation + difference spectrum (Irrigated − Rainfed)

**Optional supplementary-style diagnostics (enabled by default):**
- `FigS_<Sensor>_SNV_byDate` — SNV PCA coloured by SamplingDate (symmetry check)
- `FigS_<Sensor>_RAW_byIrrigation` — RAW PCA coloured by IrrigationRegime (symmetry check)

> Export format: **`.png` and `.fig` only** (no `.tif`).

---

## Excel outputs

A single Excel file is written to:
- `Objetivo_2/O2_PCA_Summary.xlsx`

Sheets:
- `Explained` — explained variance (%), PC1–PC5 (Table 3-ready)
- `ANOVA_PC1_PC2` — one-way ANOVA on **PC1 and PC2 scores** vs:
  - `SamplingDate`
  - `IrrigationRegime`
  with effect size **η²** (diagnostic, can support supplementary reporting)
- `TopLoadings` — top absolute loadings wavelengths (PC1/PC2), per sensor × preprocessing
- `RunInfo` — minimal run metadata (input file, options, etc.)

---

## Input requirements

### Required metadata columns (English matrix)
Your spreadsheet must include **these column names** (case-sensitive):
- `SamplingDate` (datetime or convertible to datetime)
- `IrrigationRegime` (categorical labels; supports `"irrigated/rainfed"` or `"riego/secano"`)
- `VineyardFloorManagement`

### Required spectral columns
Spectral columns must start with these prefixes:
- `FX10_nm_` for FX10 wavelengths
- `FX17_nm_` for FX17 wavelengths

Wavelengths can be encoded as:
- integer nm: `FX10_nm_720`
- decimal via `p`: `FX10_nm_720p5` → interpreted as `720.5`

The script automatically:
- extracts all columns matching each prefix,
- parses and sorts wavelengths,
- converts the table into numeric matrices `X10` and `X17`.

---

## How to run

1) Place `O2_SpectralExploration_obj2_run_v3.m` on your MATLAB path.

2) Run:

```matlab
O2_SpectralExploration_obj2_run_v3('Matriz_CHEM_HSI_MASTER_96.xlsx','Matriz');
```

3) Outputs will be written to:

```
<root_of_input_xlsx>/Objetivo_2/
```

---

## Processing details (methods traceability)

### Preprocessing
- **RAW:** spectra used as stored in the matrix.
- **SNV:** applied **row-wise** (per sample):

\[
x' = \frac{x - \mu_x}{\sigma_x}
\]

where μ and σ are computed per sample across wavelengths.

### QC steps
- Removes any sample (row) containing **NaN/Inf** in the spectra.
- Removes **near-constant wavelengths** (variance ≤ 1e−12).

### PCA computation
- PCA is computed using **SVD on mean-centred data** (`pcaSVD`), robust for **p ≫ n**.
- The script returns:
  - `score` (PC scores)
  - `coeff` (loadings)
  - `% explained variance`

### Diagnostic inference on scores
- One-way ANOVA is computed for **PC1 and PC2** scores vs:
  - `SamplingDate`
  - `IrrigationRegime`
- Effect size is reported as **η²**:

\[
\eta^2 = \frac{SSB}{SST}
\]

> Important: these ANOVAs are **diagnostic** (exploratory) and intended to support interpretation of PCA structure.

---

## Configurable options (edit near the top of the script)

Inside the script:
- `NCOMP` — number of PCs computed (default: 5)
- `DO_LOADINGS` — export loadings figures (default: `true`)
- `DO_MEAN_DIFF_SPECTRA` — export mean/difference spectra figures (default: `true`)
- `DO_EXTRA_SCORE_PLOTS` — export symmetry-check score plots (default: `true`)
- `TOPK` — number of top loadings wavelengths exported per PC (default: 15)

---

## Expected outputs checklist (quick audit)

After a successful run, `Objetivo_2/` should contain:
- `O2_log.txt`
- `O2_PCA_Summary.xlsx`
- `Fig3A_FX10_raw_byDate.(png|fig)`
- `Fig3B_FX10_SNV_byIrrigation.(png|fig)`
- `Fig3C_FX17_raw_byDate.(png|fig)`
- `Fig3D_FX17_SNV_byIrrigation.(png|fig)`
- `Loadings_FX10_RAW.(png|fig)` and `Loadings_FX10_SNV.(png|fig)` (if enabled)
- `Loadings_FX17_RAW.(png|fig)` and `Loadings_FX17_SNV.(png|fig)` (if enabled)
- `MeanDiffSpectra_*_byIrrigation.(png|fig)` for each sensor × preprocessing (if enabled)
- `FigS_*` diagnostics (if enabled)

---

## Notes for manuscript integration (Q1)

- **Table 3** should be taken from `O2_PCA_Summary.xlsx` → `Explained`.
- If you need quantitative support for statements like
  “sampling date dominates PCA structure”, use:
  - `ANOVA_PC1_PC2` (p-values + η²) as a supplementary table or text support.
- The PCA is exploratory; **do not claim anti-leakage CV** here.  
  Anti-leakage validation belongs to downstream supervised modelling sections.

---

## Provenance

Script: `O2_SpectralExploration_obj2_run_v3.m`  
Outputs folder: `Objetivo_2/` (created automatically next to the input Excel).
