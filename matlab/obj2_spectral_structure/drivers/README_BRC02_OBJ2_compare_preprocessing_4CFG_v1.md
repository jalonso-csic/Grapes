# BRC02_OBJ2_compare_preprocessing_4CFG_v1.m

Compare four spectral preprocessing configurations (RAW, SNV, SNV+Savitzky–Golay 1st derivative, SNV+Savitzky–Golay 2nd derivative) and quantify how they affect (i) descriptive spectra, (ii) PCA structure, and (iii) multivariate factorial effects (Freedman–Lane permutation MANOVA-like test) and ASCA-style variance partitioning.

## What this script does

For each preprocessing configuration (`RAW`, `SNV`, `SNV_SG1`, `SNV_SG2`), the script:

1. Reads an input table containing factor columns and spectral columns (`nm_<wavelength>`).
2. Extracts and sorts the spectral matrix `Xraw` and wavelength vector `wl`.
3. Applies preprocessing:
   - `RAW`: no transformation
   - `SNV`: row-wise SNV (per sample)
   - `SNV_SG1`: SNV + Savitzky–Golay 1st derivative
   - `SNV_SG2`: SNV + Savitzky–Golay 2nd derivative
4. Creates two versions of X:
   - `Xplot`: preprocessed but **not autoscaled** (used for meaningful spectral plots)
   - `Xmodel`: **column-autoscaled** (used for PCA, factorial test, and ASCA)
5. Exports per-configuration outputs (figures + tables).
6. Builds cross-configuration comparisons (Excel summary + comparison figures).

## Inputs

### Required input file

Edit in the **USER SETTINGS** section:

```matlab
INFILE = "Matriz_Brocoli_SUM_1nm_ASCII.xlsx";
```

### Required columns in the table

**Spectra**
- Column names must start with `nm_` and be strictly numeric wavelengths, e.g. `nm_250`, …, `nm_1800`.

**Factors (auto-detected)**
- **Part**: column name containing `parte` or `part`
- **Maturity**: column name containing `madur` or `matur`
- **N2**: column name containing `n2`, `nitro`, or `aplic`

All three are converted to `categorical`.

## Outputs

Edit in the **USER SETTINGS** section:

```matlab
OUTDIR = fullfile(pwd, "OBJ2_OUT_4CFG");
```

The script creates:

### Per-configuration folders

- `OUTDIR/CFG_RAW/`
- `OUTDIR/CFG_SNV/`
- `OUTDIR/CFG_SNV_SG1/`
- `OUTDIR/CFG_SNV_SG2/`

Each `CFG_*` folder contains:

- `OBJ2_PreprocInfo_CFG_<CFG>.mat`
- `OBJ2_Spectra_MeanSD_CFG_<CFG>.xlsx`
- `OBJ2_MeanSpectraBy_Part_CFG_<CFG>.xlsx`
- `OBJ2_MeanSpectraBy_Maturity_CFG_<CFG>.xlsx`
- `OBJ2_MeanSpectraBy_N2_CFG_<CFG>.xlsx`
- `OBJ2_GlobalMeanSD_CFG_<CFG>.fig` and `.png`
- `OBJ2_PCA_Explained_CFG_<CFG>.xlsx`
- `OBJ2_PCA_Loadings_CFG_<CFG>.xlsx`
- `OBJ2_PCA_PC1PC2_by_<Factor>_CFG_<CFG>.fig/.png`
- `OBJ2_PCA_PC1PC3_by_<Factor>_CFG_<CFG>.fig/.png` (if enabled)
- `OBJ2_SpectralFactorialTest_CFG_<CFG>.xlsx`
- `OBJ2_ASCA_VariancePartition_CFG_<CFG>.xlsx`
- `OBJ2_ASCA_LoadingPC1_<Term>_CFG_<CFG>.xlsx`
- `OBJ2_ASCA_PC1Loading_<Term>_CFG_<CFG>.fig/.png`

### Cross-configuration summary

- `OUTDIR/OBJ2_Compare_4CFG.xlsx`
  - Sheets: `PCA_Explained`, `FactorialTest`, `ASCA_Variance`, plus wide sheets `eta2p_wide` and `qBH_wide` (when available).

### Cross-configuration figures

- `OUTDIR/COMPARE_FIGS/`
  - PCA explained variance comparison (PC1–PC3)
  - Partial eta² comparison across terms
  - ASCA effect PC1 loading overlays per term

## Key parameters you may edit

In **USER SETTINGS**:

- PCA export: `N_PCS_EXPORT` (default `10`), `MAKE_PC13` (default `true`)
- Permutations: `SEED0` (default `123`), `N_PERM` (default `499`)
- Savitzky–Golay: `SG_POLY` (default `2`), `SG_WINDOW` (default `11`, must be odd)

## Runtime notes

- The permutation factorial test scales with `N_PERM` and the size of `Xmodel`.
  - Smoke test: `N_PERM = 199`
  - Final: `N_PERM >= 999` (runtime increases substantially)

## Dependencies

This script is self-contained (local functions at the end) and uses standard MATLAB functions. Toolboxes typically required:

- Statistics and Machine Learning Toolbox (`pca`, `dummyvar`)
- Signal Processing Toolbox (`sgolay`)

## How to run

1. Place the input Excel file in the working directory (or set an absolute path in `INFILE`).
2. Set MATLAB working directory to the script folder.
3. Run:

```matlab
BRC02_OBJ2_compare_preprocessing_4CFG_v1
```

## Reproducibility

- Permutation testing is seeded with `rng(SEED0)` within each CFG loop, giving reproducible permutation sequences for a given MATLAB version/platform.

## Suggested citation

If you release this script publicly, include a citation block in the repository-level README (or Zenodo record) pointing to the associated manuscript and dataset DOI.
