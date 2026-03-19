# Objective 2 — Exploratory PCA of grape hyperspectral spectra

## Overview

This folder contains the MATLAB workflow used to perform **exploratory principal component analysis (PCA)** on the sample-level hyperspectral spectra of the grape dataset.

The script:

1. reads the master Excel matrix;
2. extracts the **FX10** and **FX17** spectral blocks;
3. computes PCA for both **RAW reflectance** and **SNV-preprocessed** spectra;
4. generates **manuscript-oriented score plots**, **loading plots**, and **mean/difference spectra plots**;
5. exports an **Excel summary** and a **run log**.

This objective corresponds to the **unsupervised spectral exploration stage** of the study and supports the manuscript sections describing the dominant structure of the spectral data.

---

## Folder contents

- `obj2_exploratory_pca.m`  
  Main MATLAB workflow for exploratory PCA.

- `README.md`  
  Documentation for this objective.

---

## Folder structure

```text
obj2_exploratory_pca/
  obj2_exploratory_pca.m
  README.md
```

---

## What the script does

For each sensor and preprocessing combination, the workflow:

- extracts the relevant spectral block from the input table;
- removes rows containing **NaN** or **Inf** values;
- removes near-constant wavelengths;
- applies **SNV row-wise** when requested;
- computes PCA using **SVD on mean-centred data**;
- exports:
  - PCA score plots,
  - PCA loading plots,
  - mean and difference spectra plots,
  - an Excel summary of explained variance,
  - diagnostic ANOVA tables for PC1 and PC2,
  - top-loading wavelengths for PC1 and PC2,
  - and a run log.

The four PCA combinations are:

- **FX10 RAW**
- **FX10 SNV**
- **FX17 RAW**
- **FX17 SNV**

---

## Expected input

### File type
- `.xlsx` master table

### Structure
- **one row per sample**

### Required metadata columns
The script requires the following columns:

- `SamplingDate`
- `IrrigationRegime`
- `VineyardFloorManagement`

If any of these columns is missing, the script stops with an explicit error.

### Required spectral columns

The spectral columns must use these prefixes:

- `FX10_nm_` for FX10 wavelengths
- `FX17_nm_` for FX17 wavelengths

Examples:

- `FX10_nm_720`
- `FX10_nm_720p5`
- `FX17_nm_1450`

The script parses these names, converts them into numeric wavelengths, and sorts the spectral matrix accordingly.

---

## Preprocessing and PCA logic

### RAW
The spectra are analysed as stored in the matrix.

### SNV
Standard normal variate (SNV) is applied **row-wise**, i.e. one transformation per sample spectrum.

### PCA computation
PCA is computed by **singular value decomposition (SVD)** on mean-centred data. This implementation is robust for high-dimensional spectral matrices where the number of wavelengths exceeds the number of samples.

### Quality-control steps
Before PCA, the script:

- removes rows containing non-finite spectral values;
- removes wavelengths with variance ≤ `1e-12`.

---

## Outputs

All outputs are written next to the input Excel file under:

```text
Objetivo_2/
```

### Main output files

The script creates run-specific filenames using a timestamp tag:

- `O2_log_YYYYMMDD_HHMMSS.txt`  
  run log and diagnostics

- `O2_PCA_Summary_YYYYMMDD_HHMMSS.xlsx`  
  Excel summary of PCA outputs

### Main figure outputs

#### Manuscript-oriented score plots
- `Fig3A_FX10_raw_byDate.png/.fig`
- `Fig3B_FX10_SNV_byIrrigation.png/.fig`
- `Fig3C_FX17_raw_byDate.png/.fig`
- `Fig3D_FX17_SNV_byIrrigation.png/.fig`

#### Loading plots
- `Loadings_FX10_RAW.png/.fig`
- `Loadings_FX10_SNV.png/.fig`
- `Loadings_FX17_RAW.png/.fig`
- `Loadings_FX17_SNV.png/.fig`

#### Mean and difference spectra plots
- `MeanDiffSpectra_FX10_RAW_byIrrigation.png/.fig`
- `MeanDiffSpectra_FX10_SNV_byIrrigation.png/.fig`
- `MeanDiffSpectra_FX17_RAW_byIrrigation.png/.fig`
- `MeanDiffSpectra_FX17_SNV_byIrrigation.png/.fig`

#### Supplementary diagnostic score plots
- `FigS_FX10_RAW_byIrrigation.png/.fig`
- `FigS_FX10_SNV_byDate.png/.fig`
- `FigS_FX17_RAW_byIrrigation.png/.fig`
- `FigS_FX17_SNV_byDate.png/.fig`

---

## Excel summary contents

The PCA summary workbook includes these sheets:

- `Explained`  
  explained variance for PC1–PC5

- `ANOVA_PC1_PC2`  
  one-way ANOVA results for PC1 and PC2 scores against:
  - `SamplingDate`
  - `IrrigationRegime`

- `TopLoadings`  
  top absolute loading wavelengths for PC1 and PC2

The current script does **not** create an additional `RunInfo` sheet.

---

## Software requirements

- **MATLAB**
- a MATLAB release supporting:
  - `detectImportOptions`
  - `readtable`
  - `svd`
  - `anova1`
  - `exportgraphics`
  - `savefig`

No external helper scripts are required beyond the main MATLAB file.

---

## How to run

Example:

```matlab
obj2_exploratory_pca('Matriz_CHEM_HSI_MASTER_96.xlsx', 'Matriz')
```

If the worksheet name is omitted or left empty, the script uses:

- `Matriz`

as the default sheet name.

---

## Main script options

The following options are defined near the top of the script:

- `NCOMP = 5`  
  number of principal components to compute

- `DO_LOADINGS = true`  
  export loading plots

- `DO_MEAN_DIFF = true`  
  export mean/difference spectra plots

- `DO_SUPP_DIAGNOSTICS = true`  
  export supplementary diagnostic score plots

Top loading wavelengths are exported for:

- **PC1**
- **PC2**

with a fixed maximum of **10 wavelengths per component**.

---

## Notes for manuscript alignment

This workflow corresponds to the **exploratory PCA stage** of the study. It provides the basis for the manuscript sections describing:

- clustering by sampling date,
- separation by irrigation regime,
- loading structure across FX10 and FX17,
- and the comparison between RAW and SNV spectral representations.

It therefore acts as the bridge between the sample-level spectral matrix and the manuscript-ready exploratory figures and summaries.

---

## Reuse

If you use or adapt this workflow, please cite the associated publication once available and acknowledge the original repository.
