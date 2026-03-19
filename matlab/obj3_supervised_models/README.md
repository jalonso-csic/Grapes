# Objective 3 — Supervised PLS-R models for berry composition

## Overview

This folder contains the MATLAB workflow used to train and evaluate **supervised partial least squares regression (PLS-R)** models for grape berry composition from hyperspectral data.

The script:

1. reads the master Excel matrix and data dictionary;
2. extracts chemical endpoints from the dictionary;
3. builds spectral predictor sets for **FX10**, **FX17**, and **FUSION**;
4. applies a preprocessing grid including RAW, SNV, SG derivatives, and combined SNV+SG variants;
5. evaluates models using **repeated nested cross-validation**;
6. applies an internal outlier-control step within training folds only;
7. exports performance tables, predictions, publication-style figures, and a run log.

This objective corresponds to the **supervised modelling stage** used to support the manuscript section on prediction of berry composition by PLS-R models.

---

## Folder contents

- `obj3_plsr_berry_composition.m`  
  Main MATLAB workflow for supervised PLS-R modelling.

- `README.md`  
  Documentation for this objective.

---

## Folder structure

```text
obj3_supervised_models/
  obj3_plsr_berry_composition.m
  README.md
```

---

## What the script does

For each chemical endpoint listed in the data dictionary under **Category = "Chem endpoint"**, the workflow:

- extracts the response vector from the master matrix;
- builds spectral predictors for:
  - **FX10**
  - **FX17**
  - **FUSION** (FX10 + FX17 concatenation with a cutoff at 950 nm)
- applies the following preprocessing options:
  - `RAW`
  - `SNV`
  - `SG1`
  - `SG2`
  - `SNV_SG1`
  - `SNV_SG2`
- performs repeated nested cross-validation with:
  - outer CV for performance estimation,
  - inner CV for latent-variable selection,
  - one-standard-error LV selection rule;
- performs outlier detection **within training folds only**;
- trains the final fold-specific PLS-R models and predicts held-out samples;
- exports:
  - model performance summary,
  - sample-level predictions,
  - prediction vs reference figures for stronger models,
  - and a full run log.

---

## Input files

The script expects the following Excel file in the working directory:

```text
Matriz_CHEM_HSI_MASTER_96.xlsx
```

### Required worksheets

- `Matriz`  
  main data matrix

- `DataDictionary`  
  metadata sheet used to identify valid chemical endpoints

If either the file or the required sheets are missing, the script stops with an explicit error.

---

## Endpoint selection

Chemical endpoints are identified automatically from the `DataDictionary` sheet using:

- `Category == "Chem endpoint"`

Only endpoints present both in the dictionary and in the main matrix are modelled.

Responses with fewer than **20 valid samples** are skipped.

---

## Spectral predictors

### FX10
The script uses all columns whose names start with:

- `FX10_nm_`

### FX17
The script uses all columns whose names start with:

- `FX17_nm_`

### FUSION
The fusion block is constructed by concatenating:

- FX10 wavelengths **≤ 950 nm**
- FX17 wavelengths **> 950 nm**

This avoids overlap between the two sensors.

Wavelengths are parsed from column names and sorted numerically.

---

## Preprocessing grid

The supervised workflow evaluates the following preprocessing variants:

- `RAW`
- `SNV`
- `SG1`
- `SG2`
- `SNV_SG1`
- `SNV_SG2`

### Details

- **SNV** is applied row-wise.
- **SG1** and **SG2** correspond to first- and second-derivative Savitzky–Golay preprocessing.
- Derivatives are computed using padded convolution to reduce edge artefacts.
- In the fusion pipeline:
  - `RAW` uses concatenated FX10+FX17 spectra followed by row-wise L2 normalisation;
  - derivative-based options are applied separately to each block before concatenation.

---

## Cross-validation design

The script uses a **repeated nested CV** design:

### Outer CV
- `outerK = 10`
- `outerR = 5`

### Inner CV
- `innerK = 10`

### Latent-variable range
- forced search over `1:12`

### LV selection rule
The number of latent variables is selected using a **strict one-standard-error rule** based on inner-CV MSE.

---

## Outlier handling

Outlier control is applied **only within the training partition** of each outer fold.

### Settings
- sigma threshold: `2.5`
- maximum removal fraction: `10%`

The procedure:

1. generates fold-internal pilot predictions,
2. computes robust residual dispersion,
3. flags samples above the threshold,
4. caps removals at the allowed maximum percentage.

This design avoids leakage from the external test folds.

---

## Outputs

All outputs are written under:

```text
Objetivo_3/
  O3_PLSR_v22_YYYYMMDD_HHMMSS/
```

### Main output files

- `Polished_Performance.xlsx`  
  model-level performance summary

- `Polished_Predictions.xlsx`  
  sample-level predictions

- `run_log.txt`  
  full run log

- `Figures/`  
  prediction vs reference figures for models with `R2cv > 0.60`

---

## Performance metrics

The performance table includes:

- `R2cv`
- `RMSEcv`
- `RPD`
- `Bias`
- `LV_Mean`
- `LV_Median`
- `AvgOutliersRemoved`

These are computed from the repeated outer-CV predictions.

---

## Prediction outputs

The predictions table includes:

- endpoint
- sensor
- preprocessing
- original row index
- reference value
- cross-validated prediction

This allows downstream comparison, filtering, and figure generation.

---

## Figures

Prediction vs reference plots are exported only for models with:

- `R2cv > 0.60`

The plotting function forces:

- white background,
- black text,
- publication-safe export independent of the MATLAB UI theme.

Each plot is saved as:

- `.fig`
- `.png`

inside the `Figures/` folder.

---

## Software requirements

- **MATLAB**
- **Statistics and Machine Learning Toolbox**

The script uses functions such as:

- `plsregress`
- `cvpartition`
- `detectImportOptions`
- `readtable`
- `exportgraphics`
- `savefig`

No external helper scripts are required beyond the main MATLAB file.

---

## How to run

Open MATLAB in the objective folder and run:

```matlab
obj3_plsr_berry_composition
```

The script is configured internally to read:

```text
Matriz_CHEM_HSI_MASTER_96.xlsx
```

from the current working directory.

---

## Recommended public-facing renaming

To keep this objective consistent with the cleaned repository structure, the script should be renamed as:

- from `O3_PLSR_BerryComposition_run_v22.m`
- to `obj3_plsr_berry_composition.m`

The main function name should match the file name exactly:

```matlab
function obj3_plsr_berry_composition()
```

This is a naming cleanup only. The modelling logic does not need to be changed.

---

## Notes for manuscript alignment

This workflow corresponds to the **supervised modelling stage** of the study. It supports the manuscript section describing:

- prediction of berry composition from FX10, FX17, and fused spectra,
- the effect of preprocessing,
- the robustness of nested CV performance estimates,
- and the comparative performance of the different sensor/preprocessing configurations.

It therefore acts as the bridge between the curated spectral matrix and the manuscript-ready supervised-model outputs.

---

## Reuse

If you use or adapt this workflow, please cite the associated publication once available and acknowledge the original repository.
