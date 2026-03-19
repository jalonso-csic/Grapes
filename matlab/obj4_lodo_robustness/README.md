# Objective 4 — LODO robustness with strict evaluation and k-anchor recalibration

## Overview

This folder contains the MATLAB workflow used to evaluate the robustness of the selected PLS-R models under **leave-one-date-out (LODO)** validation.

The script compares two deployment-relevant scenarios:

1. **Strict LODO**  
   no information from the held-out sampling date is used.

2. **LODO + k anchors per date**  
   a small number of anchor samples from the held-out date are assumed to be chemically assayed, and an intercept update is applied before evaluating the remaining samples from that date.

This objective corresponds to the manuscript section on **date-aware robustness and minimal daily recalibration**.

---

## Folder contents

- `obj4_lodo_strict_kanchors.m`  
  Main MATLAB workflow for LODO robustness analysis.

- `README.md`  
  Documentation for this objective.

---

## Folder structure

```text
obj4_interpretability/
  obj4_lodo_strict_kanchors.m
  README.md
```

---

## What the script does

The workflow:

1. reads the master Excel matrix and data dictionary;
2. identifies the grouping variable for LODO:
   - `SamplingEvent`, if available;
   - otherwise `SamplingDate`;
3. reconstructs the selected models from the manuscript Table 4:
   - endpoint,
   - sensor,
   - preprocessing,
   - latent variables,
   - random-CV baseline `R2cv`;
4. evaluates each selected model under:
   - **strict LODO**,
   - **LODO + k anchors/date**;
5. exports:
   - a manuscript-style compact performance table,
   - per-fold diagnostics,
   - long-format predictions with anchor flags,
   - a slope chart comparing strict vs anchor-adjusted performance,
   - and a run log.

---

## Input files

The script expects the following Excel file in the working directory:

```text
Matriz_CHEM_HSI_MASTER_96.xlsx
```

### Required worksheets

- `Matriz`
- `DataDictionary`

If the file or any required worksheet is missing, the script stops with an explicit error.

---

## LODO design

### Group definition
The held-out group is defined as:

- `SamplingEvent`, if present;
- otherwise `SamplingDate`.

Each fold leaves out one complete sampling event/date.

### Strict LODO
For a given held-out date:

- the model is trained on all other dates;
- predictions are generated for the held-out date;
- performance is evaluated on **all test samples** from that held-out date.

### LODO + k anchors/date
For the same held-out date:

- `k` anchor samples are selected using spectral diversity;
- those anchors are assumed to be chemically assayed on that date;
- an **intercept-only bias update** is computed:

```text
bias = mean(y_true_anchor - y_pred_anchor)
```

- adjusted predictions are generated:

```text
y_pred_adj = y_pred + bias
```

- performance is evaluated **only on the non-anchor test samples**.

This avoids leakage from the anchor set into the reported test performance.

---

## Anchor strategy

### Default setting
- `kAnchors = 3`

### Anchor-selection mode
The script supports:

- `KS`  
  Kennard–Stone spectral-diversity selection on the held-out test set

- `FirstK`  
  simple selection of the first `k` samples

The default is:

- `anchorMode = 'KS'`

---

## Selected models

The workflow does not search the full model space. Instead, it evaluates a **fixed set of selected models** corresponding to the manuscript Table 4.

For each endpoint, the script stores:

- endpoint variable name,
- selected sensor,
- selected preprocessing,
- selected latent variables,
- random-CV baseline `R2cv_random`.

This makes Objective 4 a **robustness-evaluation layer**, not a model-selection stage.

---

## Spectral blocks and fusion rule

### FX10
Uses spectral columns starting with:

- `FX10_nm_`

### FX17
Uses spectral columns starting with:

- `FX17_nm_`

### FUSION
The fused predictor block is constructed as:

- FX10 wavelengths **≤ 950 nm**
- FX17 wavelengths **> 950 nm**

This avoids overlap between the two sensors.

If enabled, fusion applies row-wise L2 normalisation to each sensor block before concatenation.

---

## Preprocessing options

The script supports the preprocessing labels stored in the selected-model table:

- `RAW`
- `SNV`
- `SG1`
- `SG2`
- `SNV+SG1`
- `SNV+SG2`

### Details

- **SNV** is applied row-wise.
- **SG1** and **SG2** are Savitzky–Golay first and second derivatives.
- Derivatives are computed using padded convolution.
- For fusion, preprocessing is applied separately to FX10 and FX17 before concatenation.

---

## Outputs

All outputs are written under:

```text
Objetivo_4/
```

### Main output files

- `LODO_Performance_SelectedModels_Strict_vs_k3Anchors.xlsx`  
  performance summary

- `LODO_Predictions_SelectedModels_Strict_vs_k3Anchors.xlsx`  
  long-format predictions

- `Figure6_LODO_Strict_vs_k3Anchors.fig`
- `Figure6_LODO_Strict_vs_k3Anchors.png`  
  slope chart comparing strict vs anchor-adjusted pooled R²

- `O4_LODO_k3_runlog.txt`  
  run log

If `kAnchors` is changed in the configuration, the output filenames change accordingly.

---

## Performance workbook contents

### Sheet `Table5`
Compact manuscript-style summary including:

- endpoint
- best selected configuration
- random-CV baseline `R2cv_random`
- `R2_LODO_strict`
- `DeltaR2_strict`
- `RMSE_LODO_strict`
- `R2_LODO_k`
- `DeltaR2_k`
- `RMSE_LODO_k`
- total number of valid samples

### Sheet `PerFold`
Per-held-out-date diagnostics including:

- held-out group
- fold size
- strict-L0DO metrics
- k-anchor metrics
- number of anchors
- estimated intercept-update bias

---

## Prediction workbook contents

### Sheet `Predictions`
Long-format sample-level predictions including:

- endpoint
- selected sensor
- selected preprocessing
- Table-4 latent variables
- held-out group
- sampling date
- irrigation regime
- vineyard floor management
- block
- true response
- strict prediction
- anchor-adjusted prediction
- strict residual
- anchor-adjusted residual
- `AnchorFlag`

The `AnchorFlag` field identifies samples used for daily recalibration.

---

## Figure

The script exports a slope chart showing, for each endpoint:

- pooled `R2` under **strict LODO**
- pooled `R2` under **LODO + k anchors**

Endpoints are sorted by anchor-adjusted pooled `R2` for readability.

---

## Software requirements

- **MATLAB**
- **Statistics and Machine Learning Toolbox**

The script uses functions such as:

- `readtable`
- `plsregress`
- `sgolay`
- `pdist`
- `squareform`
- `exportgraphics`
- `savefig`

No external helper scripts are required beyond the main MATLAB file.

---

## How to run

Open MATLAB in the objective folder and run:

```matlab
obj4_lodo_strict_kanchors
```

The script is configured internally to read:

```text
Matriz_CHEM_HSI_MASTER_96.xlsx
```

from the current working directory.

---

## Recommended public-facing renaming

To keep this objective consistent with the cleaned repository structure, the script should be renamed as:

- from `O4_LODO_Strict_and_k3_run_v2.m`
- to `obj4_lodo_strict_kanchors.m`

The main function name should match the file name exactly:

```matlab
function obj4_lodo_strict_kanchors()
```

This is a naming cleanup only. The analytical logic does not need to be changed.

---

## Notes for manuscript alignment

This workflow corresponds to the **date-aware robustness evaluation stage** of the study. It supports the manuscript section describing:

- loss of performance under strict date-wise external validation,
- partial recovery through minimal daily recalibration,
- and the operational meaning of anchor-based bias correction.

It therefore acts as the bridge between the selected models from the supervised stage and the manuscript-ready robustness outputs.

---

## Reuse

If you use or adapt this workflow, please cite the associated publication once available and acknowledge the original repository.
