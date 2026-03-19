# Objective 5 — VIP interpretability for selected deployable traits

## Overview

This folder contains the MATLAB workflow used to compute and visualise **variable importance in projection (VIP)** profiles for a selected set of grape-composition endpoints.

The script focuses on traits considered **deployable** within the study, i.e. traits showing good baseline predictive performance and acceptable robustness under the date-aware validation framework. For each selected trait, the workflow:

1. reads the master Excel matrix and data dictionary;
2. reconstructs the selected modelling configuration:
   - endpoint,
   - sensor,
   - preprocessing,
   - latent variables;
3. rebuilds the spectral predictor block for **FX10**, **FX17**, or **FUSION**;
4. applies the same preprocessing logic used in the supervised pipeline;
5. fits the final PLS-R model with the fixed selected number of latent variables;
6. computes VIP scores;
7. extracts the **Top-5 VIP windows** above the threshold;
8. exports per-trait VIP profiles, manuscript-oriented tables, per-trait figures, a multi-panel summary figure, and a run log.

This objective corresponds to the manuscript section on **VIP-based interpretability and spectral regions of importance**.

---

## Folder contents

- `obj5_vip_interpretability.m`  
  Main MATLAB workflow for VIP-profile computation and figure generation.

- `README.md`  
  Documentation for this objective.

---

## Folder structure

```text
obj5_vip/
  obj5_vip_interpretability.m
  README.md
```

---

## What the script does

The workflow evaluates a fixed list of selected traits and associated best-model settings. For each trait, it:

- extracts the response vector from the master matrix;
- rebuilds the appropriate spectral block for:
  - **FX10**
  - **FX17**
  - **FUSION**
- applies the requested preprocessing:
  - `RAW`
  - `SNV`
  - `SG1`
  - `SG2`
  - `SNV+SG1`
  - `SNV+SG2`
- standardises the predictor matrix column-wise before model fitting;
- fits a PLS-R model with the selected number of latent variables;
- computes VIP scores from:
  - PLS weights,
  - score matrix,
  - and the response variance captured by each latent variable;
- exports:
  - one VIP profile sheet per trait,
  - a consolidated Top-5 VIP windows table,
  - one figure per trait,
  - a six-panel Figure 6,
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

## Selected traits

The workflow does not search traits automatically. Instead, it uses a **fixed trait list** and fixed modelling settings defined inside the script.

### Selected endpoints

- `TA_gL`
- `TSS_Brix`
- `TotalPhenolicPotential_mgkg`
- `MalicAcid_gL`
- `pH`
- `TotalAnthocyanins_mgkg`

### Stored configuration per trait

For each selected endpoint, the script specifies:

- selected sensor
- selected preprocessing
- selected number of latent variables

This makes Objective 5 an **interpretability layer** built on top of the selected supervised models.

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

If enabled, the fusion workflow also applies row-wise L2 normalisation to each sensor block before concatenation.

---

## Preprocessing

The script reproduces the same preprocessing families used in the modelling workflow:

- `RAW`
- `SNV`
- `SG1`
- `SG2`
- `SNV+SG1`
- `SNV+SG2`

### Details

- **SNV** is applied row-wise.
- **SG1** and **SG2** correspond to Savitzky–Golay first and second derivatives.
- Derivatives are computed using padded convolution.
- For fusion, preprocessing is applied separately to FX10 and FX17 before concatenation.

---

## VIP computation

VIP scores are computed from the final fitted PLS-R model using:

- PLS weights,
- latent-variable score matrix,
- and the amount of response variance explained by each latent variable.

This implementation is trait-specific and uses the selected latent-variable count defined in the script.

Before PLS fitting, predictors are standardised column-wise using a local z-score transformation.

---

## VIP windows

The script identifies VIP windows using the following logic:

- VIP threshold: `VIP >= 1.0`
- contiguous spectral regions above the threshold are treated as windows
- windows are ranked by integrated VIP area
- the **Top-5** windows are retained per trait
- windows contributing to the cumulative **75%** of total VIP area are flagged as **Primary**

Single-point windows are skipped to avoid unstable area calculations.

---

## Outputs

All outputs are written under:

```text
Objetivo_5/
```

### Main output folders

- `figures/`
- `tables/`

### Main output files

#### Tables
- `tables/O5_VIP_Profiles_SelectedTraits.xlsx`  
  one sheet per trait with:
  - wavelength,
  - VIP score,
  - spectral block

- `tables/Table6_VIP_Windows_Top5.xlsx`  
  consolidated Top-5 VIP windows table

#### Figures
- `figures/Figure6_VIP_Profiles_SelectedTraits.fig`
- `figures/Figure6_VIP_Profiles_SelectedTraits.png`

- `figures/VIP_<EndpointVar>.fig`
- `figures/VIP_<EndpointVar>.png`

#### Log
- `O5_VIP_runlog.txt`

---

## Table contents

### VIP profiles workbook
Each sheet contains:

- `wavelength_nm`
- `VIP`
- `Block`

### Top-5 windows table
The consolidated table includes:

- trait label
- endpoint variable
- best model description
- rank
- block
- start wavelength
- end wavelength
- window width
- VIP area
- VIP mean
- VIP max
- `Primary` flag

---

## Figures

### Per-trait plots
Each trait is exported as an individual VIP curve showing:

- VIP profile,
- threshold line at `VIP = 1`,
- shaded VIP windows,
- darker shading for primary windows.

### Multi-panel Figure 6
A six-panel summary figure is generated using the selected traits. It includes:

- one panel per trait,
- panel letters,
- white-background export,
- manuscript-oriented formatting.

---

## Software requirements

- **MATLAB**
- **Statistics and Machine Learning Toolbox**
- **Signal Processing Toolbox**

The script uses functions such as:

- `plsregress`
- `sgolay`
- `conv2`
- `readtable`
- `exportgraphics`
- `savefig`

No external helper scripts are required beyond the main MATLAB file.

---

## How to run

Open MATLAB in the objective folder and run:

```matlab
obj5_vip_interpretability
```

The script is configured internally to read:

```text
Matriz_CHEM_HSI_MASTER_96.xlsx
```

from the current working directory.

---

## Recommended public-facing renaming

To keep this objective consistent with the cleaned repository structure, the script should be renamed as:

- from `O5_VIP_Interpretability_LODOk3_SelectedTraits_run_v11.m`
- to `obj5_vip_interpretability.m`

The main function name should match the file name exactly:

```matlab
function obj5_vip_interpretability()
```

This is a naming cleanup only. The analytical logic does not need to be changed.

---

## Notes for manuscript alignment

This workflow corresponds to the **VIP interpretability stage** of the study. It supports the manuscript section describing:

- trait-specific VIP profiles,
- dominant spectral windows,
- primary vs secondary influential regions,
- and the translation from selected predictive models to interpretable spectral regions.

It therefore acts as the bridge between the selected supervised models and the manuscript-ready interpretability outputs.

---

## Reuse

If you use or adapt this workflow, please cite the associated publication once available and acknowledge the original repository.
