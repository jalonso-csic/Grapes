# Objective 6 — Multispectral filter-set design from VIP-derived candidate bands

## Overview

This folder contains the MATLAB workflow used to translate VIP-based interpretability into **deployable multispectral filter-set proposals**.

The script starts from the selected traits and their VIP profiles, derives candidate bands from VIP peaks, simulates multispectral band responses using an assumed optical bandwidth, and evaluates reduced band sets under a nested selection framework.

For each endpoint, the workflow:

1. reads the master Excel matrix;
2. reads the VIP-profile workbook generated in Objective 5;
3. infers whether the trait is modelled with **FX10**, **FX17**, or **FUSION**;
4. extracts a **block-aware VIP-peak candidate pool**;
5. simulates multispectral band responses using an assumed **20 nm FWHM** filter;
6. evaluates multispectral subsets across a range of band counts `K`;
7. performs greedy forward selection under repeated outer validation;
8. compares reduced multispectral subsets against the full-spectrum baseline;
9. exports candidate bands, performance tables, stability summaries, endpoint-specific K-curves, a manuscript-oriented Figure 7, and a run log.

This objective corresponds to the manuscript section on **translation from hyperspectral models to multispectral design**.

---

## Folder contents

- `obj6_multispectral_filterset_design.m`  
  Main MATLAB workflow for candidate-band generation, multispectral simulation, and subset evaluation.

- `README.md`  
  Documentation for this objective.

---

## Folder structure

```text
obj6_multispectrum/
  obj6_multispectral_filterset_design.m
  README.md
```

---

## What the script does

The workflow evaluates a fixed set of selected endpoints and, for each one:

- reads the VIP profile generated in Objective 5;
- detects local VIP maxima above a minimum threshold;
- builds a candidate band pool separately for each spectral block;
- simulates one feature per candidate band by averaging the spectral signal inside a `20 nm` bandpass window;
- compares multispectral subsets of different sizes `K`;
- performs greedy sequential forward selection inside the training loop;
- selects the number of latent variables by inner cross-validation;
- evaluates multispectral performance against the full-spectrum baseline;
- exports tables and figures for manuscript reporting.

---

## Input files

The script expects the following files in the working directory:

```text
Matriz_CHEM_HSI_MASTER_96.xlsx
O5_VIP_Profiles_SelectedTraits.xlsx
```

### Required worksheet in the matrix file

- `Matriz`

### Required sheets in the VIP workbook

The VIP workbook must contain one sheet per endpoint, with sheet names matching the endpoint variable names used in the script.

---

## Selected endpoints

The workflow evaluates the following endpoints:

- `TA_gL`
- `TSS_Brix`
- `TotalPhenolicPotential_mgkg`
- `TotalAnthocyanins_mgkg`
- `pH`
- `MalicAcid_gL`

These must match both:

- the endpoint names in the master matrix;
- and the VIP-profile sheet names in the Objective 5 workbook.

---

## Candidate-band generation

Candidate bands are derived directly from the VIP profiles.

### Rules used by the script

- only local maxima with `VIP >= 1.00` are considered;
- candidates are generated **per spectral block** (`FX10` or `FX17`);
- peaks closer than `10 nm` within the same block are merged;
- the number of retained candidates is capped at:
  - `25` per block

The candidate pool is therefore **block-aware**, which is especially important for FUSION traits.

---

## Spectral mode inference

The script infers the model mode directly from the VIP-profile sheet:

- only `FX10` bands present → `FX10`
- only `FX17` bands present → `FX17`
- both blocks present → `FUSION`

This determines how candidate bands are simulated and how the full-spectrum baseline is reconstructed.

---

## Multispectral simulation

Each candidate band is converted into a simulated multispectral response using an assumed optical filter with:

- `FWHM = 20 nm`

For a given candidate centre wavelength, the simulated band response is the mean reflectance over wavelengths inside:

```text
centre ± 10 nm
```

If no original wavelengths fall inside the target bandpass, the nearest available wavelength is used.

For FUSION traits, simulation remains block-aware and never mixes FX10 and FX17 candidates incorrectly.

---

## Full-spectrum baseline

For each endpoint, the workflow also computes a full-spectrum baseline using the spectral mode inferred from the VIP profile.

### Modes

- `FX10`
- `FX17`
- `FUSION`

### Fusion rule

For fusion traits:

- FX10 wavelengths **≤ 950 nm**
- FX17 wavelengths **> 950 nm**

This is consistent with the earlier modelling stages.

---

## Validation framework

### Outer evaluation
The default outer scheme is:

- `LODO`

with grouping variable:

- `SamplingDate`

The script also supports `KFold`, but the current configuration uses date-wise external validation.

### Repeats
- `nRepeats = 3`

### Inner selection
Inside the training folds, the workflow performs:

- greedy sequential forward selection for band-subset construction;
- latent-variable selection by inner cross-validation.

### Latent-variable limits
- multispectral models: up to `10` latent variables
- full-spectrum baseline: up to `20` latent variables

### One-standard-error rule
- enabled

---

## Number of bands evaluated

The script evaluates multispectral subsets for:

```text
K = 3:12
```

For each endpoint, this produces a K-curve relating the number of retained bands to multispectral performance.

---

## Outputs

All outputs are written under:

```text
Objetivo_6/
```

### Subfolders

- `Figures/`
- `Tables/`
- `Logs/`

### Main output files

#### Tables
- `Tables/O6_v2_CandidateBands_fromVIPPeaks.xlsx`
- `Tables/O6_v2_Multispectral_Performance_vsK.xlsx`
- `Tables/O6_v2_SelectedBandSets_Stability.xlsx`

#### Figures
- `Figures/O6_v2_<Endpoint>_<Mode>_Kcurve.fig`
- `Figures/O6_v2_<Endpoint>_<Mode>_Kcurve.png`

- `Figures/Figure7_O6_Retention_TwoBaselines.fig`
- `Figures/Figure7_O6_Retention_TwoBaselines.png`

#### Log
- `Logs/O6_v2_runlog.txt`

---

## Table contents

### Candidate-bands table
The candidate table contains, for each endpoint:

- global rank
- rank within block
- block
- peak wavelength
- peak VIP
- endpoint variable

### Multispectral-performance table
This table stores, for each endpoint and each `K`:

- endpoint
- mode
- number of bands
- full-spectrum baseline RMSE
- multispectral RMSE
- multispectral R²
- RMSE ratio relative to the full-spectrum RAW baseline

### Band-stability table
This table summarises how often candidate bands are selected across repeats and outer folds, allowing inspection of subset stability.

---

## Figures

### Endpoint-specific K-curves
For each endpoint, the script exports a K-curve showing:

- multispectral RMSE as a function of `K`;
- the full-spectrum baseline RMSE as a horizontal reference line.

### Figure 7
A summary figure is generated across endpoints, showing retention of multispectral performance relative to:

1. the full-spectrum RAW baseline;
2. the best full-spectrum LODO+k=3 baseline, when available in the configuration.

The figure annotates the selected `K` values for each endpoint.

---

## Software requirements

- **MATLAB**
- **Statistics and Machine Learning Toolbox**

The script uses functions such as:

- `plsregress`
- `cvpartition`
- `readtable`
- `writetable`
- `savefig`
- `exportgraphics`

No external helper scripts are required beyond the main MATLAB file.

---

## How to run

Open MATLAB in the objective folder and run:

```matlab
obj6_multispectral_filterset_design
```

The script is configured internally to read:

```text
Matriz_CHEM_HSI_MASTER_96.xlsx
O5_VIP_Profiles_SelectedTraits.xlsx
```

from the current working directory.

---

## Recommended public-facing renaming

To keep this objective consistent with the cleaned repository structure, the script should be renamed as:

- from `O6_Multispectral_FilterSet_Design_run_v3.m`
- to `obj6_multispectral_filterset_design.m`

The main function name should match the file name exactly:

```matlab
function obj6_multispectral_filterset_design()
```

This is a naming cleanup only. The analytical logic does not need to be changed.

---

## Notes for manuscript alignment

This workflow corresponds to the **multispectral design stage** of the study. It supports the manuscript section describing:

- candidate-band derivation from VIP peaks,
- reduction from hyperspectral to multispectral predictors,
- repeated nested selection of reduced band sets,
- retention of predictive performance with limited numbers of bands,
- and operational comparison against the full-spectrum baseline.

It therefore acts as the bridge between the interpretability stage and the final deployable filter-set proposals.

---

## Reuse

If you use or adapt this workflow, please cite the associated publication once available and acknowledge the original repository.
