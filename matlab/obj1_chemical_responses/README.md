# Objective 1 — Chemical ground-truth statistics for the grape dataset

## Overview

This folder contains the MATLAB workflow used to generate the **reference-statistics outputs** for the grape dataset, based on laboratory chemistry measured at the **experimental-unit level** (one row per experimental unit in the master Excel table).

The script:

1. reads the master chemistry table;
2. harmonises the experimental design factors;
3. fits **linear mixed-effects models (LMMs)** for each endpoint;
4. computes planned **Irrigation vs Rainfed** contrasts within each sampling date;
5. exports **manuscript-oriented Excel tables**, **publication-ready figures**, and a **run log**.

This objective corresponds to the **chemical ground-truth statistical analysis** used to support the manuscript results section on reference composition.

---

## Folder contents

- `obj1_chemical_ground_truth_statistics_lmm.m`  
  Main MATLAB workflow for chemical ground-truth statistical analysis.

- `README.md`  
  Documentation for this objective.

---

## Folder structure

```text
obj1_chemical_responses/
  obj1_chemical_ground_truth_statistics_lmm.m
  README.md
```

---

## What the script does

For each chemical endpoint, the workflow:

- resolves the relevant column in the input Excel file;
- standardises:
  - **Date**
  - **Irrigation**
  - **Management**
  - **Block**
- fits the following LMM structures:
  - **Global model:** `y ~ Date * Irrigation + (1|Block)`
  - **Irrigated subset:** `y ~ Date * Management + (1|Block)`
  - **Rainfed subset:** `y ~ Date * Management + (1|Block)`
- uses a fallback additive model when the interaction is not estimable;
- computes planned **Irrigation – Rainfed contrasts by date**;
- exports:
  - tidy ANOVA tables,
  - contrast tables,
  - means ± SE summaries,
  - manuscript-oriented p-value tables,
  - publication-quality figures,
  - and a complete run log.

---

## Expected input

### File type
- `.xlsx` master table

### Structure
- **one row per experimental unit**

### Typical required design columns
The script attempts to resolve these factors from candidate column names:

- **Sampling date**  
  e.g. `SamplingDate`, `Fecha`, `Fecha_ddmm`, `Sampling_Event_Date`

- **Irrigation regime**  
  e.g. `IrrigationRegime`, `Riego`, `Irrigation`, `WaterRegime`

- **Vineyard floor management**  
  e.g. `VineyardFloorManagement`, `Manejo_cultivo`, `Management`

- **Treatment**  
  e.g. `Treatment`, `Tratamiento`

- **Block**  
  e.g. `Block`, `BloqueRoman`, `Bloque`, `Replicate`

If any required factor is missing, the script stops with an explicit error.

---

## Endpoint coverage

The workflow is designed for the following chemical endpoints:

| Code | Endpoint | Unit |
|---|---|---|
| BW | Berry weight (50 berries) | g |
| Brix | Total soluble solids | °Brix |
| pH | pH | – |
| AT | Titratable acidity | g/L |
| Malic | Malic acid | g/L |
| Tart | Tartaric acid | g/L |
| K | Potassium | g/L |
| NOPA | YAN (alpha-amino N) | mg N/L |
| NH4 | Ammonium | mg N/L |
| PPT | Total phenolic potential | mg/kg |
| PPEx | Extractable phenolic potential | mg/kg |
| AnT | Total anthocyanins | mg/kg |
| AnEx | Extractable anthocyanins | mg/kg |

Column resolution is handled internally through candidate-name matching.

---

## Statistical design

### Global analysis
For each endpoint:

- `y ~ Date * Irrigation + (1|Block)`

### Within-irrigation management analyses
For irrigated and rainfed subsets separately:

- primary model: `y ~ Date * Management + (1|Block)`
- fallback model: `y ~ Date + Management + (1|Block)`

### Planned contrasts
The workflow computes:

- **Irrigation vs Rainfed within each sampling date**

### Inference
- ANOVA extracted from `fitlme`
- degrees of freedom based on **Satterthwaite** approximation

---

## Outputs

All outputs are written next to the input Excel file under:

```text
Objetivo_1/
  RESULTS_Obj1_YYYYMMDD_HHMMSS/
```

### Main output files

- `O1_log.txt`  
  Run log and diagnostics.

- `O1_Stats_Summary.xlsx`  
  Tidy ANOVA results, contrasts, and means.

- `O1_Table1_Manuscript.xlsx`  
  Manuscript-oriented p-value table.

- `Fig1_Irrigation_<CODE>.png/.fig`  
  Date × Irrigation means ± SE.

- `Fig2_Management_Irrigated_<CODE>.png/.fig`  
  Date × Management means ± SE within irrigated samples.

- `Fig3_Management_Rainfed_<CODE>.png/.fig`  
  Date × Management means ± SE within rainfed samples.

---

## Software requirements

- **MATLAB**
- **Statistics and Machine Learning Toolbox**
- a MATLAB release supporting:
  - `fitlme`
  - `anova`
  - `exportgraphics`
  - `savefig`

No additional helper scripts are required beyond the main MATLAB file.

---

## How to run

Example:

```matlab
obj1_chemical_ground_truth_statistics_lmm('Matriz_CHEM_HSI_MASTER_96.xlsx', 'Matriz')
```

If the worksheet name is omitted, the script attempts to use the first sheet.

---

## Notes for manuscript alignment

This workflow corresponds to the **reference-measurement statistical stage** of the study. It provides the inferential backbone for the manuscript section describing:

- temporal effects,
- irrigation effects,
- management effects within irrigation regime,
- and date-wise irrigation contrasts.

It therefore acts as the bridge between the raw chemistry table and the manuscript-ready statistical summaries.

---

## Reuse

If you use or adapt this workflow, please cite the associated publication once available and acknowledge the original repository.
