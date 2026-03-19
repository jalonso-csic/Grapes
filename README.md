# Grapes — MATLAB workflows for hyperspectral modelling, LODO robustness, VIP interpretation, and multispectral design

## Overview

This repository contains the MATLAB workflows used to support the manuscript:

**From hyperspectral modeling to multispectral design for grape maturity screening under date-wise validation**

**Zenodo DOI (v1.0.1): 10.5281/zenodo.19105661**

The code is organised as a sequence of reproducible objectives that follows the analytical logic of the study:

1. berry-level spectral extraction from hyperspectral images;
2. statistical analysis of chemical reference measurements;
3. exploratory spectral structure analysis;
4. supervised PLS-R modelling;
5. date-aware robustness analysis under strict leave-one-date-out evaluation with minimal daily recalibration;
6. VIP-based interpretability analysis;
7. translation from hyperspectral models to deployable multispectral filter-set proposals.

The repository is intended to document the **MATLAB-side analytical workflow** underlying the manuscript, from image-derived berry spectra to model simplification and multispectral design.

---

## Repository structure

```text
matlab/
  obj0_hsi_berry_extraction/
    fx10_batch_berry_extraction_surgical_precision.m
    fx17_batch_berry_extraction_hsv_sam_fullroi_compile.m
    read_envi_cube.m
    read_envi_info.m
    README.md

  obj1_chemical_responses/
    obj1_chemical_ground_truth_statistics_lmm.m
    README.md

  obj2_exploratory_pca/
    obj2_exploratory_pca.m
    README.md

  obj3_supervised_models/
    obj3_plsr_berry_composition.m
    README.md

  obj4_lodo_robustness/
    obj4_lodo_strict_kanchors.m
    README.md

  obj5_vip_interpretability/
    obj5_vip_interpretability.m
    README.md

  obj6_multispectral_design/
    obj6_multispectral_filterset_design.m
    README.md
```

Each objective folder contains:

- one principal MATLAB script;
- one folder-specific `README.md` describing purpose, inputs, outputs, and usage.

---

## Analytical workflow

### Objective 0 — HSI berry extraction
Scripts for extracting berry-level mean reflectance spectra from calibrated hyperspectral image cubes acquired with:

- **Specim FX10** (VIS–NIR)
- **Specim FX17** (SWIR)

This stage performs image-based berry segmentation, region selection, spectral extraction, and export of berry-level and sample-level spectral tables.

### Objective 1 — Chemical responses
Statistical analysis of chemical ground-truth variables using linear mixed models and manuscript-oriented summaries.

### Objective 2 — Exploratory PCA
Exploratory analysis of spectral structure across sensors and preprocessing strategies.

### Objective 3 — Supervised models
PLS-R modelling of berry-composition endpoints across sensor configurations and preprocessing pipelines.

### Objective 4 — LODO robustness
Strict date-wise validation and minimal recalibration analysis using **leave-one-date-out (LODO)** evaluation and **k-anchor** intercept updating.

### Objective 5 — VIP interpretability
Trait-specific VIP-profile computation for selected deployable endpoints, including Top-5 VIP windows and manuscript-oriented figures.

### Objective 6 — Multispectral design
Reduction of hyperspectral models to simulated multispectral filter sets using VIP-derived candidate bands, greedy subset selection, and repeated outer evaluation.

---

## Recommended execution order

The repository is modular, but the intended analytical order is:

```text
obj0 -> obj1 -> obj2 -> obj3 -> obj4 -> obj5 -> obj6
```

### Practical dependency logic

- `obj0` generates the spectral inputs used downstream.
- `obj1` analyses the chemical reference data.
- `obj2` explores spectral structure.
- `obj3` fits the supervised models.
- `obj4` stress-tests selected models under date-wise validation.
- `obj5` computes VIP profiles for selected traits.
- `obj6` uses Objective 5 outputs to derive multispectral proposals.

---

## Expected inputs

The repository combines two input layers:

### 1. Image-level inputs
Used in **Objective 0**, with sample folders containing ENVI files (`.hdr`, `.raw`) and sample-specific white/dark references.

### 2. Matrix-level inputs
Used in **Objectives 1–6**, especially:

- `Matriz_CHEM_HSI_MASTER_96.xlsx`

Depending on the objective, additional required worksheets or upstream files are described in the local `README.md` inside each folder.

### Additional upstream dependency
Objective 6 requires the VIP workbook generated in Objective 5:

- `O5_VIP_Profiles_SelectedTraits.xlsx`

---

## Software requirements

### MATLAB
The workflows were developed in MATLAB and require a standard desktop MATLAB installation.

### Toolboxes
Depending on the objective, the following MATLAB toolboxes are required:

- **Image Processing Toolbox**
- **Statistics and Machine Learning Toolbox**
- **Signal Processing Toolbox**

### Additional model support for Objective 0
Objective 0 also requires:

- **Image Processing Toolbox Model for Segment Anything Model** (`imsegsam`)

### ENVI readers
Objective 0 relies on the helper functions:

- `read_envi_info.m`
- `read_envi_cube.m`

These are included in `obj0_hsi_berry_extraction/`.

---

## How to use the repository

Each objective is designed to be run from its own folder, with the necessary input files available in the current working directory or in the expected local structure described in the corresponding `README.md`.

### General pattern

Open MATLAB in the relevant objective folder and run the main script, for example:

```matlab
obj3_plsr_berry_composition
```

or

```matlab
obj5_vip_interpretability
```

The exact input expectations, output folders, and special notes are documented in the folder-specific README of each objective.

---

## Outputs

Most scripts create their own output folders automatically. In several objectives, output folder names intentionally preserve the historical project naming used during manuscript development.

For example, some scripts write results under directories such as:

- `Objetivo_3`
- `Objetivo_4`
- `Objetivo_5`
- `Objetivo_6`

This is expected behaviour and reflects the original project workflow.

---

## Reproducibility notes

- The repository is structured around **one principal script per objective**.
- Folder-specific documentation is provided for all objectives.
- Objectives 3–6 are aligned with the manuscript’s modelling, robustness, interpretability, and multispectral-design sections.
- Objective 0 preserves the original operational extraction logic while exposing the workflow in a public-facing repository structure.

Because some workflows depend on project-specific file structures and naming conventions, the folder-level `README.md` files should be read before execution.

---

## Relationship to the manuscript

This repository supports the computational and analytical workflow behind the manuscript by linking:

- hyperspectral berry extraction,
- chemical ground-truth analysis,
- exploratory spectral analysis,
- supervised PLS-R modelling,
- date-aware robustness assessment,
- VIP-based interpretation,
- and multispectral simplification.

In that sense, the repository is not a generic hyperspectral toolbox, but a **study-specific, publication-oriented workflow package**.

---

## Suggested reading path

For readers interested in reproducing the study logic rather than executing the full pipeline immediately, the recommended order is:

1. read this main `README.md`;
2. open the `README.md` inside each objective folder;
3. inspect `obj3`, `obj4`, `obj5`, and `obj6` first if the main interest is modelling and deployment logic;
4. inspect `obj0` first if the main interest is image-to-spectrum extraction.

---

## Citation

If you use or adapt these workflows, please cite the software record:

**Alonso, J. (2026). _Grapes: MATLAB workflows for hyperspectral modelling, LODO robustness, VIP interpretability, and multispectral design_ (v1.0.1) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.19105661**

The associated publication citation can be added once formally available.

---

## Contact and maintenance

This repository was organised as a manuscript-oriented analytical package. Folder names, script names, and README files were cleaned for public-facing reproducibility while preserving the validated analytical logic used in the study.
