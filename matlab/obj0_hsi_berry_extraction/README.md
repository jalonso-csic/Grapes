# Objective 0 — Berry-level spectral extraction from hyperspectral grape images

## Overview

This folder contains the MATLAB workflows used to extract berry-level reflectance spectra from calibrated hyperspectral image cubes acquired with two Specim cameras:

- **FX10 (VIS–NIR)**: segmentation based on **SAM + strict NDVI live-mask trimming**.
- **FX17 (SWIR)**: segmentation based on **HSV pre-masking + SAM refinement**, followed by full-berry ROI extraction.

For each sample, the scripts:

1. calibrate raw hyperspectral cubes to reflectance using sample-specific **WHITE/DARK** references;
2. detect berry regions from hyperspectral-derived visual representations;
3. extract one mean reflectance spectrum per berry ROI;
4. apply a spectral quality-control selection to retain representative berry spectra;
5. export per-sample Excel tables and quality-control figures;
6. compile berry-level spectra and per-sample mean spectra across all processed samples.

These scripts generated the spectral inputs used in the downstream exploratory, supervised, interpretability, and multispectral analyses of the manuscript.

---

## Folder contents

- `fx10_batch_berry_extraction_surgical_precision.m`  
  Batch processing pipeline for **FX10** hyperspectral images.

- `fx17_batch_berry_extraction_hsv_sam_fullroi_compile.m`  
  Batch processing pipeline for **FX17** hyperspectral images.

- `read_envi_info.m`  
  Helper function to read ENVI metadata.

- `read_envi_cube.m`  
  Helper function to read ENVI hyperspectral cubes.

- `FX10_*` / `FX17_*` sample folders  
  Optional local sample folders used as the execution root during testing.

---

## Execution model

This folder is designed as a **self-contained execution root**.

The scripts use the current MATLAB working directory (`pwd`) as the root directory. Therefore, when running either script, MATLAB should be positioned inside this folder (or inside another folder with the same structure).

In practical terms, the scripts expect:

- the MATLAB scripts,
- the ENVI helper functions,
- and the sample folders (`FX10_*` or `FX17_*`)

to be available in the **same root directory**.

---

## Expected input structure

### FX10

```text
ROOT/
  FX10_.../
    capture/
      <sample>.hdr
      <sample>.raw
      WHITEREF_<sample>.hdr
      WHITEREF_<sample>.raw
      DARKREF_<sample>.hdr
      DARKREF_<sample>.raw
    <sample>.png   (optional; if absent, a pseudo-RGB image is generated)
```

### FX17

```text
ROOT/
  FX17_.../
    capture/
      <sample>.hdr
      <sample>.raw
      WHITEREF_<sample>.hdr
      WHITEREF_<sample>.raw
      DARKREF_<sample>.hdr
      DARKREF_<sample>.raw
```

The FX17 workflow also supports the alternative reference naming convention:

- `WHITEREF_FX17_<sample>.*`
- `DARKREF_FX17_<sample>.*`

---

## FX10 workflow

### Script
`fx10_batch_berry_extraction_surgical_precision.m`

### Processing logic

The FX10 pipeline implements the following sequence:

1. Reflectance calibration using sample-specific WHITE and DARK references.
2. Construction of a strict **NDVI live mask** (`NDVI > 0.38`) without hole filling.
3. **SAM** segmentation on the visual image.
4. “**Surgical cut**” of each SAM object by intersection with the NDVI live mask.
5. Object filtering using:
   - minimum area,
   - solidity threshold,
   - percentile-based size trimming.
6. Mean spectral extraction from each berry after **1-pixel erosion**.
7. Spectral quality-control selection:
   - retain the 90 spectra closest to the global median;
   - retain the 70 spectra closest to the median of the selected 90.
8. Export of per-sample spectral tables and QC figures.
9. Compilation across samples into:
   - a berry-level compiled dataset;
   - a per-sample mean spectral dataset.

### Main outputs

The script creates the following directories in `ROOT/`:

- `resultados_PASO_1_2_v34_SURGICAL/`
- `resultados_PASO_3_ROIs_full_90/`
- `resultados_PASO_3_ROIs_full_70/`

Typical outputs include:

- NDVI live-mask figures
- final berry label maps and contour overlays
- `.mat` label files
- per-sample Excel tables for the selected 90 and 70 berry spectra
- QC figures for selected ROIs and spectra
- compiled berry-level and sample-level Excel tables

### Notes

A deterministic tie-breaking rule was introduced when NDVI-based trimming splits one SAM object into multiple connected components of equal maximum area. In such cases, the retained component is the one whose centroid is closest to the centroid of the original SAM object. This prevents ambiguity during component selection and improves run-to-run reproducibility.

---

## FX17 workflow

### Script
`fx17_batch_berry_extraction_hsv_sam_fullroi_compile.m`

### Processing logic

The FX17 pipeline implements the following sequence:

1. Reflectance calibration using sample-specific WHITE and DARK references.
2. Generation of a pseudo-RGB image from the calibrated SWIR cube.
3. Construction of an **HSV pre-mask** to suppress rachis and background.
4. **SAM** segmentation on the pseudo-RGB image.
5. Restriction of SAM labels to the HSV berry mask.
6. Watershed-based splitting of oversized merged regions.
7. Selection of the **90 largest berries by area**.
8. ROI extraction using the **full berry mask eroded by 2 pixels** (fallback to the full berry mask if the eroded ROI becomes too small).
9. Spectral quality-control selection:
   - retain the 70 spectra closest to the median spectrum.
10. Export of per-sample spectral tables and QC figures.
11. Compilation across samples into:
   - a berry-level compiled dataset;
   - a per-sample mean spectral dataset.

### Main outputs

The script creates the following directories in `ROOT/`:

- `resultados_PASO_1_SAM_AntiRaquis/`
- `resultados_PASO_3_ROIs_centrales_90/`
- `resultados_PASO_3_ROIs_centrales_70/`

Typical outputs include:

- HSV base-mask figures
- raw SAM label maps
- refined berry label maps
- contour overlays on pseudo-RGB images
- `.mat` label files
- per-sample Excel tables for the selected 90 and 70 berry spectra
- QC figures for selected ROIs and spectra
- compiled berry-level and sample-level Excel tables

---

## Metadata handling

Both scripts assign sample metadata by parsing the sample folder name and matching the final numeric identifier to an internally defined design table.

The exported metadata fields may include:

- `Camara`
- `Variedad`
- `Fecha`
- `SampleIndex`
- `Tratamiento`
- `Bloque`
- `Riego`
- `Manejo_cultivo`
- `SampleName`

At present, the experimental design table is embedded directly in each script and assumes a fixed mapping for 24 sample indices.

---

## Requirements

These workflows require:

- **MATLAB**
- **Image Processing Toolbox**
- **Image Processing Toolbox Model for Segment Anything Model** (`imsegsam`)
- ENVI-reading utilities available in the same folder:
  - `read_envi_info.m`
  - `read_envi_cube.m`

---

## How to run

Open MATLAB in this folder and run:

### FX10

```matlab
fx10_batch_berry_extraction_surgical_precision
```

### FX17

```matlab
fx17_batch_berry_extraction_hsv_sam_fullroi_compile
```

---

## Scope and limitations

- These workflows were developed for the specific folder structure and sample naming conventions used in the present grape experiment.
- The design table is currently hard-coded in both scripts.
- The scripts assume that the ENVI helper functions are available in the same execution directory.
- The FX10 workflow optionally uses a pre-generated PNG image for visual segmentation; if no PNG is found, a pseudo-RGB image is generated directly from the calibrated cube.
- The FX17 workflow uses an HSV pre-mask and watershed refinement to improve berry separation in SWIR pseudo-RGB space.
- These scripts are intended for **batch spectral extraction**, not for real-time field deployment.

---

## Relation to the manuscript

This folder corresponds to the initial hyperspectral processing stage of the study. The outputs generated here were used as the spectral input matrices for the subsequent objectives:

- chemical response analysis,
- exploratory spectral analysis,
- supervised modelling,
- interpretability analysis,
- VIP-based band extraction,
- and multispectral design evaluation.

Accordingly, this objective provides the bridge between raw hyperspectral image cubes and the sample-level spectral datasets used throughout the manuscript.

---

## Citation and reuse

If you use or adapt these scripts, please cite the associated publication once available and acknowledge the original repository.
