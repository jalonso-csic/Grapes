# FX17 Berry Spectral Extraction — v4 (HSV Anti-Rachis + SAM + Central ROIs)

## Overview
`Automatizacion_bayas_por_lote_v4_meta70_compile.m` performs **batch berry-level spectral extraction** from **FX17 (SWIR)** hyperspectral data.

Core logic (audit-ready):
1. Reflectance calibration using sample-specific WHITE/DARK references.
2. HSV pre-mask (“anti-rachis”) to suppress rachis/background.
3. SAM segmentation on pseudo-RGB; labels are restricted to the HSV mask.
4. Watershed splitting of oversized blobs; berries are re-labelled.
5. Select the 90 largest berries by area; ROI = full berry eroded by 2 px (fallback to full berry).
6. Spectral QC selection: keep 70 closest to the median spectrum (Euclidean distance).
7. Per-sample exports + compilation across samples.

## Expected folder structure
Run from the experiment root:

```
ROOT/
  FX17_.../
    capture/
      <sample>.hdr/.raw
      WHITEREF_<sample>.hdr/.raw  (or WHITEREF_FX17_<sample>.*)
      DARKREF_<sample>.hdr/.raw   (or DARKREF_FX17_<sample>.*)
```

Folders starting with `FX17_` are auto-detected (excluding those containing “resultado”).

## Outputs
Created in `ROOT/`:

- `results_STEP1_2_FX17_AntiRachis_SAM/`
  - STEP1 HSV pre-mask
  - STEP2 label maps
  - `labels_SAM_<sample>.mat`

- `results_Spectra90_FX17_CentralROIs/`
  - `Spectra90_CentralROIs_<sample>.xlsx`

- `results_Spectra70_FX17_CentralROIs/`
  - `Spectra70_CentralROIs_<sample>.xlsx`
  - QC figures (`.fig` + `.png`)
  - `Spectra70_ALL_Samples_compiled.xlsx`
  - `Spectra70_MEANS_bySample.xlsx`

## Key parameters (as implemented)
- HSV: `HUE_RANGE=[0.55 0.8]`, `MIN_SATURATION_THRESHOLD=0.35`, `MAX_VALUE_THRESHOLD=0.85`, `MIN_BERRY_AREA=50`
- ROI: erosion = 2 px
- Blob refinement: `MIN_AREA_BERRY=0.25*medianArea`, `MAX_AREA_SIMPLE=1.8*medianArea`
- Selection: top 90 by area → keep 70 by spectral distance

## Dependencies
- MATLAB + Image Processing Toolbox
- Add-on: *Image Processing Toolbox Model for Segment Anything Model* (`imsegsam`)
- ENVI readers:
  - `read_envi_info.m` + `read_envi_cube.m`
  - or: `enviinfo_local.m` + `enviread_local.m`

## Run
```matlab
Automatizacion_bayas_por_lote_v4_meta70_compile
```
