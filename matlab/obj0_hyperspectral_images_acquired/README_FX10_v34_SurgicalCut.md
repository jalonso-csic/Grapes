# FX10 Berry Spectral Extraction — v34 (SAM + NDVI Surgical Cut)

## Overview
`Automatizacion_bayas_FX10_por_lote_v34_SURGICAL_PRECISION_compile.m` performs **batch berry-level spectral extraction** from **FX10 (VIS–NIR)** hyperspectral data.

Core logic (audit-ready):
1. Reflectance calibration using sample-specific WHITE/DARK references.
2. Strict NDVI “live” mask (no hole filling).
3. SAM segmentation (`imsegsam`) on the visual image.
4. “Surgical cut”: each SAM object is intersected with the NDVI mask to remove background.
5. Geometric filtering (min area + min solidity) + area trimming (P10–P99.5).
6. Mean spectra per berry ROI after 1-pixel erosion (edge purity).
7. Spectral QC selection: keep 90 closest to the global median, then 70 closest to the 90-median.
8. Per-sample exports + compilation across samples.

## Expected folder structure
Run from the experiment root:

```
ROOT/
  FX10_.../
    capture/
      <sample>.hdr/.raw
      WHITEREF_<sample>.hdr/.raw
      DARKREF_<sample>.hdr/.raw
    <sample>.png   (optional)
```

Folders starting with `FX10_` are auto-detected (excluding those containing “resultado”).

## Outputs
Created in `ROOT/`:

- `results_STEP1_2_FX10_v34_SURGICAL/`
  - STEP1 NDVI mask
  - STEP2 label maps + contours
  - `labels_SAM_<sample>.mat`

- `results_Spectra90_FX10_FullROI/`
  - `Spectra90_FullROI_<sample>.xlsx`

- `results_Spectra70_FX10_FullROI/`
  - `Spectra70_FullROI_<sample>.xlsx`
  - QC figures (`.fig` + `.png`)
  - `Spectra70_ALL_Samples_compiled.xlsx`
  - `Spectra70_MEANS_bySample.xlsx`

## Key parameters (as implemented)
- SAM: `PointGridSize = [32 32]`, `ScoreThreshold = 0.40`
- NDVI: NIR=800 nm, RED=670 nm, threshold=0.38
- Geometry: min area=150 px, min solidity=0.88
- Area trimming: P10–P99.5
- Selection: 90 → 70 by Euclidean distance to median spectra

## Dependencies
- MATLAB + Image Processing Toolbox
- Add-on: *Image Processing Toolbox Model for Segment Anything Model* (`imsegsam`)
- ENVI readers:
  - `read_envi_info.m` + `read_envi_cube.m`
  - or: `enviinfo_local.m` + `enviread_local.m`

## Run
```matlab
Automatizacion_bayas_FX10_por_lote_v34_SURGICAL_PRECISION_compile
```
