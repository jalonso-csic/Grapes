# Objective 4 — VIP stability mapping (SNV + SG 2nd derivative)

This script computes **VIP stability** across the repeated nested-CV PLSR fits from **Objective 3** and exports:

- **Figure 1**: VIP stability heatmap for **Tier 1 + selected Tier 2 + selected Tier 3S** endpoints  
  (Tier 1 displayed at the top, then Tier 2, then Tier 3S).
- **Figure 2**: VIP stability line profiles for a selected subset of Tier 3S endpoints.
- **Excel workbook**: stable wavelengths and stable contiguous regions (operational band selection output).

The implementation is **audit-friendly**: stability is computed from the raw VIP matrices saved per outer model fit.

---

## Script

- `BRC07_OBJ4_VIPStability_SelectedTiers_FigsAndBands_SNV_SG2nd_v1_public.m`

Recommended repository location:

- `matlab/obj4_interpretability/drivers/`

(Any shared helpers should live in `matlab/common/`, but this script is self-contained.)

---

## Inputs (required)

You must have run the **Objective 3** modelling scripts **with VIP raw saving enabled**, so that the following files exist (names can vary; the script searches by pattern):

### VIP raw MAT files (SNV + SG 2nd derivative)

- Tier 1: `*Tier1*SNV*SG2nd*.mat` containing `VIP_STORE`
- Tier 2: `*Tier2*SNV*SG2nd*.mat` containing `VIP_STORE`
- Tier 3S: `*Tier3S*SNV*SG2nd*.mat` containing `VIP_STORE`

Each `VIP_STORE.(endpoint)` must include at least:
- `vip`: VIP values for all outer model fits (either `[p × nModels]` or `[nModels × p]`)
- `wl`: wavelength vector (`p × 1`)

### Summary XLSX files (optional but recommended; used for ordering)

These are used **only** to order endpoints within each tier (descending by `R2_pooled`, by default):

- Tier 1 summary: `*OBJ3_Summary_Tier1*SG2nd*.xlsx`
- Tier 2 summary: `*OBJ3_Summary_Tier2*SG2nd*.xlsx`
- Tier 3S summary: `*OBJ3S_Summary_Tier3S*SG2nd*.xlsx`

If a summary file is missing, the script keeps the original endpoint order for that tier.

---

## What is “VIP stability” here?

For each model fit (outer fold across repeats), VIP is normalised so that the **mean VIP across wavelengths equals 1**:

- `VIP_norm(model, λ) = VIP(model, λ) / mean_λ(VIP(model, λ))`

Then:

- `Stability(λ) = proportion of models with VIP_norm(model, λ) > 1`

This yields a stability value in `[0, 1]` for each wavelength.

---

## Stable regions

A wavelength is considered **stable** when:

- `Stability(λ) ≥ STAB_THR`

A **stable region** is a contiguous run of stable wavelengths with:

- `n_bands ≥ MIN_BANDS`

(Assuming a 1 nm wavelength grid.)

---

## Outputs

By default, outputs are written to:

- `OBJ4_OUT/` (relative to the current working directory)

### Figures (MATLAB + PNG)

- `Fig_OBJ4_VIPStability_Combined_T1_T2sel_T3sel_SNV_SG2nd.fig`
- `Fig_OBJ4_VIPStability_Combined_T1_T2sel_T3sel_SNV_SG2nd.png` (400 dpi)
- `Fig_OBJ4_VIPStability_Tier3SelLines_SNV_SG2nd.fig`
- `Fig_OBJ4_VIPStability_Tier3SelLines_SNV_SG2nd.png` (400 dpi)

### Excel (operational band output)

- `OBJ4_StableBands_SelectedEndpoints_SNV_SG2nd.xlsx`
  - `StableRegions`: start/end nm, region width, stability statistics
  - `StableBands_nm`: per-wavelength stable list
  - `Counts`: per-endpoint summaries (n stable nm, n regions, max stability)

---

## How to run

1. Ensure Objective 3 scripts have produced the VIP raw MAT files for **SNV+SG2nd**.
2. Open MATLAB and set your working directory (e.g., `matlab/obj4_interpretability/drivers/`).
3. (Recommended) Set `OBJ3DIR_MANUAL` at the top of the script to the Objective 3 output folder.
4. Run the script.

If the Objective 3 output folder cannot be auto-detected, a folder picker will open.

---

## Customisation points (top of script)

- `Tier1_All`, `Tier2_Sel`, `Tier3_Sel`: endpoints to include
- `SORT_METRIC_COL`: metric for within-tier ordering (`R2_pooled` default)
- `STAB_THR`, `MIN_BANDS`: stability threshold and minimum region length
- Figure font/size and `PNG_DPI`

---

## Notes for publication/reproducibility

- Endpoint matching is robust to naming differences such as `.` vs `_` using `matlab.lang.makeValidName`.
- This script does not refit models; it only post-processes VIP matrices saved from Objective 3.
