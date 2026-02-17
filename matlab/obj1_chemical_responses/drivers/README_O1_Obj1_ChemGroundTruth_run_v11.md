# O1_Obj1_ChemGroundTruth_run_v11

Objective 1 (manuscript Section 3.1): **chemical ground-truth statistics** for the grape dataset (chemistry measured per experimental unit), producing manuscript-ready tables and publication-quality figures.

---

## What the script does

For each chemical endpoint, the script:

1. Reads the master chemistry Excel table (**one row per experimental unit**).
2. Harmonises design factors:
   - **Date** (ordered categorical; dd/MM)
   - **Irrigation** (Irrigation vs Rainfed)
   - **Management** (vineyard floor management)
   - **Block** (random effect)
3. Fits linear mixed-effects models (LMMs):
   - **Global model:** `y ~ Date * Irrigation + (1|Block)`
   - **Irrigated subset:** `y ~ Date * Management + (1|Block)`  
     Fallback (if the interaction is not estimable): `y ~ Date + Management + (1|Block)`
   - **Rainfed subset:** `y ~ Date * Management + (1|Block)`  
     Fallback (if the interaction is not estimable): `y ~ Date + Management + (1|Block)`
4. Computes **planned contrasts** (Irrigation vs Rainfed) **within each sampling date** using an LMM with a block random intercept.
5. Exports:
   - tidy ANOVA/contrasts/means tables (Excel),
   - a manuscript-oriented p-value table (Excel),
   - publication-quality mean ± SE figures (PNG + FIG),
   - and a run log.

---

## Requirements

- **MATLAB** with:
  - *Statistics and Machine Learning Toolbox* (for `fitlme`, `anova`, etc.)
- Functions used include `exportgraphics` and `savefig` (modern MATLAB releases).

---

## Input data

### Expected format

- **Excel file** (`.xlsx`) containing **one row per experimental unit**.
- A worksheet name (default: the first sheet in the file).

### Required columns (factor variables)

The script searches for these columns using exact matches first and then relaxed token-based matching:

- **Sampling date**: candidates include  
  `SamplingDate`, `Fecha`, `Fecha_ddmm`, `Sampling_Event_Date`
- **Irrigation regime**: candidates include  
  `IrrigationRegime`, `Riego`, `Irrigation`, `WaterRegime`
- **Vineyard floor management**: candidates include  
  `VineyardFloorManagement`, `Manejo_cultivo`, `Management`
- **Block**: candidates include  
  `Block`, `BloqueRoman`, `Bloque`, `Replicate`

If any of these are missing, the script stops with an explicit error message.

### Date parsing rules

The date column can be:

- MATLAB `datetime`
- Excel serial numbers
- Strings in formats:
  - `yyyy-MM-dd`
  - `dd/MM/yyyy`
  - `dd-MM-yyyy`
  - `dd/MM` (assumes year **2025**)

Dates are converted to an **ordered categorical** using `dd/MM` labels.

### Irrigation parsing rules

Text values are mapped as follows (case-insensitive):

- contains `"riego"` or `"irrig"` → **Irrigation**
- contains `"secano"` or `"rain"` → **Rainfed**

### Management parsing rules

Text values are mapped (case-insensitive) to:

- **open vegetation** (e.g., contains `open`, `vegetacion`)
- **tillage** (e.g., contains `tillage`, `laboreo`)
- **alternate cover crop** (e.g., contains `alternate`, `alterna`)
- **double cover crop** (e.g., contains `double`, `doble`)

Any other levels are retained as-is but placed after the preferred ordering.

### Block parsing rules

- Numeric blocks `1–4` are mapped to Roman numerals **I–IV**.
- Otherwise, block values are treated as categorical strings.

---

## Endpoints (chemical traits)

The script defines endpoint codes and resolves the corresponding Excel columns using candidate names:

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

**Important:** endpoint columns must be numeric (or importable as numeric). Non-numeric columns will be skipped.

---

## Outputs

All outputs are written under a parent folder named **`Objetivo_1`** next to the input Excel file:

```
<folder_of_inputXlsx>/
  Objetivo_1/
    RESULTS_Obj1_YYYYMMDD_HHMMSS/
      O1_log.txt
      O1_Stats_Summary.xlsx
      O1_Table1_Manuscript.xlsx
      Fig1_Irrigation_<CODE>.png / .fig
      Fig2_Management_Irrigated_<CODE>.png / .fig
      Fig3_Management_Rainfed_<CODE>.png / .fig
```

### `O1_Stats_Summary.xlsx` sheets

- `ANOVA_Global_DatexIrrig`  
  ANOVA (Satterthwaite DF) for the global LMM
- `ANOVA_Irrigated_Mgmt`  
  ANOVA for the irrigated-subset LMM (interaction or fallback)
- `ANOVA_Rainfed_Mgmt`  
  ANOVA for the rainfed-subset LMM (interaction or fallback)
- `Contr_Irrig_vs_Rainfed`  
  Planned contrast (Irrigation vs Rainfed) within each sampling date, with 95% CI
- `Means_DatexIrrig`  
  Mean, SD, n, SE for Date × Irrigation
- `Means_DatexMgmt_Irrig`  
  Mean, SD, n, SE for Date × Management (irrigated)
- `Means_DatexMgmt_Rain`  
  Mean, SD, n, SE for Date × Management (rainfed)

### `O1_Table1_Manuscript.xlsx`

- `Table1`  
  Manuscript-oriented p-values (formatted strings) for:
  - Global: Date, Irrigation, Date×Irrigation
  - Irrigated: Date, Management, Date×Management
  - Rainfed: Date, Management, Date×Management

### Figures

- **Fig1**: Date × Irrigation (mean ± SE)
- **Fig2**: Date × Management within **Irrigation**
- **Fig3**: Date × Management within **Rainfed**

Figures are exported as:
- `*.png` at **300 dpi** with white background
- `*.fig` for editable MATLAB figures

---

## How to run

1. Place `O1_Obj1_ChemGroundTruth_run_v11.m` on your MATLAB path.
2. Call:

```matlab
O1_Obj1_ChemGroundTruth_run_v11('Matriz_CHEM_HSI_MASTER_96.xlsx','Matriz');
```

If `sheetName` is omitted, the script uses the **first worksheet**.

---

## Notes for manuscript alignment (Section 3.1)

- The statistical backbone is LMMs with **Block** as a random intercept.
- p-values are obtained from **Satterthwaite** degrees of freedom (via `anova(...,'DFMethod','satterthwaite')`).
- Management effects are assessed **within irrigation regime** (irrigated and rainfed) to reflect the experimental structure.
- Planned per-date contrasts provide effect size directionality (Estimate) and uncertainty (95% CI), useful for concise reporting.

---

## Troubleshooting

- **“Missing required column …”**  
  Rename the Excel header to one of the recognised candidates or add a synonym.
- **“SamplingDate contains unparsable values …”**  
  Use a consistent date format (`yyyy-MM-dd` or `dd/MM/yyyy`) or convert to Excel serial / datetime.
- **Rank-deficient / not estimable interaction**  
  The script will attempt a **fallback model** (`Date + Management`) and logs this in `O1_log.txt`.
- **Empty or skipped endpoints**  
  Confirm the endpoint column exists and is numeric; check column naming against the endpoint candidate list.

---

## Citation / provenance

If you publish results produced by this script, cite it in your Methods as the implementation used to generate:
- Table 1 p-values (LMM ANOVA),
- per-date irrigation contrasts,
- and the mean ± SE figures used in Section 3.1.

