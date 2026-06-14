# Causal survival forest analysis of JOINT-04/05 (Anabolic-first heterogeneity)

Analysis code for a pooled causal survival forest study identifying postmenopausal
osteoporosis patients who derive the greatest benefit from choosing an anabolic agent
(teriparatide, a parathyroid hormone [PTH] analogue) over a bisphosphonate (BP) as
initial therapy, using individual-patient data from the JOINT-04 and JOINT-05 trials.

This repository contains **only the code** needed to reproduce the main analysis
(preprocessing and the causal survival forest for the fracture outcomes). **No patient
data are included** (see [Data availability](#data-availability)).

## What this code does

- Pools the minodronate (BP) arm of JOINT-04 (n = 1,623) and all of JOINT-05
  (teriparatide n = 489 vs alendronate n = 496), for a total of 2,608 patients.
- Estimates each patient's propensity score with a regression forest and fits a
  causal survival forest (`grf`) for three fracture outcomes:
  E1 (vertebral), E2 (nonvertebral), and E3 (composite = vertebral or nonvertebral).
- The effect measure is the difference in restricted mean survival time (RMST) at a
  2.0-year horizon (positive = longer fracture-free time with PTH).
- Reports the overall average treatment effect, AUTOC and Cochran's Q tests for
  heterogeneity, group average treatment effects by CATE quartile, best linear
  projection (effect modifiers), and variable importance.
- A sensitivity-analysis script reproduces the main findings after excluding prevalent
  vertebral-fracture covariates, restricting to JOINT-05 alone (the true 1:1 RCT), and
  extending the horizon to 2.5 years.

## Repository layout

```
.
├── scripts/
│   └── 02_survival_preprocessing.py   # Excel -> outputs/joint0405_survival.csv
├── R/
│   ├── 00_setup_survival.R            # covariates, constants, paths (sourced by others)
│   ├── 40_csf_main.R                  # main causal survival forest (E1/E2/E3)
│   ├── 41_csf_sensitivity.R           # sensitivity analyses
│   └── utils/
│       ├── cf_helpers.R               # GATE, Cochran's Q, RATE/AUTOC, variable importance
│       └── var_labels.R               # human-readable variable labels for figures
├── data/                              # (empty) place the raw Excel file here — not distributed
├── outputs/                           # (empty) generated CSVs, tables, and figures land here
├── renv.lock                          # pinned R package versions (R 4.6.0)
├── .Rprofile                          # renv bootstrap
├── LICENSE                            # MIT
└── README.md
```

## Requirements

- **R 4.6.0** with the packages pinned in `renv.lock` (key packages: `grf` 2.6.1,
  `survival` 3.8-6, `readr`, `dplyr`, `ggplot2`, `patchwork`, `here`).
- **Python 3.10+** with `pandas` and `openpyxl` for preprocessing.

## How to reproduce

1. **Obtain the raw data** (see [Data availability](#data-availability)) and place the
   Excel file at `data/joint0405_raw.xlsx`.

2. **Preprocess** (writes `outputs/joint0405_survival.csv`):

   ```bash
   python scripts/02_survival_preprocessing.py
   ```

3. **Restore the R environment** (first time only):

   ```r
   renv::restore()
   ```

4. **Run the analysis** from the repository root:

   ```bash
   Rscript R/40_csf_main.R          # main causal survival forest (E1/E2/E3)
   Rscript R/41_csf_sensitivity.R   # sensitivity analyses
   ```

   Tables and figures are written under `outputs/`. The random seed is fixed (seed = 42)
   and honest splitting is used, so results are reproducible.

## Data availability

The JOINT-04 and JOINT-05 data were collected by the A-TOP research group and are **not
publicly available** because of restrictions on the use of the original trial data. The
raw Excel file is therefore not included in this repository. The de-identified data that
support the findings, together with this code, are available from the corresponding
author on reasonable request and with the permission of the A-TOP research group.

The preprocessing script expects a single Excel file at `data/joint0405_raw.xlsx` whose
first sheet contains one row per patient with the covariate, treatment-assignment, and
time-to-event columns referenced in `scripts/02_survival_preprocessing.py`.

## Citation

If you use this code, please cite the associated article (citation to be added upon
publication) and this repository.

## License

Released under the [MIT License](LICENSE).
