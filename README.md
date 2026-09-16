# Meta-analysis

Following Harrer, Cuijpers, Furukawa & Ebert, *Doing Meta-Analysis in R* (<https://doing-meta.guide>)

Set working directory to code_meta/ and run the scripts in order 01 -> 07

Inputs come from `../Data/study_data/effects_long.csv`. Outputs go to `../Data/meta_outputs/` (`meta_results.csv`, `corr_results.csv`) and `../Data/figures/`, Paths are defined once in `config.R`, which creates the directories if they are missing.

R version: 4.5.0
Requires R \>= 4.3 with `meta`

## Data organization
```         
/CerebralAutoregulation
|-- /code_meta
    |-- 01_pool_odds_ratios.R
    |-- 02_pool_correlations.R
    |-- 03_forest_primary_OR.R
    |-- 04_forest_sensitivity_OR.R
    |-- 05_forest_primary_corr.R
    |-- 06_forest_sensitivity_corr.R
    |-- 07_leave_one_out.R
    |-- config.R
    |-- forest_lib.R
|-- /Data
    |-- /figures
    |-- /meta_outputs
    |-- /study_data
```

## Files

| Script | What it does |
|------------------------------------|------------------------------------|
| `config.R` | Sets paths and the model settings|
| `forest_lib.R` | Axis fitting, label columns, the shared `forest()` call |
| `01_pool_odds_ratios.R` | `metagen()` per cell |
| `02_pool_correlations.R` | `metacor()` per cell |
| `03_forest_primary_OR.R` | CPPopt and MAPopt as subgroup blocks, plus the combined estimate |
| `04_forest_sensitivity_OR.R` | Same, with overlapping cohorts re-included |
| `05_forest_primary_corr.R` | Correlation forests |
| `06_forest_sensitivity_corr.R` | Correlation forests, overlaps re-included |
| `07_leave_one_out.R` | Sensitivity analysis to check robustness of results |


``` bash
Rscript 01_pool_odds_ratios.R
Rscript 02_pool_correlations.R
Rscript 03_forest_primary_OR.R
Rscript 04_forest_sensitivity_OR.R
Rscript 05_forest_primary_corr.R
Rscript 06_forest_sensitivity_corr.R
Rscript 07_leave_one_out.R
```

## Figure inventory

Written to `../Data/figures/`:

| Figure | Content |
|------------------------------------|------------------------------------|
| `forest_OR_<outcome>_<metric>.png` | CPPopt and MAPopt as separate subgroup blocks, each with its own subtotal and the combined estimate on the bottom line |
| `forest_sens_OR_<outcome>_<metric>.png` | Overlapping cohorts re-included, with a Dataset column and daggers on re-entered records |
| `figS_corr_<metric>.png` | Correlations, primary set |
| `figS_corr_sens_<metric>.png` | Correlations, overlapping cohorts re-included |

Note that **every correlation in the dataset is CPPopt-based**: no MAPopt study reported a correlation coefficient, and none reported a correlation with mortality. The correlation forests therefore have no target split to show.
