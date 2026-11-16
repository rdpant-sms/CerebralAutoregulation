# Meta-analysis
Deviation From Autoregulation-Derived Optimal Pressure Targets and Outcome After Acute Brain Injury: A Systematic Review and Meta-Analysis

Following Harrer, Cuijpers, Furukawa & Ebert, *Doing Meta-Analysis in R* (<https://doing-meta.guide>)

Set working directory to code_meta/ and run `Rscript run_all.R`, which runs the scripts in order 01 -> 09 and writes the R and package versions to `../Data/meta_outputs/session_info.txt`. Each script can also be run on its own, in order.

Inputs come from `../Data/study_data/effects_long.csv`. Outputs go to `../Data/meta_outputs/` (`meta_results.csv`, `corr_results.csv`) and `../Data/figures/`, Paths are defined once in `config.R`, which creates the directories if they are missing.

R version: 4.5.0, `meta` 8.3.0

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
    |-- 08_forest_etiology.R
    |-- 09_forest_unpooled.R
    |-- run_all.R
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
| `config.R` | Paths, the model settings |
| `forest_lib.R` | Axis fitting, label columns, the shared `forest()` call |
| `01_pool_odds_ratios.R` | `metagen()` per cell |
| `02_pool_correlations.R` | `metacor()` per cell |
| `03_forest_primary_OR.R` | CPPopt and MAPopt as subgroup blocks, plus the combined estimate |
| `04_forest_sensitivity_OR.R` | Same, with overlapping cohorts re-included |
| `05_forest_primary_corr.R` | Correlation forests |
| `06_forest_sensitivity_corr.R` | Correlation forests, overlaps re-included |
| `07_leave_one_out.R` | Leave one out analysis to check robustness of results |
| `08_forest_etiology.R` | Primary cells grouped by diagnosis, tests for subgroup differences by diagnosis and by target, and pooling within a diagnosis where it reaches 3 studies |
| `09_forest_unpooled.R` | Odds ratios from cells below the 3-study floor, shown study by study without pooling |
| `run_all.R` | Runs 01 to 09 in order and records the session information |

Requires R \>= 4.3 with `meta`. Neither `metafor` nor the guide's companion package `dmetar` is used.

``` bash
Rscript 01_pool_odds_ratios.R
Rscript 02_pool_correlations.R
Rscript 03_forest_primary_OR.R
Rscript 04_forest_sensitivity_OR.R
Rscript 05_forest_primary_corr.R
Rscript 06_forest_sensitivity_corr.R
Rscript 07_leave_one_out.R
Rscript 08_forest_etiology.R
Rscript 09_forest_unpooled.R
```

or, equivalently, `Rscript run_all.R`.

## Figure inventory

Written to `../Data/figures/`:

| Figure | Content |
|------------------------------------|------------------------------------|
| `forest_OR_<outcome>_<metric>.png` | CPPopt and MAPopt as separate subgroup blocks, each with its own subtotal and the combined estimate on the bottom line |
| `forest_sens_OR_<outcome>_<metric>.png` | Overlapping cohorts re-included, with a Dataset column and daggers on re-entered records |
| `figS_corr_<metric>.png` | Correlations forest plots |
| `figS_corr_sens_<metric>.png` | Correlations forest plots, overlapping cohorts re-included |
| `forest_etiology_OR_<outcome>_<metric>.png` | Primary cells with diagnoses as subgroups, with the test for subgroup differences |
| `forest_unpooled_<outcome>_<cell>.png` | Cells below the 3-study floor, one estimate per study, no pooled estimate; `forest_unpooled_panels.csv` lists the panel order and titles |

## Results files

Written to `../Data/meta_outputs/`:

| File | Content |
|------------------------------------|------------------------------------|
| `meta_results.csv` | Pooled odds ratios: combined targets, separate targets, and overlapping cohorts re-included |
| `corr_results.csv` | Pooled correlations |
| `leave_one_out.csv` | Each primary cell refitted with each study omitted in turn |
| `subgroup_tests.csv` | Tests for subgroup differences by diagnosis and by target |
| `etiology_results.csv` | Pooled odds ratios within a single diagnosis |
| `session_info.txt` | R and package versions, written by `run_all.R` |

Odds ratios for percentage-of-monitoring-time metrics are reported per 10% of monitoring time (`OR_UNIT` in `config.R`); estimates are stored per 1% in `effects_long.csv`.

