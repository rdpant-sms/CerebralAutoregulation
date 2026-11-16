# Run the whole meta-analysis in order, from code_meta/:
#
#   Rscript run_all.R
#
# Each script is run in a fresh environment, so none depends on objects left
# behind by another. The R and package versions used are written to
# ../Data/meta_outputs/session_info.txt for the methods section.

scripts <- c(
  "01_pool_odds_ratios.R",       # pooled odds ratios          -> meta_results.csv
  "02_pool_correlations.R",      # pooled correlations         -> corr_results.csv
  "03_forest_primary_OR.R",      # primary OR forests          (Figures 3, 5, 6)
  "04_forest_sensitivity_OR.R",  # overlapping cohorts re-included (eFigure 1)
  "05_forest_primary_corr.R",    # correlation forests         (eFigure 4)
  "06_forest_sensitivity_corr.R",# correlation forests, overlaps re-included (eFigure 5)
  "07_leave_one_out.R",          # leave-one-out               -> leave_one_out.csv
  "08_forest_etiology.R",        # by diagnosis (eFigure 2)    -> subgroup_tests.csv, etiology_results.csv
  "09_forest_unpooled.R"         # cells below MIN_K, unpooled (eFigure 3)
)

for (script in scripts) {
  cat("\n#####", script, "#####\n")
  sys.source(script, envir = new.env(parent = globalenv()), chdir = TRUE)
}

source("config.R")
writeLines(capture.output(sessionInfo()), paste0(META_OUTPUTS, "session_info.txt"))
cat("\nALL DONE. Versions written to", paste0(META_OUTPUTS, "session_info.txt"), "\n")
