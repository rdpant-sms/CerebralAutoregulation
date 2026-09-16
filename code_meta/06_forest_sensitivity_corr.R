# Supplementary correlation forests with overlapping cohorts re-included.

source("config.R")
source("forest_lib.R")

data <- add_label_columns(read_cells("correlation", c("primary", "sensitivity")))
figures_made <- 0

for (outcome in c("mortality", "poor_outcome")) {
  for (metric in METRICS) {

    rows <- data[data$outcome == outcome & data$exposure_metric == metric, ]
    if (nrow(rows) < MIN_K) next

    if (!any(rows$analysis_set == "sensitivity")) {
      cat("skipped", outcome, metric, "- no overlapping records in this cell\n")
      next
    }

    rows  <- sort_rows(rows, by_dataset = TRUE)
    model <- fit_model(rows, "correlation")
    file  <- paste0("figS_corr_sens_", metric, ".png")

    draw_forest(model, rows, file, COLUMNS_CORR_SENS, "correlation", width = 13.2)
    figures_made <- figures_made + 1
  }
}
cat("SENSITIVITY CORRELATION FORESTS DONE:", figures_made, "figures\n")
