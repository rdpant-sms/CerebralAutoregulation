# Sensitivity odds-ratio forests: overlapping cohorts re-included.

source("config.R")
source("forest_lib.R")

data <- add_label_columns(read_cells("OR", c("primary", "sensitivity")))
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
    model <- fit_model(rows, "OR")
    file  <- paste0("forest_sens_OR_", outcome, "_", metric, ".png")

    draw_forest(model, rows, file, COLUMNS_OR_SENS, "OR", width = 13.0)
    figures_made <- figures_made + 1
  }
}
cat("SENSITIVITY OR FORESTS DONE:", figures_made, "figures\n")
