# Supplementary correlation forest plots, primary set.

source("config.R")
source("forest_lib.R")

data <- add_label_columns(read_cells("correlation", "primary"))

for (outcome in c("mortality", "poor_outcome")) {
  for (metric in METRICS) {

    rows <- data[data$outcome == outcome & data$exposure_metric == metric, ]
    if (nrow(rows) < MIN_K) next

    rows  <- sort_rows(rows)
    model <- fit_model(rows, "correlation")
    file  <- paste0("figS_corr_", metric, ".png")

    draw_forest(model, rows, file, COLUMNS_CORR, "correlation", width = 12.6)
  }
}
cat("CORRELATION FORESTS DONE\n")
