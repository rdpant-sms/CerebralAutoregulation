# Primary odds-ratio forest plots: one record per institutional cohort.

source("config.R") 
source("forest_lib.R")

data <- add_label_columns(read_cells("OR", "primary"))

for (outcome in c("mortality", "poor_outcome")) {
  for (metric in METRICS) {

    rows <- data[data$outcome == outcome & data$exposure_metric == metric, ]
    if (nrow(rows) < MIN_K) next

    rows  <- sort_rows(rows)
    model <- fit_model(rows, "OR")
    file  <- paste0("forest_OR_", outcome, "_", metric, ".png")

    draw_forest(model, rows, file, COLUMNS_OR, "OR")
  }
}
cat("PRIMARY OR FORESTS DONE\n")
