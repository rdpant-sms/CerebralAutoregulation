# ---------------------------------------------------------------------------
# Pooling the correlations, with the {meta} package.
# ---------------------------------------------------------------------------
suppressMessages(library(meta))
source("config.R")

pool <- function(rows, label) {
  metacor(cor     = as.numeric(rows$estimate),
          n       = rows$n_total,
          studlab = rows$study_id,
          common  = FALSE,
          random  = TRUE,
          method.tau       = TAU_METHOD,
          method.random.ci = CI_METHOD,
          method.I2        = I2_METHOD,
          prediction       = TRUE,
          title            = label)
}

results <- NULL

for (analysis in c("primary", "sensitivity_incl_overlap")) {

  if (analysis == "primary") {
    data <- read_cells("correlation", "primary")
  } else {
    data <- read_cells("correlation", c("primary", "sensitivity"))
  }

  for (target in c("CPPopt", "MAPopt")) {
    for (outcome in c("mortality", "poor_outcome")) {
      for (metric in METRICS) {

        rows <- data[data$target == target &
                     data$outcome == outcome &
                     data$exposure_metric == metric, ]
        if (nrow(rows) < MIN_K) next

        model <- pool(rows, paste(analysis, target, outcome, metric))
        results <- rbind(results,
          tidy_meta(model, "correlation", analysis, target, outcome, metric))

        if (analysis == "primary") {
          cat("\n=====", target, outcome, metric, "=====\n")
          print(summary(model))
        }
      }
    }
  }
}

for (column in NUMCOLS) {
  results[[column]] <- round(as.numeric(results[[column]]), 4)
}
write.csv(results, CORR_RESULTS, row.names = FALSE)

cat("\nwrote corr_results.csv:", nrow(results), "rows\n")
print(results[, c("analysis", "target", "outcome", "metric", "k",
                  "est", "ci_low", "ci_high", "I2_Q", "I2_tau2")],
      row.names = FALSE)

