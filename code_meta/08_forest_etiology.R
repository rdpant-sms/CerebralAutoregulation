# Subgroup analysis by etiology (diagnosis), for the supplement.
#
# Two outputs:
#   1. a forest plot of each primary cell subgrouped by diagnosis rather than by
#      target, which keeps every study and gives a test for subgroup differences;
#   2. etiology_results.csv, the same cells pooled separately within each diagnosis,
#      for the cells that reach MIN_K on their own.
#
# Separate within-diagnosis pooling is under-powered by construction: splitting a
# cell three ways leaves k = 3, where the Hartung-Knapp interval uses t on 2 df
# (critical value 4.30). The subgroup figure is the more informative of the two.

source("config.R")
source("forest_lib.R")

ETIOLOGY_NAMES <- c(TBI = "Traumatic brain injury",
                    SAH = "Subarachnoid hemorrhage",
                    ICH = "Intracerebral hemorrhage",
                    IS  = "Ischemic stroke")

fit_by_etiology <- function(data) {
  present <- ETIOLOGY_ORDER[ETIOLOGY_ORDER %in% data$Etiology]
  groups  <- factor(unname(ETIOLOGY_NAMES[data$Etiology]),
                    levels = unname(ETIOLOGY_NAMES[present]))
  metagen(TE = data$yi, seTE = data$sei, studlab = data$slab, data = data,
          sm = "OR", common = FALSE, random = TRUE,
          method.tau = TAU_METHOD, method.random.ci = CI_METHOD,
          method.I2 = I2_METHOD, prediction = TRUE,
          subgroup = groups, print.subgroup.name = FALSE)
}

# rows are ordered by diagnosis first, then by target, then study
sort_by_etiology <- function(data) {
  data[order(match(data$Etiology, ETIOLOGY_ORDER),
             match(data$target, c("CPPopt", "MAPopt")), data$study_id), ]
}

COLUMNS_ETIOLOGY <- list(
  labels  = c("Target", "N", "Outcome def.", "Thr. (mmHg)"),
  columns = c("target", "N", "OutcomeDef",   "Thr"))

data <- add_label_columns(read_cells("OR", "primary"))
figures <- 0
# tests for subgroup differences, by diagnosis here and by target (the grouping of
# the primary forest plots, fit_model() in forest_lib.R), for the manuscript text
subgroup_tests <- NULL
add_test <- function(by, model) data.frame(
  outcome = outcome, metric = metric, by = by, Q = round(model$Q.b.random, 4),
  df = model$df.Q.b, pval = round(model$pval.Q.b.random, 4))

for (outcome in c("mortality", "poor_outcome")) {
  for (metric in METRICS) {
    rows <- data[data$outcome == outcome & data$exposure_metric == metric, ]
    if (nrow(rows) < MIN_K) next
    if (length(unique(rows$Etiology)) < 2) next     # nothing to compare

    rows  <- sort_by_etiology(rows)
    model <- fit_by_etiology(rows)
    file  <- paste0("forest_etiology_OR_", outcome, "_", metric, ".png")
    draw_forest(model, rows, file, COLUMNS_ETIOLOGY, "OR")
    figures <- figures + 1

    cat(sprintf("   test for subgroup differences: Q = %.2f, df = %d, p = %.4f\n",
                model$Q.b.random, model$df.Q.b, model$pval.Q.b.random))
    subgroup_tests <- rbind(subgroup_tests, add_test("etiology", model))
    if (length(unique(rows$target)) > 1)
      subgroup_tests <- rbind(subgroup_tests, add_test("target", fit_model(rows, "OR")))
  }
}
write.csv(subgroup_tests, paste0(META_OUTPUTS, "subgroup_tests.csv"), row.names = FALSE)

# --- separately pooled within each diagnosis ---------------------------------
results <- NULL
for (sets in list("primary", c("primary", "sensitivity"))) {
  d <- read_cells("OR", sets)
  label <- if (length(sets) == 1) "primary" else "sensitivity_incl_overlap"
  for (et in ETIOLOGY_ORDER) {
    for (outcome in c("mortality", "poor_outcome")) {
      for (metric in METRICS) {
        rows <- d[d$etiology_short == et & d$outcome == outcome &
                  d$exposure_metric == metric, ]
        if (nrow(rows) < MIN_K) next
        model <- metagen(TE = rows$yi, seTE = rows$sei, studlab = rows$study_id,
                         sm = "OR", common = FALSE, random = TRUE,
                         method.tau = TAU_METHOD, method.random.ci = CI_METHOD,
                         method.I2 = I2_METHOD, prediction = TRUE)
        row <- tidy_meta(model, "OR", label, et, outcome, metric)
        row$targets <- paste(sort(unique(rows$target)), collapse = " + ")
        row$studies <- paste(sort(rows$study_id), collapse = "; ")
        results <- rbind(results, row)
      }
    }
  }
}
for (column in NUMCOLS) results[[column]] <- round(as.numeric(results[[column]]), 4)
write.csv(results, paste0(META_OUTPUTS, "etiology_results.csv"), row.names = FALSE)

cat("\nETIOLOGY SUBGROUPS DONE:", figures, "figures,", nrow(results), "pooled cells\n")
print(results[, c("analysis", "target", "outcome", "metric", "k", "est",
                  "ci_low", "ci_high", "pval", "targets")], row.names = FALSE)
