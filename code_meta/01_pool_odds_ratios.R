# ---------------------------------------------------------------------------
# Pooling the odds ratios, with the {meta} package.
#
# Three analyses are run:
#   1. combined targets             CPPopt and MAPopt pooled together, as subgroups
#   2. separate targets             one record per institutional cohort, each target separately
#   3. overlapping included         the same as no. 1, with overlapping cohorts re-included

# ---------------------------------------------------------------------------
suppressMessages(library(meta))
source("config.R")

# Pool one set of studies. 
# TAU_METHOD, CI_METHOD, I2_METHOD set in Config.R
pool <- function(rows, label) {
  metagen(TE      = rows$yi,
          seTE    = rows$sei,
          studlab = rows$study_id,
          data    = rows,
          sm      = "OR",
          common  = FALSE,
          random  = TRUE,
          method.tau       = TAU_METHOD,
          method.random.ci = CI_METHOD,
          method.I2        = I2_METHOD,
          prediction       = TRUE,
          title            = label)
}

results <- NULL 

# --- 1. CPPopt and MAPopt pooled together ------------------------
# Pooling CPPopt and MAPopt targets together
# Also gives the test for a difference between the subgroups 
# Only runs where both targets actually contribute.
for (outcome in c("mortality", "poor_outcome")) {
  for (metric in METRICS) {
    
    rows <- primary[primary$outcome == outcome & primary$exposure_metric == metric, ]
    if (nrow(rows) < MIN_K) next
    if (length(unique(rows$target)) < 2) next
    
    model <- metagen(TE      = rows$yi,
                     seTE    = rows$sei,
                     studlab = rows$study_id,
                     data    = rows,
                     sm      = "OR",
                     common  = FALSE,
                     random  = TRUE,
                     method.tau       = TAU_METHOD,
                     method.random.ci = CI_METHOD,
                     method.I2        = I2_METHOD,
                     prediction       = TRUE,
                     subgroup         = rows$target,
                     title            = paste("Any target:", outcome, metric))
    
    cat("\n===== ANY TARGET:", outcome, metric, "=====\n")
    print(summary(model))
    
    results <- rbind(results,
                     tidy_meta(model, "OR", "any_target", "Any", outcome, metric))
  }
}

# --- 2. separate targets (CPPopt and MAPopt separate) ----------------------------------------------------
primary <- read_cells("OR", "primary")

for (target in c("CPPopt", "MAPopt")) {
  for (outcome in c("mortality", "poor_outcome")) {
    for (metric in METRICS) {

      rows <- primary[primary$target == target &
                      primary$outcome == outcome &
                      primary$exposure_metric == metric, ]
      if (nrow(rows) < MIN_K) next

      model <- pool(rows, paste(target, outcome, metric))
      cat("\n=====", target, outcome, metric, "=====\n")
      print(summary(model))

      results <- rbind(results,
        tidy_meta(model, "OR", "separate targets", target, outcome, metric))
    }
  }
}


# --- 3. overlapping cohorts re-included ------------------------
sensitivity <- read_cells("OR", c("primary", "sensitivity"))

for (target in c("CPPopt", "MAPopt")) {
  for (outcome in c("mortality", "poor_outcome")) {
    for (metric in METRICS) {

      rows <- sensitivity[sensitivity$target == target &
                          sensitivity$outcome == outcome &
                          sensitivity$exposure_metric == metric, ]
      if (nrow(rows) < MIN_K) next

      model <- pool(rows, paste(target, outcome, metric, "(overlap re-included)"))
      cat("\n=====", target, outcome, metric, "(overlap re-included) =====\n")
      print(summary(model))
      
      results <- rbind(results,
        tidy_meta(model, "OR", "sensitivity_incl_overlap", target, outcome, metric))
    }
  }
}


# --- write out ---
for (column in NUMCOLS) {
  results[[column]] <- round(as.numeric(results[[column]]), 4)
}
write.csv(results, META_RESULTS, row.names = FALSE)

cat("\nwrote meta_results.csv:", nrow(results), "rows\n")
print(results[, c("analysis", "target", "outcome", "metric", "k",
                  "est", "ci_low", "ci_high", "I2_Q", "I2_tau2")],
      row.names = FALSE)


