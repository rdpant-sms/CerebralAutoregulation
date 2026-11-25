# ---------------------------------------------------------------------------
# Pooling the differences of medians, with the {metamedian} and {meta} packages.
#
# For studies that reported each outcome group's median (IQR) or mean (SD)
# instead of an odds ratio. The effect is worse-outcome group minus
# better-outcome group, in percentage points of monitoring time or mm Hg.
#
# Writes ../Data/meta_outputs/median_diff_results.csv,
#        ../Data/meta_outputs/median_diff_study_values.csv
#        ../Data/figures/forest_mediandiff_<outcome>_<metric>.png
# ---------------------------------------------------------------------------
suppressMessages(library(meta))
suppressMessages(library(metamedian))
source("config.R")
source("forest_lib.R")

md  <- read.csv(paste0(STUDY_DATA, "meandiff_long.csv"), check.names = FALSE, stringsAsFactors = FALSE)
cov <- read.csv(paste0(STUDY_DATA, "covidence_consensus.csv"), check.names = FALSE, stringsAsFactors = FALSE)
cov$`Study ID` <- trimws(cov$`Study ID`)

# Cells to try, in the order they are written out. deltaP is not pooled because
# studies define it as CPP - CPPopt, CPPopt - CPP or |CPP - CPPopt|.
# MAPopt is left out because each MAPopt cell has a single study.
CELLS <- list(
  list(outcome = "poor_outcome", metric = "pct_below"),
  list(outcome = "poor_outcome", metric = "pct_above"),
  list(outcome = "poor_outcome", metric = "pct_within"),
  list(outcome = "poor_outcome", metric = "Popt"),
  list(outcome = "mortality",    metric = "Popt"),
  list(outcome = "mortality",    metric = "pct_below"),
  list(outcome = "mortality",    metric = "pct_above")
)

# Outcomes measured at discharge are not pooled (Kevci 2025).
discharge <- cov$`Study ID`[cov$`Follow-up duration` == "Other: Discharge"]

UNIT <- c(pct_below = "percentage points", pct_above = "percentage points", Popt = "mm Hg")
N_WORSE_COLUMN <- c(mortality = "number of non-survivors N", poor_outcome = "number with poor outcome N")

num <- function(x) suppressWarnings(as.numeric(trimws(x)))

# Klein 2018 gives a 95% CI of each group mean; turn it into an SD.
sd_from_ci <- function(text, n) {
  m <- regmatches(text, regexec("CI\\s*(-?[0-9.]+)\\s*-\\s*(-?[0-9.]+)", text))[[1]]
  if (length(m) < 3 || is.na(n)) return(NA_real_)
  sqrt(n) * (as.numeric(m[3]) - as.numeric(m[2])) / (2 * qt(0.975, n - 1))
}

# One group's summary as median/q1/q3 or mean/sd.
group_values <- function(stat, value, disp, low, high, n) {
  out <- c(med = NA, q1 = NA, q3 = NA, mean = NA, sd = NA)
  if (startsWith(stat, "median")) {
    out[c("med", "q1", "q3")] <- c(num(value), num(low), num(high))
  } else if (grepl("^CI", trimws(disp))) {
    out[c("mean", "sd")] <- c(num(value), sd_from_ci(disp, n))
  } else {
    out[c("mean", "sd")] <- c(num(value), num(disp))
  }
  out
}

group_text <- function(v) {
  if (!is.na(v["med"])) sprintf("%.1f [%.1f, %.1f]", v["med"], v["q1"], v["q3"])
  else sprintf("%.1f (%.1f)*", v["mean"], v["sd"])
}

# One row per study. Group 1 = worse outcome, group 2 = better outcome.
cell_rows <- function(cell) {
  d <- md[md$target == "CPPopt" & md$outcome == cell$outcome &
          md$exposure_metric == cell$metric & !(md$study_id %in% discharge), ]

  rows <- NULL
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    n_total  <- num(r$n_total)
    n_worse  <- num(cov[cov$`Study ID` == r$study_id, N_WORSE_COLUMN[cell$outcome]][1])
    n_better <- n_total - n_worse
    w <- group_values(r$stat_type, r$worse_value,  r$worse_disp,  r$worse_disp_low,  r$worse_disp_high,  n_worse)
    b <- group_values(r$stat_type, r$better_value, r$better_disp, r$better_disp_low, r$better_disp_high, n_better)

    rows <- rbind(rows, data.frame(
      study_id = r$study_id, target = r$target, etiology = r$etiology,
      outcome_def = r$outcome_def, stat_type = r$stat_type, threshold = r$threshold,
      n.g1 = n_worse,  med.g1 = w["med"], q1.g1 = w["q1"], q3.g1 = w["q3"], mean.g1 = w["mean"], sd.g1 = w["sd"],
      n.g2 = n_better, med.g2 = b["med"], q1.g2 = b["q1"], q3.g2 = b["q3"], mean.g2 = b["mean"], sd.g2 = b["sd"],
      worse_text = group_text(w), better_text = group_text(b),
      row.names = NULL))
  }

  # Drop studies that give no spread (IQR or SD), e.g. Lilja-Cyron 2021 for time below target
  has_spread <- ifelse(is.na(rows$med.g1),
                       !is.na(rows$sd.g1) & !is.na(rows$sd.g2),
                       !is.na(rows$q1.g1) & !is.na(rows$q3.g1) & !is.na(rows$q1.g2) & !is.na(rows$q3.g2))
  rows <- rows[has_spread, ]

  # One study per cohort: stop if two studies in a cell come from the same dataset
  datasets <- trimws(cov$`Dataset Cluster`[match(rows$study_id, cov$`Study ID`)])
  if (any(duplicated(datasets))) stop("overlapping cohorts in ", cell$outcome, " ", cell$metric)

  rows$Target     <- rows$target
  rows$Etiology   <- rows$etiology
  rows$OutcomeDef <- tidy_outcome_def(rows$outcome_def)
  rows$Thr        <- ifelse(trimws(rows$threshold) == "", "-", trimws(rows$threshold))
  rows[order(match(rows$Etiology, ETIOLOGY_ORDER), rows$study_id), ]
}

# Study-level difference of medians and its SE (quantile estimation, McGrath 2020).
# Median and mean studies go to metamedian separately; a mean is used as the median.
study_estimates <- function(rows) {
  rows$yi <- NA_real_
  rows$sei <- NA_real_
  is_med <- !is.na(rows$med.g1)
  if (any(is_med)) {
    e <- metamedian(rows[is_med, c("n.g1", "med.g1", "q1.g1", "q3.g1", "n.g2", "med.g2", "q1.g2", "q3.g2")],
                    median_method = "qe", pool_studies = FALSE)
    rows$yi[is_med] <- e$yi; rows$sei[is_med] <- e$sei
  }
  if (any(!is_med)) {
    e <- metamedian(rows[!is_med, c("n.g1", "mean.g1", "sd.g1", "n.g2", "mean.g2", "sd.g2")],
                    median_method = "qe", pool_studies = FALSE)
    rows$yi[!is_med] <- e$yi; rows$sei[!is_med] <- e$sei
  }
  rows
}

results <- NULL
study_values <- NULL

for (cell in CELLS) {
  rows <- study_estimates(cell_rows(cell))
  if (nrow(rows) < MIN_K) next

  model <- metagen(TE      = yi,
                   seTE    = sei,
                   studlab = study_id,
                   data    = rows,
                   sm      = "MD",
                   common  = FALSE,
                   random  = TRUE,
                   method.tau       = TAU_METHOD,
                   method.random.ci = CI_METHOD,
                   method.I2        = I2_METHOD,
                   prediction       = TRUE,
                   title            = paste(cell$outcome, cell$metric))

  results <- rbind(results, data.frame(
    outcome = cell$outcome, metric = cell$metric, unit = UNIT[cell$metric],
    k = model$k, n_patients = sum(rows$n.g1 + rows$n.g2),
    studies = paste(rows$study_id, collapse = "; "),
    n_reporting_means = sum(is.na(rows$med.g1)),
    diff_medians = model$TE.random, ci_low = model$lower.random, ci_high = model$upper.random,
    p = model$pval.random,
    I2 = 100 * model$I2, I2_low = 100 * model$lower.I2, I2_high = 100 * model$upper.I2,
    tau2 = model$tau2, pi_low = model$lower.predict, pi_high = model$upper.predict,
    row.names = NULL))

  study_values <- rbind(study_values, data.frame(
    outcome = cell$outcome, metric = cell$metric, study_id = rows$study_id,
    target = rows$target, etiology = rows$etiology, outcome_def = rows$OutcomeDef,
    stat_type = rows$stat_type, threshold = rows$threshold,
    n_worse = rows$n.g1, worse = rows$worse_text,
    n_better = rows$n.g2, better = rows$better_text,
    diff = rows$yi, se = rows$sei, ci_low = model$lower, ci_high = model$upper,
    weight_pct = 100 * model$w.random / sum(model$w.random), row.names = NULL))

  # --- forest plot ---
  shown <- drop_empty_columns(rows, list(labels  = c("Target", "Etiology", "Outcome def.", "Thr. (mmHg)"),
                                         columns = c("Target", "Etiology", "OutcomeDef", "Thr")))
  file <- paste0(FIG, "forest_mediandiff_", cell$outcome, "_", cell$metric, ".png")
  png(file, width = 11.5 + 1.1 * length(shown$columns), height = 2.8 + 0.32 * nrow(rows),
      units = "in", res = 200)
  forest(model, sortvar = NULL,
         leftcols  = c("studlab", shown$columns, "n.g1", "worse_text", "n.g2", "better_text"),
         leftlabs  = c("Study", shown$labels, "N", "Worse outcome", "N", "Better outcome"),
         rightcols = c("effect", "ci", "w.random"),
         rightlabs = c("Difference", "95% CI", "Weight"),
         label.left = "higher in better outcome", label.right = "higher in worse outcome",
         smlab = paste0("Difference of medians\n(", UNIT[cell$metric], ")"),
         digits = 1, print.tau2 = TRUE, print.I2.ci = TRUE,
         text.addline1 = "Median [Q1, Q3]; * mean (SD), median taken as the mean")
  invisible(dev.off())
  cat("wrote", basename(file), "-", nrow(rows), "studies\n")
}

write.csv(results,      paste0(META_OUTPUTS, "median_diff_results.csv"),      row.names = FALSE)
write.csv(study_values, paste0(META_OUTPUTS, "median_diff_study_values.csv"), row.names = FALSE)
cat("\nwrote median_diff_results.csv:", nrow(results), "cells\n")
