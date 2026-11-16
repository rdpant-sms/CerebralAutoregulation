# Paths and shared settings for the {meta}-package implementation.

BASE <- "../"
DATA <- paste0(BASE, "Data/")

# Data/ is organised into study_data/ (one row per study), meta_outputs/ (pooled
# results) and figures/. Word documents and the analysis log stay at the top level.
FIG          <- paste0(DATA, "figures/")
STUDY_DATA   <- paste0(DATA, "study_data/")
META_OUTPUTS <- paste0(DATA, "meta_outputs/")
for (d in c(FIG, STUDY_DATA, META_OUTPUTS))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

EFFECTS      <- paste0(STUDY_DATA, "effects_long.csv")
META_RESULTS <- paste0(META_OUTPUTS, "meta_results.csv")
CORR_RESULTS <- paste0(META_OUTPUTS, "corr_results.csv")

# Analysis settings, fixed in one place so every script uses the same model.
# These follow the recommendations in Doing Meta-Analysis in R (Harrer et al.):
#   - always fit a random-effects model in clinical/observational work
#   - REML for the between-study variance with continuous-type effect sizes
#   - Knapp-Hartung adjustment, which matters most when k is small (ours are 3-7)
#   - always report a prediction interval alongside the pooled effect
TAU_METHOD <- "REML"
CI_METHOD  <- "HK"

# Which I^2 the reported I2 column holds: "tau2" (from the fitted REML tau^2) or
# "Q" (Higgins & Thompson). Both are always written out; see tidy_meta() below.
I2_METHOD  <- "tau2"
MIN_K      <- 3

# Reporting increment for the odds ratios. The estimates are stored per 1% of
# monitoring time, which is too small an increment to read: a real effect looks
# like "OR 1.02". Everything is therefore reported per 10%. log(OR) is linear in
# the increment, so multiplying yi and sei by 10 is an exact re-expression. It
# leaves I^2, the weights, Q and every p-value unchanged, and multiplies tau^2
# by 100. Set OR_UNIT to 1 to go back to per-1% reporting.
OR_UNIT       <- 10
OR_UNIT_LABEL <- "per 10% of monitoring time"

METRICS    <- c("pct_below", "pct_above", "pct_within", "pct_outside")

ETIOLOGY <- c("Traumatic Brain Injury" = "TBI", "Subarachnoid Hemorrhage" = "SAH",
              "Intracerebral Hemorrhage" = "ICH", "Ischemic Stroke" = "IS")

METRIC_LABEL <- c(pct_below = "% of monitoring time below the target",
                  pct_above = "% of monitoring time above the target",
                  pct_within = "% of monitoring time within the optimal range",
                  pct_outside = "% of monitoring time outside the limits of autoregulation")

# Read the analysis file and keep the poolable rows of one effect type.
# `sets` selects the analysis set: "primary" is cohort-overlap de-duplicated,
# adding "sensitivity" re-includes overlapping records.
read_cells <- function(effect, sets) {
  d <- read.csv(EFFECTS, check.names = FALSE, stringsAsFactors = FALSE)
  d <- d[d$effect_type == effect & d$poolable == "yes" &
         d$analysis_set %in% sets & d$exposure_metric %in% METRICS, ]
  d$yi <- as.numeric(d$yi); d$sei <- as.numeric(d$sei)
  d$n_total <- suppressWarnings(as.numeric(d$n_total))
  d <- d[!is.na(d$yi) & !is.na(d$sei) & d$sei > 0, ]
  # re-express the odds ratios on the reporting increment (see OR_UNIT)
  if (effect == "OR") {
    d$yi  <- d$yi  * OR_UNIT
    d$sei <- d$sei * OR_UNIT
  }
  d$etiology_short <- ifelse(d$etiology %in% names(ETIOLOGY), ETIOLOGY[d$etiology],
                             substr(d$etiology, 1, 3))
  d
}

# Convert a number from the scale the model works on to the scale we report.
# Odds ratios are modelled as log odds ratios, correlations as Fisher z.
to_plot_scale <- function(x, effect) {
  if (effect == "OR") exp(x) else tanh(x)
}

# One row of results from a fitted model, so every script reports the same
# numbers in the same order.
tidy_meta <- function(model, effect, analysis, target, outcome, metric) {
  i2_q  <- 100 * update(model, method.I2 = "Q")$I2
  i2_t2 <- 100 * update(model, method.I2 = "tau2")$I2

  data.frame(
    analysis = analysis,
    target   = target,
    outcome  = outcome,
    metric   = metric,
    k        = model$k,
    est      = to_plot_scale(model$TE.random,     effect),
    ci_low   = to_plot_scale(model$lower.random,  effect),
    ci_high  = to_plot_scale(model$upper.random,  effect),
    pi_low   = to_plot_scale(model$lower.predict, effect),
    pi_high  = to_plot_scale(model$upper.predict, effect),
    I2       = if (I2_METHOD == "Q") i2_q else i2_t2,
    I2_Q     = i2_q,
    I2_tau2  = i2_t2,
    tau2     = model$tau2,
    pval     = model$pval.random,
    Q        = model$Q,
    Q_pval   = model$pval.Q,
    stringsAsFactors = FALSE)
}

# Columns that get rounded before the results are written out.
NUMCOLS <- c("est", "ci_low", "ci_high", "pi_low", "pi_high",
             "I2", "I2_Q", "I2_tau2", "tau2", "pval", "Q", "Q_pval")
