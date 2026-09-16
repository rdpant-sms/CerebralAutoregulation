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

I2_METHOD  <- "tau2"
MIN_K      <- 3
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
#
# Two I^2 values are stored on purpose. {meta} uses the Higgins & Thompson (2002)
# formula I2 = (Q - df) / Q, which is also the Cochrane Handbook one. The other
# definition works from the fitted between-study variance, I2 = tau2/(tau2 + s2).
# They agree under a DerSimonian-Laird tau^2 but can differ sharply under REML,
# so both are written out and I2_METHOD picks which one the `I2` column holds.
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
