# ---------------------------------------------------------------------------
# Forest plots, drawn with meta::forest().
# ---------------------------------------------------------------------------
suppressMessages(library(meta))

# Headings for the two target types, used as the subgroup labels.
TARGET_NAMES <- c(CPPopt = "CPPopt",
                  MAPopt = "MAPopt")

# Diagnoses are listed in this order inside each target block.
ETIOLOGY_ORDER <- c("TBI", "SAH", "ICH", "IS")


# Cleaning up text
tidy_outcome_def <- function(text) {
  text <- sub("^ *Other: *", "", text)
  text <- gsub(" *>= *", " ≥ ", text)
  text <- gsub(" *<= *", " ≤ ", text)
  text <- trimws(text)
  ifelse(text == "", "-", text)
}

# Which threshold to print depends on the metric: "below" uses the lower bound,
# "above" uses the upper bound, "within" and "outside" use both.
threshold_text <- function(metric, lower, upper) {
  lower <- trimws(ifelse(is.na(lower), "", as.character(lower)))
  upper <- trimws(ifelse(is.na(upper), "", as.character(upper)))
  out <- rep("-", length(metric))
  for (i in seq_along(metric)) {
    if (metric[i] == "pct_below") {
      if (lower[i] != "") out[i] <- lower[i]
    } else if (metric[i] == "pct_above") {
      if (upper[i] != "") out[i] <- upper[i]
    } else {
      if (lower[i] != "" || upper[i] != "") out[i] <- paste0(lower[i], " to ", upper[i])
    }
  }
  out
}

# Add every column a figure might print. 
add_label_columns <- function(data) {
  data$Etiology   <- data$etiology_short
  data$Dataset    <- data$dataset_cluster
  data$N          <- data$n_total
  data$OutcomeDef <- tidy_outcome_def(data$outcome_def)
  data$Scale      <- data$outcome_scale
  data$Method     <- data$corr_method
  data$Thr        <- threshold_text(data$exposure_metric,
                                    data$threshold_lower, data$threshold_upper)
  # records re-entered from the overlap de-duplication are marked with a dagger
  data$slab <- ifelse(data$analysis_set == "sensitivity",
                      paste0(data$study_id, " †"), data$study_id)
  data
}

# Which label columns each kind of figure shows.
COLUMNS_OR <- list(
  labels  = c("Etiology", "N", "Outcome def.", "Thr. (mmHg)"),
  columns = c("Etiology", "N", "OutcomeDef",   "Thr"))

COLUMNS_OR_SENS <- list(
  labels  = c("Etiology", "Dataset", "N", "Outcome def.", "Thr. (mmHg)"),
  columns = c("Etiology", "Dataset", "N", "OutcomeDef",   "Thr"))

COLUMNS_CORR <- list(
  labels  = c("Etiology", "N", "Scale", "Method", "Thr. (mmHg)"),
  columns = c("Etiology", "N", "Scale", "Method", "Thr"))

COLUMNS_CORR_SENS <- list(
  labels  = c("Etiology", "Dataset", "N", "Scale", "Method", "Thr. (mmHg)"),
  columns = c("Etiology", "Dataset", "N", "Scale", "Method", "Thr"))

# A column printed only if at least one row has something in it.
# Mortality rows have no outcome definition, and some metrics have no threshold.
column_has_content <- function(values) {
  values <- trimws(as.character(values))
  any(!(values %in% c("", "-", "—", "NA", "NR")))
}

drop_empty_columns <- function(data, choice) {
  keep <- c()
  for (i in seq_along(choice$columns)) {
    if (column_has_content(data[[choice$columns[i]]])) keep <- c(keep, i)
  }
  list(labels = choice$labels[keep], columns = choice$columns[keep])
}


# ---- ordering and fitting --------------------------------------------------

# Rows are sorted before fitting by target first, then diagnosis, then study name.
sort_rows <- function(data, by_dataset = FALSE) {
  target_rank   <- match(data$target, c("CPPopt", "MAPopt"))
  etiology_rank <- match(data$Etiology, ETIOLOGY_ORDER)
  if (by_dataset) {
    data[order(target_rank, etiology_rank, data$Dataset, data$study_id), ]
  } else {
    data[order(target_rank, etiology_rank, data$study_id), ]
  }
}

# Fit one cell, with CPPopt and MAPopt as subgroups.

fit_model <- function(data, effect) {
  present <- c("CPPopt", "MAPopt")
  present <- present[present %in% data$target]
  groups <- factor(unname(TARGET_NAMES[data$target]),
                   levels = unname(TARGET_NAMES[present]))

  metagen(TE = data$yi, seTE = data$sei, studlab = data$slab, data = data,
          sm = if (effect == "OR") "OR" else "ZCOR",
          common = FALSE, random = TRUE,
          method.tau = TAU_METHOD,
          method.random.ci = CI_METHOD,
          method.I2 = I2_METHOD,
          prediction = TRUE,
          subgroup = groups,
          print.subgroup.name = FALSE)
}


# ---- axes clean-up -------------------------------------------------

axis_limits <- function(model) {
  study_low   <- model$TE - 1.96 * model$seTE
  study_high  <- model$TE + 1.96 * model$seTE
  all_limits  <- c(study_low, study_high)

  trimmed <- as.numeric(quantile(all_limits, c(0.05, 0.95)))

  must_show <- c(model$TE, model$lower.random, model$upper.random,
                 model$lower.predict, model$upper.predict)
  must_show <- must_show[is.finite(must_show)]

  low  <- min(trimmed[1], must_show)
  high <- max(trimmed[2], must_show)

  slack     <- 0.15 * (high - low)
  near_miss <- all_limits[all_limits >= low - slack & all_limits <= high + slack]
  low  <- min(low,  near_miss)
  high <- max(high, near_miss)

  padding <- 0.08 * (high - low)
  c(low - padding, high + padding)
}

nice_ticks <- function(limits, effect) {
  low  <- to_plot_scale(limits[1], effect)
  high <- to_plot_scale(limits[2], effect)

  if (effect == "OR") {
    steps <- c(0.05, 0.1, 0.25, 0.5)
  } else {
    steps <- c(0.1, 0.2, 0.25, 0.5)
  }

  for (step in steps) {
    ticks <- seq(floor(low / step) * step, ceiling(high / step) * step, by = step)
    ticks <- round(ticks, 3)
    # an odds ratio cannot be zero or less; a correlation cannot reach 1
    if (effect == "OR") {
      ticks <- ticks[ticks > 0]
    } else {
      ticks <- ticks[ticks > -0.99 & ticks < 0.99]
    }
    if (length(ticks) >= 3 && length(ticks) <= 8) return(ticks)
  }
  signif(seq(low, high, length.out = 5), 2)
}


# ---- drawing forest plots -------------------------------------------------------

draw_forest <- function(model, data, file, choice, effect, width = 12.4) {
  shown <- drop_empty_columns(data, choice)

  ticks     <- nice_ticks(axis_limits(model), effect)
  axis_low  <- min(ticks)
  axis_high <- max(ticks)

  # A subgroup of two studies uses a t distribution with 1 degree of freedom, 
  # so its interval can be enormous. Those are drawn as squares, which get 
  # arrowheads, rather than as diamonds, which do not. The true interval is 
  # still printed on the right.
  sub_low  <- to_plot_scale(model$lower.random.w, effect)
  sub_high <- to_plot_scale(model$upper.random.w, effect)
  sub_low  <- sub_low[is.finite(sub_low)]
  sub_high <- sub_high[is.finite(sub_high)]
  clipped  <- any(sub_low < axis_low) || any(sub_high > axis_high)

  # With one target the subgroup block is the whole analysis, so the separate
  # overall block is switched off and the subgroup block carries the diamond,
  # the prediction interval and the heterogeneity line. That keeps the layout
  # the same as a two-target figure, heading included, with nothing repeated.
  n_targets     <- nlevels(factor(model$subgroup))
  single_target <- n_targets == 1

  if (effect == "OR") {
    effect_label <- "OR"
    x_label      <- "Odds ratio per 1% of monitoring time"
  } else {
    effect_label <- "r"
    x_label      <- "Correlation with worse functional outcome (r)"
  }

  rows_per_target <- if (single_target) 2.6 else 2.4
  height <- 2.3 + 0.30 * (model$k + rows_per_target * n_targets)

  png(paste0(FIG, file), width = width, height = height, units = "in", res = 200)
  forest(model,
         sortvar   = NULL,                  # rows were already sorted by sort_rows()
         leftcols  = c("studlab", shown$columns),
         leftlabs  = c("Study",   shown$labels),
         rightcols = c("effect", "ci", "w.random"),
         rightlabs = c(effect_label, "95% CI", "Weight"),
         xlim = c(axis_low, axis_high),
         at   = ticks,
         type.subgroup.random = if (clipped) "square" else "diamond",
         overall             = !single_target,
         random.subgroup     = TRUE,
         prediction.subgroup = single_target,
         subgroup.hetstat    = !single_target,
         test.subgroup       = !single_target,
         arrow.type = "closed", arrow.length = 0.10,
         prediction = TRUE, print.tau2 = TRUE, print.I2 = TRUE,
         label.left = "better outcome", label.right = "worse outcome",
         smlab = "", xlab = x_label,
         # study squares keep the package default grey; summaries are black
         col.diamond = "black", col.predict = "black",
         digits = 2, digits.se = 3)
  invisible(dev.off())

  if (clipped) {
    cat("wrote", file, "[subtotal clipped]\n")
  } else {
    cat("wrote", file, "\n")
  }
}
