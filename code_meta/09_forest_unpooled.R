# Odds ratios that were not pooled because the cell has fewer than MIN_K
# non-overlapping (primary) studies. Each study's estimate is shown as reported
# on the review's increment, with no pooled estimate, weight or heterogeneity.
#
# Only odds ratios are shown; correlation coefficients below the floor are left out.
#
# One panel per cell. Pressure-level exposures (the absolute optimal pressure and the
# limits of autoregulation, per 1 mm Hg) share one panel per outcome, because they
# come from a single study and are on a different scale from the time metrics.
#
# Writes figures/forest_unpooled_<outcome>_<cell>.png and
# figures/forest_unpooled_panels.csv (letter, title, file), which the supplement
# builder reads for the panel order and titles.

source("config.R")
source("forest_lib.R")

PRESSURE <- c("Popt", "LLA", "ULA")
PRESSURE_LABEL <- c(Popt = "Optimal pressure", LLA = "Lower limit (LLA)",
                    ULA = "Upper limit (ULA)")
OUTCOME_TEXT <- c(mortality = "mortality", poor_outcome = "poor functional outcome")

d <- read.csv(EFFECTS, check.names = FALSE, stringsAsFactors = FALSE)
d <- d[d$effect_type == "OR" & d$poolable == "yes" & d$analysis_set == "primary" &
       d$exposure_metric %in% c(METRICS, PRESSURE), ]
d$yi  <- as.numeric(d$yi)
d$sei <- as.numeric(d$sei)
d$n_total <- suppressWarnings(as.numeric(d$n_total))
d$etiology_short <- ifelse(d$etiology %in% names(ETIOLOGY), ETIOLOGY[d$etiology],
                           substr(d$etiology, 1, 3))
d <- add_label_columns(d)

# time metrics are re-expressed per OR_UNIT % as in every other figure;
# pressure metrics stay per 1 mm Hg
is_time <- d$exposure_metric %in% METRICS
d$yi[is_time]  <- d$yi[is_time]  * OR_UNIT
d$sei[is_time] <- d$sei[is_time] * OR_UNIT
d$cell <- ifelse(is_time, d$exposure_metric, "pressure")
d$Exposure <- ifelse(is_time, "", PRESSURE_LABEL[d$exposure_metric])

panels <- NULL
letter <- 1
for (outcome in c("mortality", "poor_outcome")) {
  for (cell in c(METRICS, "pressure")) {
    rows <- d[d$outcome == outcome & d$cell == cell, ]
    if (nrow(rows) == 0) next
    # a pressure panel is shown when every metric in it is below the floor
    per_metric <- table(rows$exposure_metric)
    if (any(per_metric >= MIN_K)) next

    rows <- rows[order(match(rows$target, c("CPPopt", "MAPopt")),
                       match(rows$exposure_metric, c(METRICS, PRESSURE)),
                       rows$study_id), ]
    model <- metagen(TE = rows$yi, seTE = rows$sei, studlab = rows$slab, data = rows,
                     sm = "OR", common = FALSE, random = FALSE)

    lo <- min(rows$yi - 1.96 * rows$sei, 0)
    hi <- max(rows$yi + 1.96 * rows$sei, 0)
    pad <- 0.08 * (hi - lo)
    ticks <- nice_ticks(c(lo - pad, hi + pad), "OR")

    if (cell == "pressure") {
      choice <- list(labels  = c("Exposure", "Target", "Etiology", "N", "Outcome def."),
                     columns = c("Exposure", "target", "Etiology", "N", "OutcomeDef"))
      x_label <- "Odds ratio per 1 mm Hg"
      what <- paste("Absolute optimal pressure and limits of autoregulation and",
                    OUTCOME_TEXT[outcome], "(odds ratio per 1 mm Hg)")
    } else {
      choice <- list(labels  = c("Target", "Etiology", "N", "Outcome def.", "Thr. (mmHg)"),
                     columns = c("target", "Etiology", "N", "OutcomeDef", "Thr"))
      x_label <- paste("Odds ratio", OR_UNIT_LABEL)
      what <- paste(sub("the target", "the optimal target",
                        sub("% of monitoring time", "Percentage of monitoring time",
                            METRIC_LABEL[cell])), "and", OUTCOME_TEXT[outcome],
                    paste0("(odds ratio ", OR_UNIT_LABEL, ")"))
    }
    shown <- drop_empty_columns(rows, choice)
    file <- paste0("forest_unpooled_", outcome, "_", cell, ".png")

    png(paste0(FIG, file), width = 11.4, height = 1.9 + 0.30 * nrow(rows),
        units = "in", res = 200)
    forest(model,
           sortvar   = NULL,
           leftcols  = c("studlab", shown$columns),
           leftlabs  = c("Study",   shown$labels),
           rightcols = c("effect", "ci"),
           rightlabs = c("OR", "95% CI"),
           xlim = c(min(ticks), max(ticks)), at = ticks,
           overall = FALSE, hetstat = FALSE, overall.hetstat = FALSE,
           test.overall = FALSE, weight.study = "same",
           arrow.type = "closed", arrow.length = 0.10,
           label.left = "better outcome", label.right = "worse outcome",
           smlab = "", xlab = x_label, digits = 2)
    invisible(dev.off())
    cat("wrote", file, "-", nrow(rows), "estimates\n")

    panels <- rbind(panels, data.frame(letter = LETTERS[letter], title = what,
                                       file = file, k = nrow(rows)))
    letter <- letter + 1
  }
}
write.csv(panels, paste0(FIG, "forest_unpooled_panels.csv"), row.names = FALSE)
cat("UNPOOLED OR FORESTS DONE:", nrow(panels), "panels\n")
