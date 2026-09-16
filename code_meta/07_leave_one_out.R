# Leave-one-out sensitivity analysis for the primary odds-ratio forest plots.
#
# Each of the four primary figures is refitted k+1 times: once with every study,
# then once with each study omitted in turn. 
#
suppressMessages(library(meta))
source("config.R")

pool <- function(rows) {
  metagen(TE = rows$yi, seTE = rows$sei, studlab = rows$study_id, data = rows,
          sm = "OR", common = FALSE, random = TRUE,
          method.tau = TAU_METHOD, method.random.ci = CI_METHOD,
          method.I2 = I2_METHOD, prediction = TRUE)
}

d <- read_cells("OR", "primary")
out <- NULL

for (outcome in c("mortality", "poor_outcome")) {
  for (metric in METRICS) {

    rows <- d[d$outcome == outcome & d$exposure_metric == metric, ]
    if (nrow(rows) < MIN_K) next

    add <- function(m, omitted, k) data.frame(
      outcome = outcome, metric = metric, omitted = omitted, k = k,
      est = exp(m$TE.random), ci_low = exp(m$lower.random),
      ci_high = exp(m$upper.random), pval = m$pval.random,
      I2 = 100 * m$I2, tau2 = m$tau2, stringsAsFactors = FALSE)

    full <- pool(rows)
    out <- rbind(out, add(full, "(all studies)", nrow(rows)))
    cat(sprintf("\n%s / %s  full: OR %.3f (%.3f-%.3f) k=%d I2=%.0f%%\n",
                outcome, metric, exp(full$TE.random), exp(full$lower.random),
                exp(full$upper.random), full$k, 100 * full$I2))

    for (i in seq_len(nrow(rows))) {
      sub <- rows[-i, , drop = FALSE]
      m <- pool(sub)
      out <- rbind(out, add(m, rows$study_id[i], nrow(sub)))
      cat(sprintf("   omit %-24s OR %.3f (%.3f-%.3f) p=%.3f I2=%.0f%%\n",
                  rows$study_id[i], exp(m$TE.random), exp(m$lower.random),
                  exp(m$upper.random), m$pval.random, 100 * m$I2))
    }
  }
}

num <- c("est","ci_low","ci_high","pval","I2","tau2")
out[num] <- lapply(out[num], function(x) round(as.numeric(x), 4))
write.csv(out, paste0(META_OUTPUTS, "leave_one_out.csv"), row.names = FALSE)
cat("\nwrote leave_one_out.csv:", nrow(out), "rows\n")
