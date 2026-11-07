# MICE-based temporal-validation metrics for the frozen M1 model, C1-unflagged headline stratum,
# with design-based 95% CIs (Protocol §9.4-9.6, Amendment v1.1 §4-5).
#
# Completes the items that 08/11/13 left complete-case only: calibration slope/intercept and the
# survey-weighted C-statistic under MICE; CIs for O/E, every subgroup O/E and NB, and the full
# 2-20% decision curve for all strategies; calibration slopes under the >=50-event rule and the
# pre-specified collapsing hierarchy; and S1/S3/S9 re-run on the MI predictions.
#
# Inputs: data/processed/val_headline_mi_preds.rds (from 11), the validation design from 05.
# Outputs: results/tables/mice_validation_metrics.csv, decision_curve_mice.csv,
#          subgroup_oe_mice_ci.csv, sensitivity_mice.csv
#
# NOTE: this script was drafted alongside a Python re-implementation of the same logic that was
# used to fill the manuscript (the Python version reproduced 11's O/E = 1.040 and NB@5% = 33.3/1000
# exactly). It has not been executed in R; run it and compare against the manuscript tables before
# submission. Expected: O/E 1.04 (0.92-1.16), slope 1.06 (0.84-1.29), weighted C 0.70 (0.64-0.75),
# model NB@5% 33.3 (24.9-43.3)/1000; small differences in CI limits are expected from bootstrap
# randomness (seed below) but point estimates should match to the reported precision.

source("R/00_config.R")
source("R/05_survey_design.R")
library(survey); library(dplyr)

set.seed(20260921)
R_BOOT <- 500   # pre-specified (§7.2); not SCALE$boot_R, so that the MI metrics use the protocol number

mi <- readRDS(file.path(PROCESSED_DIR, "val_headline_mi_preds.rds"))   # n = 3,064 with m1_mi_pred
stopifnot(all(mi$SEQN %in% val_design$variables$SEQN))

# Attach the MI prediction to the FULL validation design, then subset (domain estimation, §7.2)
val_design <- update(val_design,
                     m1_mi_pred = mi$m1_mi_pred[match(SEQN, mi$SEQN)],
                     ckd_s1 = (egfr_race_free < EGFR_THRESHOLD) | (acr_uncorrected >= ACR_THRESHOLD))
boot_full <- as.svrepdesign(val_design, type = "mrbbootstrap", replicates = R_BOOT)
hl <- subset(boot_full, in_analytic_domain & !c1_flag & !is.na(m1_mi_pred))
cat(sprintf("[15] headline MI rows: n=%d, events=%d\n", nrow(hl$variables), sum(hl$variables$ckd_case)))

# ---- statistic functions (all take data + weights so they work inside withReplicates) ----------
wnb <- function(treat, y, w, pt) (sum(w[treat & y]) - sum(w[treat & !y]) * pt / (1 - pt)) / sum(w)
wlogit_slope <- function(y, lp, w) {
  f <- suppressWarnings(glm(y ~ lp, weights = w, family = quasibinomial()))
  unname(coef(f)[2])
}
wlogit_int <- function(y, lp, w) {
  f <- suppressWarnings(glm(y ~ 1 + offset(lp), weights = w, family = quasibinomial()))
  unname(coef(f)[1])
}
wcstat <- function(y, p, w) {   # weighted concordance, ties count 1/2
  o <- order(p); p <- p[o]; y <- y[o]; w <- w[o]
  u <- match(p, unique(p))
  wc <- tapply(w * y, u, sum); wn <- tapply(w * !y, u, sum)
  cum_wn <- cumsum(wn) - wn
  sum(wc * (cum_wn + 0.5 * wn)) / (sum(wc) * sum(wn))
}
all_stats <- function(d, w, thresholds = DCA_THRESHOLDS) {
  y <- as.logical(d$ckd_case); p <- d$m1_mi_pred
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  out <- c(obs = sum(w * y) / sum(w), exp = sum(w * p) / sum(w))
  out["oe"] <- out["obs"] / out["exp"]
  out["cal_int"] <- wlogit_int(y, lp, w)
  out["cal_slope"] <- wlogit_slope(y, lp, w)
  out["c"] <- wcstat(y, p, w)
  for (pt in thresholds) {
    k <- sprintf("%.2f", pt)
    out[paste0("nb_model_", k)] <- wnb(p >= pt, y, w, pt)
    out[paste0("nb_all_", k)]   <- wnb(rep(TRUE, length(y)), y, w, pt)
    out[paste0("nb_c2_", k)]    <- wnb(d$c2_flag, y, w, pt)
  }
  out["nb_diff_model_minus_all_0.05"] <- out["nb_model_0.05"] - out["nb_all_0.05"]
  out["nb_diff_model_minus_c2_0.05"]  <- out["nb_model_0.05"] - out["nb_c2_0.05"]
  out
}
summarise_rep <- function(est) {
  reps <- attr(est, "replicates")
  data.frame(metric = names(coef(est)), est = as.numeric(coef(est)),
             lo = apply(reps, 2, quantile, 0.025), hi = apply(reps, 2, quantile, 0.975),
             se = as.numeric(SE(est)), row.names = NULL)
}

# ---- primary ------------------------------------------------------------------------------------
prim <- withReplicates(hl, function(w, d) all_stats(d, w), return.replicates = TRUE)
prim_tab <- summarise_rep(prim)
cat("\n=== PRIMARY (MICE, WTPH2YR, C1-unflagged) ===\n")
print(prim_tab[prim_tab$metric %in% c("obs", "exp", "oe", "cal_int", "cal_slope", "c",
                                      "nb_model_0.05", "nb_all_0.05", "nb_c2_0.05", "nb_model_0.10",
                                      "nb_diff_model_minus_all_0.05", "nb_diff_model_minus_c2_0.05"), ], digits = 4)
write.csv(prim_tab, file.path(TABLES_DIR, "mice_validation_metrics.csv"), row.names = FALSE)

dca <- data.frame(threshold = DCA_THRESHOLDS)
for (s in c("model", "all", "c2")) {
  k <- sprintf("nb_%s_%.2f", s, DCA_THRESHOLDS)
  dca[[paste0("nb_", s)]] <- prim_tab$est[match(k, prim_tab$metric)]
  dca[[paste0("nb_", s, "_lo")]] <- prim_tab$lo[match(k, prim_tab$metric)]
  dca[[paste0("nb_", s, "_hi")]] <- prim_tab$hi[match(k, prim_tab$metric)]
}
write.csv(dca, file.path(TABLES_DIR, "decision_curve_mice.csv"), row.names = FALSE)

# ---- subgroups: O/E, NB@5%, slope where >=50 events (full granularity + collapsing hierarchy) ---
RIDRETH3_LABELS <- c(`1`="Mexican American", `2`="Other Hispanic", `3`="Non-Hispanic White",
                      `4`="Non-Hispanic Black", `6`="Non-Hispanic Asian", `7`="Other/Multiracial")
v <- hl$variables
groups <- list()
for (code in names(RIDRETH3_LABELS)) groups[[paste0("race|", RIDRETH3_LABELS[code])]] <- v$race_eth == as.numeric(code)
groups[["race_collapsed|Hispanic (combined)"]] <- v$race_eth %in% c(1, 2)
inc <- cut(v$income_pir, breaks = INCOME_BREAKS, labels = INCOME_LABELS)
for (lv in INCOME_LABELS) groups[[paste0("income|", lv)]] <- !is.na(inc) & inc == lv
groups[["income_collapsed|<2.0x FPL"]]  <- !is.na(inc) & inc %in% INCOME_LABELS[1:2]
groups[["income_collapsed|>=2.0x FPL"]] <- !is.na(inc) & inc %in% INCOME_LABELS[3:4]
groups[["insurance|Insured"]]   <- !is.na(v$hiq011) & v$hiq011 == 1
groups[["insurance|Uninsured"]] <- !is.na(v$hiq011) & v$hiq011 == 2

sub_rows <- list()
for (g in names(groups)) {
  m <- groups[[g]]
  ev <- sum(v$ckd_case[m]); nn <- sum(m)
  f <- function(w, d) {
    y <- as.logical(d$ckd_case[m]); p <- d$m1_mi_pred[m]; ww <- w[m]
    out <- c(oe = (sum(ww * y) / sum(ww)) / (sum(ww * p) / sum(ww)),
             nb5 = wnb(p >= DCA_PRIMARY, y, ww, DCA_PRIMARY))
    if (ev >= MIN_EVENTS_FOR_SLOPE) out["slope"] <- wlogit_slope(y, qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6)), ww)
    out
  }
  est <- withReplicates(hl, f, return.replicates = TRUE)
  t <- summarise_rep(est)
  parts <- strsplit(g, "|", fixed = TRUE)[[1]]
  sub_rows[[g]] <- data.frame(dimension = parts[1], level = parts[2], n = nn, events = ev,
                              oe = t$est[t$metric == "oe"], oe_lo = t$lo[t$metric == "oe"], oe_hi = t$hi[t$metric == "oe"],
                              nb_5pct = t$est[t$metric == "nb5"], nb_lo = t$lo[t$metric == "nb5"], nb_hi = t$hi[t$metric == "nb5"],
                              cal_slope = if (ev >= MIN_EVENTS_FOR_SLOPE) t$est[t$metric == "slope"] else NA,
                              slope_lo = if (ev >= MIN_EVENTS_FOR_SLOPE) t$lo[t$metric == "slope"] else NA,
                              slope_hi = if (ev >= MIN_EVENTS_FOR_SLOPE) t$hi[t$metric == "slope"] else NA)
}
sub_tab <- do.call(rbind, sub_rows); rownames(sub_tab) <- NULL
cat("\n=== SUBGROUPS (MICE) ===\n"); print(sub_tab, digits = 3)
write.csv(sub_tab, file.path(TABLES_DIR, "subgroup_oe_mice_ci.csv"), row.names = FALSE)

# weighted SD of log(O/E) across race/ethnicity (weights = unweighted n, as in 13) + NHB/NHW ratio, with intervals
race_masks <- lapply(names(RIDRETH3_LABELS), function(code) v$race_eth == as.numeric(code))
race_summary <- withReplicates(hl, function(w, d) {
  oes <- sapply(race_masks, function(m) { ww <- w[m]; y <- as.logical(d$ckd_case[m]); p <- d$m1_mi_pred[m]
                                            (sum(ww * y) / sum(ww)) / (sum(ww * p) / sum(ww)) })
  ns <- sapply(race_masks, sum); logs <- log(oes); mu <- sum(ns * logs) / sum(ns)
  c(sd_log_oe = sqrt(sum(ns * (logs - mu)^2) / sum(ns)), max_abs_diff = max(oes) - min(oes),
    nhb_over_nhw = oes[4] / oes[3])
}, return.replicates = TRUE)
cat("\n=== RACE/ETHNICITY SUMMARY ===\n"); print(summarise_rep(race_summary), digits = 3)

# ---- sensitivity analyses on the MI predictions: S1 (no bridging), S3 (WTMEC2YR), S9 (C2-unflagged) ----
sens_stats <- function(dsn, label, outcome = "ckd_case") {
  est <- withReplicates(dsn, function(w, d) {
    d$ckd_case <- d[[outcome]]
    all_stats(d, w, thresholds = DCA_PRIMARY)[c("obs", "exp", "oe", "cal_int", "cal_slope", "c", "nb_model_0.05", "nb_all_0.05", "nb_c2_0.05")]
  }, return.replicates = TRUE)
  t <- summarise_rep(est); t$analysis <- label; t$n <- nrow(dsn$variables); t$events <- sum(dsn$variables[[outcome]])
  t
}
s1 <- sens_stats(hl, "S1_no_bridging", outcome = "ckd_s1")
val_design_mec2 <- update(val_design_mec, m1_mi_pred = mi$m1_mi_pred[match(SEQN, mi$SEQN)])
boot_mec <- as.svrepdesign(val_design_mec2, type = "mrbbootstrap", replicates = R_BOOT)
s3 <- sens_stats(subset(boot_mec, in_analytic_domain & !c1_flag & !is.na(m1_mi_pred)), "S3_wtmec2yr")
s9 <- sens_stats(subset(hl, !c2_flag), "S9_c2_unflagged")
sens <- rbind(s1, s3, s9)
cat("\n=== SENSITIVITY (MICE) ===\n"); print(sens[sens$metric %in% c("oe", "cal_slope", "c", "nb_model_0.05", "nb_all_0.05", "nb_c2_0.05"), ], digits = 4)
write.csv(sens, file.path(TABLES_DIR, "sensitivity_mice.csv"), row.names = FALSE)
cat("[15] done.\n")
