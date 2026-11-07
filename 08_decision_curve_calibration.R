# Applies the FROZEN M1 model to the temporal validation set (no recalibration -- §9.4 "frozen
# model is the primary result"), then computes calibration and decision-curve analysis restricted
# to the C1-unflagged headline stratum (§6.3, §9.5), all as weighted, design-based quantities.
#
# Scope note (draft run): the full curve (2%-20%) is reported as weighted point estimates for every
# threshold; bootstrap CIs (Rao-Wu rescaled bootstrap, §7.2) are computed at the primary (5%) and
# secondary (10%) thresholds only, to keep this pass tractable. At RUN_SCALE=full, extend the CI
# loop to every threshold in DCA_THRESHOLDS.

source("R/00_config.R")
source("R/05_survey_design.R")
library(survey)
library(splines)

m1_fit <- readRDS(file.path(PROCESSED_DIR, "m1_fit.rds"))

val_df <- val_design$variables
val_df <- val_df[val_df$in_analytic_domain & !val_df$c1_flag, ]  # headline stratum, §6.3
cat(sprintf("[08] validation C1-unflagged headline stratum: n=%d\n", nrow(val_df)))

model_predictor_cols <- c("age","female","income_pir","education","hiq011","smq020",
                           "bmxbmi","bmxwaist","diq010","bpq020","cvd_any","race_eth")
complete_rows <- stats::complete.cases(val_df[, model_predictor_cols])
val_df$m1_pred <- NA_real_
val_df$m1_pred[complete_rows] <- as.numeric(predict(m1_fit, newdata = val_df[complete_rows, ], type = "response"))
cat(sprintf("[08] predicted for %d/%d (complete predictors)\n", sum(complete_rows), nrow(val_df)))

val_df <- val_df[complete_rows, ]  # analysis set for this draft pass; missingness handled by MICE at RUN_SCALE=full

# Re-attach to a survey design restricted to this exact headline analytic set, preserving design info
headline_design <- subset(val_design, in_analytic_domain & !c1_flag)
headline_design <- headline_design[complete_rows, ]

# --- Calibration (Protocol §9.4) --------------------------------------------------------------
calib_df <- headline_design$variables
stopifnot(nrow(calib_df) == nrow(val_df))
calib_df$m1_pred <- val_df$m1_pred

w <- weights(headline_design)
obs_w   <- sum(w * calib_df$ckd_case) / sum(w)
exp_w   <- sum(w * calib_df$m1_pred)  / sum(w)
oe_ratio <- obs_w / exp_w

# Calibration slope/intercept via survey-weighted logistic regression of outcome on logit(pred) (§9.4)
calib_df$logit_pred <- qlogis(pmin(pmax(calib_df$m1_pred, 1e-6), 1 - 1e-6))
cal_design <- update(headline_design, logit_pred = calib_df$logit_pred, ckd_case_num = as.numeric(calib_df$ckd_case))
cal_slope_fit <- svyglm(ckd_case_num ~ logit_pred, design = cal_design, family = quasibinomial())
cal_intercept_fit <- svyglm(ckd_case_num ~ offset(logit_pred), design = cal_design, family = quasibinomial())

cat("\n=== CALIBRATION (M1, frozen, temporal validation, C1-unflagged headline stratum) ===\n")
cat(sprintf("Weighted observed prevalence: %.4f\n", obs_w))
cat(sprintf("Weighted mean predicted risk: %.4f\n", exp_w))
cat(sprintf("O/E ratio: %.3f\n", oe_ratio))
cat("Calibration slope (target 1.0):\n"); print(summary(cal_slope_fit)$coefficients)
cat("Calibration intercept (calibration-in-the-large, target 0.0):\n"); print(summary(cal_intercept_fit)$coefficients)

# Weighted C-statistic (secondary, §9.4)
c_stat <- tryCatch(svyauc <- {
  library(pROC)
  # design-unweighted AUC as a fallback if svy-native AUC isn't available; weighted approx via
  # case/non-case weighted rank comparison (Newcombe-style) would be more exact -- flagged for full run.
  as.numeric(pROC::auc(pROC::roc(calib_df$ckd_case, calib_df$m1_pred, quiet = TRUE)))
}, error = function(e) NA)
cat(sprintf("C-statistic (unweighted approximation this draft pass; design-weighted version at RUN_SCALE=full): %.3f\n", c_stat))

# --- Decision curve analysis (Protocol §9.5, Amendment §1) -----------------------------------
# Strategies in the C1-unflagged headline stratum: model (M1), test-all, test-none, C2.
# C1 == test-none by stratum construction (Amendment §1); SCORED unavailable (see 04 script).
weighted_net_benefit <- function(treat_indicator, outcome, w, pt) {
  tp <- sum(w[treat_indicator & outcome])
  fp <- sum(w[treat_indicator & !outcome])
  n  <- sum(w)
  (tp/n) - (fp/n) * (pt/(1-pt))
}

dca_rows <- list()
for (pt in DCA_THRESHOLDS) {
  nb_model <- weighted_net_benefit(calib_df$m1_pred >= pt, calib_df$ckd_case, w, pt)
  nb_all   <- weighted_net_benefit(rep(TRUE, nrow(calib_df)), calib_df$ckd_case, w, pt)
  nb_none  <- 0
  nb_c2    <- weighted_net_benefit(calib_df$c2_flag, calib_df$ckd_case, w, pt)
  dca_rows[[length(dca_rows)+1]] <- data.frame(threshold = pt, model = nb_model, test_all = nb_all,
                                                test_none = nb_none, c1 = 0, c2 = nb_c2)
}
dca <- do.call(rbind, dca_rows)
cat("\n=== DECISION CURVE (net benefit per person), C1-unflagged headline stratum, validation ===\n")
print(round(dca, 5))
write.csv(dca, file.path(TABLES_DIR, "decision_curve.csv"), row.names = FALSE)

# --- Net benefit in clinically legible units at p_t=5% (Protocol §6.4) -----------------------
pt5 <- DCA_PRIMARY
nb5 <- dca[dca$threshold == pt5, ]
cases_per_1000 <- function(nb) nb * 1000
cat(sprintf("\nAt p_t=5%%: model NB=%.5f -> %.1f cases detected per 1,000 screened (net of harm-weighted false positives)\n",
            nb5$model, cases_per_1000(nb5$model)))
cat(sprintf("Test-all NB=%.5f -> %.1f/1000; C2 NB=%.5f -> %.1f/1000; test-none/C1 NB=0\n",
            nb5$test_all, cases_per_1000(nb5$test_all), nb5$c2, cases_per_1000(nb5$c2)))

# --- Bootstrap CI at primary (5%) and secondary (10%) thresholds only, this draft pass ---------
cat(sprintf("\n[08] Building Rao-Wu rescaled bootstrap replicate design (R=%d) for CI at 5%% and 10%%...\n", BOOT_REPLICATES))
boot_design <- tryCatch(build_boot_design(headline_design, R = BOOT_REPLICATES), error = function(e) { cat("BOOT DESIGN FAILED:", conditionMessage(e), "\n"); NULL })

if (!is.null(boot_design)) {
  nb_statistic <- function(pt) {
    function(data, w) {
      treat <- data$m1_pred >= pt
      tp <- sum(w[treat & data$ckd_case]); fp <- sum(w[treat & !data$ckd_case]); n <- sum(w)
      (tp/n) - (fp/n) * (pt/(1-pt))
    }
  }
  for (pt in c(DCA_PRIMARY, DCA_SECONDARY)) {
    est <- withReplicates(boot_design, function(w, data) nb_statistic(pt)(calib_df, w))
    cat(sprintf("Model NB at p_t=%.0f%%: %.5f (bootstrap SE %.5f, 95%% CI [%.5f, %.5f])\n",
                pt*100, coef(est), SE(est), coef(est) - 1.96*SE(est), coef(est) + 1.96*SE(est)))
  }
}

save(oe_ratio, obs_w, exp_w, cal_slope_fit, cal_intercept_fit, c_stat, dca,
     file = file.path(TABLES_DIR, "calibration_decision_curve.RData"))
