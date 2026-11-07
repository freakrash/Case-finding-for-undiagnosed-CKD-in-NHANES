# Sensitivity analyses S1, S3, S9 (Protocol §10) -- cheap to add given 08's infrastructure.
# S1: validation WITHOUT the urinary albumin bridging correction (naive analysis).
# S3: validation weighted by WTMEC2YR instead of WTPH2YR.
# S9: C2 (instead of C1) as the "unflagged" stratum definition, testing whether the headline
#     result depends on which rule defines the deployment population.
# S4-S8 (persistence adjustment, widened diagnosed definition, alternate thresholds, secondary BP
# model, formal complete-case label) are NOT implemented here -- S8 is simply what 07/08 already are
# (complete-case, pre-MI); S4/S5/S6/S7 need additional data/modelling choices not yet made and are
# left for the study team, not silently faked.

source("R/00_config.R")
source("R/05_survey_design.R")
library(survey); library(splines)

m1_fit <- readRDS(file.path(PROCESSED_DIR, "m1_fit.rds"))
model_predictor_cols <- c("age","female","income_pir","education","hiq011","smq020",
                           "bmxbmi","bmxwaist","diq010","bpq020","cvd_any","race_eth")

weighted_net_benefit <- function(treat, outcome, w, pt) {
  tp <- sum(w[treat & outcome]); fp <- sum(w[treat & !outcome]); n <- sum(w)
  (tp/n) - (fp/n) * (pt/(1-pt))
}
calib_and_dca <- function(df, w, label) {
  cc <- stats::complete.cases(df[, model_predictor_cols])
  df <- df[cc, ]; w <- w[cc]
  pred <- as.numeric(predict(m1_fit, newdata = df, type = "response"))
  obs <- sum(w * df$ckd_case) / sum(w); exp_ <- sum(w * pred) / sum(w)
  nb5 <- weighted_net_benefit(pred >= DCA_PRIMARY, df$ckd_case, w, DCA_PRIMARY)
  cat(sprintf("[%s] n=%d, O/E=%.3f (obs=%.4f, exp=%.4f), model NB@5%%=%.5f (%.1f/1000)\n",
              label, nrow(df), obs/exp_, obs, exp_, nb5, nb5*1000))
  list(n = nrow(df), oe = obs/exp_, nb5 = nb5)
}

val_df_full <- val_design$variables
headline_mask <- val_df_full$in_analytic_domain & !val_df_full$c1_flag

cat("=== S1: validation WITHOUT bridging correction (naive uncorrected albumin) ===\n")
s1_df <- val_df_full[headline_mask, ]
s1_df$ckd_case <- (s1_df$egfr_race_free < EGFR_THRESHOLD) | (s1_df$acr_uncorrected >= ACR_THRESHOLD)
s1_result <- calib_and_dca(s1_df, s1_df$WTPH2YR, "S1 uncorrected")
cat("   (primary/corrected result for comparison, from 08: O/E=1.08, NB@5%=32.6/1000)\n\n")

cat("=== S3: validation weighted by WTMEC2YR instead of WTPH2YR ===\n")
s3_df <- val_df_full[headline_mask, ]
s3_result <- calib_and_dca(s3_df, s3_df$wt_mec, "S3 WTMEC2YR-weighted")
cat("\n")

cat("=== S9: C2 (not C1) defines the 'unflagged' deployment stratum ===\n")
s9_mask <- val_df_full$in_analytic_domain & !val_df_full$c2_flag
s9_df <- val_df_full[s9_mask, ]
s9_result <- calib_and_dca(s9_df, s9_df$WTPH2YR, "S9 C2-unflagged stratum")
cat(sprintf("   (C2-unflagged stratum is %d of %d analytic-domain validation rows, vs %d for C1-unflagged)\n",
            sum(s9_mask), sum(val_df_full$in_analytic_domain), sum(headline_mask)))

save(s1_result, s3_result, s9_result, file = file.path(TABLES_DIR, "sensitivity_s1_s3_s9.RData"))
