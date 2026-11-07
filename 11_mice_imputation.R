# Multiple imputation by chained equations (Protocol §9.9): m imputations, imputation model
# includes the outcome, design strata, PSU identifier, and sampling weight; performed SEPARATELY
# within development and validation sets so no information crosses the temporal boundary.
# Participants missing the outcome components or KIQ022 are excluded rather than imputed (already
# enforced by in_analytic_domain upstream). Model refit on each imputed development set and pooled
# via Rubin's rules (mitools::MIcombine, survey-aware). Validation-side predictions: the pooled
# coefficient vector is applied to each imputed validation set and the m_val predicted risks per
# participant are averaged -- a standard, defensible way to combine an external model's predictions
# across imputations; this does not additionally propagate imputation uncertainty into the
# decision-curve bootstrap CI (documented simplification, not hidden).

source("R/00_config.R")
source("R/05_survey_design.R")
library(survey); library(splines); library(mice); library(mitools)

dev_all <- dev_design$variables
dev_dom <- dev_all[dev_all$in_analytic_domain, ]
val_all <- val_design$variables
val_headline <- val_all[val_all$in_analytic_domain & !val_all$c1_flag, ]

impute_vars <- c("ckd_case","age","female","income_pir","education","hiq011","smq020",
                  "bmxbmi","bmxwaist","diq010","bpq020","cvd_any","race_eth",
                  "combined_stratum","combined_psu")

# Canonical factor levels from the FULL analytic domain (dev + validation), so a level absent from
# the small validation headline stratum (n=3,066) doesn't silently produce a different model.matrix
# structure than development's -- this is what caused a real "contrasts can be applied only to
# factors with >=2 levels" failure on the first attempt (some levels only appear in one set).
factor_cols <- c("education","hiq011","smq020","diq010","bpq020","race_eth","combined_stratum","combined_psu")
canonical_levels <- lapply(factor_cols, function(v) sort(unique(c(dev_dom[[v]], val_headline[[v]]))))
names(canonical_levels) <- factor_cols

impute_one <- function(df, m, label) {
  d <- df[, impute_vars]
  d$ckd_case <- as.numeric(d$ckd_case)
  for (v in c("education","hiq011","smq020","diq010","bpq020","race_eth","combined_stratum","combined_psu"))
    d[[v]] <- factor(d[[v]], levels = canonical_levels[[v]])
  cat(sprintf("[11] imputing %s: n=%d, m=%d, vars with missingness: %s\n",
              label, nrow(d), m, paste(names(which(colSums(is.na(d)) > 0)), collapse=", ")))
  meth <- mice::make.method(d)
  # design strata/PSU/weight are auxiliary predictors, not themselves imputed (weight has no NAs here)
  meth[c("combined_stratum","combined_psu")] <- ""
  imp <- mice::mice(d, m = m, method = meth, printFlag = FALSE, seed = 20260921)
  imp
}

t0 <- Sys.time()
dev_imp <- impute_one(dev_dom, SCALE$mi_m, "development")
cat(sprintf("[11] development imputation done in %.1f sec\n", as.numeric(Sys.time()-t0, units="secs")))

t0 <- Sys.time()
val_imp <- impute_one(val_headline, SCALE$mi_m, "validation (C1-unflagged headline stratum)")
cat(sprintf("[11] validation imputation done in %.1f sec\n", as.numeric(Sys.time()-t0, units="secs")))

m1_formula <- ckd_case ~ ns(age, df = 3) + female + income_pir + education +
  hiq011 + smq020 + ns(bmxbmi, df = 3) + bmxwaist + diq010 + bpq020 + cvd_any + race_eth

# --- Fit M1 on each imputed development set, using the ORIGINAL survey weights (design vars are
# complete, only re-attach the imputed predictor columns) -----------------------------------------
dev_w <- dev_dom$wt_dev
fits <- vector("list", SCALE$mi_m)
for (i in seq_len(SCALE$mi_m)) {
  di <- mice::complete(dev_imp, i)
  di$wt_dev <- dev_w
  di$combined_stratum <- dev_dom$combined_stratum
  di$combined_psu <- dev_dom$combined_psu
  des_i <- svydesign(id = ~combined_psu, strata = ~combined_stratum, weights = ~wt_dev, nest = TRUE, data = di)
  fits[[i]] <- svyglm(m1_formula, design = des_i, family = quasibinomial())
}
pooled_fit <- mitools::MIcombine(fits)
cat("\n=== M1, pooled across m=", SCALE$mi_m, " imputations (Rubin's rules) ===\n", sep="")
print(summary(pooled_fit))
saveRDS(pooled_fit, file.path(PROCESSED_DIR, "m1_mi_pooled.rds"))
saveRDS(fits, file.path(PROCESSED_DIR, "m1_mi_fits_checkpoint.rds"))
saveRDS(val_imp, file.path(PROCESSED_DIR, "val_imp_checkpoint.rds"))
saveRDS(val_headline, file.path(PROCESSED_DIR, "val_headline_checkpoint.rds"))
saveRDS(m1_formula, file.path(PROCESSED_DIR, "m1_formula_checkpoint.rds"))

# --- Predict on each imputed validation set: pair imputation i's dev-fitted model (fits[[i]]) with
# imputation i's validation data, average the m predicted probabilities per participant. This is the
# standard "predict-then-pool" MI approach, and -- critically -- using predict() on the ORIGINAL
# fitted svyglm object (rather than rebuilding X with model.matrix() independently on validation data)
# is required for correctness here: ns() computes spline knots from whatever data is passed to it, so
# an independently-rebuilt model.matrix() on validation data would silently use DIFFERENT knots than
# the ones the model was fit with, producing wrong (not just missing) predictions. predict.svyglm
# correctly reuses the training-data knots stored in the fitted model's terms object.
# mice also occasionally leaves a handful of cells unimputed in a rare category (observed: 2/3066
# rows on bpq020 here) -- handled by complete.cases + index-based assignment, not raw recycling.
val_preds <- matrix(NA_real_, nrow = nrow(val_headline), ncol = SCALE$mi_m)
for (i in seq_len(SCALE$mi_m)) {
  vi <- mice::complete(val_imp, i)
  resid_na <- !stats::complete.cases(vi[, all.vars(m1_formula)])
  if (any(resid_na)) cat(sprintf("[11] imputation %d: %d row(s) with a residual NA after MICE, excluded from this imputation's prediction\n", i, sum(resid_na)))
  val_preds[!resid_na, i] <- as.numeric(predict(fits[[i]], newdata = vi[!resid_na, , drop = FALSE], type = "response"))
}
val_headline$m1_mi_pred <- rowMeans(val_preds, na.rm = TRUE)
# rowMeans(na.rm=TRUE) on an all-NA row returns NaN (not NA) -- a handful of rows missing a
# prediction in every one of the m imputations (the same ~2 rows mice couldn't fully impute).
# Downstream summaries (min/max/mean/sum below) do NOT use na.rm, deliberately: dropping these rows
# explicitly here, once, is safer than silently ignoring NaN/NA at every later step (which bit us
# once already -- min()/max() without na.rm on 2 leftover NA rows had silently NA'd every subsequent
# summary statistic, even though 3064/3066 predictions were actually fine).
n_dropped <- sum(is.nan(val_headline$m1_mi_pred))
if (n_dropped > 0) cat(sprintf("[11] dropping %d row(s) with no valid prediction in any imputation\n", n_dropped))
val_headline <- val_headline[!is.nan(val_headline$m1_mi_pred), ]

cat(sprintf("\n[11] MI-averaged validation predicted risk: range [%.4f, %.4f], mean %.4f (n=%d)\n",
            min(val_headline$m1_mi_pred), max(val_headline$m1_mi_pred), mean(val_headline$m1_mi_pred),
            nrow(val_headline)))

# --- Recompute calibration + primary-threshold decision curve with the MI-based predictions -------
w <- val_headline$WTPH2YR
obs_w <- sum(w * val_headline$ckd_case) / sum(w)
exp_w <- sum(w * val_headline$m1_mi_pred) / sum(w)
cat(sprintf("[11-MI] Weighted observed=%.4f, predicted=%.4f, O/E=%.3f (complete-case version was O/E=1.08)\n",
            obs_w, exp_w, obs_w/exp_w))

weighted_net_benefit <- function(treat, outcome, w, pt) {
  tp <- sum(w[treat & outcome]); fp <- sum(w[treat & !outcome]); n <- sum(w)
  (tp/n) - (fp/n) * (pt/(1-pt))
}
nb5_mi <- weighted_net_benefit(val_headline$m1_mi_pred >= DCA_PRIMARY, val_headline$ckd_case, w, DCA_PRIMARY)
cat(sprintf("[11-MI] Model net benefit at p_t=5%% (MI-based, n=%d, no longer complete-case-restricted): %.5f -> %.1f/1000\n",
            nrow(val_headline), nb5_mi, nb5_mi*1000))

saveRDS(val_headline, file.path(PROCESSED_DIR, "val_headline_mi_preds.rds"))
save(obs_w, exp_w, nb5_mi, file = file.path(TABLES_DIR, "mi_results.RData"))
