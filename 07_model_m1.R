# Fits M1 (Protocol §9.2): survey-weighted logistic regression, restricted cubic splines (4 knots)
# for age and BMI, on the DEVELOPMENT design. M1 is the pre-specified fallback under the parsimony
# rule and the model reported in full even if M2/M3 are selected (§9.2), so it is fit first.
# Predictors restricted to non-laboratory information only (§9.1) -- BP excluded from the primary
# set per §7.4 (device-change discontinuity; BP only enters the secondary S7 model).

source("R/00_config.R")
source("R/05_survey_design.R")
library(survey)
library(splines)

# --- Predictor set (Protocol §9.1), limited to items harmonized and available across all cycles ---
# income_pir, education, hiq011 (insurance), smq020 (smoking), bmxbmi, bmxwaist, diabetes/HTN/CVD
# self-report, race_eth. (Some §9.1 candidates -- usual source of care, alcohol, physical activity,
# family history of kidney disease -- were not in the files pulled for this draft-scale run; noted
# as a gap for the full run, not silently dropped.)
model_vars <- c("age","female","income_pir","education","hiq011","smq020","bmxbmi","bmxwaist",
                 "diq010","bpq020","cvd_any","race_eth","ckd_case")

dev_df <- dev_design$variables
dev_model_df <- dev_df[dev_df$in_analytic_domain, ]
cat(sprintf("[07_model_m1] development analytic domain for model fitting: n=%d\n", nrow(dev_model_df)))
cat("Missingness in candidate predictors (development analytic domain):\n")
print(sapply(dev_model_df[model_vars], function(x) sum(is.na(x))))

# --- Complete-case fit for this draft pass (S8 in the protocol is complete-case as a robustness
# check; full MICE m=20 imputation, §9.9, is the primary pre-specified analysis and is implemented
# separately -- see 07b_model_m1_mi.R for the imputed version at RUN_SCALE=full). ------------------
m1_design <- subset(dev_design, in_analytic_domain)

m1_formula <- ckd_case ~ ns(age, df = 3) + female + income_pir + factor(education) +
  factor(hiq011) + factor(smq020) + ns(bmxbmi, df = 3) + bmxwaist +
  factor(diq010) + factor(bpq020) + cvd_any + factor(race_eth)

m1_fit <- tryCatch(
  svyglm(m1_formula, design = m1_design, family = quasibinomial()),
  error = function(e) { cat("[07_model_m1] FIT ERROR:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(m1_fit)) {
  cat("\n[07_model_m1] M1 fit summary:\n")
  print(summary(m1_fit))
  saveRDS(m1_fit, file.path(PROCESSED_DIR, "m1_fit.rds"))

  # predict.svyglm does not reliably honour na.action=na.pass across survey package versions, so
  # predict on the complete-case rows explicitly and assign back by row index (rather than relying
  # on vector recycling, which silently misaligns if any rows were dropped).
  complete_rows <- stats::complete.cases(dev_model_df[, c("age","female","income_pir","education",
                                                            "hiq011","smq020","bmxbmi","bmxwaist",
                                                            "diq010","bpq020","cvd_any","race_eth")])
  dev_model_df$m1_pred <- NA_real_
  dev_model_df$m1_pred[complete_rows] <- as.numeric(predict(m1_fit, newdata = dev_model_df[complete_rows, ], type = "response"))
  cat(sprintf("[07_model_m1] predicted for %d/%d complete-case rows (%d missing a predictor)\n",
              sum(complete_rows), nrow(dev_model_df), sum(!complete_rows)))
  cat(sprintf("\n[07_model_m1] Development predicted-risk range: [%.4f, %.4f], mean=%.4f\n",
              min(dev_model_df$m1_pred, na.rm=TRUE), max(dev_model_df$m1_pred, na.rm=TRUE),
              mean(dev_model_df$m1_pred, na.rm=TRUE)))
  saveRDS(dev_model_df, file.path(PROCESSED_DIR, "dev_model_df_with_preds.rds"))
}
