# Runs the full pipeline end to end, in order. Set CKD_RUN_SCALE=full in the environment before
# running to reproduce the pre-specified protocol settings (m=20 imputations, 500 bootstrap
# replicates, 5x10-fold CV) instead of this draft pass's reduced settings (see 00_config.R).
#
# Usage:  Rscript R/run_pipeline.R
#         CKD_RUN_SCALE=full Rscript R/run_pipeline.R

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
steps <- c(
  "R/00_config.R",
  "R/01_download_data.R",       # idempotent: skips files already cached in data/raw/
  "R/02_bridging_correction.R",
  "R/03_outcome_comparators.R",
  "R/04_harmonize_dataset.R",
  "R/05_survey_design.R",
  "R/06_table1_flow.R",
  "R/07_model_m1.R",
  "R/08_decision_curve_calibration.R",
  "R/09_subgroup_equity.R",
  "R/10_model_m2_m3_parsimony.R",
  "R/11_mice_imputation.R",
  "R/12_sensitivity_s1_s3_s9.R",
  "R/13_subgroup_equity_mi.R"
  # Not yet built: SHAP, S4-S8, MICE-based decision-curve/calibration-slope/CV, race-included eGFR
  # secondary comparison -- see manuscript/manuscript_draft.md "Status and scope of this draft".
)
for (s in steps) { cat("\n", strrep("=", 70), "\n[run_pipeline] ", s, "\n", strrep("=", 70), "\n", sep = ""); source(s) }
cat("\n[run_pipeline] Done. See results/tables/ and manuscript/manuscript_draft.md.\n")
