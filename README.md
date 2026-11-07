# Case-finding-for-undiagnosed-CKD-in-NHANES
Case-finding for undiagnosed CKD in NHANES
Analysis code for the manuscript Non-laboratory case-finding for undiagnosed chronic kidney disease in US adults without diagnosed diabetes or hypertension: development and temporal external validation in NHANES (Reza Omid).

The pipeline develops a non-laboratory prediction model for undiagnosed chronic kidney disease (eGFR <60 mL/min/1.73 m² by the race-free 2021 CKD-EPI equation and/or urinary albumin-to-creatinine ratio ≥30 mg/g, among adults reporting no prior kidney-disease diagnosis) in pooled NHANES 2011–March 2020 cycles, and evaluates the frozen model, without recalibration, in the August 2021–August 2023 cycle. Reporting follows the TRIPOD+AI statement.

Data

All data are public-use NHANES files from the National Center for Health Statistics (https://wwwn.cdc.gov/nchs/nhanes/). No data are stored in this repository: R/01_download_data.R downloads the required files into data/raw/.

Requirements

R with the packages survey, mitools, mice, glmnet, xgboost, pROC, haven, dplyr, ggplot2 and scales:

r
install.packages(c("survey", "mitools", "mice", "glmnet", "xgboost",
                   "pROC", "haven", "dplyr", "ggplot2", "scales"))
Running the pipeline

Run from the repository root (all paths are relative to it):

sh
Rscript R/run_pipeline.R

The default settings are the ones reported in the manuscript: 5 imputations, 2 × 10-fold PSU-clustered cross-validation, an 8-combination XGBoost grid, and 500 bootstrap replicates for the multiple-imputation validation metrics (200 for the complete-case pass). To use the protocol's pre-specified computational settings (20 imputations, 500 bootstrap replicates, 5 × 10-fold cross-validation):

sh
CKD_RUN_SCALE=full Rscript R/run_pipeline.R

Deviations from the registered protocol are listed in Table S3 of the manuscript. Tables are written to results/tables/, figures to results/figures/, and intermediate files to data/processed/.

Scripts
Script	Purpose
00_config.R	Locked protocol settings: cycles and weights, outcome and comparator definitions, decision thresholds, bridging coefficients, run scale
01_download_data.R	Downloads and caches the NHANES component files
02_bridging_correction.R	NCHS urinary albumin bridging correction for the validation cycle
03_outcome_comparators.R	eGFR, ACR, the undiagnosed-CKD outcome, and comparator rules C1, C2 and SCORED
04_harmonize_dataset.R	Merges cycles into one harmonized analytic dataset
05_survey_design.R	Complex survey design objects (strata, PSU, pooled weights)
06_table1_flow.R	Participant flow and weighted Table 1
07_model_m1.R	M1: survey-weighted logistic regression with natural splines for age and BMI
08_decision_curve_calibration.R	Frozen-model temporal validation, complete-case: calibration, discrimination, decision curve
09_subgroup_equity.R	Complete-case subgroup calibration and net benefit
10_model_m2_m3_parsimony.R	M2 (elastic net) and M3 (XGBoost), PSU-clustered cross-validation, parsimony rule
11_mice_imputation.R	Multiple imputation within development and validation sets, pooled with Rubin's rules
12_sensitivity_s1_s3_s9.R	Sensitivity analyses S1 (no bridging correction), S3 (WTMEC2YR weight), S9 (C2-defined stratum)
13_subgroup_equity_mi.R	Subgroup calibration with multiple-imputation predictions
14_figures.R	Figures
15_mice_validation_metrics.R	Multiple-imputation validation metrics with design-based bootstrap CIs (headline results)
run_pipeline.R	Runs all steps in order
Citation

If you use this code, please cite the article; the full reference will be added on publication. An archived release will be deposited on Zenodo with a DOI on acceptance.
