# Config: locked specifications from Milestone1_Protocol_Undiagnosed_CKD_NHANES.pdf (v1.0, 7 Aug 2026)
# and Protocol_Amendment_v1.1_1.pdf (2 Sep 2026). Every constant below cites the section that locks it.
# Do not edit thresholds/definitions here without recording a dated protocol amendment (per the study's own rule).

suppressPackageStartupMessages({
  library(dplyr)
})

# --- Reduced-scale run flags (this session) ------------------------------
# The pre-specified protocol settings (§9.3, §7.2, §9.9) are computationally heavy
# (5x10-fold PSU-clustered CV over 3 models, 500-replicate Rao-Wu bootstrap, m=20 MI).
# RUN_SCALE = "draft" runs a faithful but lighter pass so we get real, honest numbers
# without a multi-hour blocking run. Switch to "full" to reproduce the exact protocol
# settings for the submission-ready analysis.
RUN_SCALE <- Sys.getenv("CKD_RUN_SCALE", unset = "draft")  # "draft" | "full"

SCALE <- if (RUN_SCALE == "full") {
  list(mi_m = 20, boot_R = 500, cv_folds = 10, cv_repeats = 5,
       xgb_grid_n = 20, xgb_nrounds_max = 1000)
} else {
  list(mi_m = 5, boot_R = 200, cv_folds = 10, cv_repeats = 2,
       xgb_grid_n = 8, xgb_nrounds_max = 400)
}

# --- Paths -----------------------------------------------------------------
# All paths are relative to the project root; scripts are always run/sourced from there
# (see R/run_pipeline.R), so no working-directory detection is needed here.
RAW_DIR       <- "data/raw"
PROCESSED_DIR <- "data/processed"
RESULTS_DIR   <- "results"
TABLES_DIR    <- file.path(RESULTS_DIR, "tables")
FIGURES_DIR   <- file.path(RESULTS_DIR, "figures")
for (d in c(RAW_DIR, PROCESSED_DIR, TABLES_DIR, FIGURES_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --- NHANES cycle manifest (Protocol §7.1, §8.1-8.2) ------------------------
# role: "development" (pooled, MEC-weighted) or "validation" (temporal external, Aug2021-Aug2023)
CYCLES <- list(
  G = list(role = "development", years = "2011-2012", suffix = "G", url_dir = "2011",
           weight_var = "WTMEC2YR", int_weight_var = "WTINT2YR", weight_frac = 2/9.2),
  H = list(role = "development", years = "2013-2014", suffix = "H", url_dir = "2013",
           weight_var = "WTMEC2YR", int_weight_var = "WTINT2YR", weight_frac = 2/9.2),
  I = list(role = "development", years = "2015-2016", suffix = "I", url_dir = "2015",
           weight_var = "WTMEC2YR", int_weight_var = "WTINT2YR", weight_frac = 2/9.2),
  P = list(role = "development", years = "2017-March2020", suffix = "P_", url_dir = "2017",
           weight_var = "WTMECPRP", int_weight_var = "WTINTPRP", weight_frac = 3.2/9.2, prefix_style = TRUE),
  L = list(role = "validation", years = "Aug2021-Aug2023", suffix = "L", url_dir = "2021",
           weight_var = "WTMEC2YR", int_weight_var = "WTINT2YR", weight_frac = NA)  # primary WTPH2YR joined separately from BIOPRO_L
)

# File-name stem -> which NHANES component file (with contents used, protocol §8)
FILE_STEMS <- c(
  DEMO    = "DEMO",     # demographics, design vars, weights, INDFMPIR, RIDRETH3, RIDEXPRG
  BIOPRO  = "BIOPRO",   # LBXSCR serum creatinine
  ALB_CR  = "ALB_CR",   # URXUMA, URXUCR, URDACT
  KIQ_U   = "KIQ_U",    # KIQ022 (undiagnosed restriction), KIQ025 (dialysis)
  DIQ     = "DIQ",      # DIQ010 self-reported diabetes (comparator C1)
  BPQ     = "BPQ",      # BPQ020 self-reported hypertension (comparator C1)
  MCQ     = "MCQ",      # self-reported CVD (comparator C2)
  HIQ     = "HIQ",      # insurance (subgroup analysis)
  SMQ     = "SMQ",      # smoking (predictor)
  BMX     = "BMX"       # BMI, waist circumference (predictor)
)

nhanes_file_name <- function(stem, cyc) {
  info <- CYCLES[[cyc]]
  if (isTRUE(info$prefix_style)) paste0("P_", stem) else paste0(stem, "_", info$suffix)
}

nhanes_url <- function(stem, cyc) {
  info <- CYCLES[[cyc]]
  fname <- nhanes_file_name(stem, cyc)
  sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%s/DataFiles/%s.XPT", info$url_dir, fname)
}

# --- Outcome definition (Protocol §5) ---------------------------------------
EGFR_THRESHOLD <- 60      # mL/min/1.73m^2, KDIGO G3a-G5
ACR_THRESHOLD  <- 30      # mg/g, KDIGO A2-A3

# --- Undiagnosed restriction (Protocol §5.4) --------------------------------
# Eligible: KIQ022 == 2 (No). Excluded: KIQ022 in c(1,7,9,NA) or KIQ025 == 1 (dialysis, past 12mo).

# --- Comparator rules (Protocol §6) -----------------------------------------
# C1: DIQ010==1 OR BPQ020==1  (diabetes or hypertension, self-report)
# C2: C1 OR age>=60 OR self-reported CVD (MI/CHD/heart failure/stroke)
# C3 / SCORED >= 4 (Amendment names this "SCORED" comparator; 9-item score, see 03 script)

# --- Decision thresholds (Protocol §6.4) ------------------------------------
DCA_THRESHOLDS <- seq(0.02, 0.20, by = 0.01)   # full curve 2%-20%
DCA_PRIMARY    <- 0.05                          # primary reporting threshold
DCA_SECONDARY  <- 0.10                          # secondary threshold reported in text

# --- Subgroup definitions (Protocol §9.6) ------------------------------------
INCOME_BREAKS  <- c(-Inf, 1.0, 2.0, 4.0, Inf)
INCOME_LABELS  <- c("<1.0 (below FPL)", "1.0-1.99", "2.0-3.99", ">=4.0")
MIN_EVENTS_FOR_SLOPE <- 50   # §9.6: levels with <50 events suppress calibration slope, report O/E only

# Amendment §... collapsing hierarchy for the C1-unflagged (restricted) stratum ONLY,
# applied to the calibration slope only, when the full-granularity race/ethnicity level
# does not meet MIN_EVENTS_FOR_SLOPE. Order matters: apply top to bottom until met.
RACE_COLLAPSE_HIERARCHY <- list(
  full        = c("Mexican American","Other Hispanic","Non-Hispanic White","Non-Hispanic Black","Non-Hispanic Asian","Other/Multiracial"),
  hispanic_combined = c("Non-Hispanic White","Non-Hispanic Black","Hispanic","Non-Hispanic Asian","Other/Multiracial")
)
INCOME_COLLAPSE_BINARY   <- c("<2.0x FPL", ">=2.0x FPL")
INSURANCE_COLLAPSE_BINARY <- c("Insured", "Uninsured")

# --- Survey design (Protocol §7.2) -------------------------------------------
# Design object built ONCE on the pooled/full file; every restriction (undiagnosed, C1-unflagged,
# subgroup) is applied via subset()/domain estimation on that object, never by pre-filtering rows.
BOOT_REPLICATES <- SCALE$boot_R   # Rao-Wu rescaled bootstrap, resampling PSU within strata

# --- Urinary albumin bridging correction (Protocol §7.3, Amendment §3) ------
# Backward (new->old) weighted Deming equation. Amendment §3.3 made this a pre-registration
# gate ("not yet confirmed against source documentation ... blocked on data access"). CONFIRMED
# 21 Sep 2026 directly against the CDC/NCHS ALB_CR_L documentation page
# (https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/ALB_CR_L.htm), which states:
#   URXUMA_old = -1.643 + 1.189 x URXUMA_new (mg/L); 95% CI intercept -1.930 to -1.355, slope 1.158-1.221.
# This matches Protocol v1.0 exactly. The Amendment v1.1 §3.1 figure of "-1.642" for the intercept
# is a transcription rounding error relative to the confirmed source value of -1.643; the confirmed
# value is used throughout (this changes no downstream conclusion -- see Amendment §3.2 re: the floor
# window -- but should be corrected in the amendment text before OSF posting).
BRIDGE_INTERCEPT <- -1.643
BRIDGE_SLOPE     <- 1.189
BRIDGE_INTERCEPT_CI <- c(-1.930, -1.355)
BRIDGE_SLOPE_CI     <- c(1.158, 1.221)
OLD_METHOD_LLOD_UMA <- 0.30    # mg/L, old method (solid-phase fluorescent immunoassay)
NEW_METHOD_LLOD_UMA <- 0.02    # mg/L, LC-MS/MS
FLOOR_VALUE         <- OLD_METHOD_LLOD_UMA / sqrt(2)   # 0.212 mg/L, Amendment §3.2

# --- Models (Protocol §9.2) ---------------------------------------------------
# M0: Rule C1 (not fitted) | M1: survey-weighted logistic + rcs(age,4)+rcs(bmi,4)
# M2: elastic-net penalised logistic, survey weights as obs weights
# M3: XGBoost, survey weights as instance weights
# Parsimony rule: simplest of M1-M3 within 1 SE of best CV weighted log loss (fixed pre-hoc).

cat(sprintf("[config] RUN_SCALE=%s  mi_m=%d boot_R=%d cv=%dx%d\n",
            RUN_SCALE, SCALE$mi_m, SCALE$boot_R, SCALE$cv_repeats, SCALE$cv_folds))
