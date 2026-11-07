# Derives, per participant: eGFR (race-free CKD-EPI 2021, primary; race-included, secondary per §9.7),
# ACR, the CKD case definition, the undiagnosed restriction, and the comparator rules C1/C2/SCORED.
# Protocol §5 (endpoint, locked), §6 (comparators, locked). All functions are pure / vectorised so
# they can be unit-tested against known values before being run on real NHANES rows.

source("R/00_config.R")
source("R/02_bridging_correction.R")

# --- eGFR: CKD-EPI 2021 creatinine equation, race-free (Protocol §5.1; Inker et al., NEJM 2021) ---
# sex: 1 = male, 2 = female (NHANES RIAGENDR coding)
egfr_ckdepi_2021 <- function(scr, age, sex) {
  female <- sex == 2
  kappa  <- ifelse(female, 0.7, 0.9)
  alpha  <- ifelse(female, -0.241, -0.302)
  sex_mult <- ifelse(female, 1.012, 1.0)
  ratio <- scr / kappa
  egfr <- 142 * pmin(ratio, 1)^alpha * pmax(ratio, 1)^(-1.200) * 0.9938^age * sex_mult
  egfr
}

# Race-included CKD-EPI 2021 variant (Protocol §9.7, secondary comparison only -- the headline
# outcome/model per D3 is race-free). NHANES RIDRETH3-derived "Black" indicator required by caller;
# coefficients per Inker et al. 2021 Table for the race-included refit.
egfr_ckdepi_2021_race <- function(scr, age, sex, black) {
  female <- sex == 2
  kappa  <- ifelse(female, 0.7, 0.9)
  alpha  <- ifelse(female, -0.241, -0.302)
  sex_mult <- ifelse(female, 1.018, 1.0)
  race_mult <- ifelse(black, 1.159, 1.0)
  ratio <- scr / kappa
  142 * pmin(ratio, 1)^alpha * pmax(ratio, 1)^(-1.200) * 0.9938^age * sex_mult * race_mult
}

.egfr_selftest <- function() {
  # 50-year-old woman, Scr 0.9 mg/dL: ratio=0.9/0.7>1 -> use max branch exponent -1.200
  e <- egfr_ckdepi_2021(scr = 0.9, age = 50, sex = 2)
  stopifnot(e > 60 && e < 100)  # sanity: normal-ish creatinine, mid-age -> not floridly abnormal
  # Higher creatinine -> lower eGFR (monotonicity)
  stopifnot(egfr_ckdepi_2021(1.5, 50, 1) < egfr_ckdepi_2021(0.9, 50, 1))
  # Older age -> lower eGFR at fixed creatinine
  stopifnot(egfr_ckdepi_2021(0.9, 70, 1) < egfr_ckdepi_2021(0.9, 30, 1))
  cat("[03_outcome_comparators] eGFR self-test OK (example: 50F, Scr=0.9 -> eGFR=", round(e,1), ")\n")
}
.egfr_selftest()

# --- Outcome construction (Protocol §5.3) -----------------------------------
#' @param urxuma urinary albumin (mg/L), ALREADY bridging-corrected if validation cycle (see caller)
#' @param urxucr urine creatinine (mg/dL)
derive_acr <- function(urxuma, urxucr) (urxuma / urxucr) * 100  # URDACT formula, Protocol §5.2

derive_ckd_case <- function(egfr, acr) (egfr < EGFR_THRESHOLD) | (acr >= ACR_THRESHOLD)

# --- Undiagnosed restriction (Protocol §5.4) --------------------------------
derive_undiagnosed_eligible <- function(kiq022, kiq025) {
  eligible <- kiq022 == 2
  eligible[is.na(kiq022) | !(kiq022 %in% c(1, 2))] <- NA   # refused/DK/missing -> indeterminate, excluded
  eligible[!is.na(kiq025) & kiq025 == 1] <- FALSE            # dialysis in past 12mo -> excluded regardless
  eligible
}

# --- Comparator C1 (Protocol §6.1): diabetes OR hypertension, self-report ---
derive_c1 <- function(diq010, bpq020) {
  (!is.na(diq010) & diq010 == 1) | (!is.na(bpq020) & bpq020 == 1)
}

# --- Comparator C2 (Protocol §6.2): C1 OR age>=60 OR self-reported CVD ------
derive_c2 <- function(c1, age, cvd_any) {
  c1 | (age >= 60) | (!is.na(cvd_any) & cvd_any)
}

# --- SCORED (Protocol §6.2/Table3; Bang et al., Arch Intern Med 2007;167(4):374-81) --------------
# Nine items, integer points, cut-point >=4 flags. Points confirmed against the UNC Kidney Center's
# patient-facing SCORED card, which cites the Bang et al. 2007 reference directly:
#   Age: 50-59y=2, 60-69y=3, >=70y=4 (one tier only, <50y=0 -- NOT cumulative across tiers)
#   Female sex=1; anemia hx=1; hypertension=1; diabetes=1; history of MI or stroke=1;
#   history of CHF/heart failure=1; peripheral vascular ("circulation disease in legs")=1;
#   protein in urine (self-report)=1.
# (Note: MI/stroke and CHF are two SEPARATE 1-point items here, distinct from Protocol §6.2's C2
# definition, which deliberately combines MI/CHD/HF/stroke into a single flag -- that is a genuine
# difference between the two comparators, not an inconsistency.)
# Availability check (Amendment §2): computed only for participants with all nine items present in
# a given cycle; missing-item share reported by cycle rather than imputed across this comparison.
SCORED_ITEMS <- c("female","anemia_hx","hypertension","diabetes","mi_stroke_hx","chf_hx","pad_hx","proteinuria_hx","age")

scored_age_points <- function(age) {
  ifelse(is.na(age), NA_real_, ifelse(age >= 70, 4, ifelse(age >= 60, 3, ifelse(age >= 50, 2, 0))))
}

derive_scored <- function(df) {
  # df must contain: female, anemia_hx, hypertension, diabetes, mi_stroke_hx, chf_hx, pad_hx,
  # proteinuria_hx (all logical/0-1) and age (numeric). NA in any -> that item unavailable.
  item_cols <- c("female","anemia_hx","hypertension","diabetes","mi_stroke_hx","chf_hx","pad_hx","proteinuria_hx")
  all_present <- Reduce(`&`, lapply(c(item_cols, "age"), function(cn) !is.na(df[[cn]])))

  age_pts <- scored_age_points(df$age)

  pts <- as.numeric(df$female) + as.numeric(df$anemia_hx) + as.numeric(df$hypertension) +
    as.numeric(df$diabetes) + as.numeric(df$mi_stroke_hx) + as.numeric(df$chf_hx) +
    as.numeric(df$pad_hx) + as.numeric(df$proteinuria_hx) + age_pts

  list(score = ifelse(all_present, pts, NA_real_),
       flag  = ifelse(all_present, pts >= 4, NA),
       available = all_present)
}

.scored_selftest <- function() {
  # 62-year-old woman, hypertensive, diabetic, no anemia/CVD/PAD/proteinuria:
  # age(60-69)=3 + female=1 + hypertension=1 + diabetes=1 = 6 -> flag TRUE
  df <- data.frame(female = TRUE, anemia_hx = FALSE, hypertension = TRUE, diabetes = TRUE,
                    mi_stroke_hx = FALSE, chf_hx = FALSE, pad_hx = FALSE, proteinuria_hx = FALSE,
                    age = 62)
  s <- derive_scored(df)
  stopifnot(s$score == 6, s$flag == TRUE, s$available == TRUE)

  # 30-year-old man, no risk factors: age=0, all else 0 -> score 0, flag FALSE
  df2 <- data.frame(female = FALSE, anemia_hx = FALSE, hypertension = FALSE, diabetes = FALSE,
                     mi_stroke_hx = FALSE, chf_hx = FALSE, pad_hx = FALSE, proteinuria_hx = FALSE,
                     age = 30)
  s2 <- derive_scored(df2)
  stopifnot(s2$score == 0, s2$flag == FALSE)

  # Missing item (anemia_hx = NA) -> unavailable, not imputed as 0
  df3 <- df; df3$anemia_hx <- NA
  s3 <- derive_scored(df3)
  stopifnot(is.na(s3$score), is.na(s3$flag), s3$available == FALSE)

  cat("[03_outcome_comparators] SCORED self-test OK (62F+HTN+DM example -> score=6, flag=TRUE)\n")
}
.scored_selftest()

.comparator_selftest <- function() {
  stopifnot(derive_c1(1, 2) == TRUE)      # diabetes yes -> flagged
  stopifnot(derive_c1(2, 2) == FALSE)     # neither -> unflagged
  stopifnot(derive_ckd_case(55, 10) == TRUE)   # low eGFR alone -> case
  stopifnot(derive_ckd_case(90, 45) == TRUE)   # high ACR alone -> case
  stopifnot(derive_ckd_case(90, 10) == FALSE)  # neither -> not a case
  stopifnot(is.na(derive_undiagnosed_eligible(9, 2)))   # "don't know" -> indeterminate
  stopifnot(derive_undiagnosed_eligible(2, 1) == FALSE) # No prior dx but on dialysis -> excluded
  stopifnot(derive_undiagnosed_eligible(2, 2) == TRUE)  # No prior dx, no dialysis -> eligible
  cat("[03_outcome_comparators] comparator self-test OK\n")
}
.comparator_selftest()
