# Merges the per-cycle component files into one analytic dataset with a harmonized column schema.
# Protocol §5 (endpoint), §6 (comparators), §9.1 (candidate predictors), §7.1-7.2 (split & weights).
# Variable codes confirmed directly against CDC/NCHS data-file documentation (2026-09-21):
#   DEMO: RIAGENDR, RIDAGEYR, RIDRETH3, INDFMPIR, DMDEDUC2, RIDEXPRG, SDMVSTRA, SDMVPSU, WTINT2YR,
#         WTMEC2YR (WTMECPRP for the P_ pre-pandemic file), WTPH2YR (L cycle only, in BIOPRO_L)
#   BIOPRO: LBXSCR | ALB_CR: URXUMA, URXUCR, URDACT | KIQ_U: KIQ022, KIQ025
#   DIQ: DIQ010 | BPQ: BPQ020 | MCQ: MCQ160B(CHF) MCQ160C(CHD) MCQ160E(MI) MCQ160F(stroke)
#   HIQ: HIQ011 (any coverage) | SMQ: SMQ020, SMQ040 | BMX: BMXBMI, BMXWAIST

source("R/00_config.R")
source("R/02_bridging_correction.R")
source("R/03_outcome_comparators.R")
library(dplyr)

read_component <- function(stem, cyc) {
  fname <- nhanes_file_name(stem, cyc)
  path <- file.path(RAW_DIR, paste0(fname, ".rds"))
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

harmonize_cycle <- function(cyc) {
  info <- CYCLES[[cyc]]
  demo   <- read_component("DEMO", cyc)
  biopro <- read_component("BIOPRO", cyc)
  albcr  <- read_component("ALB_CR", cyc)
  kiq    <- read_component("KIQ_U", cyc)
  diq    <- read_component("DIQ", cyc)
  bpq    <- read_component("BPQ", cyc)
  mcq    <- read_component("MCQ", cyc)
  hiq    <- read_component("HIQ", cyc)
  smq    <- read_component("SMQ", cyc)
  bmx    <- read_component("BMX", cyc)
  if (is.null(demo) || is.null(biopro) || is.null(albcr) || is.null(kiq)) {
    warning(sprintf("cycle %s: missing a required file, skipping", cyc)); return(NULL)
  }

  d <- demo %>%
    transmute(SEQN, cycle = cyc, role = info$role,
              age = RIDAGEYR, sex = RIAGENDR,
              race_eth = RIDRETH3,
              income_pir = INDFMPIR,
              education = if ("DMDEDUC2" %in% names(demo)) DMDEDUC2 else NA_real_,
              pregnant  = if ("RIDEXPRG" %in% names(demo)) RIDEXPRG else NA_real_,
              sdmvstra = SDMVSTRA, sdmvpsu = SDMVPSU,
              wt_int = .data[[info$int_weight_var]],
              wt_mec = .data[[info$weight_var]])

  d <- d %>% left_join(biopro %>% select(SEQN, LBXSCR), by = "SEQN")
  d <- d %>% left_join(albcr  %>% select(SEQN, URXUMA, URXUCR), by = "SEQN")
  d <- d %>% left_join(kiq    %>% select(SEQN, KIQ022, KIQ025), by = "SEQN")

  if (cyc == "L" && "WTPH2YR" %in% names(biopro)) {
    d <- d %>% left_join(biopro %>% select(SEQN, WTPH2YR), by = "SEQN")
  } else {
    d$WTPH2YR <- NA_real_
  }

  pick <- function(df, col) if (!is.null(df) && col %in% names(df)) df[[col]][match(d$SEQN, df$SEQN)] else NA

  d$diq010 <- if (!is.null(diq)) pick(diq, "DIQ010") else NA
  d$bpq020 <- if (!is.null(bpq)) pick(bpq, "BPQ020") else NA
  d$mcq_mi   <- if (!is.null(mcq)) pick(mcq, "MCQ160E") else NA
  d$mcq_chd  <- if (!is.null(mcq)) pick(mcq, "MCQ160C") else NA
  d$mcq_chf  <- if (!is.null(mcq)) pick(mcq, "MCQ160B") else NA
  d$mcq_stroke <- if (!is.null(mcq)) pick(mcq, "MCQ160F") else NA
  d$hiq011  <- if (!is.null(hiq)) pick(hiq, "HIQ011") else NA
  d$smq020  <- if (!is.null(smq)) pick(smq, "SMQ020") else NA
  d$smq040  <- if (!is.null(smq)) pick(smq, "SMQ040") else NA
  d$bmxbmi  <- if (!is.null(bmx)) pick(bmx, "BMXBMI") else NA
  d$bmxwaist<- if (!is.null(bmx)) pick(bmx, "BMXWAIST") else NA

  d
}

cat("Harmonizing cycles...\n")
all_cycles <- lapply(names(CYCLES), harmonize_cycle)
names(all_cycles) <- names(CYCLES)
pooled <- bind_rows(all_cycles)
cat(sprintf("Pooled raw rows (pre-eligibility): %d across %d cycles\n", nrow(pooled), length(all_cycles)))

# --- Apply protocol-locked inclusion/exclusion (§5.5) -----------------------
pooled <- pooled %>%
  mutate(
    female = sex == 2,
    cvd_any = (!is.na(mcq_mi) & mcq_mi == 1) | (!is.na(mcq_chd) & mcq_chd == 1) |
              (!is.na(mcq_chf) & mcq_chf == 1) | (!is.na(mcq_stroke) & mcq_stroke == 1),
    mi_stroke_hx = (!is.na(mcq_mi) & mcq_mi == 1) | (!is.na(mcq_stroke) & mcq_stroke == 1),
    chf_hx = !is.na(mcq_chf) & mcq_chf == 1,
    # PAD is not collected as a discrete self-report item in the MCQ files reviewed; SCORED's PAD
    # item is therefore unavailable in NHANES, which will make SCORED "available" = FALSE for
    # everyone unless a proxy is substituted -- flagged explicitly rather than silently assumed 0.
    pad_hx = NA,
    anemia_hx = NA,      # NHANES has no single "ever told you have anemia" item in files reviewed here
    proteinuria_hx = NA, # not self-reported in NHANES; would require a lab-based proxy, which would
                          # break SCORED's "non-laboratory, self-report" design -- left NA, not imputed
    age_eligible = age >= 20,
    egfr_race_free = egfr_ckdepi_2021(LBXSCR, age, sex),
    black = race_eth == 4,  # RIDRETH3: 4 = Non-Hispanic Black
    egfr_race_included = egfr_ckdepi_2021_race(LBXSCR, age, sex, black)
  )

# Bridging correction (§7.3) applied ONLY to validation-cycle (L) urinary albumin, before ACR.
pooled <- pooled %>%
  mutate(urxuma_corrected = ifelse(cycle == "L" & !is.na(URXUMA),
                                    apply_bridging_correction(URXUMA), URXUMA),
         acr = derive_acr(urxuma_corrected, URXUCR),
         acr_uncorrected = derive_acr(URXUMA, URXUCR))  # Sensitivity A (§7.3), L cycle only meaningful

pooled <- pooled %>%
  mutate(
    ckd_case = derive_ckd_case(egfr_race_free, acr),
    undiagnosed_eligible = derive_undiagnosed_eligible(KIQ022, KIQ025),
    c1_flag = derive_c1(diq010, bpq020)
  ) %>%
  mutate(c2_flag = derive_c2(c1_flag, age, cvd_any))

# NHANES codes 7 ("Refused") and 9 ("Don't know") on these items are missing-data codes, not
# substantive categories -- recoding them to NA before they reach any model as factor levels.
# (This was caught by a real CV failure: education=7 appeared in a held-out fold but not in that
# fold's training data, so predict() errored on an unseen factor level -- the fix is not just
# technical, treating "Refused" as a legitimate, predictive education category would be wrong
# regardless.) diq010=3 ("Borderline") is a genuine informative category and is kept.
recode_missing <- function(x, codes) { x[x %in% codes] <- NA; x }
pooled <- pooled %>%
  mutate(
    education = recode_missing(education, c(7, 9)),
    hiq011    = recode_missing(hiq011, c(7, 9)),
    smq020    = recode_missing(smq020, c(7, 9)),
    diq010    = recode_missing(diq010, c(7, 9)),
    bpq020    = recode_missing(bpq020, c(7, 9))
  )

scored_input <- pooled %>%
  transmute(female, anemia_hx, pad_hx, proteinuria_hx, mi_stroke_hx, chf_hx, age,
            hypertension = ifelse(bpq020 %in% c(1, 2), bpq020 == 1, NA),
            diabetes     = ifelse(diq010 %in% c(1, 2), diq010 == 1, NA))
scored_out <- derive_scored(scored_input)
pooled$scored_score <- scored_out$score
pooled$scored_flag  <- scored_out$flag
pooled$scored_available <- scored_out$available

# Analytic-set inclusion (§5.5): non-pregnant adults 20+, valid creatinine, valid ACR, valid KIQ022,
# undiagnosed (eligible==TRUE). Exclusions are recorded as a domain flag, NOT by pre-filtering rows,
# so the full survey design object retains every row for correct variance estimation (§5.5, §11).
pooled <- pooled %>%
  mutate(
    excl_pregnant = !is.na(pregnant) & pregnant == 1,
    excl_missing_outcome = is.na(LBXSCR) | is.na(URXUMA) | is.na(URXUCR),
    excl_missing_kiq = is.na(KIQ022),
    in_analytic_domain = age_eligible & !excl_pregnant & !excl_missing_outcome &
                          !excl_missing_kiq & !is.na(undiagnosed_eligible) & undiagnosed_eligible
  )

# --- Development pooled examination weight (§7.2): WTMEC2YR x (2/9.2) for G/H/I; WTMECPRP x (3.2/9.2) for P
pooled <- pooled %>%
  mutate(wt_dev = ifelse(role == "development", wt_mec * sapply(cycle, function(c) CYCLES[[c]]$weight_frac), NA_real_))

saveRDS(pooled, file.path(PROCESSED_DIR, "analytic_pooled.rds"))

cat(sprintf("\nAnalytic domain (undiagnosed, complete outcome, adult, non-pregnant): %d of %d rows\n",
            sum(pooled$in_analytic_domain, na.rm = TRUE), nrow(pooled)))
cat("By role:\n")
print(pooled %>% filter(in_analytic_domain) %>% count(role))
cat("By cycle:\n")
print(pooled %>% filter(in_analytic_domain) %>% count(cycle))
cat(sprintf("Weighted CKD case count (unweighted n, analytic domain): %d cases / %d\n",
            sum(pooled$ckd_case & pooled$in_analytic_domain, na.rm = TRUE),
            sum(pooled$in_analytic_domain, na.rm = TRUE)))
cat(sprintf("C1-unflagged share of analytic domain: %.1f%%\n",
            100 * mean(!pooled$c1_flag[pooled$in_analytic_domain], na.rm = TRUE)))
cat("SCORED availability (all 9 items present):\n")
print(table(pooled$scored_available[pooled$in_analytic_domain], useNA = "ifany"))
