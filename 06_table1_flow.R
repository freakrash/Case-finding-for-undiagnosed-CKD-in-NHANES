# Real participant flow + weighted Table 1 (TRIPOD+AI items 20a/20b), computed from the actual
# harmonized NHANES data -- not a projection. Protocol §8.3 projected dev ~20-22k/2.2-2.6k events,
# validation ~5.5-6.5k/600-750 events (whole undiagnosed set); Amendment §1 revised the validation
# C1-unflagged headline-stratum projection to ~180-240 events.

source("R/00_config.R")
source("R/05_survey_design.R")
library(survey)

flow <- function(df, label) {
  n_total <- nrow(df)
  n_adult <- sum(df$age_eligible, na.rm = TRUE)
  n_notpreg <- sum(df$age_eligible & !df$excl_pregnant, na.rm = TRUE)
  n_complete_outcome <- sum(df$age_eligible & !df$excl_pregnant & !df$excl_missing_outcome, na.rm = TRUE)
  n_valid_kiq <- sum(df$age_eligible & !df$excl_pregnant & !df$excl_missing_outcome & !df$excl_missing_kiq, na.rm = TRUE)
  n_undiagnosed <- sum(df$in_analytic_domain, na.rm = TRUE)
  data.frame(set = label, n_total_file = n_total, n_adult_20plus = n_adult,
             n_nonpregnant = n_notpreg, n_complete_outcome = n_complete_outcome,
             n_valid_kiq022 = n_valid_kiq, n_undiagnosed_analytic = n_undiagnosed)
}

pooled_all <- readRDS(file.path(PROCESSED_DIR, "analytic_pooled.rds"))
flow_table <- rbind(
  flow(pooled_all %>% dplyr::filter(role == "development"), "Development (2011 - March 2020)"),
  flow(pooled_all %>% dplyr::filter(role == "validation"), "Temporal validation (Aug 2021 - Aug 2023)")
)
print(flow_table)
write.csv(flow_table, file.path(TABLES_DIR, "participant_flow.csv"), row.names = FALSE)

# --- Weighted, design-based counts within the analytic (undiagnosed) domain -----------------------
weighted_summary <- function(design, label) {
  dom <- subset(design, in_analytic_domain)
  n_unw <- nrow(dom$variables)
  cases_unw <- sum(dom$variables$ckd_case, na.rm = TRUE)
  prev <- svymean(~ckd_case, dom, na.rm = TRUE)
  p <- coef(prev)["ckd_caseTRUE"]; se <- SE(prev)["ckd_caseTRUE"]
  n_pop_est <- sum(weights(dom))
  cat(sprintf("\n[%s] unweighted n=%d, unweighted cases=%d (%.1f%%)\n", label, n_unw, cases_unw, 100*cases_unw/n_unw))
  cat(sprintf("[%s] weighted prevalence: %.2f%% (SE %.2f), estimated population size ~%.0f\n",
              label, 100*p, 100*se, n_pop_est))
  invisible(list(n_unw = n_unw, cases_unw = cases_unw, weighted_prev = p, weighted_se = se))
}

dev_sum <- weighted_summary(dev_design, "Development, full undiagnosed set")
val_sum <- weighted_summary(val_design, "Validation, full undiagnosed set (WTPH2YR)")

# Headline C1-unflagged stratum (Protocol §6.3) -- domain estimate on top of the same design object
dev_headline  <- subset(dev_design, in_analytic_domain & !c1_flag)
val_headline  <- subset(val_design, in_analytic_domain & !c1_flag)
dev_h_sum <- (function() {
  n_unw <- nrow(dev_headline$variables); cases_unw <- sum(dev_headline$variables$ckd_case, na.rm=TRUE)
  prev <- svymean(~ckd_case, dev_headline, na.rm = TRUE)
  p <- coef(prev)["ckd_caseTRUE"]; se <- SE(prev)["ckd_caseTRUE"]
  cat(sprintf("\n[Development, C1-unflagged headline stratum] unweighted n=%d, cases=%d (%.1f%%), weighted prev=%.2f%% (SE %.2f)\n",
              n_unw, cases_unw, 100*cases_unw/n_unw, 100*p, 100*se))
  list(n_unw=n_unw, cases_unw=cases_unw, weighted_prev=p)
})()
val_h_sum <- (function() {
  n_unw <- nrow(val_headline$variables); cases_unw <- sum(val_headline$variables$ckd_case, na.rm=TRUE)
  prev <- svymean(~ckd_case, val_headline, na.rm = TRUE)
  p <- coef(prev)["ckd_caseTRUE"]; se <- SE(prev)["ckd_caseTRUE"]
  cat(sprintf("[Validation, C1-unflagged headline stratum] unweighted n=%d, cases=%d (%.1f%%), weighted prev=%.2f%% (SE %.2f)\n",
              n_unw, cases_unw, 100*cases_unw/n_unw, 100*p, 100*se))
  list(n_unw=n_unw, cases_unw=cases_unw, weighted_prev=p)
})()

# --- Bridging correction impact on validation-cycle case status (Amendment §3.1, §S2 primary result) --
val_dat <- val_headline$variables
flip_summary <- val_dat %>%
  dplyr::mutate(case_corrected = ckd_case,
                case_uncorrected = derive_ckd_case_uncorr <- (egfr_race_free < EGFR_THRESHOLD) | (acr_uncorrected >= ACR_THRESHOLD)) %>%
  dplyr::summarise(
    n = dplyr::n(),
    newly_case = sum(case_corrected & !case_uncorrected, na.rm = TRUE),
    newly_noncase = sum(!case_corrected & case_uncorrected, na.rm = TRUE),
    unchanged = sum(case_corrected == case_uncorrected, na.rm = TRUE)
  )
cat("\n[Bridging correction impact, validation C1-unflagged stratum, unweighted counts]\n")
print(flip_summary)

save(flow_table, dev_sum, val_sum, dev_h_sum, val_h_sum, flip_summary,
     file = file.path(TABLES_DIR, "flow_and_prevalence.RData"))
