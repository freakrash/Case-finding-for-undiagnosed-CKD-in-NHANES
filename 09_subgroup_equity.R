# Subgroup calibration (Protocol §9.6, "primary, not supplementary"; Amendment §4). Weighted O/E
# ratio reported at every level regardless of event count (Amendment: "never post hoc" collapsing).
# Draft-pass scope: O/E ratio + subgroup net benefit at p_t=5% for every level; calibration slope
# with the >=50-event rule and collapsing hierarchy is deferred to the full run (needs the full
# imputed dataset to have stable per-level event counts worth reporting).

source("R/00_config.R")
source("R/05_survey_design.R")
source("R/08_decision_curve_calibration.R")  # reuses calib_df / headline_design / w / DCA_PRIMARY

RIDRETH3_LABELS <- c(`1`="Mexican American", `2`="Other Hispanic", `3`="Non-Hispanic White",
                      `4`="Non-Hispanic Black", `6`="Non-Hispanic Asian", `7`="Other/Multiracial")

subgroup_oe <- function(mask, label) {
  ww <- w[mask]; oo <- calib_df$ckd_case[mask]; pp <- calib_df$m1_pred[mask]
  n_unw <- sum(mask); events_unw <- sum(oo)
  if (n_unw < 5) return(data.frame(level=label, n_unw=n_unw, events_unw=events_unw, oe=NA, oe_se=NA, nb_5pct=NA))
  obs <- sum(ww*oo)/sum(ww); exp_ <- sum(ww*pp)/sum(ww)
  nb5 <- weighted_net_benefit(pp >= DCA_PRIMARY, oo, ww, DCA_PRIMARY)
  data.frame(level=label, n_unw=n_unw, events_unw=events_unw, oe=obs/exp_, nb_5pct=nb5)
}

cat("\n=== SUBGROUP O/E + net benefit @5%, C1-unflagged headline stratum, validation (M1, frozen) ===\n")

race_levels <- unique(calib_df$race_eth)
race_tab <- do.call(rbind, lapply(sort(race_levels), function(lv)
  subgroup_oe(calib_df$race_eth == lv, RIDRETH3_LABELS[as.character(lv)])))
cat("\n-- Race/ethnicity --\n"); print(race_tab)

income_cat <- cut(calib_df$income_pir, breaks = INCOME_BREAKS, labels = INCOME_LABELS)
income_tab <- do.call(rbind, lapply(levels(income_cat), function(lv)
  subgroup_oe(income_cat == lv & !is.na(income_cat), lv)))
cat("\n-- Family income-to-poverty ratio --\n"); print(income_tab)

ins_tab <- do.call(rbind, list(
  subgroup_oe(calib_df$hiq011 == 1, "Insured"),
  subgroup_oe(calib_df$hiq011 == 2, "Uninsured")
))
cat("\n-- Insurance --\n"); print(ins_tab)

# Amendment §5: primary cross-subgroup summary = weighted SD of log(O/E) across levels (replaces
# max-abs-difference, which is a biased-upward max over noisy estimates).
log_oe_summary <- function(tab) {
  ok <- !is.na(tab$oe) & tab$oe > 0
  if (sum(ok) < 2) return(NA)
  logs <- log(tab$oe[ok]); ws <- tab$n_unw[ok]
  wm <- weighted.mean(logs, ws)
  sqrt(sum(ws * (logs - wm)^2) / sum(ws))
}
cat(sprintf("\nWeighted SD of log(O/E) across race/ethnicity levels (Amendment §5 primary summary): %.3f\n", log_oe_summary(race_tab)))
cat(sprintf("Max absolute O/E difference across race/ethnicity levels (secondary descriptor): %.3f\n",
            diff(range(race_tab$oe, na.rm=TRUE))))

write.csv(race_tab, file.path(TABLES_DIR, "subgroup_oe_race.csv"), row.names = FALSE)
write.csv(income_tab, file.path(TABLES_DIR, "subgroup_oe_income.csv"), row.names = FALSE)
write.csv(ins_tab, file.path(TABLES_DIR, "subgroup_oe_insurance.csv"), row.names = FALSE)
