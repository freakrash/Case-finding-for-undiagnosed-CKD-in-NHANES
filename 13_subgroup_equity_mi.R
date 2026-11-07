# Re-runs the subgroup O/E analysis (Protocol §9.6) using the MI-based predictions from 11 instead
# of the complete-case M1 predictions from 08 -- this is the version that should be reported as
# primary, since MICE (not complete-case) is the protocol's pre-specified missing-data handling.

source("R/00_config.R")
val_headline <- readRDS(file.path(PROCESSED_DIR, "val_headline_mi_preds.rds"))
w <- val_headline$WTPH2YR

RIDRETH3_LABELS <- c(`1`="Mexican American", `2`="Other Hispanic", `3`="Non-Hispanic White",
                      `4`="Non-Hispanic Black", `6`="Non-Hispanic Asian", `7`="Other/Multiracial")

weighted_net_benefit <- function(treat, outcome, w, pt) {
  tp <- sum(w[treat & outcome]); fp <- sum(w[treat & !outcome]); n <- sum(w)
  (tp/n) - (fp/n) * (pt/(1-pt))
}
subgroup_oe <- function(mask, label) {
  ww <- w[mask]; oo <- val_headline$ckd_case[mask]; pp <- val_headline$m1_mi_pred[mask]
  n_unw <- sum(mask); events_unw <- sum(oo)
  if (n_unw < 5) return(data.frame(level=label, n_unw=n_unw, events_unw=events_unw, oe=NA, nb_5pct=NA))
  obs <- sum(ww*oo)/sum(ww); exp_ <- sum(ww*pp)/sum(ww)
  nb5 <- weighted_net_benefit(pp >= DCA_PRIMARY, oo, ww, DCA_PRIMARY)
  data.frame(level=label, n_unw=n_unw, events_unw=events_unw, oe=obs/exp_, nb_5pct=nb5)
}

cat("=== SUBGROUP O/E (MICE-based predictions), C1-unflagged headline stratum, validation ===\n")
race_tab <- do.call(rbind, lapply(sort(unique(val_headline$race_eth)), function(lv)
  subgroup_oe(val_headline$race_eth == lv, RIDRETH3_LABELS[as.character(lv)])))
cat("\n-- Race/ethnicity --\n"); print(race_tab)

income_cat <- cut(val_headline$income_pir, breaks = INCOME_BREAKS, labels = INCOME_LABELS)
income_tab <- do.call(rbind, lapply(levels(income_cat), function(lv)
  subgroup_oe(income_cat == lv & !is.na(income_cat), lv)))
cat("\n-- Income --\n"); print(income_tab)

ins_tab <- do.call(rbind, list(
  subgroup_oe(!is.na(val_headline$hiq011) & val_headline$hiq011 == 1, "Insured"),
  subgroup_oe(!is.na(val_headline$hiq011) & val_headline$hiq011 == 2, "Uninsured")
))
cat("\n-- Insurance --\n"); print(ins_tab)

log_oe_summary <- function(tab) {
  ok <- !is.na(tab$oe) & tab$oe > 0
  logs <- log(tab$oe[ok]); ws <- tab$n_unw[ok]
  sqrt(sum(ws * (logs - weighted.mean(logs, ws))^2) / sum(ws))
}
cat(sprintf("\nWeighted SD log(O/E) across race/ethnicity (MI-based): %.3f\n", log_oe_summary(race_tab)))

write.csv(race_tab, file.path(TABLES_DIR, "subgroup_oe_race_mi.csv"), row.names = FALSE)
write.csv(income_tab, file.path(TABLES_DIR, "subgroup_oe_income_mi.csv"), row.names = FALSE)
write.csv(ins_tab, file.path(TABLES_DIR, "subgroup_oe_insurance_mi.csv"), row.names = FALSE)
