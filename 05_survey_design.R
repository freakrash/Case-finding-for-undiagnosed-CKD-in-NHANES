# Builds the NHANES complex-survey design objects (Protocol §7.2). Built ONCE per role on every
# row of that role's pooled file (not just the analytic domain) -- every downstream restriction
# (undiagnosed, C1-unflagged stratum, subgroup) is applied via subset()/domain estimation on top of
# this object, per the protocol's explicit warning: "this is the single most common survey-analysis
# error in this literature."
#
# NHANES SDMVSTRA/SDMVPSU codes are only unique WITHIN a release cycle, so pooling development
# cycles (G/H/I/P) requires combining cycle+stratum (and cycle+stratum+psu) into unique ids before
# building the design object -- standard NHANES multi-cycle pooling practice.

source("R/00_config.R")
library(survey)
library(dplyr)

pooled <- readRDS(file.path(PROCESSED_DIR, "analytic_pooled.rds")) %>%
  mutate(
    combined_stratum = paste(cycle, sdmvstra, sep = "_"),
    combined_psu      = paste(cycle, sdmvstra, sdmvpsu, sep = "_")
  )

dev_df <- pooled %>% filter(role == "development")
val_df_all <- pooled %>% filter(role == "validation")

# WTPH2YR is undefined (NA), not zero, for participants outside the phlebotomy-weighting frame
# (e.g. under the MEC phlebotomy age cutoff). This is a data-existence restriction on the weight
# variable itself, not a pre-filter of the scientific analytic domain: every in_analytic_domain row
# already has a valid, non-missing WTPH2YR (verified below), so this does not touch the undiagnosed/
# C1-unflagged/subgroup restrictions the protocol's domain-subsetting rule protects (§5.5, §11).
stopifnot(all(!is.na(val_df_all$WTPH2YR[val_df_all$in_analytic_domain])))
val_df <- val_df_all %>% filter(!is.na(WTPH2YR))
n_dropped_no_ph_weight <- nrow(val_df_all) - nrow(val_df)
cat(sprintf("[05_survey_design] validation rows outside the WTPH2YR frame (dropped from val design, none analytic-domain): %d\n",
            n_dropped_no_ph_weight))

stopifnot(nrow(dev_df) > 0, nrow(val_df) > 0)
stopifnot(!any(is.na(dev_df$wt_dev)))
stopifnot(!any(is.na(val_df$WTPH2YR)))

# --- Development design (Taylor linearisation; pooled MEC/pre-pandemic weight) --------------------
dev_design <- svydesign(
  id = ~combined_psu, strata = ~combined_stratum, weights = ~wt_dev,
  nest = TRUE, data = dev_df
)

# --- Validation design, primary weight WTPH2YR (Protocol §7.2 "Validation weight (primary)") ------
val_design <- svydesign(
  id = ~combined_psu, strata = ~combined_stratum, weights = ~WTPH2YR,
  nest = TRUE, data = val_df
)

# --- Validation design, sensitivity weight WTMEC2YR (Protocol §7.2 sensitivity; Amendment S3) -----
val_design_mec <- svydesign(
  id = ~combined_psu, strata = ~combined_stratum, weights = ~wt_mec,
  nest = TRUE, data = val_df
)

# --- Rao-Wu-type rescaled bootstrap replicate designs (Protocol §7.2: "Rao-Wu rescaled bootstrap
# resampling PSUs within strata"; BOOT_REPLICATES set by RUN_SCALE in 00_config.R) -----------------
# survey::as.svrepdesign(type="mrbbootstrap") implements the Rao-Wu-Yue multistage rescaling
# bootstrap for stratified multistage designs, which is what "Rao-Wu" refers to here.
build_boot_design <- function(design, R = BOOT_REPLICATES) {
  as.svrepdesign(design, type = "mrbbootstrap", replicates = R)
}

cat(sprintf("[05_survey_design] development design: n=%d, %d strata, %d PSUs\n",
            nrow(dev_df), length(unique(dev_df$combined_stratum)), length(unique(dev_df$combined_psu))))
cat(sprintf("[05_survey_design] validation design: n=%d, %d strata, %d PSUs\n",
            nrow(val_df), length(unique(val_df$combined_stratum)), length(unique(val_df$combined_psu))))

# Sanity: undiagnosed, C1-unflagged headline subpopulation size (domain estimate, not pre-filtered)
headline_domain <- subset(dev_design, in_analytic_domain & !c1_flag)
cat(sprintf("[05_survey_design] dev C1-unflagged analytic domain (unweighted n): %d\n", nrow(headline_domain$variables)))
headline_domain_val <- subset(val_design, in_analytic_domain & !c1_flag)
cat(sprintf("[05_survey_design] validation C1-unflagged analytic domain (unweighted n): %d\n", nrow(headline_domain_val$variables)))
