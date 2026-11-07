# Urinary albumin bridging correction (Protocol §7.3; Amendment §3.1-3.3, §S2).
# Applies ONLY to validation-cycle (L) URXUMA before ACR is recomputed, so the outcome in the
# validation set is defined on the same measurement scale as development (Protocol §7.3 "Primary
# validation"). Sensitivity A = uncorrected; Sensitivity B = coefficient uncertainty propagated (§S2).

source("R/00_config.R")

#' Backward Deming correction: new-method (LC-MS/MS) -> old-method-equivalent scale, with the
#' floor rule (Amendment §3.2): any corrected value below 0.212 mg/L (including negative values
#' and small positive values under the floor) is set to 0.212 mg/L.
#' @param urxuma_new numeric vector, new-method urinary albumin (mg/L)
#' @param intercept,slope Deming coefficients (default: confirmed values, see 00_config.R)
apply_bridging_correction <- function(urxuma_new,
                                       intercept = BRIDGE_INTERCEPT,
                                       slope = BRIDGE_SLOPE,
                                       floor_value = FLOOR_VALUE) {
  corrected <- intercept + slope * urxuma_new
  pmax(corrected, floor_value)
}

#' Crossover point (Amendment §3.1): reported (new-method) albumin value at which the correction
#' is self-consistent (corrected == reported). Below it the correction lowers the value; above it,
#' the correction raises it.
bridging_crossover_albumin <- function(intercept = BRIDGE_INTERCEPT, slope = BRIDGE_SLOPE) {
  intercept / (1 - slope)
}

#' Equivalent crossover expressed in urine creatinine (mg/dL), since the case-defining threshold
#' is on ACR = albumin/creatinine*100 >= 30 mg/g, i.e. albumin = 0.3 x creatinine at the threshold.
bridging_crossover_creatinine <- function(intercept = BRIDGE_INTERCEPT, slope = BRIDGE_SLOPE) {
  bridging_crossover_albumin(intercept, slope) / 0.3
}

#' Sensitivity B (§S2 / Amendment §6): propagate Deming coefficient CI through the correction,
#' returning corrected values at each corner of the intercept/slope confidence box, for flip-count
#' analysis (case status changes within the CI on the coefficients).
apply_bridging_correction_ci <- function(urxuma_new,
                                          intercept_ci = BRIDGE_INTERCEPT_CI,
                                          slope_ci = BRIDGE_SLOPE_CI,
                                          floor_value = FLOOR_VALUE) {
  corners <- expand.grid(intercept = intercept_ci, slope = slope_ci)
  out <- lapply(seq_len(nrow(corners)), function(i) {
    apply_bridging_correction(urxuma_new, corners$intercept[i], corners$slope[i], floor_value)
  })
  do.call(cbind, out)  # n x 4 matrix, one column per (intercept,slope) corner
}

# --- Self-test against the worked examples in Amendment §3.1-3.2 -----------
.bridging_selftest <- function() {
  stopifnot(abs(bridging_crossover_albumin() - 8.69) < 0.02)
  stopifnot(abs(bridging_crossover_creatinine() - 29.0) < 0.5)

  # 1.5 mg/L reported -> raw corrected ~0.141 (positive, but still under the 0.212 floor) -> floors to 0.212
  raw_15 <- BRIDGE_INTERCEPT + BRIDGE_SLOPE * 1.5
  stopifnot(abs(raw_15 - 0.141) < 0.01)
  stopifnot(abs(apply_bridging_correction(1.5) - FLOOR_VALUE) < 1e-9)
  # 1.0 mg/L reported -> raw corrected negative (~-0.454) -> floors to 0.212
  raw_10 <- BRIDGE_INTERCEPT + BRIDGE_SLOPE * 1.0
  stopifnot(raw_10 < 0)
  stopifnot(abs(apply_bridging_correction(1.0) - FLOOR_VALUE) < 1e-9)
  # The raw correction is monotonic (slope > 1: larger input -> larger raw output), so no inversion
  # exists in the correction itself. The inversion Amendment §3.2 warns about comes from the OLD,
  # naive floor rule ("floor only negative values"): that rule would leave 1.5's positive-but-
  # sub-floor output (0.141) UNfloored while flooring 1.0's negative output (-0.454) to 0.212,
  # inverting the ordering (0.141 < 0.212 even though input 1.5 > input 1.0). Our corrected rule
  # floors anything below FLOOR_VALUE regardless of sign, so both collapse identically instead.
  stopifnot(raw_15 > raw_10)                              # raw correction itself IS monotonic
  naive_floor <- function(x) ifelse(x < 0, FLOOR_VALUE, x)  # the old, buggy rule
  stopifnot(naive_floor(raw_15) < naive_floor(raw_10))       # ...and it really does invert
  stopifnot(apply_bridging_correction(1.5) == apply_bridging_correction(1.0))  # corrected rule: no inversion

  # New detection limit (0.02 mg/L) -> corrected ~ -1.62 -> floors to 0.212
  stopifnot(abs(apply_bridging_correction(NEW_METHOD_LLOD_UMA) - FLOOR_VALUE) < 1e-9)

  # No case-status determination is affected in the 1.382-1.560 window (Amendment §3.2): floored
  # value 0.212 mg/L must correspond to ACR < 30 mg/g even at the lowest plausible urine creatinine.
  # ACR = (0.212 / URXUCR) * 100; solve URXUCR at which this equals 30: URXUCR = 0.212*100/30 = 0.707 mg/dL.
  # Any physiologically plausible urine creatinine exceeds that, so ACR stays far below 30 at the floor.
  stopifnot((FLOOR_VALUE / 1.1) * 100 < 30)  # even at the L-cycle's own creatinine LLOD (1.1 mg/dL)

  cat("[02_bridging_correction] self-test OK. Crossover: albumin=", round(bridging_crossover_albumin(),3),
      "mg/L, urine creatinine=", round(bridging_crossover_creatinine(),1), "mg/dL\n")
}
.bridging_selftest()
