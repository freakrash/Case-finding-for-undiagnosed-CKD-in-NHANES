# Downloads every NHANES file required by the protocol (§8.1-8.2) for every cycle in the
# development/validation split (§7.1), and reads each .XPT into an RDS cache under data/raw/.
# Idempotent: skips files already downloaded/parsed.

source("R/00_config.R")
library(haven)

download_one <- function(stem, cyc) {
  fname <- nhanes_file_name(stem, cyc)
  xpt_path <- file.path(RAW_DIR, paste0(fname, ".XPT"))
  rds_path <- file.path(RAW_DIR, paste0(fname, ".rds"))
  if (file.exists(rds_path)) return(invisible(rds_path))

  url <- nhanes_url(stem, cyc)
  if (!file.exists(xpt_path)) {
    ok <- tryCatch({
      utils::download.file(url, xpt_path, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) { message(sprintf("  [FAIL] %s: %s", url, conditionMessage(e))); FALSE })
    if (!ok) return(invisible(NULL))
  }

  df <- tryCatch(haven::read_xpt(xpt_path), error = function(e) {
    message(sprintf("  [PARSE FAIL] %s: %s", xpt_path, conditionMessage(e))); NULL
  })
  if (is.null(df)) return(invisible(NULL))
  saveRDS(df, rds_path)
  invisible(rds_path)
}

manifest <- expand.grid(stem = names(FILE_STEMS), cyc = names(CYCLES), stringsAsFactors = FALSE)
# BPXO (oscillometric BP) only exists for validation cycle L; BPX (auscultatory) for development.
bp_rows <- data.frame(
  stem = c(rep("BPX", 4), "BPXO"),
  cyc  = c("G","H","I","P","L"),
  stringsAsFactors = FALSE
)
manifest <- rbind(manifest, bp_rows)

results <- character(0)
for (i in seq_len(nrow(manifest))) {
  stem <- manifest$stem[i]; cyc <- manifest$cyc[i]
  fname <- nhanes_file_name(stem, cyc)
  cat(sprintf("[%3d/%3d] %s ... ", i, nrow(manifest), fname))
  path <- download_one(stem, cyc)
  status <- if (is.null(path)) "MISSING" else "ok"
  cat(status, "\n")
  results <- c(results, sprintf("%s\t%s\t%s", fname, cyc, status))
}

writeLines(results, file.path(RAW_DIR, "_download_manifest.tsv"))
missing <- grep("MISSING$", results, value = TRUE)
if (length(missing)) {
  cat("\n=== FILES NOT FOUND (check file/variable names for this cycle) ===\n")
  cat(paste(missing, collapse = "\n"), "\n")
} else {
  cat("\nAll files downloaded/parsed successfully.\n")
}
