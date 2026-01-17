# scoring/03_sanity_check_norms.R
# ------------------------------------------------------------
# Sanity checks for exported norms tables
# VALIDATES THE FORMAL CONTRACT OF 02
# FAIL-FAST, NO HEURISTICS, NO GUESSING
# ------------------------------------------------------------

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
})

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------
MEASURE_SPECS <- here("scoring/measure_specs.csv")
NORMS_LOOKUP  <- here("scoring/norms_lookup.csv")
SS_INTERVALS  <- here("scoring/ss_to_raw_intervals.csv")

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
fail <- function(msg) stop(paste0("[FAIL] ", msg), call. = FALSE)

# ------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------
if (!file.exists(MEASURE_SPECS)) fail("measure_specs.csv not found")
if (!file.exists(NORMS_LOOKUP))  fail("norms_lookup.csv not found")
if (!file.exists(SS_INTERVALS))  fail("ss_to_raw_intervals.csv not found")

specs <- read_csv(MEASURE_SPECS, show_col_types = FALSE)
norms <- read_csv(NORMS_LOOKUP,  show_col_types = FALSE)
ssint <- read_csv(SS_INTERVALS,  show_col_types = FALSE)

names(specs) <- tolower(names(specs))
names(norms) <- tolower(names(norms))
names(ssint) <- tolower(names(ssint))

# ------------------------------------------------------------
# Validate measure_specs
# ------------------------------------------------------------
needed_specs <- c(
  "measure", "raw_type", "direction_flag",
  "raw_min_theoretical", "raw_max_theoretical",
  "abs_min", "abs_max"
)
miss <- setdiff(needed_specs, names(specs))
if (length(miss) > 0)
  fail(paste("measure_specs.csv missing columns:", paste(miss, collapse = ", ")))

specs <- specs %>%
  mutate(
    abs_min = as.numeric(abs_min),
    abs_max = ifelse(abs_max == "Inf", Inf, as.numeric(abs_max))
  )

# ------------------------------------------------------------
# 1) Core checks on norms_lookup.csv
# ------------------------------------------------------------
needed_norms <- c("measure", "age", "raw", "p", "z_raw", "z_aligned", "ss_1_19")
miss <- setdiff(needed_norms, names(norms))
if (length(miss) > 0)
  fail(paste("norms_lookup.csv missing columns:", paste(miss, collapse = ", ")))

if (anyNA(norms$z_aligned))
  fail("NA found in z_aligned")

if (anyNA(norms$ss_1_19))
  fail("NA found in ss_1_19")

if (any(norms$ss_1_19 < 1 | norms$ss_1_19 > 19))
  fail("SS values outside [1,19] detected")

# IMPORTANT:
# No monotonicity checks on norms_lookup.csv.
# norms_lookup is percentile-based and discretized.
# Local inversions are allowed and expected.

# ------------------------------------------------------------
# 2) Validate ss_to_raw_intervals.csv structure
# ------------------------------------------------------------
needed_ssint <- c("measure", "age", "ss_1_19", "raw_min", "raw_max")
miss <- setdiff(needed_ssint, names(ssint))
if (length(miss) > 0)
  fail(paste("ss_to_raw_intervals.csv missing columns:", paste(miss, collapse = ", ")))

if (anyNA(ssint$raw_min) || anyNA(ssint$raw_max))
  fail("NA found in SS->raw intervals")

if (any(ssint$raw_min > ssint$raw_max))
  fail("Found SS->raw intervals with raw_min > raw_max")

# ------------------------------------------------------------
# 3) Domain coverage and monotonicity on SS -> raw intervals
# ------------------------------------------------------------
for (m in unique(ssint$measure)) {
  
  sp <- specs %>% filter(measure == m)
  if (nrow(sp) == 0) next
  
  for (a in unique(ssint$age[ssint$measure == m])) {
    
    sub <- ssint %>%
      filter(measure == m, age == a) %>%
      arrange(ss_1_19)
    
    if (nrow(sub) == 0) next
    
    covered_min <- min(sub$raw_min, na.rm = TRUE)
    covered_max <- max(sub$raw_max, na.rm = TRUE)
    
    # -------------------------
    # Domain coverage
    # -------------------------
    if (sp$raw_type == "continuous") {
      
      if (covered_min > sp$abs_min + 1e-8)
        fail(paste("SS->raw lower gap for", m, "age", a))
      
      if (is.finite(sp$abs_max) &&
          covered_max < sp$abs_max - 1e-8)
        fail(paste("SS->raw upper gap for", m, "age", a))
    }
    
    if (sp$raw_type == "integer") {
      
      expected <- seq(
        sp$raw_min_theoretical,
        sp$raw_max_theoretical,
        by = 1
      )
      
      covered <- seq(
        floor(covered_min),
        ceiling(covered_max),
        by = 1
      )
      
      if (!all(expected %in% covered))
        fail(paste(
          "SS->raw intervals do not cover full integer domain for",
          m, "age", a
        ))
    }
    
    # -------------------------
    # Monotonicity of intervals (ONLY HERE)
    # -------------------------
    if (sp$raw_type == "continuous") {
      
      if (sp$direction_flag == "lower_better") {
        # higher SS = better -> raw must decrease
        if (any(diff(sub$raw_min) > 1e-8))
          fail(paste("Non-monotone SS->raw (raw_min) for", m, "age", a))
      }
      
      if (sp$direction_flag == "higher_better") {
        # higher SS = better -> raw must increase
        if (any(diff(sub$raw_max) < -1e-8))
          fail(paste("Non-monotone SS->raw (raw_max) for", m, "age", a))
      }
    }
  }
}

# ------------------------------------------------------------
# Summary (informative, does not change pass/fail logic)
# ------------------------------------------------------------
cat("\n--- SUMMARY ---\n")

cat("Measures in specs:", n_distinct(specs$measure), "\n")
cat("Measures in norms_lookup:", n_distinct(norms$measure), "\n")
cat("Measures in ss_to_raw_intervals:", n_distinct(ssint$measure), "\n")

cat("\nAges (norms_lookup):", length(unique(norms$age)),
    "range:", min(norms$age), "-", max(norms$age), "\n")

cat("Ages (ss_to_raw_intervals):", length(unique(ssint$age)),
    "range:", min(ssint$age), "-", max(ssint$age), "\n")

cat("\nRows:\n")
cat("  norms_lookup:", nrow(norms), "\n")
cat("  ss_to_raw_intervals:", nrow(ssint), "\n")

# quick check: do all measures in specs appear in both outputs?
missing_lookup <- setdiff(specs$measure, unique(norms$measure))
missing_ssint  <- setdiff(specs$measure, unique(ssint$measure))

if (length(missing_lookup) > 0) {
  cat("\n[WARN] Measures in specs missing from norms_lookup:\n  ",
      paste(missing_lookup, collapse = ", "), "\n")
} else {
  cat("\n[OK] All measures in specs are present in norms_lookup\n")
}

if (length(missing_ssint) > 0) {
  cat("[WARN] Measures in specs missing from ss_to_raw_intervals:\n  ",
      paste(missing_ssint, collapse = ", "), "\n")
} else {
  cat("[OK] All measures in specs are present in ss_to_raw_intervals\n")
}

cat("\nSS coverage per measure (interval table):\n")
ss_cov <- ssint %>%
  group_by(measure) %>%
  summarise(
    n_age = n_distinct(age),
    ss_min = min(ss_1_19),
    ss_max = max(ss_1_19),
    .groups = "drop"
  )
print(ss_cov, n = Inf)

cat("\n--- END SUMMARY ---\n\n")

message("03_sanity_check_norms.R completed successfully.")
