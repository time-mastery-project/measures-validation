# scoring/02_export_norms_to_tables.R
# ------------------------------------------------------------
# Build explicit, auditable conversion tables for the Sense of Time battery.
#
# OUTPUTS
# A) scoring/norms_lookup.csv
#    measure | age | raw | p | z_raw | z_aligned | ss_1_19
#
# B) scoring/ss_to_raw_intervals.csv
#    measure | age | ss_1_19 | raw_min | raw_max
#
# Notes for integer (discrete) measures:
# - Some SS may be unattainable at a given age, in that case raw_min/raw_max are NA.
# - Each integer raw is mapped to exactly one SS to avoid human confusion.
# ------------------------------------------------------------

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------
RDS_FULL <- here("scoring/normPercTabs_full.RDS")
RDS_CORE <- here("scoring/normPercTabs.RDS")
MEASURE_SPECS <- here("scoring/measure_specs.csv")

OUT_LOOKUP <- here("scoring/norms_lookup.csv")
OUT_SSINT  <- here("scoring/ss_to_raw_intervals.csv")

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
fail <- function(msg) stop(paste0("[FAIL] ", msg), call. = FALSE)
clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

pp_from_z <- function(z) {
  clamp(round(10 + 3 * z), 1, 19)
}

approx_rule2 <- function(x, y, xout) {
  approx(x = x, y = y, xout = xout, rule = 2)$y
}

# SS boundaries in aligned-z space for rounding rule:
# SS = round(10 + 3*z), clamped to [1, 19]
ss_to_z_bounds <- function(ss) {
  if (ss < 1 || ss > 19) stop("ss must be 1..19")
  z_lo <- if (ss == 1)  -Inf else (ss - 0.5 - 10) / 3
  z_hi <- if (ss == 19)  Inf else (ss + 0.5 - 10) / 3
  c(z_lo = z_lo, z_hi = z_hi)
}

# Convert aligned-z bounds to percentile bounds p (CDF of raw, before direction alignment).
zs_to_p_bounds <- function(z_lo, z_hi, direction_flag) {
  if (direction_flag == "higher_better") {
    p_lo <- pnorm(z_lo)
    p_hi <- pnorm(z_hi)
  } else if (direction_flag == "lower_better") {
    p_lo <- pnorm(-z_hi)
    p_hi <- pnorm(-z_lo)
  } else {
    fail(paste0("Unknown direction_flag: ", direction_flag))
  }
  c(p_lo = p_lo, p_hi = p_hi)
}

# ------------------------------------------------------------
# Load measure specs
# ------------------------------------------------------------
if (!file.exists(MEASURE_SPECS)) fail("measure_specs.csv not found")

specs <- read_csv(MEASURE_SPECS, show_col_types = FALSE)
names(specs) <- tolower(names(specs))

needed_specs <- c(
  "measure", "raw_type", "direction_flag",
  "raw_min_theoretical", "raw_max_theoretical",
  "abs_min", "abs_max", "raw_grid_n"
)
miss <- setdiff(needed_specs, names(specs))
if (length(miss) > 0) {
  fail(paste("measure_specs.csv missing columns:", paste(miss, collapse = ", ")))
}

specs <- specs %>%
  mutate(
    abs_min = as.numeric(abs_min),
    abs_max = ifelse(abs_max == "Inf", Inf, as.numeric(abs_max)),
    raw_grid_n = ifelse(is.na(raw_grid_n), 201L, as.integer(raw_grid_n))
  )

# ------------------------------------------------------------
# Load norms RDS
# ------------------------------------------------------------
rds_path <- if (file.exists(RDS_FULL)) RDS_FULL else RDS_CORE
if (!file.exists(rds_path)) fail("No norms RDS found (normPercTabs_full.RDS or normPercTabs.RDS)")

tab <- readRDS(rds_path)

# ------------------------------------------------------------
# Extract matrices to long format (robust, order-safe)
# ------------------------------------------------------------
extract_matrix <- function(mat, measure_name) {
  if (!is.matrix(mat)) return(NULL)
  
  p <- suppressWarnings(as.numeric(rownames(mat)))
  age <- suppressWarnings(as.numeric(gsub("^age-", "", colnames(mat))))
  
  if (anyNA(p) || anyNA(age)) {
    fail(paste("Invalid row/col names (percentiles/ages) in", measure_name))
  }
  if (!all(diff(p) > 0)) {
    fail(paste("Percentiles not strictly increasing in", measure_name))
  }
  if (!all(diff(age) > 0)) {
    fail(paste("Ages not strictly increasing in", measure_name))
  }
  
  age_chr <- format(age, scientific = FALSE, trim = TRUE)
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(df) <- paste0("age_", age_chr)
  df$p <- p
  
  df %>%
    pivot_longer(
      cols = starts_with("age_"),
      names_to = "age",
      values_to = "raw"
    ) %>%
    mutate(
      age = as.numeric(str_replace(age, "^age_", "")),
      measure = measure_name
    ) %>%
    select(measure, p, age, raw)
}

extract_any <- function(obj, name) {
  if (is.matrix(obj)) {
    return(extract_matrix(obj, name))
  }
  if (is.list(obj)) {
    nms <- names(obj)
    if (is.null(nms)) nms <- rep("", length(obj))
    
    map2_dfr(obj, seq_along(obj), function(el, i) {
      el_nm <- nms[[i]]
      if (is.null(el_nm) || el_nm == "") el_nm <- as.character(i)
      
      candidate <- el_nm
      if (!(candidate %in% specs$measure)) {
        candidate <- paste0(name, "_", el_nm)
      }
      extract_any(el, candidate)
    })
  } else {
    NULL
  }
}

raw_long <- imap_dfr(tab, ~ extract_any(.x, .y))
if (nrow(raw_long) == 0) fail("No valid matrices extracted from RDS")

# ------------------------------------------------------------
# Build norms_lookup.csv (percentile-based)
# ------------------------------------------------------------
lookup_list <- list()

for (m in unique(raw_long$measure)) {
  sp <- specs %>% filter(measure == m)
  if (nrow(sp) == 0) next
  sp <- sp[1, ]
  
  dat <- raw_long %>% filter(measure == m)
  
  p_clip <- clamp(dat$p, 1e-6, 1 - 1e-6)
  z_raw <- qnorm(p_clip)
  
  z_aligned <- if (sp$direction_flag == "lower_better") -z_raw else z_raw
  ss <- pp_from_z(z_aligned)
  
  lookup_list[[m]] <- tibble(
    measure = m,
    age = dat$age,
    raw = dat$raw,
    p = dat$p,
    z_raw = z_raw,
    z_aligned = z_aligned,
    ss_1_19 = ss
  )
}

norms_lookup <- bind_rows(lookup_list)
if (nrow(norms_lookup) == 0) fail("norms_lookup ended up empty (measure names mismatch?)")

# ------------------------------------------------------------
# Build ss_to_raw_intervals.csv
# - continuous: SS-specific percentile bounds -> raw quantiles
# - integer: deterministic raw->SS mapping (unique), then invert, unattainable SS -> NA
# ------------------------------------------------------------

build_ss_intervals_one <- function(m, a, sp, dat_pa) {
  dat_pa <- dat_pa %>%
    filter(is.finite(p), is.finite(raw)) %>%
    arrange(p)
  
  if (nrow(dat_pa) < 2) return(NULL)
  
  # Monotonicity of raw(p) is required
  dr <- diff(dat_pa$raw)
  if (any(dr < -1e-8, na.rm = TRUE)) {
    fail(paste0(
      "Non-monotone raw(p) for measure=", m,
      ", age=", a,
      ". This usually indicates extraction misalignment."
    ))
  }
  
  p_min <- min(dat_pa$p, na.rm = TRUE)
  p_max <- max(dat_pa$p, na.rm = TRUE)
  
  # Quantile function raw(p), linear is fine for continuous measures
  qraw <- function(pout) {
    pout2 <- clamp(pout, p_min, p_max)
    approx_rule2(dat_pa$p, dat_pa$raw, pout2)
  }
  
  # -------------------------
  # INTEGER MEASURES: unique raw -> SS, then invert
  # -------------------------
  if (sp$raw_type == "integer") {
    if (is.na(sp$raw_min_theoretical) || is.na(sp$raw_max_theoretical)) {
      fail(paste0("Integer measure ", m, " has NA theoretical bounds in measure_specs.csv"))
    }
    lo <- as.integer(sp$raw_min_theoretical)
    hi <- as.integer(sp$raw_max_theoretical)
    if (hi < lo) fail(paste0("Invalid theoretical bounds for ", m))
    
    # Use norms_lookup z_aligned to build a deterministic function z_aligned(raw)
    dz <- norms_lookup %>%
      filter(measure == m, age == a) %>%
      select(raw, z_aligned) %>%
      filter(is.finite(raw), is.finite(z_aligned)) %>%
      group_by(raw) %>%
      summarise(z_aligned = median(z_aligned), .groups = "drop") %>%
      arrange(raw)
    
    if (nrow(dz) < 2) {
      fail(paste0("Not enough points to build z(raw) for integer measure=", m, ", age=", a))
    }
    
    # Enforce monotonicity (z_aligned should follow performance direction)
    if (sp$direction_flag == "higher_better") {
      dz <- dz %>% mutate(z_aligned = cummax(z_aligned))
    } else {
      dz <- dz %>% mutate(z_aligned = cummin(z_aligned))
    }
    
    z_fun <- function(r) approx_rule2(dz$raw, dz$z_aligned, r)
    
    raw_vals <- lo:hi
    ss_vals <- pp_from_z(z_fun(raw_vals))
    
    raw_to_ss <- tibble(raw = raw_vals, ss_1_19 = ss_vals)
    
    # Invert to intervals, unattainable SS -> NA
    inv <- raw_to_ss %>%
      group_by(ss_1_19) %>%
      summarise(
        raw_min = min(raw),
        raw_max = max(raw),
        .groups = "drop"
      )
    
    out <- tibble(ss_1_19 = 1:19) %>%
      left_join(inv, by = "ss_1_19") %>%
      mutate(measure = m, age = a) %>%
      select(measure, age, ss_1_19, raw_min, raw_max)
    
    return(out)
  }
  
  # -------------------------
  # CONTINUOUS MEASURES: SS percentile bounds -> raw bounds
  # -------------------------
  if (sp$raw_type == "continuous") {
    if (!is.finite(sp$abs_min)) {
      fail(paste0("abs_min must be finite for continuous measure ", m))
    }
  }
  
  out <- map_dfr(1:19, function(ss) {
    zb <- ss_to_z_bounds(ss)
    pb <- zs_to_p_bounds(zb[["z_lo"]], zb[["z_hi"]], sp$direction_flag)
    
    p_lo <- pb[["p_lo"]]
    p_hi <- pb[["p_hi"]]
    
    raw_lo <- qraw(p_lo)
    raw_hi <- qraw(p_hi)
    
    tibble(
      measure = m,
      age = a,
      ss_1_19 = ss,
      raw_min = raw_lo,
      raw_max = raw_hi
    )
  })
  
  # Enforce absolute domain endpoints where requested
  if (sp$raw_type == "continuous") {
    if (sp$direction_flag == "lower_better") {
      out <- out %>%
        mutate(
          raw_min = ifelse(ss_1_19 == 19, sp$abs_min, raw_min),
          raw_max = ifelse(ss_1_19 == 1,  sp$abs_max, raw_max)
        )
    } else {
      out <- out %>%
        mutate(
          raw_min = ifelse(ss_1_19 == 1,  sp$abs_min, raw_min),
          raw_max = ifelse(ss_1_19 == 19, sp$abs_max, raw_max)
        )
    }
  }
  
  out %>% select(measure, age, ss_1_19, raw_min, raw_max)
}

ssint_list <- list()

for (m in unique(norms_lookup$measure)) {
  sp <- specs %>% filter(measure == m)
  if (nrow(sp) == 0) next
  sp <- sp[1, ]
  
  ages_m <- sort(unique(raw_long$age[raw_long$measure == m]))
  
  for (a in ages_m) {
    dat_pa <- raw_long %>%
      filter(measure == m, age == a) %>%
      select(p, raw)
    
    tib <- build_ss_intervals_one(m, a, sp, dat_pa)
    if (is.null(tib)) {
      fail(paste0(
        "Could not build SS intervals for measure=", m,
        ", age=", a,
        ". Check that the norms table contains at least two finite quantile points."
      ))
    }
    
    ssint_list[[paste(m, a, sep = "_")]] <- tib
  }
}

ss_to_raw_intervals <- bind_rows(ssint_list)
if (nrow(ss_to_raw_intervals) == 0) fail("ss_to_raw_intervals ended up empty")

# ------------------------------------------------------------
# Minimal internal checks before writing
# ------------------------------------------------------------
if (anyNA(norms_lookup$z_aligned) || anyNA(norms_lookup$ss_1_19)) {
  fail("NA in norms_lookup (z_aligned or ss_1_19)")
}

# Guarantee 19 SS rows for each (measure, age), but allow NA intervals for integer measures
coverage <- ss_to_raw_intervals %>%
  group_by(measure, age) %>%
  summarise(n_rows = n(), .groups = "drop")

if (any(coverage$n_rows != 19)) {
  bad <- coverage %>% filter(n_rows != 19)
  fail(paste0(
    "Expected exactly 19 rows for each (measure, age). First offenders: ",
    paste0(head(paste(bad$measure, bad$age, bad$n_rows, sep = ":"), 5), collapse = ", ")
  ))
}

# Check: for continuous measures, no NA allowed in intervals
cont_measures <- specs %>% filter(raw_type == "continuous") %>% pull(measure)

bad_cont <- ss_to_raw_intervals %>%
  filter(measure %in% cont_measures) %>%
  filter(is.na(raw_min) | is.na(raw_max))

if (nrow(bad_cont) > 0) {
  fail(paste0(
    "Continuous measures have NA intervals. First offenders: ",
    paste0(head(paste(bad_cont$measure, bad_cont$age, bad_cont$ss_1_19, sep = ":"), 5), collapse = ", ")
  ))
}

# Basic validity: where both finite, raw_min <= raw_max
bad_order <- ss_to_raw_intervals %>%
  filter(is.finite(raw_min), is.finite(raw_max)) %>%
  filter(raw_min > raw_max)

if (nrow(bad_order) > 0) {
  fail("Invalid SS->raw intervals: raw_min > raw_max for some rows")
}

# ------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------
write_csv(norms_lookup, OUT_LOOKUP)
write_csv(ss_to_raw_intervals, OUT_SSINT)

# ------------------------------------------------------------
# Transfer all scoring files also to shiny app folder
# ------------------------------------------------------------
dir.create("shiny/", recursive = TRUE, showWarnings = FALSE)
files <- list.files("scoring/", full.names = TRUE, recursive = FALSE)
files <- files[file.info(files)$isdir == FALSE]
file.copy(files, to = "shiny/", overwrite = TRUE)

# ------------------------------------------------------------
# Final message
# ------------------------------------------------------------
message("02_export_norms_to_tables.R completed successfully.")

