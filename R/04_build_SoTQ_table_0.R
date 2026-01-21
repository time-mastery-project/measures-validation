# R/04_build_SoTQ_table.R
# ------------------------------------------------------------
# Build a table-based Sense of Time Quotient (SoTQ), replicable manually.
#
# Outputs:
#  - scoring/SoT_totalConversionTable.csv
#      sumSS (sum of 6 SS indicators) -> SoTQ (M=100, SD=15)
#  - scoring/SoTQ_profileSpreadThreshold.csv
#      deltaSS = max(SS) - min(SS), 95th percentile (heterogeneity flag)
#
# Paper-aligned 6 indicators for SoTQ:
#   1) TE composite = mean of aligned z from TE_Barca_dd and TE_Ladro_dd
#   2) TR = AverageDevAbs_TR  (temporary, item-level 2..12 s to be added later)
#   3) TD = RatioTD mapped via measure "TD"
#   4) Child questionnaire total = OTm_total
#   5) Parent questionnaire total = QSTp_Total_parent mapped via "QST_parent_total"
#   6) Teacher questionnaire total = QSTp_Total_teacher mapped via "QST_teacher_total"
#
# IMPORTANT:
# - Uses scoring/norms_lookup.csv (new pipeline), not the old exported grid.
# - SoTQ table is derived from the empirical correlation matrix of the 6 SS indicators.
# ------------------------------------------------------------

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(readxl)
  library(tibble)
})

# -----------------------------
# Config
# -----------------------------
NORMS_LOOKUP <- here("scoring/norms_lookup.csv")
DATA_XLSX    <- here("data/Database_Time_16_05_25_PD_MI.xlsx")

OUT_TABLE    <- here("scoring/SoT_totalConversionTable.csv")
OUT_SPREAD   <- here("scoring/SoTQ_profileSpreadThreshold.csv")

STOP_IF_MISSING_TR_TOTAL <- TRUE
SOTQ_CLAMP <- c(40, 160)  # set to c(-Inf, Inf) to disable clamping

# -----------------------------
# Helpers
# -----------------------------
need_file <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
}

to_num <- function(x) {
  if (is.null(x)) return(NA_real_)
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x <- trimws(x)
    x <- gsub(",", ".", x)
  }
  suppressWarnings(as.numeric(x))
}

clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

# linear interpolation with boundary capping
lin_interp <- function(x, y, x0) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  x0 <- suppressWarnings(as.numeric(x0))
  
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (!length(x) || !is.finite(x0)) return(NA_real_)
  
  if (any(x == x0)) return(y[which(x == x0)[1]])
  if (x0 <= min(x)) return(y[which.min(x)])
  if (x0 >= max(x)) return(y[which.max(x)])
  
  xs <- sort(unique(x))
  hi <- which(xs > x0)[1]
  lo <- hi - 1
  x1 <- xs[lo]; x2 <- xs[hi]
  
  y1 <- y[which(x == x1)[1]]
  y2 <- y[which(x == x2)[1]]
  
  y1 + (y2 - y1) * (x0 - x1) / (x2 - x1)
}

nearest_age <- function(age_value, age_grid) {
  age_value <- suppressWarnings(as.numeric(age_value))
  age_grid <- suppressWarnings(as.numeric(age_grid))
  age_grid <- age_grid[is.finite(age_grid)]
  if (!is.finite(age_value) || !length(age_grid)) return(NA_real_)
  age_grid[which.min(abs(age_grid - age_value))]
}

# z -> SS (1-19), vectorized and safe
z_to_ss19 <- function(z) {
  z <- suppressWarnings(as.numeric(z))
  ss <- round(10 + 3 * z)
  ss[!is.finite(ss)] <- NA_real_
  clamp(ss, 1, 19)
}

lookup_z_aligned <- function(norms, measure, age, raw) {
  if (is.null(measure) || !nzchar(measure)) return(NA_real_)
  age <- suppressWarnings(as.numeric(age))
  raw <- suppressWarnings(as.numeric(raw))
  if (!is.finite(age) || !is.finite(raw)) return(NA_real_)
  
  n_m <- norms %>% filter(.data$measure == !!measure)
  if (nrow(n_m) == 0) return(NA_real_)
  
  a0 <- nearest_age(age, unique(n_m$age))
  if (!is.finite(a0)) return(NA_real_)
  
  g <- n_m %>%
    filter(.data$age == !!a0) %>%
    transmute(raw = to_num(.data$raw), z_aligned = to_num(.data$z_aligned)) %>%
    filter(is.finite(.data$raw), is.finite(.data$z_aligned))
  
  if (nrow(g) < 2) return(NA_real_)
  
  # ensure strictly increasing x for interpolation (collapse duplicates)
  g2 <- g %>%
    group_by(.data$raw) %>%
    summarise(z_aligned = mean(.data$z_aligned), .groups = "drop") %>%
    arrange(.data$raw)
  
  if (nrow(g2) < 2) return(NA_real_)
  lin_interp(g2$raw, g2$z_aligned, raw)
}

# -----------------------------
# Load inputs
# -----------------------------
need_file(NORMS_LOOKUP)
need_file(DATA_XLSX)

norms <- suppressMessages(read_csv(NORMS_LOOKUP, show_col_types = FALSE))
names(norms) <- tolower(names(norms))

required_norms_cols <- c("measure", "age", "raw", "z_aligned", "ss_1_19")
miss_norms <- setdiff(required_norms_cols, names(norms))
if (length(miss_norms)) {
  stop("norms_lookup.csv missing columns: ", paste(miss_norms, collapse = ", "), call. = FALSE)
}

norms <- norms %>%
  mutate(
    age = to_num(.data$age),
    raw = to_num(.data$raw),
    z_aligned = to_num(.data$z_aligned),
    ss_1_19 = to_num(.data$ss_1_19)
  ) %>%
  filter(is.finite(.data$age), is.finite(.data$raw), is.finite(.data$z_aligned))

dat <- as.data.frame(read_excel(DATA_XLSX))

# optional inclusion filter used in your pipeline
if ("Class" %in% names(dat)) dat <- dat[!is.na(dat$Class), , drop = FALSE]

if (!("Age" %in% names(dat))) stop("Dataset missing required column: Age", call. = FALSE)
dat$Age <- to_num(dat$Age)

# -----------------------------
# Build raw variables needed for the 6 SoTQ indicators
# -----------------------------
# Child questionnaire total (OTm_total)
if (all(paste0("OTm_", 1:16) %in% names(dat))) {
  dat$OTm_total <- rowSums(dat[, paste0("OTm_", 1:16)], na.rm = FALSE)
} else if (!("OTm_total" %in% names(dat))) {
  stop("Dataset missing OTm_total and OTm_1..OTm_16, cannot compute child total.", call. = FALSE)
} else {
  dat$OTm_total <- to_num(dat$OTm_total)
}

# TE raws must exist (absolute deviation for each video)
if (!("TE_Barca_dd" %in% names(dat))) stop("Dataset missing TE_Barca_dd", call. = FALSE)
if (!("TE_Ladro_dd" %in% names(dat))) stop("Dataset missing TE_Ladro_dd", call. = FALSE)
dat$TE_Barca_dd <- to_num(dat$TE_Barca_dd)
dat$TE_Ladro_dd <- to_num(dat$TE_Ladro_dd)

# TD raw
if (!("RatioTD" %in% names(dat))) stop("Dataset missing RatioTD", call. = FALSE)
dat$RatioTD <- to_num(dat$RatioTD)

# Parent/Teacher totals in dataset
if (!("QSTp_Total_parent" %in% names(dat))) stop("Dataset missing QSTp_Total_parent", call. = FALSE)
if (!("QSTp_Total_teacher" %in% names(dat))) stop("Dataset missing QSTp_Total_teacher", call. = FALSE)
dat$QSTp_Total_parent  <- to_num(dat$QSTp_Total_parent)
dat$QSTp_Total_teacher <- to_num(dat$QSTp_Total_teacher)

# TR total (temporary)
has_TR_total <- "AverageDevAbs_TR" %in% names(dat)
if (has_TR_total) dat$AverageDevAbs_TR <- to_num(dat$AverageDevAbs_TR)

if (!has_TR_total && isTRUE(STOP_IF_MISSING_TR_TOTAL)) {
  stop("Dataset missing AverageDevAbs_TR. Add it, or set STOP_IF_MISSING_TR_TOTAL = FALSE.", call. = FALSE)
}

# -----------------------------
# Check that measures exist in norms_lookup
# -----------------------------
needed_measures <- c(
  "TE_Barca_dd", "TE_Ladro_dd",
  "AverageDevAbs_TR",
  "TD",
  "OTm_total",
  "QST_parent_total",
  "QST_teacher_total"
)

present_measures <- unique(norms$measure)
missing_measures <- setdiff(needed_measures, present_measures)

# If TR total not in dataset, do not require its norms here
if (!has_TR_total) missing_measures <- setdiff(missing_measures, "AverageDevAbs_TR")

if (length(missing_measures)) {
  stop(
    "These measures are missing in norms_lookup (norms$measure): ",
    paste(missing_measures, collapse = ", "),
    call. = FALSE
  )
}

# -----------------------------
# Compute aligned z for the 6 indicators, row-wise
# -----------------------------
get_row_z <- function(i) {
  age <- dat$Age[i]
  
  z_te_b <- lookup_z_aligned(norms, "TE_Barca_dd", age, dat$TE_Barca_dd[i])
  z_te_l <- lookup_z_aligned(norms, "TE_Ladro_dd", age, dat$TE_Ladro_dd[i])
  z_te   <- if (is.finite(z_te_b) && is.finite(z_te_l)) mean(c(z_te_b, z_te_l)) else NA_real_
  
  z_tr <- if (has_TR_total) lookup_z_aligned(norms, "AverageDevAbs_TR", age, dat$AverageDevAbs_TR[i]) else NA_real_
  z_td <- lookup_z_aligned(norms, "TD", age, dat$RatioTD[i])
  
  z_child   <- lookup_z_aligned(norms, "OTm_total",         age, dat$OTm_total[i])
  z_parent  <- lookup_z_aligned(norms, "QST_parent_total",  age, dat$QSTp_Total_parent[i])
  z_teacher <- lookup_z_aligned(norms, "QST_teacher_total", age, dat$QSTp_Total_teacher[i])
  
  c(
    z_TE = z_te,
    z_TR = z_tr,
    z_TD = z_td,
    z_TMOQ_child_total   = z_child,
    z_TMOQ_parent_total  = z_parent,
    z_TMOQ_teacher_total = z_teacher
  )
}

Z <- t(vapply(seq_len(nrow(dat)), get_row_z, FUN.VALUE = numeric(6)))
Z <- as.data.frame(Z)

# Convert to SS (1-19)
SS <- Z %>% mutate(across(everything(), z_to_ss19))
names(SS) <- sub("^z_", "SS_", names(SS))

# complete cases for correlation estimation
ss_cols <- names(SS)
complete_idx <- complete.cases(SS[, ss_cols, drop = FALSE])
df_ss <- SS[complete_idx, ss_cols, drop = FALSE]

k <- ncol(df_ss)
if (k != 6) stop("Unexpected number of indicators: expected 6, got ", k, call. = FALSE)

if (nrow(df_ss) < 50) {
  warning("Few complete cases for correlation estimation (N = ", nrow(df_ss), ").", call. = FALSE)
}

# -----------------------------
# Build SoTQ conversion table (theoretical distribution of sumSS)
# -----------------------------
R <- suppressWarnings(cor(df_ss, use = "pairwise.complete.obs"))
if (any(!is.finite(R))) stop("Correlation matrix contains NA/Inf, cannot build SoTQ table.", call. = FALSE)
diag(R) <- 1

mu_sum <- k * 10
var_sum <- (3^2) * sum(R)  # includes diagonals
sd_sum  <- sqrt(var_sum)

if (!is.finite(sd_sum) || sd_sum <= 0) stop("Invalid sd_sum computed for SoTQ table.", call. = FALSE)

sum_min <- k * 1
sum_max <- k * 19
sumSS_vals <- sum_min:sum_max

z_sum   <- (sumSS_vals - mu_sum) / sd_sum
SoTQ_raw <- round(100 + 15 * z_sum)
SoTQ <- clamp(SoTQ_raw, SOTQ_CLAMP[1], SOTQ_CLAMP[2])

out <- tibble(
  sumSS = sumSS_vals,
  SoTQ  = SoTQ
)

write_csv(out, OUT_TABLE)

# -----------------------------
# Profile heterogeneity threshold (empirical)
# -----------------------------
deltaSS <- apply(df_ss, 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
thr95 <- as.numeric(quantile(deltaSS, probs = 0.95, na.rm = TRUE, names = FALSE))

spread_out <- tibble(
  definition = "deltaSS = max(SS_6_indicators) - min(SS_6_indicators)",
  percentile = 0.95,
  threshold = thr95,
  n_complete = nrow(df_ss),
  k_indicators = k
)

write_csv(spread_out, OUT_SPREAD)

# -----------------------------
# Console summary
# -----------------------------
message("Saved:")
message("  - ", OUT_TABLE)
message("  - ", OUT_SPREAD)
message("")
message("Sanity info:")
message("  k = ", k)
message("  N complete = ", nrow(df_ss))
message("  mu_sum = ", sprintf("%.2f", mu_sum))
message("  sd_sum = ", sprintf("%.2f", sd_sum))
message("  thr95(deltaSS) = ", sprintf("%.2f", thr95))
message("")
message("04_build_SoTQ_table.R completed successfully.")
