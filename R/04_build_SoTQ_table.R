# R/04_build_SoTQ_table.R
# ------------------------------------------------------------
# Build:
#  1) scoring/SoT_totalConversionTable.csv (sumSS -> SoTQ)
#  2) scoring/SoTQ_profileSpreadThreshold.csv (profile spread threshold)
#  3) scoring/TR_total_norm_params.csv (age-specific mu/sigma for TR total)
#
# TR total is defined as:
#  - compute z_aligned for PercDevAbs_TR_2..12 (each normed separately)
#  - m = mean(z_aligned_2..12) if all 11 are present
#  - z_TR_total = (m - mu_age) / sigma_age using exported TR_total_norm_params.csv
#
# Run from project root.
# ------------------------------------------------------------

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(readr)
  library(dplyr)
  library(tibble)
})

# -----------------------------
# Config
# -----------------------------
NORMS_LOOKUP <- here("scoring", "norms_lookup.csv")

# Pick the latest Database_Time_*.xlsx in data/
db_files <- list.files(here("data"), pattern = "^Database_Time_.*\\.xlsx$", full.names = TRUE)
if (!length(db_files)) stop("No Database_Time_*.xlsx found in data/.", call. = FALSE)
DATA_XLSX <- db_files[which.max(file.info(db_files)$mtime)]

OUT_TABLE    <- here("scoring", "SoT_totalConversionTable.csv")
OUT_SPREAD   <- here("scoring", "SoTQ_profileSpreadThreshold.csv")
OUT_TRPARAMS <- here("scoring", "TR_total_norm_params.csv")

SOTQ_CLAMP <- c(40, 160)

TR_SECONDS <- 2:12
TR_ITEMS   <- paste0("PercDevAbs_TR_", TR_SECONDS)

# If an age cell has small N, shrink mu/sigma toward global values
MIN_N_PER_AGE <- 5

# -----------------------------
# Helpers
# -----------------------------
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

normalize_names <- function(x) {
  x <- enc2utf8(x)
  x <- gsub("^\ufeff", "", x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x <- gsub("_+", "_", x)
  make.unique(x, sep = "_")
}

lin_interp <- function(x, y, x0) {
  x <- to_num(x); y <- to_num(y); x0 <- to_num(x0)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (!is.finite(x0) || length(x) < 2) return(NA_real_)
  
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
  age_value <- to_num(age_value)
  age_grid <- to_num(age_grid)
  age_grid <- age_grid[is.finite(age_grid)]
  if (!is.finite(age_value) || !length(age_grid)) return(NA_real_)
  age_grid[which.min(abs(age_grid - age_value))]
}

lookup_z_aligned <- function(norms, measure, age, raw) {
  age <- to_num(age); raw <- to_num(raw)
  if (!is.finite(age) || !is.finite(raw)) return(NA_real_)
  
  dfm <- norms %>% filter(.data$measure == !!measure)
  if (nrow(dfm) < 2) return(NA_real_)
  
  ages <- sort(unique(dfm$age))
  if (!length(ages)) return(NA_real_)
  
  age_cap <- clamp(age, min(ages), max(ages))
  a_near <- nearest_age(age_cap, ages)
  if (!is.finite(a_near)) return(NA_real_)
  
  dfa <- dfm %>%
    filter(.data$age == !!a_near) %>%
    transmute(raw = .data$raw, z_aligned = .data$z_aligned) %>%
    filter(is.finite(.data$raw), is.finite(.data$z_aligned))
  
  if (nrow(dfa) < 2) return(NA_real_)
  
  # collapse duplicates for stable interpolation
  dfa2 <- dfa %>%
    group_by(.data$raw) %>%
    summarise(z_aligned = mean(.data$z_aligned), .groups = "drop") %>%
    arrange(.data$raw)
  
  if (nrow(dfa2) < 2) return(NA_real_)
  lin_interp(dfa2$raw, dfa2$z_aligned, raw)
}

# vectorized, safe for dplyr::across
z_to_ss19 <- function(z) {
  z <- to_num(z)
  ss <- rep(NA_real_, length(z))
  ok <- is.finite(z)
  ss[ok] <- round(10 + 3 * z[ok])
  ss <- clamp(ss, 1, 19)
  as.integer(ss)
}

# -----------------------------
# Load inputs
# -----------------------------
if (!file.exists(DATA_XLSX)) stop("DATA_XLSX not found: ", DATA_XLSX, call. = FALSE)
if (!file.exists(NORMS_LOOKUP)) stop("Missing norms lookup: ", NORMS_LOOKUP, call. = FALSE)

dat_raw <- read_excel(DATA_XLSX)
dat <- as.data.frame(dat_raw)

norms <- read_csv(NORMS_LOOKUP, show_col_types = FALSE)
names(norms) <- normalize_names(names(norms))

req_norms <- c("measure", "age", "raw", "z_aligned")
miss_norms <- setdiff(req_norms, names(norms))
if (length(miss_norms) > 0) stop("norms_lookup.csv missing columns: ", paste(miss_norms, collapse = ", "), call. = FALSE)

norms <- norms %>%
  mutate(
    measure = as.character(.data$measure),
    age = to_num(.data$age),
    raw = to_num(.data$raw),
    z_aligned = to_num(.data$z_aligned)
  ) %>%
  filter(is.finite(.data$age), is.finite(.data$raw), is.finite(.data$z_aligned), !is.na(.data$measure))

age_grid <- sort(unique(norms$age))

# -----------------------------
# Required columns in dataset
# -----------------------------
if (!("Age" %in% names(dat))) stop("Dataset must contain column 'Age' (years).", call. = FALSE)
dat$Age <- to_num(dat$Age)

if (!("TE_Barca_dd" %in% names(dat))) stop("Dataset missing TE_Barca_dd", call. = FALSE)
if (!("TE_Ladro_dd" %in% names(dat))) stop("Dataset missing TE_Ladro_dd", call. = FALSE)
if (!("RatioTD" %in% names(dat))) stop("Dataset missing RatioTD", call. = FALSE)
if (!("QSTp_Total_parent" %in% names(dat))) stop("Dataset missing QSTp_Total_parent", call. = FALSE)
if (!("QSTp_Total_teacher" %in% names(dat))) stop("Dataset missing QSTp_Total_teacher", call. = FALSE)

dat$TE_Barca_dd <- to_num(dat$TE_Barca_dd)
dat$TE_Ladro_dd <- to_num(dat$TE_Ladro_dd)
dat$RatioTD <- to_num(dat$RatioTD)
dat$QSTp_Total_parent <- to_num(dat$QSTp_Total_parent)
dat$QSTp_Total_teacher <- to_num(dat$QSTp_Total_teacher)

# Child total: either OTm_total exists, or compute from OTm_1..OTm_16
if (all(paste0("OTm_", 1:16) %in% names(dat))) {
  dat$OTm_total <- rowSums(dat[, paste0("OTm_", 1:16)], na.rm = FALSE)
} else if ("OTm_total" %in% names(dat)) {
  dat$OTm_total <- to_num(dat$OTm_total)
} else {
  stop("Dataset missing OTm_total and OTm_1..OTm_16, cannot compute child total.", call. = FALSE)
}

# TR item columns
dat_norm_names <- dat
names(dat_norm_names) <- normalize_names(names(dat_norm_names))
tr_cols_norm <- normalize_names(TR_ITEMS)
missing_tr <- setdiff(tr_cols_norm, names(dat_norm_names))
if (length(missing_tr) > 0) {
  stop(
    paste0(
      "Dataset is missing TR item columns required for the new TR total. Missing: ",
      paste(missing_tr, collapse = ", "),
      ". Expected (original names): ",
      paste(TR_ITEMS, collapse = ", ")
    ),
    call. = FALSE
  )
}

# Check norms contain all measures needed
needed_measures <- c(
  "TE_Barca_dd", "TE_Ladro_dd",
  TR_ITEMS,
  "TD",
  "OTm_total",
  "QST_parent_total",
  "QST_teacher_total"
)
missing_measures <- setdiff(needed_measures, unique(norms$measure))
if (length(missing_measures)) {
  stop(
    "These measures are missing in norms_lookup (norms$measure): ",
    paste(missing_measures, collapse = ", "),
    call. = FALSE
  )
}

# -----------------------------
# Compute z-aligned indicators (row-wise)
# -----------------------------
n <- nrow(dat)

# TE total from 2 z's
z_te_b <- vapply(seq_len(n), function(i) lookup_z_aligned(norms, "TE_Barca_dd", dat$Age[i], dat$TE_Barca_dd[i]), numeric(1))
z_te_l <- vapply(seq_len(n), function(i) lookup_z_aligned(norms, "TE_Ladro_dd", dat$Age[i], dat$TE_Ladro_dd[i]), numeric(1))
z_te   <- ifelse(is.finite(z_te_b) & is.finite(z_te_l), (z_te_b + z_te_l) / 2, NA_real_)

# TD
z_td <- vapply(seq_len(n), function(i) lookup_z_aligned(norms, "TD", dat$Age[i], dat$RatioTD[i]), numeric(1))

# Questionnaires
z_child_total <- vapply(seq_len(n), function(i) lookup_z_aligned(norms, "OTm_total", dat$Age[i], dat$OTm_total[i]), numeric(1))
z_parent      <- vapply(seq_len(n), function(i) lookup_z_aligned(norms, "QST_parent_total", dat$Age[i], dat$QSTp_Total_parent[i]), numeric(1))
z_teacher     <- vapply(seq_len(n), function(i) lookup_z_aligned(norms, "QST_teacher_total", dat$Age[i], dat$QSTp_Total_teacher[i]), numeric(1))

# TR items: z_aligned per duration, then m = mean(z) if complete
tr_raw_mat <- as.data.frame(dat_norm_names[, tr_cols_norm, drop = FALSE])
for (j in seq_along(tr_raw_mat)) tr_raw_mat[[j]] <- to_num(tr_raw_mat[[j]])

z_tr_items <- matrix(NA_real_, nrow = n, ncol = length(TR_ITEMS))
colnames(z_tr_items) <- TR_ITEMS

for (j in seq_along(TR_ITEMS)) {
  meas <- TR_ITEMS[j]
  rawv <- tr_raw_mat[[j]]
  z_tr_items[, j] <- vapply(seq_len(n), function(i) lookup_z_aligned(norms, meas, dat$Age[i], rawv[i]), numeric(1))
}

tr_complete <- apply(z_tr_items, 1, function(x) all(is.finite(x)))
m_tr <- rep(NA_real_, n)
m_tr[tr_complete] <- rowMeans(z_tr_items[tr_complete, , drop = FALSE])

# -----------------------------
# Build TR total norm parameters (mu/sigma of m_tr by age cell)
# -----------------------------
age_cell <- vapply(dat$Age, function(a) nearest_age(a, age_grid), numeric(1))
age_cell <- clamp(age_cell, min(age_grid), max(age_grid))

m_ok <- is.finite(m_tr) & is.finite(age_cell)
m_global <- m_tr[m_ok]

mu_global <- mean(m_global)
sigma_global <- sd(m_global)
n_global <- length(m_global)

if (!is.finite(sigma_global) || sigma_global <= 0) sigma_global <- 1

by_age <- tibble(age = age_cell[m_ok], m = m_tr[m_ok]) %>%
  group_by(.data$age) %>%
  summarise(
    n_complete = dplyr::n(),
    mu_raw = mean(.data$m),
    sigma_raw = sd(.data$m),
    .groups = "drop"
  )

params <- tibble(age = age_grid) %>%
  left_join(by_age, by = "age") %>%
  mutate(
    n_complete = ifelse(is.na(.data$n_complete), 0L, as.integer(.data$n_complete)),
    mu_raw = ifelse(is.finite(.data$mu_raw), .data$mu_raw, NA_real_),
    sigma_raw = ifelse(is.finite(.data$sigma_raw) & .data$sigma_raw > 0, .data$sigma_raw, NA_real_),
    w = pmin(1, .data$n_complete / MIN_N_PER_AGE),
    mu = ifelse(is.finite(.data$mu_raw), .data$w * .data$mu_raw + (1 - .data$w) * mu_global, mu_global),
    sigma = ifelse(is.finite(.data$sigma_raw), .data$w * .data$sigma_raw + (1 - .data$w) * sigma_global, sigma_global),
    sigma = ifelse(!is.finite(.data$sigma) | .data$sigma <= 0, sigma_global, .data$sigma),
    mu_global = mu_global,
    sigma_global = sigma_global,
    n_global = n_global
  ) %>%
  select(age, n_complete, mu, sigma, mu_global, sigma_global, n_global)

write_csv(params, OUT_TRPARAMS)

# z_TR_total: re-standardize m_tr using params (nearest age cell)
mu_by_cell <- setNames(params$mu, params$age)
sd_by_cell <- setNames(params$sigma, params$age)

z_tr_total <- rep(NA_real_, n)
for (i in seq_len(n)) {
  if (!is.finite(m_tr[i]) || !is.finite(age_cell[i])) next
  a <- as.character(age_cell[i])
  mu <- mu_by_cell[[a]]
  sd <- sd_by_cell[[a]]
  if (!is.finite(mu) || !is.finite(sd) || sd <= 0) next
  z_tr_total[i] <- (m_tr[i] - mu) / sd
}

# -----------------------------
# Assemble the 6 SoTQ indicators (aligned z)
# -----------------------------
Z <- tibble(
  z_TE = z_te,
  z_TR = z_tr_total,
  z_TD = z_td,
  z_Child = z_child_total,
  z_Parent = z_parent,
  z_Teacher = z_teacher
)

# Convert to scaled scores (SS 1..19)
SS <- Z %>% mutate(across(everything(), z_to_ss19))

# complete cases for correlation estimation
complete_idx <- complete.cases(SS)
df_ss <- SS[complete_idx, , drop = FALSE]

k <- ncol(df_ss)
if (k != 6) stop("Unexpected number of indicators: expected 6, got ", k, call. = FALSE)

if (nrow(df_ss) < 50) warning("Few complete cases for correlation estimation (N = ", nrow(df_ss), ").", call. = FALSE)

# -----------------------------
# Build SoTQ conversion table (theoretical distribution of sumSS)
# -----------------------------
R <- suppressWarnings(cor(df_ss, use = "pairwise.complete.obs"))
if (any(!is.finite(R))) stop("Correlation matrix contains NA/Inf, cannot build SoTQ table.", call. = FALSE)
diag(R) <- 1

mu_sum <- k * 10
var_sum <- (3^2) * sum(R)  # includes diagonals
sd_sum <- sqrt(var_sum)

if (!is.finite(sd_sum) || sd_sum <= 0) stop("Invalid sd_sum computed for SoTQ table.", call. = FALSE)

sum_min <- k * 1
sum_max <- k * 19
sumSS_vals <- sum_min:sum_max

z_sum <- (sumSS_vals - mu_sum) / sd_sum
SoTQ_raw <- round(100 + 15 * z_sum)
SoTQ <- clamp(SoTQ_raw, SOTQ_CLAMP[1], SOTQ_CLAMP[2])

out <- tibble(sumSS = sumSS_vals, SoTQ = SoTQ)
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

message("Saved:")
message("  - ", OUT_TRPARAMS)
message("  - ", OUT_TABLE)
message("  - ", OUT_SPREAD)
message("Using DATA_XLSX: ", DATA_XLSX)
message("k = ", k, ", N complete = ", nrow(df_ss))
message("04_build_SoTQ_table.R completed successfully.")
