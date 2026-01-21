# R/05_estimate_SoTQ_CI_params.R
# ------------------------------------------------------------
# Estimate reliability for the SoTQ composite (6 indicators) and derive
# a 90% CI half-width on the SoTQ scale (M=100, SD=15).
#
# Decisions implemented (per our discussion):
#  1) Compute reliability on aligned z indicators (continuous), not on SS 1–19.
#  2) Avoid low N by using pairwise-complete correlations (no full-row complete-cases requirement).
#  3) Prefer omega_total, fall back to alpha, clamp rho_used to [0, 1] for SEM.
#  4) Keep output compatible with the Shiny app, and add extra transparency columns.
#
# TR total logic (preferred):
#  - standardize each TR duration (2..12 s) separately via norms_lookup.csv
#  - mean of the 11 z_aligned values (per subject), only if all 11 present
#  - re-standardize that mean by age using scoring/TR_total_norm_params.csv
#
# Output:
#   scoring/SoTQ_CI_params.csv with (at least):
#     data_file, N_used, rho_alpha, omega_total, SD_SoTQ, SEM, zcrit90, halfwidth90
#   plus additional diagnostics columns.
# ------------------------------------------------------------

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(readr)
  library(dplyr)
  library(tibble)
  library(psych)
  library(stringr)
})

# -------------------------
# Config
# -------------------------

# If you want to hardcode a specific file, set DATA_XLSX explicitly.
# Otherwise, the script will automatically pick the most recent Database_Time_*.xlsx in data/.
DATA_XLSX <- NA_character_

DATA_DIR     <- here("data")
NORMS_LOOKUP  <- here("scoring", "norms_lookup.csv")
TR_PARAMS     <- here("scoring", "TR_total_norm_params.csv")
OUT_CSV       <- here("scoring", "SoTQ_CI_params.csv")

TR_SECONDS <- 2:12
TR_ITEMS   <- paste0("PercDevAbs_TR_", TR_SECONDS)

ZCRIT_90 <- qnorm(0.95)  # 90% CI two-sided
SD_SoTQ  <- 15

MIN_N_FOR_RELIABILITY <- 20  # used only for warnings, pairwise can work with less but becomes unstable

# -------------------------
# Helpers
# -------------------------

normalize_names <- function(x) {
  x <- enc2utf8(x)
  x <- gsub("^\ufeff", "", x)              # BOM
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x <- gsub("_+", "_", x)
  make.unique(x, sep = "_")
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

infer_age_unit <- function(age_vec) {
  age_vec <- to_num(age_vec)
  age_vec <- age_vec[is.finite(age_vec)]
  if (!length(age_vec)) return("years")
  if (max(age_vec) <= 20) "years" else "months"
}

# Vector-safe: years and months can be vectors (months can also be scalar)
age_to_grid_unit <- function(years, months = 0, unit = c("years", "months")) {
  unit <- match.arg(unit)
  years <- to_num(years)
  months <- to_num(months)
  
  if (length(months) == 1 && length(years) > 1) months <- rep(months, length(years))
  if (length(months) != length(years) && length(years) != 1) {
    stop("age_to_grid_unit: 'months' must be length 1 or same length as 'years'.", call. = FALSE)
  }
  
  out <- rep(NA_real_, length(years))
  ok <- is.finite(years)
  if (!any(ok)) return(out)
  
  m2 <- months
  m2[!is.finite(m2)] <- 0
  
  if (unit == "years") {
    out[ok] <- years[ok] + m2[ok] / 12
  } else {
    out[ok] <- years[ok] * 12 + m2[ok]
  }
  out
}

nearest_age <- function(age_value, age_grid) {
  age_value <- to_num(age_value)
  age_grid <- to_num(age_grid)
  age_grid <- age_grid[is.finite(age_grid)]
  if (!is.finite(age_value) || !length(age_grid)) return(NA_real_)
  age_grid[which.min(abs(age_grid - age_value))]
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

pick_existing <- function(nm_vec, candidates) {
  hit <- intersect(candidates, nm_vec)
  if (length(hit) >= 1) hit[1] else NA_character_
}

find_tr_colname <- function(nm_vec, tsec) {
  candidates <- c(
    paste0("percdevabs_tr_", tsec),
    paste0("percdevabs_tr", tsec),
    paste0("tr_percdevabs_", tsec),
    paste0("tr_percdevabs", tsec)
  )
  hit <- pick_existing(nm_vec, candidates)
  if (!is.na(hit)) return(hit)
  
  rx <- paste0("^percdevabs_tr_?0*", tsec, "$")
  m <- nm_vec[str_detect(nm_vec, rx)]
  if (length(m) >= 1) return(m[1])
  
  NA_character_
}

# Pairwise N matrix for a numeric matrix
pairwise_n <- function(M) {
  k <- ncol(M)
  out <- matrix(NA_integer_, k, k)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      out[i, j] <- sum(is.finite(M[, i]) & is.finite(M[, j]))
    }
  }
  colnames(out) <- colnames(M)
  rownames(out) <- colnames(M)
  out
}

alpha_from_cov <- function(S) {
  # General Cronbach alpha from covariance/correlation matrix S:
  # alpha = k/(k-1) * (1 - trace(S)/sum(S))
  k <- ncol(S)
  if (k < 2) return(NA_real_)
  if (!all(is.finite(S))) return(NA_real_)
  k/(k - 1) * (1 - sum(diag(S)) / sum(S))
}

# -------------------------
# Resolve input dataset
# -------------------------

if (is.na(DATA_XLSX) || !nzchar(DATA_XLSX)) {
  cand <- list.files(DATA_DIR, pattern = "^Database_Time_.*\\.xlsx$", full.names = TRUE)
  if (length(cand) < 1) stop("No Database_Time_*.xlsx found in data/.", call. = FALSE)
  info <- file.info(cand)
  DATA_XLSX <- rownames(info)[which.max(info$mtime)]
}

if (!file.exists(DATA_XLSX)) stop("DATA_XLSX not found: ", DATA_XLSX, call. = FALSE)
if (!file.exists(NORMS_LOOKUP)) stop("Missing: ", NORMS_LOOKUP, call. = FALSE)

cat("Using dataset:\n  ", DATA_XLSX, "\n\n", sep = "")

# -------------------------
# Load data and tables
# -------------------------

dat_raw <- read_excel(DATA_XLSX)
dat <- as.data.frame(dat_raw)
names(dat) <- normalize_names(names(dat))

norms <- read_csv(NORMS_LOOKUP, show_col_types = FALSE)
names(norms) <- normalize_names(names(norms))

req_norms <- c("measure", "age", "raw", "z_aligned")
miss_norms <- setdiff(req_norms, names(norms))
if (length(miss_norms) > 0) {
  stop("norms_lookup.csv missing required columns: ", paste(miss_norms, collapse = ", "), call. = FALSE)
}

norms <- norms %>%
  mutate(
    measure    = as.character(.data$measure),
    age        = to_num(.data$age),
    raw        = to_num(.data$raw),
    z_aligned  = to_num(.data$z_aligned)
  ) %>%
  filter(!is.na(.data$measure), is.finite(.data$age), is.finite(.data$raw), is.finite(.data$z_aligned))

if (nrow(norms) < 10) stop("norms_lookup.csv looks empty after filtering.", call. = FALSE)

GRID_AGE_UNIT     <- infer_age_unit(norms$age)
NORMS_BY_MEASURE  <- split(norms, norms$measure)
MEASURES_AVAILABLE <- sort(unique(norms$measure))

# TR mode: prefer per-duration pipeline, else stop (keeps consistency with app scoring)
TR_MODE <- "per_duration"
if (!all(TR_ITEMS %in% MEASURES_AVAILABLE)) {
  missing_meas <- TR_ITEMS[!TR_ITEMS %in% MEASURES_AVAILABLE]
  stop(
    "norms_lookup.csv does not contain all TR per-duration measures (PercDevAbs_TR_2..12).\n",
    "Missing: ", paste(missing_meas, collapse = ", "), "\n",
    "To proceed, regenerate norms_lookup.csv to include these measures (recommended).",
    call. = FALSE
  )
}

# TR params needed for per-duration TR total
if (!file.exists(TR_PARAMS)) stop("Missing: ", TR_PARAMS, call. = FALSE)

trp <- read_csv(TR_PARAMS, show_col_types = FALSE)
names(trp) <- normalize_names(names(trp))
if (!all(c("age", "mu", "sigma") %in% names(trp))) {
  stop("TR_total_norm_params.csv must contain columns: age, mu, sigma.", call. = FALSE)
}

trp <- trp %>%
  mutate(age = to_num(.data$age), mu = to_num(.data$mu), sigma = to_num(.data$sigma)) %>%
  filter(is.finite(.data$age), is.finite(.data$mu), is.finite(.data$sigma), .data$sigma > 0) %>%
  arrange(.data$age)

if (nrow(trp) < 2) stop("TR_total_norm_params.csv looks invalid after filtering.", call. = FALSE)

TRP_AGE_UNIT <- infer_age_unit(trp$age)

# -------------------------
# Lookup z_aligned (raw interpolation within age, then age interpolation)
# -------------------------

lookup_z_from_grid <- function(measure, age_grid_value, raw_value) {
  dfm <- NORMS_BY_MEASURE[[measure]]
  if (is.null(dfm)) return(NA_real_)
  
  ages <- sort(unique(dfm$age))
  if (!length(ages)) return(NA_real_)
  
  age_in <- to_num(age_grid_value)
  raw_in <- to_num(raw_value)
  if (!is.finite(age_in) || !is.finite(raw_in)) return(NA_real_)
  
  age_cap <- clamp(age_in, min(ages), max(ages))
  a1 <- max(ages[ages <= age_cap])
  a2 <- min(ages[ages >= age_cap])
  
  z_at_age <- function(a) {
    dfa <- dfm %>%
      filter(.data$age == a) %>%
      select(.data$raw, .data$z_aligned) %>%
      filter(is.finite(.data$raw), is.finite(.data$z_aligned))
    
    if (nrow(dfa) < 2) return(NA_real_)
    
    dfa2 <- dfa %>%
      group_by(.data$raw) %>%
      summarise(z_aligned = mean(.data$z_aligned), .groups = "drop") %>%
      arrange(.data$raw)
    
    if (nrow(dfa2) < 2) return(NA_real_)
    lin_interp(dfa2$raw, dfa2$z_aligned, raw_in)
  }
  
  z1 <- z_at_age(a1)
  z2 <- z_at_age(a2)
  
  if (!is.finite(z1) && !is.finite(z2)) return(NA_real_)
  if (isTRUE(all.equal(a1, a2))) return(z1)
  
  w <- (age_cap - a1) / (a2 - a1)
  z1 + w * (z2 - z1)
}

lookup_vec <- function(measure, age_vec, raw_vec) {
  n <- length(age_vec)
  if (length(raw_vec) != n) stop("lookup_vec: age and raw vectors have different lengths.", call. = FALSE)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) out[i] <- lookup_z_from_grid(measure, age_vec[i], raw_vec[i])
  out
}

get_mu_sigma_for_age <- function(trp_df, age_grid_value) {
  age_in <- to_num(age_grid_value)
  if (!is.finite(age_in)) return(list(mu = NA_real_, sigma = NA_real_))
  ages <- trp_df$age
  age_cap <- clamp(age_in, min(ages), max(ages))
  a0 <- nearest_age(age_cap, ages)
  row <- trp_df[which.min(abs(ages - a0)), , drop = FALSE]
  list(mu = row$mu[1], sigma = row$sigma[1])
}

# -------------------------
# Column mapping (robust)
# -------------------------

age_name <- pick_existing(names(dat), c("age", "age_years", "eta", "eta_anni", "chronological_age", "age_in_years"))
if (is.na(age_name)) stop("Could not find an age column in the dataset.", call. = FALSE)
age_years <- to_num(dat[[age_name]])

# TE
te_boat_name  <- pick_existing(names(dat), c("te_barca_dd", "te_boat_dd", "boat_dd", "te_barca"))
te_thief_name <- pick_existing(names(dat), c("te_ladro_dd", "te_thief_dd", "thief_dd", "te_ladro"))
if (is.na(te_boat_name) || is.na(te_thief_name)) {
  stop("Missing TE columns (expected TE_Barca_dd and TE_Ladro_dd raw deviations) in the dataset.", call. = FALSE)
}

# TD
td_name <- pick_existing(names(dat), c("ratiotd", "ratio_td", "td_ratio", "ratio", "td"))
if (is.na(td_name)) stop("Missing TD ratio column in the dataset.", call. = FALSE)

# Questionnaires
child_total_name   <- pick_existing(names(dat), c("otm_total", "otm_tot", "tmoq_c_total"))
parent_total_name  <- pick_existing(names(dat), c("qst_parent_total", "qstp_total_parent", "tmoq_p_total"))
teacher_total_name <- pick_existing(names(dat), c("qst_teacher_total", "qstp_total_teacher", "tmoq_t_total"))
if (is.na(child_total_name) || is.na(parent_total_name) || is.na(teacher_total_name)) {
  stop("Missing questionnaire total columns (child, parent, teacher) in the dataset.", call. = FALSE)
}

# TR items (2..12)
tr_item_cols <- vapply(TR_SECONDS, function(t) find_tr_colname(names(dat), t), character(1))
if (any(is.na(tr_item_cols))) {
  miss_t <- TR_SECONDS[is.na(tr_item_cols)]
  stop(
    "Missing TR per-duration columns in the dataset for targets: ",
    paste(miss_t, collapse = ", "),
    ". Expected normalized names like percdevabs_tr_2 .. percdevabs_tr_12.",
    call. = FALSE
  )
}

cat("Matched columns:\n")
cat("  Age:", age_name, "\n")
cat("  TE boat:", te_boat_name, "\n")
cat("  TE thief:", te_thief_name, "\n")
cat("  TD:", td_name, "\n")
cat("  Child total:", child_total_name, "\n")
cat("  Parent total:", parent_total_name, "\n")
cat("  Teacher total:", teacher_total_name, "\n")
cat("  TR items:\n")
for (i in seq_along(TR_SECONDS)) cat("    ", TR_SECONDS[i], "s: ", tr_item_cols[i], "\n", sep = "")
cat("\n")

# -------------------------
# Compute age in the norms grid units
# -------------------------

age_grid_norms <- if (GRID_AGE_UNIT == "years") {
  age_to_grid_unit(age_years, 0, "years")
} else {
  age_to_grid_unit(age_years, 0, "months")
}

age_grid_trp <- if (TRP_AGE_UNIT == "years") {
  age_to_grid_unit(age_years, 0, "years")
} else {
  age_to_grid_unit(age_years, 0, "months")
}

# -------------------------
# Compute aligned z indicators (6 indicators)
# -------------------------

# TE composite (mean of the two aligned z values, only if both present)
z_te_b <- lookup_vec("TE_Barca_dd", age_grid_norms, to_num(dat[[te_boat_name]]))
z_te_l <- lookup_vec("TE_Ladro_dd", age_grid_norms, to_num(dat[[te_thief_name]]))
z_te   <- ifelse(is.finite(z_te_b) & is.finite(z_te_l), (z_te_b + z_te_l) / 2, NA_real_)

# TD
z_td <- lookup_vec("TD", age_grid_norms, to_num(dat[[td_name]]))

# Questionnaires
z_child   <- lookup_vec("OTm_total",          age_grid_norms, to_num(dat[[child_total_name]]))
z_parent  <- lookup_vec("QST_parent_total",   age_grid_norms, to_num(dat[[parent_total_name]]))
z_teacher <- lookup_vec("QST_teacher_total",  age_grid_norms, to_num(dat[[teacher_total_name]]))

# TR items (2..12), mean, then re-standardize using TR_total_norm_params.csv
z_tr_items <- sapply(seq_along(TR_SECONDS), function(j) {
  tsec <- TR_SECONDS[j]
  raw_vec <- to_num(dat[[tr_item_cols[j]]])
  lookup_vec(paste0("PercDevAbs_TR_", tsec), age_grid_norms, raw_vec)
})
colnames(z_tr_items) <- paste0("z_tr_", TR_SECONDS)

z_tr_mean <- apply(z_tr_items, 1, function(row) {
  if (any(!is.finite(row))) return(NA_real_)
  mean(row)
})

pars <- lapply(age_grid_trp, function(a) get_mu_sigma_for_age(trp, a))
mu_by_row  <- vapply(pars, `[[`, numeric(1), "mu")
sig_by_row <- vapply(pars, `[[`, numeric(1), "sigma")

z_tr_total <- ifelse(
  is.finite(z_tr_mean) & is.finite(mu_by_row) & is.finite(sig_by_row) & sig_by_row > 0,
  (z_tr_mean - mu_by_row) / sig_by_row,
  NA_real_
)

Z <- tibble(
  TE      = z_te,
  TR      = z_tr_total,
  TD      = z_td,
  Child   = z_child,
  Parent  = z_parent,
  Teacher = z_teacher
)

# -------------------------
# Diagnostics: availability
# -------------------------

n_total <- nrow(Z)
n_by_indicator <- sapply(Z, function(x) sum(is.finite(x)))
n_complete6 <- sum(complete.cases(Z))

cat("Data availability:\n")
cat("  N total rows:", n_total, "\n")
cat("  N complete (all 6 indicators):", n_complete6, "\n")
cat("  N finite by indicator:\n")
print(n_by_indicator)

# -------------------------
# Reliability and CI half-width (on aligned z, using pairwise correlations)
# -------------------------

Z_mat <- as.matrix(Z)

# Pairwise N matrix and pairwise correlation matrix
N_pair <- pairwise_n(Z_mat)
Rz <- suppressWarnings(cor(Z_mat, use = "pairwise.complete.obs"))

# Basic sanity
diag(Rz) <- 1
Rz <- (Rz + t(Rz)) / 2

# Effective N for reporting: median and min pairwise N (off-diagonal)
off <- lower.tri(N_pair, diag = FALSE)
N_pair_min <- suppressWarnings(min(N_pair[off], na.rm = TRUE))
N_pair_med <- suppressWarnings(as.integer(round(median(N_pair[off], na.rm = TRUE))))

cat("\nPairwise N (off-diagonal):\n")
cat("  min =", N_pair_min, "\n")
cat("  median =", N_pair_med, "\n")

if (!is.finite(N_pair_min) || N_pair_min < MIN_N_FOR_RELIABILITY) {
  warning("Low effective pairwise N for at least one indicator pair, reliability may be unstable.")
}

# If R has NA entries (no overlap for some pair), fall back to complete cases correlation
R_source <- "pairwise"
if (any(!is.finite(Rz))) {
  cat("\nSome pairwise correlations are NA/Inf, falling back to complete-case correlations.\n")
  Z_cc <- Z[complete.cases(Z), , drop = FALSE]
  if (nrow(Z_cc) < 3) {
    warning("Too few complete cases to estimate reliability, setting reliability to 0.")
    Rz <- diag(ncol(Z_mat))
    colnames(Rz) <- colnames(Z_mat)
    rownames(Rz) <- colnames(Z_mat)
    R_source <- "none"
  } else {
    Rz <- suppressWarnings(cor(as.matrix(Z_cc)))
    diag(Rz) <- 1
    Rz <- (Rz + t(Rz)) / 2
    R_source <- "complete_case"
  }
}

cat("\nCorrelation matrix used for reliability (rounded):\n")
print(round(Rz, 3))

# Smooth to nearest PSD correlation matrix if needed (omega can fail on non-PD matrices)
Rz_s <- tryCatch(psych::cor.smooth(Rz), error = function(e) Rz)

# Alpha from correlation matrix (as covariance matrix since diag=1)
rho_alpha <- alpha_from_cov(Rz_s)

# Omega total (1-factor) from correlation matrix
omega_total <- tryCatch({
  om <- suppressWarnings(psych::omega(Rz_s, nfactors = 1, n.obs = max(3, N_pair_med), plot = FALSE))
  as.numeric(om$omega.tot)
}, error = function(e) NA_real_)

cat("\nReliability diagnostics (aligned z, ", R_source, " correlations):\n", sep = "")
cat("  alpha =", ifelse(is.finite(rho_alpha), round(rho_alpha, 3), NA), "\n")
cat("  omega_total =", ifelse(is.finite(omega_total), round(omega_total, 3), NA), "\n")

# Prefer omega_total, else alpha, else 0. Clamp for SEM.
rho_used <- omega_total
rho_source <- "omega_total"
if (!is.finite(rho_used)) {
  rho_used <- rho_alpha
  rho_source <- "alpha"
}
if (!is.finite(rho_used)) {
  rho_used <- 0
  rho_source <- "none"
}
rho_used <- max(0, min(1, rho_used))

SEM <- SD_SoTQ * sqrt(1 - rho_used)
halfwidth90 <- ZCRIT_90 * SEM

# -------------------------
# Save output
# -------------------------

out <- tibble(
  data_file = basename(DATA_XLSX),
  
  # Backward-compatible field for the app:
  # with pairwise reliability there is no single N, we store the median pairwise N as "N_used"
  N_used = N_pair_med,
  
  # Original columns expected by the app:
  rho_alpha = rho_alpha,
  omega_total = omega_total,
  SD_SoTQ = SD_SoTQ,
  SEM = SEM,
  zcrit90 = ZCRIT_90,
  halfwidth90 = halfwidth90,
  
  # Extra transparency columns (harmless for the app):
  rho_used = rho_used,
  rho_source = rho_source,
  R_source = R_source,
  N_total = n_total,
  N_complete6 = n_complete6,
  N_pairwise_min = N_pair_min,
  N_pairwise_median = N_pair_med
)

write_csv(out, OUT_CSV)

cat("\nSaved CI parameters to:\n  ", OUT_CSV, "\n\n", sep = "")
print(out)
cat("\n05_estimate_SoTQ_CI_params.R completed successfully.\n")
