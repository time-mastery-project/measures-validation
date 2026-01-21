# R/05_estimate_SoTQ_CI_params.R
# ------------------------------------------------------------
# Estimate reliability for the observed SoTQ composite (6 indicators) and derive
# a 90% CI half-width on the SoTQ scale (M=100, SD=15).
#
# Uses table-based age norming via scoring/norms_lookup.csv (z_aligned lookup).
#
# Output:
#   scoring/SoTQ_CI_params.csv with:
#     rho_alpha, omega_total (optional), SEM, halfwidth90, N_used
# ------------------------------------------------------------

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(readr)
  library(dplyr)
  library(tibble)
  library(psych)
})

# -------------------------
# Paths
# -------------------------
DATA_XLSX    <- here("data", "Database_Time_16_05_25_PD_MI.xlsx")
NORMS_LOOKUP <- here("scoring", "norms_lookup.csv")
OUT_CSV      <- here("scoring", "SoTQ_CI_params.csv")

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

pick_col <- function(df, candidates) {
  n <- nrow(df)
  hit <- intersect(candidates, names(df))
  if (length(hit) < 1) return(list(vec = rep(NA_real_, n), name = NA_character_))
  v <- to_num(df[[hit[1]]])
  if (length(v) != n) v <- rep(NA_real_, n)
  list(vec = v, name = hit[1])
}

safe_rowmean2 <- function(a, b) {
  if (length(a) != length(b)) stop("Vectors have different lengths in safe_rowmean2().", call. = FALSE)
  out <- (a + b) / 2
  out[is.na(a) | is.na(b)] <- NA_real_
  out
}

lin_interp <- function(x, y, x0) {
  x <- to_num(x); y <- to_num(y); x0 <- to_num(x0)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (!is.finite(x0) || length(x) < 2) return(NA_real_)
  
  # exact hit
  if (any(x == x0)) return(y[which(x == x0)[1]])
  
  # cap
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

# -------------------------
# Load data + norms_lookup
# -------------------------
stopifnot(file.exists(DATA_XLSX))
stopifnot(file.exists(NORMS_LOOKUP))

dat_raw <- read_excel(DATA_XLSX)
dat <- as.data.frame(dat_raw)
names(dat) <- normalize_names(names(dat))

norms <- read_csv(NORMS_LOOKUP, show_col_types = FALSE)
names(norms) <- normalize_names(names(norms))

req_norms <- c("measure", "age", "raw", "z_aligned")
miss <- setdiff(req_norms, names(norms))
if (length(miss) > 0) stop("norms_lookup.csv missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)

norms$age <- to_num(norms$age)
norms$raw <- to_num(norms$raw)
norms$z_aligned <- to_num(norms$z_aligned)
norms <- norms %>% filter(is.finite(age), is.finite(raw), is.finite(z_aligned), !is.na(measure))

# -------------------------
# Lookup: nearest age, interpolate raw -> z_aligned (rule = 2 semantics)
# -------------------------
lookup_z_aligned <- function(measure, age, raw) {
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
  
  # collapse duplicate raw values (needed for stable interpolation)
  dfa2 <- dfa %>%
    group_by(.data$raw) %>%
    summarise(z_aligned = mean(.data$z_aligned), .groups = "drop") %>%
    arrange(.data$raw)
  
  if (nrow(dfa2) < 2) return(NA_real_)
  lin_interp(dfa2$raw, dfa2$z_aligned, raw)
}

lookup_vec <- function(measure, age_vec, raw_vec) {
  n <- length(age_vec)
  if (length(raw_vec) != n) stop("lookup_vec: age and raw vectors have different lengths.", call. = FALSE)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) out[i] <- lookup_z_aligned(measure, age_vec[i], raw_vec[i])
  out
}

# -------------------------
# Column mapping (robust)
# -------------------------
age_col <- pick_col(dat, c("age", "age_years", "eta", "eta_anni", "chronological_age", "age_in_years"))
age_years <- age_col$vec

te_boat_col  <- pick_col(dat, c("te_barca_dd", "te_barca", "barca_dd", "te_boat_dd", "boat_dd"))
te_thief_col <- pick_col(dat, c("te_ladro_dd", "te_ladro", "ladro_dd", "te_thief_dd", "thief_dd"))

tr_col <- pick_col(dat, c("averagedevabs_tr", "average_dev_abs_tr", "tr_averagedevabs", "tr_devabs", "tr"))
td_col <- pick_col(dat, c("ratiotd", "ratio_td", "td_ratio", "td", "ratio"))

child_col <- pick_col(dat, c("otm_total", "tmoq_c_total", "tmoq_children_total", "tmoq_c", "tmoqc_total"))
parent_col <- pick_col(dat, c("qstp_total_parent", "qst_parent_total", "tmoq_p_total", "tmoq_parent_total", "tmoq_p"))
teacher_col <- pick_col(dat, c("qstp_total_teacher", "qst_teacher_total", "tmoq_t_total", "tmoq_teacher_total", "tmoq_t"))

cat("Matched columns:\n")
cat("  Age:", age_col$name, "\n")
cat("  TE boat:", te_boat_col$name, "\n")
cat("  TE thief:", te_thief_col$name, "\n")
cat("  TR:", tr_col$name, "\n")
cat("  TD:", td_col$name, "\n")
cat("  Child total:", child_col$name, "\n")
cat("  Parent total:", parent_col$name, "\n")
cat("  Teacher total:", teacher_col$name, "\n\n")

if (all(is.na(age_years))) stop("Could not find an age column (in years). Check DATA_XLSX column names.", call. = FALSE)

# -------------------------
# Compute z_aligned indicators (6 indicators)
# -------------------------
z_te_b <- lookup_vec("TE_Barca_dd", age_years, te_boat_col$vec)
z_te_l <- lookup_vec("TE_Ladro_dd", age_years, te_thief_col$vec)
z_te   <- safe_rowmean2(z_te_b, z_te_l)

z_tr    <- lookup_vec("AverageDevAbs_TR",  age_years, tr_col$vec)
z_td    <- lookup_vec("TD",                age_years, td_col$vec)
z_child <- lookup_vec("OTm_total",         age_years, child_col$vec)
z_par   <- lookup_vec("QST_parent_total",  age_years, parent_col$vec)
z_teach <- lookup_vec("QST_teacher_total", age_years, teacher_col$vec)

Z <- tibble(
  TE      = z_te,
  TR      = z_tr,
  TD      = z_td,
  Child   = z_child,
  Parent  = z_par,
  Teacher = z_teach
) %>% filter(if_all(everything(), ~ !is.na(.x)))

N_used <- nrow(Z)
if (N_used < 30) {
  stop(
    paste0("Too few complete cases for SoTQ reliability (N=", N_used, "). ",
           "Likely missing columns or too many missing values."),
    call. = FALSE
  )
}

# -------------------------
# Reliability: alpha (primary), omega (optional)
# -------------------------
alpha_res <- psych::alpha(Z, check.keys = FALSE)
rho_alpha <- unname(alpha_res$total$raw_alpha)

omega_res <- tryCatch(psych::omega(Z, nfactors = 1, plot = FALSE), error = function(e) NULL)
omega_total <- if (!is.null(omega_res) && !is.null(omega_res$omega.tot)) unname(omega_res$omega.tot) else NA_real_

# -------------------------
# CI constants: 90% CI on SoTQ scale
# -------------------------
SD_SoTQ <- 15
SEM <- SD_SoTQ * sqrt(1 - rho_alpha)
zcrit90 <- 1.645
halfwidth90 <- zcrit90 * SEM

out <- tibble(
  data_file = basename(DATA_XLSX),
  N_used = N_used,
  rho_alpha = rho_alpha,
  omega_total = omega_total,
  SD_SoTQ = SD_SoTQ,
  SEM = SEM,
  zcrit90 = zcrit90,
  halfwidth90 = halfwidth90
)

write_csv(out, OUT_CSV)
print(out)

cat("\nSaved CI parameters to:\n  ", OUT_CSV, "\n", sep = "")

# ------------------------------------------------------------
# Again trransfer everything to shiny app folder
# ------------------------------------------------------------
dir.create("shiny/", recursive = TRUE, showWarnings = FALSE)
files <- list.files("scoring/", full.names = TRUE, recursive = FALSE)
files <- files[file.info(files)$isdir == FALSE]
file.copy(files, to = "shiny/", overwrite = TRUE)
