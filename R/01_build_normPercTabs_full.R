# R/build_normPercTabs_full.R
rm(list = ls())

library(readxl)
library(dplyr)
library(here)
library(gamlss)
library(gamlss.dist)
library(gamlss.add)

normalize = function(x, naIgnore=T){
  q = rep(NA,length(x))
  if(naIgnore==F & sum(is.na(x))>0) return(q)
  p = rank(x[!is.na(x)])/(length(x[!is.na(x)])+1)
  q[!is.na(x)] = qnorm(p)
  return(q)
}

# -------------------------
# Data
# -------------------------
df <- data.frame(read_excel(here("data/Database_Time_16_05_25_PD_MI.xlsx")))
dd <- data.frame(read_excel(here("data/Data_Dictionary.xlsx")))
dd$Var_name <- gsub("[ °-]", ".", dd$Var_name)

df <- df[!is.na(df$Class), ]  # consistent with your report pipeline
dat <- df

# -------------------------
# Grids
# -------------------------
ageBrack <- seq(from = 6, to = 12, by = 0.25)

# 19 probabilities matching standard-score cutoffs (1–19) based on the normal CDF grid
perc_list_19 <- c(
  0.00135, 0.00383, 0.00982,
  0.02275, 0.04779, 0.09121,
  0.15866, 0.25249, 0.36944,
  0.50,
  0.63056, 0.74751, 0.84134,
  0.90879, 0.95221, 0.97725,
  0.99018, 0.99617, 0.99865
)

# OPTIONAL: if you want the paper’s “5th–95th” available too, keep it separate
perc_list_5_95 <- seq(0.05, 0.95, by = 0.05)

new_norm_mat <- function(perc, ages) {
  m <- matrix(NA_real_, nrow = length(perc), ncol = length(ages))
  rownames(m) <- perc
  colnames(m) <- paste0("age-", ages)
  m
}

fill_quantiles <- function(fit, ages, perc, qfun) {
  m <- new_norm_mat(perc, ages)
  for (i in seq_along(ages)) {
    pa <- predictAll(fit, newdata = data.frame(Age = ages[i]), type = "response")
    m[, i] <- qfun(pa)
  }
  m
}

# -------------------------
# Child questionnaire OTm: total + subscales
# -------------------------
dat$OTm_total       <- rowSums(dat[, paste0("OTm_", 1:16)], na.rm = FALSE)
dat$OTm_orientation <- rowSums(dat[, paste0("OTm_", 1:10)], na.rm = FALSE)  # “Time Orientation/Knowledge” split in your docs
dat$OTm_management  <- rowSums(dat[, paste0("OTm_", 11:16)], na.rm = FALSE)

fit_OTm <- list()
norms_OTm_19 <- list()
norms_OTm_5_95 <- list()

otm_specs <- list(OTm_total = 16, OTm_orientation = 10, OTm_management = 6)

for (nm in names(otm_specs)) {
  bd <- otm_specs[[nm]]
  d <- dat %>% select(Age, all_of(nm)) %>% filter(!is.na(.data[[nm]]))
  
  fit_OTm[[nm]] <- gamlss(
    as.formula(paste0("cbind(", nm, ", ", bd, " - ", nm, ") ~ pbm(Age)")),
    sigma.fo = ~ pbm(Age),
    family = BB,
    data = d,
    trace = FALSE
  )
  
  norms_OTm_19[[nm]] <- fill_quantiles(
    fit_OTm[[nm]], ageBrack, perc_list_19,
    qfun = function(pa) qBB(p = perc_list_19, mu = pa$mu, sigma = pa$sigma, bd = bd)
  )
  
  norms_OTm_5_95[[nm]] <- fill_quantiles(
    fit_OTm[[nm]], ageBrack, perc_list_5_95,
    qfun = function(pa) qBB(p = perc_list_5_95, mu = pa$mu, sigma = pa$sigma, bd = bd)
  )
}

# -------------------------
# Parent / Teacher totals (use totals already in dataset)
# -------------------------
stopifnot("QSTp_Total_parent" %in% names(dat))
stopifnot("QSTp_Total_teacher" %in% names(dat))

fit_QSTp <- gamlss(
  cbind(QSTp_Total_parent, 30 - QSTp_Total_parent) ~ pbm(Age),
  sigma.fo = ~ pbm(Age),
  family = BB,
  data = dat[!is.na(dat$QSTp_Total_parent), c("Age", "QSTp_Total_parent")],
  trace = FALSE
)

fit_QSTt <- gamlss(
  cbind(QSTp_Total_teacher, 30 - QSTp_Total_teacher) ~ pbm(Age),
  sigma.fo = ~ pbm(Age),
  family = BB,
  data = dat[!is.na(dat$QSTp_Total_teacher), c("Age", "QSTp_Total_teacher")],
  trace = FALSE
)

norms_QSTp_19 <- fill_quantiles(
  fit_QSTp, ageBrack, perc_list_19,
  qfun = function(pa) qBB(p = perc_list_19, mu = pa$mu, sigma = pa$sigma, bd = 30)
)

norms_QSTt_19 <- fill_quantiles(
  fit_QSTt, ageBrack, perc_list_19,
  qfun = function(pa) qBB(p = perc_list_19, mu = pa$mu, sigma = pa$sigma, bd = 30)
)

# -------------------------
# Time Discrimination (TD): RatioTD, lower is better
# -------------------------
stopifnot("RatioTD" %in% names(dat))
dat_td <- dat %>% select(Age, RatioTD) %>% filter(!is.na(RatioTD))

# rescale ratio (~1..3) to 0..1 for modeling, then back-transform
td_min <- 1
td_max <- 3.004
dat_td$TD_prop <- (dat_td$RatioTD - td_min) / (td_max - td_min)

fit_TD <- gamlss(
  TD_prop ~ pb(Age),
  sigma.fo = ~ pb(Age),
  nu.fo = ~ 1,
  tau.fo = ~ 1,
  family = BCTo,
  data = dat_td,
  trace = FALSE
)

norms_TD_prop_19 <- fill_quantiles(
  fit_TD, ageBrack, perc_list_19,
  qfun = function(pa) qBCTo(p = perc_list_19, mu = pa$mu, sigma = pa$sigma, nu = pa$nu, tau = pa$tau)
)

norms_TD_19 <- td_min + norms_TD_prop_19 * (td_max - td_min)

# -------------------------
# Time Estimation (TE): two items
# -------------------------
stopifnot(all(c("TE_Barca_dd", "TE_Ladro_dd") %in% names(dat)))

fit_TE <- list()
norms_TE_19 <- list()

for (v in c("TE_Barca_dd", "TE_Ladro_dd")) {
  d <- dat %>% select(Age, all_of(v)) %>% filter(!is.na(.data[[v]]))
  
  fit_TE[[v]] <- gamlss(
    as.formula(paste0(v, " ~ pbm(-Age)")),
    sigma.fo = ~ pbm(Age),
    family = ZAGA,
    data = d,
    trace = FALSE
  )
  
  norms_TE_19[[v]] <- fill_quantiles(
    fit_TE[[v]], ageBrack, perc_list_19,
    qfun = function(pa) qZAGA(p = perc_list_19, mu = pa$mu, sigma = pa$sigma, nu = pa$nu)
  )
}

# optional aggregate TE table if TEdd_Mean exists (does not replace item scoring)
norms_TEdd_Mean_19 <- NULL
if ("TEdd_Mean" %in% names(dat)) {
  d <- dat %>% select(Age, TEdd_Mean) %>% filter(!is.na(TEdd_Mean))
  fit_TEdd <- gamlss(TEdd_Mean ~ pbm(-Age), sigma.fo = ~ pbm(Age), family = ZAGA, data = d, trace = FALSE)
  norms_TEdd_Mean_19 <- fill_quantiles(
    fit_TEdd, ageBrack, perc_list_19,
    qfun = function(pa) qZAGA(p = perc_list_19, mu = pa$mu, sigma = pa$sigma, nu = pa$nu)
  )
}

# -------------------------
# Time Reproduction (TR): item-level norms
# -------------------------
tr_items <- dd$Var_name[dd$ScaleName == "TR" & dd$Level == "item"]
tr_items <- tr_items[!is.na(tr_items) & grepl("^PercDevAbs_TR_", tr_items)]

stopifnot(length(tr_items) > 0)

dat_tr <- dat %>% select(Age, all_of(tr_items))
dat_tr <- dat_tr[complete.cases(dat_tr), ]

fit_TR_items <- list()
norms_TR_items_19 <- list()

for (v in tr_items) {
  fit_TR_items[[v]] <- gamlss(
    as.formula(paste0(v, " ~ pb(Age)")),
    sigma.fo = ~ pb(Age),
    nu.fo = ~ 1,
    tau.fo = ~ 1,
    family = BCTo,
    data = dat_tr,
    trace = FALSE
  )
  
  norms_TR_items_19[[v]] <- fill_quantiles(
    fit_TR_items[[v]], ageBrack, perc_list_19,
    qfun = function(pa) qBCTo(p = perc_list_19, mu = pa$mu, sigma = pa$sigma, nu = pa$nu, tau = pa$tau)
  )
}

# optional aggregate TR table if AverageDevAbs_TR exists (again, does not replace item scoring)
norms_AverageDevAbs_TR_19 <- NULL
if ("AverageDevAbs_TR" %in% names(dat)) {
  d <- dat %>% select(Age, AverageDevAbs_TR) %>% filter(!is.na(AverageDevAbs_TR))
  fit_TR_total <- gamlss(
    AverageDevAbs_TR ~ pb(Age),
    sigma.fo = ~ pb(Age),
    nu.fo = ~ 1,
    tau.fo = ~ 1,
    family = BCTo,
    data = d,
    trace = FALSE
  )
  norms_AverageDevAbs_TR_19 <- fill_quantiles(
    fit_TR_total, ageBrack, perc_list_19,
    qfun = function(pa) qBCTo(p = perc_list_19, mu = pa$mu, sigma = pa$sigma, nu = pa$nu, tau = pa$tau)
  )
}

# -------------------------
# Save outputs
# -------------------------

# (1) Paper-compatible object (minimal, stable naming)
tabList_core <- list(
  OTm  = norms_OTm_19$OTm_total,
  QSTp = norms_QSTp_19,
  QSTt = norms_QSTt_19,
  TD   = norms_TD_19,
  teB  = norms_TE_19$TE_Barca_dd,
  teL  = norms_TE_19$TE_Ladro_dd,
  TR   = norms_TR_items_19
)

saveRDS(tabList_core, here("scoring/normPercTabs.RDS"))

# (2) Expanded object (for Shiny + manual conversion needs)
direction <- c(
  OTm_total = "higher_better",
  OTm_orientation = "higher_better",
  OTm_management = "higher_better",
  QSTp_Total_parent = "higher_better",
  QSTp_Total_teacher = "higher_better",
  TD = "lower_better",
  TE_Barca_dd = "lower_better",
  TE_Ladro_dd = "lower_better",
  TEdd_Mean = "lower_better",
  TR_items = "lower_better",
  AverageDevAbs_TR = "lower_better"
)

tabList_full <- list(
  meta = list(
    age_grid = ageBrack,
    perc_list_19 = perc_list_19,
    perc_list_5_95 = perc_list_5_95,
    direction = direction
  ),
  # questionnaires
  OTm_total       = norms_OTm_19$OTm_total,
  OTm_orientation = norms_OTm_19$OTm_orientation,
  OTm_management  = norms_OTm_19$OTm_management,
  QST_parent_total  = norms_QSTp_19,
  QST_teacher_total = norms_QSTt_19,
  
  # tasks
  TD = norms_TD_19,
  TE_Barca_dd = norms_TE_19$TE_Barca_dd,
  TE_Ladro_dd = norms_TE_19$TE_Ladro_dd,
  TEdd_Mean = norms_TEdd_Mean_19,
  TR_items = norms_TR_items_19,
  AverageDevAbs_TR = norms_AverageDevAbs_TR_19,
  
  # optional: descriptive percentiles (paper sentence) for questionnaires if you want
  descriptive_5_95 = list(
    OTm_total = norms_OTm_5_95$OTm_total,
    OTm_orientation = norms_OTm_5_95$OTm_orientation,
    OTm_management = norms_OTm_5_95$OTm_management
  )
)

saveRDS(tabList_full, here("scoring/normPercTabs_full.RDS"))

message("Saved:\n  scoring/normPercTabs.RDS (core)\n  scoring/normPercTabs_full.RDS (expanded)")

#######################
