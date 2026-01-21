# app.R
# Sense of Time Battery, scoring and standardization (table-based, no GAMLSS)
#
# Tables used:
# - scoring/norms_lookup.csv  (main norms lookup; columns: measure, age, raw, z_aligned, ...)
# - scoring/SoT_totalConversionTable.csv
# - scoring/SoTQ_CI_params.csv (optional)
# - scoring/SoTQ_profileSpreadThreshold.csv (optional)
#
# Plot tuning:
# - Change PLOT_TEXT_SIZE to quickly tune all plot labels.

suppressPackageStartupMessages({
  library(shiny)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(tibble)
})

# -------------------------
# Global tuning parameters
# -------------------------

PLOT_TEXT_SIZE <- 18  # single knob to tune plot label sizes

# -------------------------
# Helpers
# -------------------------

normalize_names <- function(x) {
  x <- enc2utf8(x)
  x <- gsub("^\ufeff", "", x) # BOM
  x <- trimws(x)
  x[is.na(x) | x == ""] <- paste0("col_", seq_along(x))[is.na(x) | x == ""]
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

# Scaled score on the 1-19 scale (mean 10, SD 3), derived from aligned z
scaled_from_z <- function(z_aligned) {
  ss <- round(10 + 3 * z_aligned)
  clamp(ss, 1, 19)
}

find_file <- function(relpath) {
  candidates <- unique(c(
    relpath,
    file.path(getwd(), relpath),
    file.path("..", relpath),
    file.path("/mnt/data", basename(relpath))
  ))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit) || length(hit) == 0) return(NA_character_)
  hit
}

infer_age_unit <- function(age_vec) {
  age_vec <- suppressWarnings(as.numeric(age_vec))
  age_vec <- age_vec[is.finite(age_vec)]
  if (!length(age_vec)) return("years")
  if (max(age_vec) <= 20) "years" else "months"
}

age_to_grid_unit <- function(years, months, unit = c("years", "months")) {
  unit <- match.arg(unit)
  years <- suppressWarnings(as.numeric(years))
  months <- suppressWarnings(as.numeric(months))
  if (!is.finite(years) || !is.finite(months)) return(NA_real_)
  if (unit == "years") return(years + months / 12)
  years * 12 + months
}

grid_age_to_ym <- function(age_grid, unit = c("years", "months")) {
  unit <- match.arg(unit)
  age_grid <- suppressWarnings(as.numeric(age_grid))
  if (!is.finite(age_grid)) return(list(years = NA_integer_, months = NA_integer_))
  if (unit == "years") {
    y <- floor(age_grid)
    m <- round((age_grid - y) * 12)
    if (m == 12) { y <- y + 1; m <- 0 }
    return(list(years = as.integer(y), months = as.integer(m)))
  } else {
    total_m <- round(age_grid)
    y <- floor(total_m / 12)
    m <- total_m - 12 * y
    return(list(years = as.integer(y), months = as.integer(m)))
  }
}

is_whole_number <- function(x, tol = 1e-9) {
  x <- to_num(x)
  if (!is.finite(x)) return(FALSE)
  abs(x - round(x)) < tol
}

# -------------------------
# Measures mapping (names must match the norms table)
# -------------------------

MEASURES <- list(
  TE_BOAT = "TE_Barca_dd",
  TE_THIEF = "TE_Ladro_dd",
  TR = "AverageDevAbs_TR",
  TD = "TD",
  CHILD_ORIENTATION = "OTm_orientation",
  CHILD_MANAGEMENT = "OTm_management",
  CHILD_TOTAL = "OTm_total",
  PARENT_TOTAL = "QST_parent_total",
  TEACHER_TOTAL = "QST_teacher_total"
)

# These 6 scaled scores are summed for the SoTQ conversion table
SOTQ_COMPONENTS <- c(
  "Time Estimation (composite)",
  "Time Reproduction",
  "Time Discrimination",
  "Child questionnaire, Total",
  "Parent questionnaire, Total",
  "Teacher questionnaire, Total"
)

# -------------------------
# Load tables (once at startup)
# -------------------------

.NORMS_PATH <- find_file(file.path("scoring", "norms_lookup.csv"))
.NORMS <- NULL
.NORMS_BY_MEASURE <- NULL
.GRID_AGE_UNIT <- "years"
.NORMS_AGE_MIN <- NA_real_
.NORMS_AGE_MAX <- NA_real_

if (!is.na(.NORMS_PATH)) {
  .NORMS <- suppressWarnings(read_csv(.NORMS_PATH, show_col_types = FALSE))
  names(.NORMS) <- normalize_names(names(.NORMS))
  
  # Required columns for lookup
  req_cols <- c("measure", "age", "raw", "z_aligned")
  if (all(req_cols %in% names(.NORMS))) {
    .NORMS <- .NORMS %>%
      mutate(
        measure = as.character(.data$measure),
        age = to_num(.data$age),
        raw = to_num(.data$raw),
        z_aligned = to_num(.data$z_aligned)
      ) %>%
      filter(!is.na(.data$measure), is.finite(.data$age), is.finite(.data$raw), is.finite(.data$z_aligned))
    
    if (nrow(.NORMS) > 0) {
      .GRID_AGE_UNIT <- infer_age_unit(.NORMS$age)
      .NORMS_BY_MEASURE <- split(.NORMS, .NORMS$measure)
      .NORMS_AGE_MIN <- min(.NORMS$age, na.rm = TRUE)
      .NORMS_AGE_MAX <- max(.NORMS$age, na.rm = TRUE)
    } else {
      .NORMS <- NULL
    }
  } else {
    .NORMS <- NULL
  }
}

.SOTQ_TABLE_PATH <- find_file(file.path("scoring", "SoT_totalConversionTable.csv"))
.SOTQ_TABLE <- NULL
if (!is.na(.SOTQ_TABLE_PATH)) {
  .SOTQ_TABLE <- suppressWarnings(read_csv(.SOTQ_TABLE_PATH, show_col_types = FALSE))
  names(.SOTQ_TABLE) <- normalize_names(names(.SOTQ_TABLE))
  if (all(c("sumss", "sotq") %in% names(.SOTQ_TABLE))) {
    .SOTQ_TABLE <- .SOTQ_TABLE %>%
      mutate(sumss = to_num(.data$sumss), sotq = to_num(.data$sotq)) %>%
      filter(is.finite(.data$sumss), is.finite(.data$sotq)) %>%
      arrange(.data$sumss)
    if (nrow(.SOTQ_TABLE) == 0) .SOTQ_TABLE <- NULL
  } else {
    .SOTQ_TABLE <- NULL
  }
}

.SOTQ_CI_PATH <- find_file(file.path("scoring", "SoTQ_CI_params.csv"))
.SOTQ_HALF90 <- NA_real_
if (!is.na(.SOTQ_CI_PATH)) {
  ci <- suppressWarnings(read_csv(.SOTQ_CI_PATH, show_col_types = FALSE))
  names(ci) <- normalize_names(names(ci))
  if ("halfwidth90" %in% names(ci)) {
    hw <- to_num(ci$halfwidth90[1])
    if (is.finite(hw)) .SOTQ_HALF90 <- hw
  }
}

.SOTQ_SPREAD_PATH <- find_file(file.path("scoring", "SoTQ_profileSpreadThreshold.csv"))
.SOTQ_SPREAD_THRESH <- NA_real_
if (!is.na(.SOTQ_SPREAD_PATH)) {
  sp <- suppressWarnings(read_csv(.SOTQ_SPREAD_PATH, show_col_types = FALSE))
  names(sp) <- normalize_names(names(sp))
  if ("threshold" %in% names(sp)) {
    th <- to_num(sp$threshold[1])
    if (is.finite(th)) .SOTQ_SPREAD_THRESH <- th
  }
}

# -------------------------
# Norm lookup (interpolate on raw within nearest ages, then interpolate in age, on aligned z)
# -------------------------

lookup_z_from_grid <- function(measure, age_grid_value, raw_value) {
  if (is.null(.NORMS_BY_MEASURE) || is.null(.NORMS_BY_MEASURE[[measure]])) {
    return(list(z = NA_real_, note = "missing norms"))
  }
  
  dfm <- .NORMS_BY_MEASURE[[measure]]
  ages <- sort(unique(dfm$age))
  if (!length(ages)) return(list(z = NA_real_, note = "missing norms"))
  
  age_in <- suppressWarnings(as.numeric(age_grid_value))
  raw_in <- suppressWarnings(as.numeric(raw_value))
  if (!is.finite(age_in) || !is.finite(raw_in)) return(list(z = NA_real_, note = NA_character_))
  
  age_min <- min(ages); age_max <- max(ages)
  age_cap <- clamp(age_in, age_min, age_max)
  age_clamped <- !isTRUE(all.equal(age_in, age_cap))
  
  a1 <- max(ages[ages <= age_cap])
  a2 <- min(ages[ages >= age_cap])
  
  z_at_age <- function(a) {
    dfa <- dfm %>%
      filter(.data$age == a) %>%
      arrange(.data$raw) %>%
      select(.data$raw, .data$z_aligned)
    
    if (nrow(dfa) < 2) return(list(z = NA_real_, raw_clamped = FALSE))
    
    raw_min <- min(dfa$raw, na.rm = TRUE)
    raw_max <- max(dfa$raw, na.rm = TRUE)
    raw_cap <- clamp(raw_in, raw_min, raw_max)
    raw_clamped <- !isTRUE(all.equal(raw_in, raw_cap))
    
    list(
      z = lin_interp(dfa$raw, dfa$z_aligned, raw_in),
      raw_clamped = raw_clamped
    )
  }
  
  r1 <- z_at_age(a1)
  r2 <- z_at_age(a2)
  
  if (!is.finite(r1$z) && !is.finite(r2$z)) return(list(z = NA_real_, note = NA_character_))
  
  if (isTRUE(all.equal(a1, a2))) {
    z <- r1$z
    raw_clamped <- isTRUE(r1$raw_clamped)
  } else {
    w <- (age_cap - a1) / (a2 - a1)
    z <- r1$z + w * (r2$z - r1$z)
    raw_clamped <- isTRUE(r1$raw_clamped) || isTRUE(r2$raw_clamped)
  }
  
  outside <- character(0)
  if (age_clamped) outside <- c(outside, "age outside norms range")
  if (raw_clamped) outside <- c(outside, "raw outside norms range")
  note <- if (length(outside)) paste(outside, collapse = ", ") else NA_character_
  
  list(z = z, note = note)
}

# -------------------------
# PsychoPy scoring (TR and TD) - aligned with TimeTasks.R rules
# -------------------------

score_td_from_csv <- function(path) {
  d <- tryCatch(read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(d)) d <- tryCatch(read.csv(path, header = TRUE, sep = ","), error = function(e) NULL)
  if (is.null(d)) stop("Could not read the TD file as CSV.")
  
  names(d) <- normalize_names(names(d))
  if (!"ratio" %in% names(d)) stop("TD file does not contain a 'ratio' column.")
  
  if ("key_resp_corr" %in% names(d)) {
    d$correct <- d$key_resp_corr
  } else if ("correct_response_target" %in% names(d)) {
    d$correct <- d$correct_response_target
  } else if ("correct" %in% names(d)) {
    d$correct <- d$correct
  } else {
    stop("TD file does not contain a recognized correctness column.")
  }
  
  d$ratio <- to_num(d$ratio)
  d$correct <- to_num(d$correct)
  
  d <- d[!is.na(d$ratio), , drop = FALSE]
  d <- d[d$ratio != 3, , drop = FALSE]
  if (nrow(d) < 3) stop("TD file has too few usable trials after filtering.")
  
  rev_flag <- rep(0L, nrow(d))
  for (i in 2:nrow(d)) {
    if (is.finite(d$correct[i]) && is.finite(d$correct[i - 1]) && (d$correct[i] != d$correct[i - 1])) {
      rev_flag[i] <- 1L
    }
  }
  
  rev_rows_all <- which(rev_flag == 1L)
  if (length(rev_rows_all) < 3) stop("TD has fewer than 3 reversals, cannot compute the score.")
  
  rev_rows <- rev_rows_all[-c(1, 2)] # drop first 2 reversals
  ratios <- d$ratio[rev_rows]
  ratios <- ratios[is.finite(ratios)]
  if (length(ratios) < 1) stop("TD has no reversals left after dropping the first two.")
  
  last_n <- min(6, length(ratios))
  score <- mean(tail(ratios, last_n))
  
  list(
    raw = score,
    message = sprintf("TD, the mean of the last %d reversals is Ratio = %.3f.", last_n, score)
  )
}

score_tr_from_csv <- function(path) {
  d <- tryCatch(read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(d)) d <- tryCatch(read.table(path, header = TRUE, sep = ","), error = function(e) NULL)
  if (is.null(d) || ncol(d) < 3) d <- tryCatch(read.table(path, header = TRUE), error = function(e) NULL)
  if (is.null(d)) stop("Could not read the TR file.")
  
  names(d) <- normalize_names(names(d))
  task_cols <- names(d)[startsWith(names(d), "task")]
  
  if (length(task_cols) > 0) {
    if (!"task1" %in% names(d)) stop("TR, task columns detected but 'Task1' is missing.")
    d <- d[!is.na(d$task1), , drop = FALSE]
    if (nrow(d) < 1) stop("TR, no task rows found.")
    
    target <- rep(NA_real_, nrow(d))
    mapping <- c(
      task1 = 8, task2 = 5, task3 = 2, task4 = 6, task5 = 12,
      task6 = 4, task7 = 9, task8 = 10, task9 = 3, task10 = 11, task11 = 7
    )
    for (nm in names(mapping)) {
      if (nm %in% names(d)) target[d[[nm]] == 1] <- mapping[[nm]]
    }
    
    resp_col <- "key_resp_10_duration"
    if (!resp_col %in% names(d)) {
      alt <- names(d)[str_detect(names(d), "duration")]
      if (length(alt) < 1) stop("TR, could not find a response duration column.")
      resp_col <- alt[1]
    }
    
    response <- suppressWarnings(as.numeric(gsub("\\[|\\]", "", d[[resp_col]])))
    out <- tibble(target = target, response = response)
    
  } else {
    if (!all(c("type", "duration", "reproducedtime") %in% names(d))) {
      stop("TR, unrecognized file format, expected Task* columns or (type, duration, reproducedTime).")
    }
    d <- d[d$type != "practice", , drop = FALSE]
    out <- tibble(
      target = to_num(d$duration) / 1000,
      response = to_num(d$reproducedtime) / 1000
    )
  }
  
  out <- out %>%
    mutate(devabs = 100 * abs(.data$response - .data$target) / .data$target) %>%
    filter(.data$response > 0.100, .data$response < 36)
  
  if (nrow(out) < 1) stop("TR, no valid trials after filtering (<100 ms or >36 s removed).")
  
  score <- mean(out$devabs, na.rm = TRUE)
  list(
    raw = score,
    message = sprintf("TR, absolute average deviation is %.2f%%.", score)
  )
}

# -------------------------
# SoTQ conversion and interpretation
# -------------------------

lookup_sotq_from_sumss <- function(sumSS) {
  if (is.null(.SOTQ_TABLE)) return(list(sotq = NA_real_, note = "SoTQ table missing"))
  sumSS <- as.integer(round(to_num(sumSS)))
  if (!is.finite(sumSS)) return(list(sotq = NA_real_, note = NA_character_))
  
  idx <- which(.SOTQ_TABLE$sumss == sumSS)
  if (length(idx) >= 1) {
    sotq <- .SOTQ_TABLE$sotq[idx[1]]
    return(list(sotq = round(sotq), note = NA_character_))
  }
  
  # Should not happen if table covers all integer sums, keep robust:
  list(sotq = NA_real_, note = "sumSS not found in SoTQ conversion table")
}

interpret_sotq <- function(sotq) {
  if (!is.finite(sotq)) return(list(label = "Not computed", color = "grey60"))
  if (sotq < 85)   return(list(label = "Low",     color = "#c2182b"))
  if (sotq <= 115) return(list(label = "Typical", color = "#1b9e77"))
  list(label = "High", color = "#2166ac")
}

# -------------------------
# UI
# -------------------------

ui <- fluidPage(
  titlePanel("Sense of Time Battery, scoring and standardization"),
  sidebarLayout(
    sidebarPanel(
      h4("Age"),
      fluidRow(
        column(6, numericInput("age_years", "Years (integer)", value = 8, min = 0, max = 25, step = 1)),
        column(6, numericInput("age_months", "Months (0 to 11)", value = 0, min = 0, max = 11, step = 1))
      ),
      tags$hr(),
      
      h4("Time Estimation, 30 s videos"),
      helpText("Enter absolute deviation from 30 seconds, in seconds."),
      numericInput("te_boat", "Boat, |estimate - 30|", value = NA, min = 0),
      numericInput("te_thief", "Thief, |estimate - 30|", value = NA, min = 0),
      
      tags$hr(),
      h4("Time Reproduction"),
      fileInput("tr_file", "Upload TR PsychoPy or OpenSesame CSV", accept = c(".csv")),
      numericInput("tr_manual", "Manual entry, absolute average percent deviation", value = NA, min = 0),
      
      tags$hr(),
      h4("Time Discrimination"),
      fileInput("td_file", "Upload TD PsychoPy CSV", accept = c(".csv")),
      numericInput("td_manual", "Manual entry, Ratio", value = NA, min = 0),
      
      tags$hr(),
      h4("Child questionnaire"),
      numericInput("child_orientation", "Time orientation, raw", value = NA, min = 0),
      numericInput("child_management", "Time management, raw", value = NA, min = 0),
      numericInput("child_total", "Total, raw, auto computed if both domains are provided", value = NA, min = 0),
      
      tags$hr(),
      h4("Parent questionnaire"),
      numericInput("parent_total", "Total, raw", value = NA, min = 0),
      
      tags$hr(),
      h4("Teacher questionnaire"),
      numericInput("teacher_total", "Total, raw", value = NA, min = 0),
      
      tags$hr(),
      actionButton("compute", "Compute", class = "btn-primary")
    ),
    
    mainPanel(
      uiOutput("tables_status"),
      tags$hr(),
      uiOutput("warnings_top"),
      tags$hr(),
      uiOutput("sotq_box"),
      tags$hr(),
      h4("Profile, scaled scores (1-19)"),
      plotOutput("profile_plot", height = "420px"),
      tags$hr(),
      h4("Results"),
      tableOutput("results_table"),
      tags$hr(),
      h4("Processing details"),
      verbatimTextOutput("task_messages")
    )
  )
)

# -------------------------
# Server
# -------------------------

server <- function(input, output, session) {
  
  # If norms are loaded, tighten the age "Years" UI range to something reasonable
  observe({
    if (is.null(.NORMS) || !is.finite(.NORMS_AGE_MAX)) return()
    
    max_years <- if (.GRID_AGE_UNIT == "years") {
      ceiling(.NORMS_AGE_MAX)
    } else {
      ceiling(.NORMS_AGE_MAX / 12)
    }
    
    # Keep some safety, but make it hard to type extreme ages by mistake
    max_years <- clamp(max_years, 1, 30)
    
    updateNumericInput(session, "age_years", max = max_years)
  })
  
  output$tables_status <- renderUI({
    items <- list()
    
    if (is.null(.NORMS)) {
      items <- c(items, list(div(
        style = "color:#b30000; font-weight:600;",
        "ERROR: Norms file not found or invalid, expected scoring/norms_lookup.csv."
      )))
    } else {
      age_min_ym <- grid_age_to_ym(.NORMS_AGE_MIN, .GRID_AGE_UNIT)
      age_max_ym <- grid_age_to_ym(.NORMS_AGE_MAX, .GRID_AGE_UNIT)
      items <- c(items, list(div(
        style = "color:#1f7a1f; font-weight:600;",
        sprintf(
          "Norms loaded (%s), %s age scale. Normative age range is approximately %d years %d months to %d years %d months.",
          .NORMS_PATH,
          .GRID_AGE_UNIT,
          age_min_ym$years, age_min_ym$months,
          age_max_ym$years, age_max_ym$months
        )
      )))
    }
    
    if (is.null(.SOTQ_TABLE)) {
      items <- c(items, list(div(
        style = "color:#b30000; font-weight:600;",
        "WARNING: SoTQ conversion table missing, SoTQ cannot be computed."
      )))
    } else {
      items <- c(items, list(div(
        style = "color:#1f7a1f; font-weight:600;",
        sprintf("SoTQ conversion table loaded, %s.", .SOTQ_TABLE_PATH)
      )))
    }
    
    if (!is.finite(.SOTQ_HALF90)) {
      items <- c(items, list(div(
        style = "color:#8a6d3b; font-weight:600;",
        "NOTE: SoTQ CI parameters missing or invalid, CI will not be shown."
      )))
    }
    
    if (!is.finite(.SOTQ_SPREAD_THRESH)) {
      items <- c(items, list(div(
        style = "color:#8a6d3b; font-weight:600;",
        "NOTE: Profile spread threshold missing, heterogeneity warning disabled."
      )))
    }
    
    do.call(tagList, items)
  })
  
  computed <- eventReactive(input$compute, {
    
    validate(
      need(!is.null(.NORMS), "Norms file missing, cannot compute."),
      need(is.finite(input$age_years) && input$age_years >= 0, "Age (years) must be a non-negative integer."),
      need(is_whole_number(input$age_years), "Age (years) must be an integer."),
      need(is.finite(input$age_months) && input$age_months >= 0 && input$age_months <= 11, "Age (months) must be between 0 and 11.")
    )
    
    # Compute age in the norms grid units
    age_grid_entered <- age_to_grid_unit(input$age_years, input$age_months, .GRID_AGE_UNIT)
    validate(need(is.finite(age_grid_entered), "Age is invalid."))
    
    info_msgs <- character(0)
    warn_msgs <- character(0)
    
    # SERIOUS AGE GUARD: clamp to normative extremes and raise a very prominent warning
    age_grid_used <- age_grid_entered
    if (is.finite(.NORMS_AGE_MIN) && is.finite(.NORMS_AGE_MAX)) {
      age_grid_used <- clamp(age_grid_entered, .NORMS_AGE_MIN, .NORMS_AGE_MAX)
      if (!isTRUE(all.equal(age_grid_entered, age_grid_used))) {
        entered_ym <- grid_age_to_ym(age_grid_entered, .GRID_AGE_UNIT)
        used_ym <- grid_age_to_ym(age_grid_used, .GRID_AGE_UNIT)
        warn_msgs <- c(
          warn_msgs,
          sprintf(
            "VERY IMPORTANT: entered age (%d years %d months) is outside the normative range. The app uses the nearest extreme (%d years %d months).",
            entered_ym$years, entered_ym$months, used_ym$years, used_ym$months
          )
        )
      }
    }
    
    # TD: upload has priority over manual, but manual+upload is a top warning
    td_raw <- NA_real_
    td_upload_ok <- FALSE
    td_manual_present <- is.finite(to_num(input$td_manual))
    td_upload_present <- !is.null(input$td_file) && nzchar(input$td_file$datapath)
    
    if (td_upload_present) {
      td_res <- tryCatch(score_td_from_csv(input$td_file$datapath), error = function(e) e)
      if (inherits(td_res, "error")) {
        warn_msgs <- c(warn_msgs, paste0("TD upload error, ", td_res$message))
      } else {
        td_raw <- td_res$raw
        td_upload_ok <- TRUE
        info_msgs <- c(info_msgs, td_res$message)
      }
    }
    
    if (td_upload_present && td_manual_present) {
      warn_msgs <- c(
        warn_msgs,
        "VERY IMPORTANT: TD, both a file upload and a manual value were provided. The app uses the uploaded file score and ignores the manual value."
      )
    }
    
    if (!td_upload_ok && td_manual_present) {
      td_raw <- to_num(input$td_manual)
      info_msgs <- c(info_msgs, sprintf("TD, manual entry used, Ratio = %.3f.", td_raw))
    }
    
    # TR: upload has priority over manual, but manual+upload is a top warning
    tr_raw <- NA_real_
    tr_upload_ok <- FALSE
    tr_manual_present <- is.finite(to_num(input$tr_manual))
    tr_upload_present <- !is.null(input$tr_file) && nzchar(input$tr_file$datapath)
    
    if (tr_upload_present) {
      tr_res <- tryCatch(score_tr_from_csv(input$tr_file$datapath), error = function(e) e)
      if (inherits(tr_res, "error")) {
        warn_msgs <- c(warn_msgs, paste0("TR upload error, ", tr_res$message))
      } else {
        tr_raw <- tr_res$raw
        tr_upload_ok <- TRUE
        info_msgs <- c(info_msgs, tr_res$message)
      }
    }
    
    if (tr_upload_present && tr_manual_present) {
      warn_msgs <- c(
        warn_msgs,
        "VERY IMPORTANT: TR, both a file upload and a manual value were provided. The app uses the uploaded file score and ignores the manual value."
      )
    }
    
    if (!tr_upload_ok && tr_manual_present) {
      tr_raw <- to_num(input$tr_manual)
      info_msgs <- c(info_msgs, sprintf("TR, manual entry used, deviation = %.2f%%.", tr_raw))
    }
    
    # TE raw inputs
    te_boat <- to_num(input$te_boat)
    te_thief <- to_num(input$te_thief)
    
    # Child questionnaire, smart resolution of Total vs domains
    child_orientation <- to_num(input$child_orientation)
    child_management <- to_num(input$child_management)
    child_total_entered <- to_num(input$child_total)
    
    have_o <- is.finite(child_orientation)
    have_m <- is.finite(child_management)
    have_t <- is.finite(child_total_entered)
    
    child_total_used <- NA_real_
    child_total_note <- NA_character_
    
    if (have_o && have_m) {
      total_calc <- child_orientation + child_management
      child_total_used <- total_calc
      
      if (!have_t) {
        info_msgs <- c(info_msgs, sprintf("Child questionnaire, total computed as the sum of domains, Total = %.0f.", total_calc))
        child_total_note <- "Total computed from domains"
      } else {
        if (!isTRUE(all.equal(child_total_entered, total_calc))) {
          warn_msgs <- c(
            warn_msgs,
            sprintf(
              "Child questionnaire, entered Total (%.0f) differs from the sum of domains (%.0f). The app uses the sum of domains for internal consistency.",
              child_total_entered, total_calc
            )
          )
          child_total_note <- "Total computed from domains, manual total ignored"
        } else {
          child_total_note <- "Total confirmed as sum of domains"
        }
      }
      
    } else if (have_t) {
      child_total_used <- child_total_entered
      child_total_note <- "Total entered by user"
      
      if (have_o || have_m) {
        miss <- character(0)
        if (!have_o) miss <- c(miss, "Time orientation")
        if (!have_m) miss <- c(miss, "Time management")
        warn_msgs <- c(
          warn_msgs,
          paste0(
            "Child questionnaire, domains are incomplete, the app uses the entered total for scoring. Missing: ",
            paste(miss, collapse = ", "),
            "."
          )
        )
      }
      
    } else {
      if (have_o || have_m) {
        miss <- character(0)
        if (!have_o) miss <- c(miss, "Time orientation")
        if (!have_m) miss <- c(miss, "Time management")
        warn_msgs <- c(
          warn_msgs,
          paste0(
            "Child questionnaire, total cannot be computed because a domain score is missing. Missing: ",
            paste(miss, collapse = ", "),
            "."
          )
        )
      }
    }
    
    # Parent and teacher
    parent_total <- to_num(input$parent_total)
    teacher_total <- to_num(input$teacher_total)
    
    # TE composite, compute only if both videos are present
    te_z <- NA_real_
    te_ss <- NA_real_
    te_note <- NA_character_
    if (is.finite(te_boat) && is.finite(te_thief)) {
      z1 <- lookup_z_from_grid(MEASURES$TE_BOAT, age_grid_used, te_boat)
      z2 <- lookup_z_from_grid(MEASURES$TE_THIEF, age_grid_used, te_thief)
      te_z <- mean(c(z1$z, z2$z))
      te_ss <- scaled_from_z(te_z)
      
      te_note <- paste(na.omit(c(z1$note, z2$note)), collapse = "; ")
      if (!nzchar(te_note)) te_note <- NA_character_
      
      info_msgs <- c(
        info_msgs,
        sprintf(
          "Time Estimation details, Boat deviation = %.2f s, Thief deviation = %.2f s. Composite is the mean of the two age-standardized aligned z scores.",
          te_boat, te_thief
        )
      )
    }
    
    add_row <- function(label, raw_val, measure_name, extra_note = NA_character_) {
      if (!is.finite(raw_val)) {
        return(tibble(subtest = label, raw = NA_real_, ss = NA_real_, z = NA_real_, note = NA_character_))
      }
      st <- lookup_z_from_grid(measure_name, age_grid_used, raw_val)
      tibble(
        subtest = label,
        raw = raw_val,
        ss = scaled_from_z(st$z),
        z = st$z,
        note = if (!is.na(extra_note)) extra_note else st$note
      )
    }
    
    results <- tibble(
      subtest = character(0),
      raw = numeric(0),
      ss = numeric(0),
      z = numeric(0),
      note = character(0)
    )
    
    results <- bind_rows(results, tibble(
      subtest = "Time Estimation (composite)",
      raw = if (is.finite(te_boat) && is.finite(te_thief)) mean(c(te_boat, te_thief)) else NA_real_,
      ss = te_ss,
      z = te_z,
      note = te_note
    ))
    results <- bind_rows(results, add_row("Time Reproduction", tr_raw, MEASURES$TR))
    results <- bind_rows(results, add_row("Time Discrimination", td_raw, MEASURES$TD))
    
    results <- bind_rows(results, add_row("Child questionnaire, Time orientation", child_orientation, MEASURES$CHILD_ORIENTATION))
    results <- bind_rows(results, add_row("Child questionnaire, Time management", child_management, MEASURES$CHILD_MANAGEMENT))
    results <- bind_rows(results, add_row("Child questionnaire, Total", child_total_used, MEASURES$CHILD_TOTAL, extra_note = child_total_note))
    
    results <- bind_rows(results, add_row("Parent questionnaire, Total", parent_total, MEASURES$PARENT_TOTAL))
    results <- bind_rows(results, add_row("Teacher questionnaire, Total", teacher_total, MEASURES$TEACHER_TOTAL))
    
    # SoTQ computation
    sotq <- NA_real_
    sumSS <- NA_integer_
    ci_low <- NA_real_
    ci_high <- NA_real_
    sotq_note <- NA_character_
    deltaSS <- NA_real_
    spread_warn <- FALSE
    missing_components <- character(0)
    
    if (!is.null(.SOTQ_TABLE)) {
      comp <- results %>% filter(.data$subtest %in% SOTQ_COMPONENTS)
      
      for (nm in SOTQ_COMPONENTS) {
        row_nm <- comp %>% filter(.data$subtest == nm)
        if (nrow(row_nm) == 0 || !is.finite(row_nm$ss[1])) missing_components <- c(missing_components, nm)
      }
      
      if (length(missing_components) == 0) {
        sumSS <- as.integer(sum(comp$ss))
        sotq_res <- lookup_sotq_from_sumss(sumSS)
        sotq <- sotq_res$sotq
        sotq_note <- sotq_res$note
        
        if (is.finite(sotq) && is.finite(.SOTQ_HALF90)) {
          ci_low <- round(sotq - .SOTQ_HALF90)
          ci_high <- round(sotq + .SOTQ_HALF90)
        }
        
        deltaSS <- max(comp$ss) - min(comp$ss)
        if (is.finite(.SOTQ_SPREAD_THRESH) && is.finite(deltaSS) && deltaSS > .SOTQ_SPREAD_THRESH) {
          spread_warn <- TRUE
        }
      } else {
        warn_msgs <- c(
          warn_msgs,
          paste0(
            "VERY IMPORTANT: SoTQ not computed because required components are missing: ",
            paste(missing_components, collapse = ", "),
            "."
          )
        )
      }
    }
    
    list(
      info = info_msgs,
      warnings = warn_msgs,
      results = results,
      sotq = sotq,
      sumSS = sumSS,
      ci_low = ci_low,
      ci_high = ci_high,
      sotq_note = sotq_note,
      deltaSS = deltaSS,
      spread_warn = spread_warn,
      missing_components = missing_components
    )
  })
  
  # Top warnings (must stay on top)
  output$warnings_top <- renderUI({
    req(input$compute)
    res <- computed()
    req(is.list(res))
    
    msgs <- res$warnings
    if (length(msgs) == 0) return(NULL)
    
    div(
      style = "border-left:8px solid #b2182b; padding:14px; background:#fff5f5; box-shadow:0 2px 8px rgba(0,0,0,0.08);",
      h4(style = "margin-top:0; margin-bottom:10px; color:#b2182b;", "Warnings"),
      tagList(lapply(msgs, function(m) {
        div(style = "color:#b2182b; font-weight:600; margin-bottom:6px;", HTML(paste0("&#9888; ", m)))
      }))
    )
  })
  
  # Bottom details (informational messages)
  output$task_messages <- renderText({
    req(input$compute)
    res <- computed()
    req(is.list(res))
    
    msgs <- res$info
    if (length(msgs) == 0) return("No processing details.")
    paste(msgs, collapse = "\n")
  })
  
  output$sotq_box <- renderUI({
    req(input$compute)
    res <- computed()
    req(is.list(res))
    
    if (!is.finite(res$sotq)) {
      missing_txt <- if (length(res$missing_components) > 0) {
        paste0("Missing components: ", paste(res$missing_components, collapse = ", "), ".")
      } else if (is.null(.SOTQ_TABLE)) {
        "SoTQ table missing, SoTQ cannot be computed."
      } else {
        "SoTQ not computed."
      }
      
      div(
        style = "border-left:10px solid #888888; padding:22px; background:#fafafa; box-shadow:0 2px 10px rgba(0,0,0,0.10);",
        h2(HTML("Sense of Time Quotient (SoTQ)")),
        p("SoTQ is computed only when all 6 components are available: Time Estimation composite, Time Reproduction, Time Discrimination, Child Total, Parent Total, Teacher Total."),
        p(style = "color:#555555;", missing_txt)
      )
    } else {
      it <- interpret_sotq(res$sotq)
      
      ci_txt <- if (is.finite(res$ci_low) && is.finite(res$ci_high)) {
        sprintf("90%% confidence interval: [%d, %d].", as.integer(res$ci_low), as.integer(res$ci_high))
      } else {
        "90% confidence interval not available."
      }
      
      spread_txt <- if (isTRUE(res$spread_warn)) {
        sprintf("Caution: profile spread is high, deltaSS = %d, exceeds the threshold.", as.integer(res$deltaSS))
      } else if (is.finite(res$deltaSS)) {
        sprintf("Profile spread: deltaSS = %d.", as.integer(res$deltaSS))
      } else {
        NULL
      }
      
      note_txt <- if (!is.na(res$sotq_note)) paste0("Note: ", res$sotq_note, ".") else NULL
      
      div(
        style = sprintf("border-left:10px solid %s; padding:22px; background:#ffffff; box-shadow:0 2px 10px rgba(0,0,0,0.10);", it$color),
        h2(HTML(sprintf("Sense of Time Quotient (SoTQ): <b>%d</b>", as.integer(res$sotq)))),
        p(HTML(sprintf("Interpretation: <b style='color:%s'>%s</b>.", it$color, it$label))),
        p(sprintf("Sum of scaled scores (sumSS): %d.", as.integer(res$sumSS))),
        p(ci_txt),
        if (!is.null(spread_txt) && isTRUE(res$spread_warn)) div(style = "color:#b2182b; font-weight:700; margin-top:6px;", HTML(paste0("&#9888; ", spread_txt))),
        if (!is.null(spread_txt) && !isTRUE(res$spread_warn)) div(style = "color:#444444; margin-top:6px;", spread_txt),
        if (!is.null(note_txt)) p(note_txt)
      )
    }
  })
  
  output$results_table <- renderTable({
    req(input$compute)
    res <- computed()
    req(is.list(res))
    
    res$results %>%
      mutate(
        raw = ifelse(is.na(.data$raw), NA, round(.data$raw, 3)),
        ss = ifelse(is.na(.data$ss), NA, as.integer(.data$ss)),
        z = ifelse(is.na(.data$z), NA, round(.data$z, 3)),
        note = ifelse(is.na(.data$note), "", .data$note)
      ) %>%
      transmute(
        `Subtest` = .data$subtest,
        `Raw score` = .data$raw,
        `Scaled score (SS, 1-19)` = .data$ss,
        `Aligned z` = .data$z,
        `Note` = .data$note
      )
  }, striped = TRUE, hover = TRUE, spacing = "s", na = "")
  
  output$profile_plot <- renderPlot({
    req(input$compute)
    res <- computed()
    req(is.list(res))
    
    df_plot <- res$results %>%
      filter(is.finite(.data$ss)) %>%
      mutate(x = seq_len(dplyr::n()))
    
    validate(need(nrow(df_plot) > 0, "Enter at least one score, then click Compute."))
    
    n <- nrow(df_plot)
    
    p <- ggplot() +
      # Typical range band (robust, simple)
      geom_rect(
        aes(xmin = 0.5, xmax = n + 0.5, ymin = 7, ymax = 13),
        alpha = 0.20,
        fill = "grey70",
        color = NA
      ) +
      geom_point(data = df_plot, aes(x = .data$x, y = .data$ss), size = 5)
    
    if (n >= 2) {
      p <- p + geom_line(data = df_plot, aes(x = .data$x, y = .data$ss, group = 1), linewidth = 0.6)
    }
    
    p +
      scale_y_continuous(limits = c(1, 19), breaks = 1:19) +
      scale_x_continuous(
        limits = c(0.5, n + 0.5),
        breaks = df_plot$x,
        labels = df_plot$subtest
      ) +
      labs(
        x = NULL,
        y = "Scaled score (SS, 1-19)",
        title = "Subtest profile, scaled scores"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        text = element_text(size = PLOT_TEXT_SIZE),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 15, hjust = 1)
      )
  })
}

shinyApp(ui, server)
