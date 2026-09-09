# Reviewer 1: Time Discrimination reversals, matching checks,
# and within-session stability

library(tidyverse)
library(readxl)

raw_dir <- "data/TimeDiscrimination-raw"
db_file <- "data/Database_Time_17.07.25_PD&MI.xlsx"

clean_id <- function(x) {
  x |>
    as.character() |>
    str_to_upper() |>
    str_remove_all("[^A-Z0-9]")
}

raw <- map_dfr(list.files(raw_dir, "\\.csv$", full.names = TRUE), \(f) {
  d <- read_csv(f, show_col_types = FALSE, name_repair = "unique_quiet")
  names(d) <- str_to_lower(str_trim(names(d)))

  if (!"ratio" %in% names(d)) return(tibble())

  d |>
    transmute(
      file = f,
      id = clean_id(str_extract(basename(f), "^[^_]+")),
      ratio = suppressWarnings(as.numeric(ratio))
    ) |>
    filter(!is.na(ratio))
})

chosen_files <- raw |>
  count(id, file, name = "n_trials") |>
  arrange(id, desc(n_trials), desc(file)) |>
  slice_head(n = 1, by = id)

raw <- raw |>
  semi_join(chosen_files, by = c("id", "file")) |>
  group_by(id, file) |>
  mutate(trial = row_number()) |>
  ungroup()

reversals <- raw |>
  group_by(id, file) |>
  arrange(trial, .by_group = TRUE) |>
  mutate(direction = sign(lead(ratio) - ratio)) |>
  filter(!is.na(direction), direction != 0) |>
  mutate(reversal = replace_na(direction != lag(direction), FALSE)) |>
  filter(reversal) |>
  summarise(
    reversal_values = list(ratio),
    n_reversals = n(),
    .groups = "drop"
  )

ages <- read_excel(db_file) |>
  transmute(
    id = clean_id(Code),
    age = as.numeric(Age)
  ) |>
  distinct(id, .keep_all = TRUE)

results <- chosen_files |>
  select(id, file, n_trials) |>
  left_join(reversals, by = c("id", "file")) |>
  left_join(ages, by = "id") |>
  mutate(
    n_reversals = replace_na(n_reversals, 0L),
    reversal_values = map(reversal_values, \(x) if (is.null(x)) numeric() else x),
    last6 = map(
      reversal_values,
      \(x) if (length(x) >= 8) tail(x[-c(1, 2)], 6) else numeric()
    ),
    first3_mean = map_dbl(
      last6,
      \(x) if (length(x) == 6) mean(x[1:3]) else NA_real_
    ),
    last3_mean = map_dbl(
      last6,
      \(x) if (length(x) == 6) mean(x[4:6]) else NA_real_
    ),
    abs_change = abs(last3_mean - first3_mean),
    age_band = case_when(
      age < 8 ~ "6-7 years",
      age >= 8 ~ "8-11 years",
      TRUE ~ NA_character_
    )
  )

cat("\nNumber of valid trials in selected files:\n")
print(table(results$n_trials))

cat("\nBreakdown of inclusion/exclusion:\n")
print(
  results |>
    mutate(
      age_match = !is.na(age),
      complete = n_trials == 36,
      status = case_when(
        age_match & complete ~ "Included",
        !age_match & complete ~ "No database match, complete raw file",
        age_match & !complete ~ "Database match, incomplete raw file",
        !age_match & !complete ~ "No database match, incomplete raw file"
      )
    ) |>
    count(status),
  n = Inf
)

cat("\nDatabase-matched but incomplete protocols:\n")
print(
  results |>
    filter(!is.na(age), n_trials != 36) |>
    select(id, age, n_trials, file) |>
    arrange(n_trials, id),
  n = Inf
)

cat("\nRaw IDs not matched to the main database:\n")
print(
  results |>
    filter(is.na(age)) |>
    select(id, n_trials, file) |>
    arrange(id),
  n = Inf
)

results_valid <- results |>
  filter(!is.na(age), n_trials == 36)

cat("\nValid analysis sample size:\n")
print(nrow(results_valid))

cat("\nOverall reversal distribution in valid cases:\n")
print(
  results_valid |>
    summarise(
      n = n(),
      median_reversals = median(n_reversals),
      IQR_reversals = IQR(n_reversals),
      min_reversals = min(n_reversals),
      max_reversals = max(n_reversals),
      pct_with_8plus = mean(n_reversals >= 8) * 100
    )
)

cat("\nResults by age band in valid cases:\n")
print(
  results_valid |>
    group_by(age_band) |>
    summarise(
      n = n(),
      median_reversals = median(n_reversals),
      IQR_reversals = IQR(n_reversals),
      pct_with_8plus = mean(n_reversals >= 8) * 100,
      stability_r = cor(first3_mean, last3_mean, use = "complete.obs"),
      median_abs_change = median(abs_change, na.rm = TRUE),
      .groups = "drop"
    )
)

cat("\nOverall within-session stability in valid cases:\n")
print(
  results_valid |>
    summarise(
      n = sum(!is.na(first3_mean)),
      stability_r = cor(first3_mean, last3_mean, use = "complete.obs"),
      median_abs_change = median(abs_change, na.rm = TRUE)
    )
)
