# VERSION v2: sources DR_01_functions_stepwise_by_arm_v2.R.
#######################################################################################################
############################ DR_02_data_prep_stepwise_by_arm.R ########################################
#######################################################################################################

# This script loads all trial CSV files, excludes placebo, creates one dataset per trial x active arm,
# imputes active-arm zero doses, groups sparse doses, and prepares the longitudinal variables.

source("01_functions.R")
load_required_packages()

###################################### User settings #######################################

DATA_DIR <- "/Users/kchalkou/Desktop/Projects_ongoing/Dose-response/Saved Data"
CSV_PATTERN <- "\\.csv$"

PLACEBO_LABELS <- c("PLACEBO")
EXCLUDE_TREATMENTS <- c()

VISIT_MAX <- 9
MAX_FOLLOWUP_VISIT <- VISIT_MAX
ASSUME_VISIT_IS_DAYS <- NULL

MIN_N_PER_DOSE <- 10
ZERO_DOSE_IMPUTE <- "nearest"   # options: "nearest", "previous", "next"

SIMULATE_SIDE_EFFECTS_IF_MISSING <- TRUE
SIDE_EFFECT_SEED <- 2025

###################################### Read and standardise all trials #######################################

csv_files <- list.files(
  DATA_DIR,
  pattern = CSV_PATTERN,
  full.names = TRUE
)

if (length(csv_files) == 0) {
  stop("No CSV files found in DATA_DIR: ", DATA_DIR)
}

trial_names <- tools::file_path_sans_ext(basename(csv_files))

raw_trials <- purrr::map2(
  csv_files,
  trial_names,
  function(path, trial_name) {
    read.csv(path, stringsAsFactors = FALSE) |>
      standardise_trial_columns(trial_name = trial_name) |>
      convert_visit_to_week(
        visit_max = VISIT_MAX,
        assume_visit_is_days = ASSUME_VISIT_IS_DAYS
      ) |>
      keep_last_record_per_patient_visit()
  }
)

names(raw_trials) <- trial_names

raw_all <- dplyr::bind_rows(raw_trials, .id = "trial_file")

################################################################################
######################## OBSERVED PLACEBO IMPROVEMENT ###########################
################################################################################

PLACEBO_MAX_VISIT <- 9

placebo_patient_visits <- raw_all |>
  dplyr::filter(
    toupper(as.character(treat)) %in% toupper(PLACEBO_LABELS),
    visit >= 0,
    visit <= PLACEBO_MAX_VISIT
  ) |>
  dplyr::arrange(trial_name, pid, visit) |>
  dplyr::group_by(trial_name, pid) |>
  dplyr::mutate(
    placebo_outcome_0 = outcome[which.min(visit)],
    placebo_delta_outcome = placebo_outcome_0 - outcome
  ) |>
  dplyr::ungroup()

placebo_counts_by_trial <- placebo_patient_visits |>
  dplyr::filter(
    visit > 0,
    visit <= 8,
    !is.na(placebo_delta_outcome)
  ) |>
  dplyr::group_by(trial_name, visit) |>
  dplyr::summarise(
    placebo_n_patients = dplyr::n_distinct(pid),
    placebo_n_visit_rows = dplyr::n(),
    .groups = "drop"
  )

placebo_counts_by_trial

placebo_observed_by_trial <- placebo_patient_visits |>
  dplyr::filter(
    visit > 0,
    visit <= PLACEBO_MAX_VISIT,
    !is.na(placebo_delta_outcome)
  ) |>
  dplyr::group_by(trial_name, visit) |>
  dplyr::summarise(
    n_placebo = dplyr::n_distinct(pid),
    placebo_mean_delta = mean(placebo_delta_outcome, na.rm = TRUE),
    placebo_sd_delta = stats::sd(placebo_delta_outcome, na.rm = TRUE),
    placebo_se_delta = placebo_sd_delta / sqrt(n_placebo),
    placebo_lower_95 = placebo_mean_delta - stats::qnorm(0.975) * placebo_se_delta,
    placebo_upper_95 = placebo_mean_delta + stats::qnorm(0.975) * placebo_se_delta,
    .groups = "drop"
  )

placebo_observed_by_trial
###################################### Define active arms #######################################

active_all <- raw_all |>
  dplyr::filter(
    !(toupper(treat) %in% toupper(PLACEBO_LABELS)),
    !(toupper(treat) %in% toupper(EXCLUDE_TREATMENTS))
  ) |>
  dplyr::mutate(
    arm_name = paste0(trial_name, "__", treat),
    arm_name = safe_name(arm_name)
  )

active_arm_raw <- split(active_all, active_all$arm_name)

###################################### Preprocess each active arm #######################################

active_arm_datasets <- purrr::imap(
  active_arm_raw,
  function(dat, arm_name) {
    dat |>
      impute_active_zero_dose(method = ZERO_DOSE_IMPUTE) |>
      group_sparse_doses(min_n_per_dose = MIN_N_PER_DOSE) |>
      simulate_side_effects_if_needed(seed = SIDE_EFFECT_SEED)
  }
)

###################################### Prepare analysis variables per active arm #######################################

analysis_arms <- purrr::imap(
  active_arm_datasets,
  function(dat, arm_name) {
    dose_levels <- sort(unique(dat$dose[dat$visit > 0 & dat$dose > 0 & !is.na(dat$dose)]))
    dose_history_levels <- sort(unique(c(0, dose_levels)))

    prepare_dose_response_arm_data(
      data = dat,
      visit_min = 0,
      visit_max = VISIT_MAX,
      max_followup_visit = MAX_FOLLOWUP_VISIT,
      dose_levels = dose_levels,
      dose_history_levels = dose_history_levels
    )
  }
)

###################################### Basic summaries #######################################

arm_summary <- purrr::imap_dfr(
  analysis_arms,
  function(dat, arm_name) {
    summarise_dose_response_data(dat) |>
      dplyr::mutate(arm_name = arm_name, .before = 1)
  }
)

# Separate sparse-dose grouping table per arm.
# Do not bind these with dose factors; arms may have different dose levels.
dose_grouping_tables <- purrr::imap(
  analysis_arms,
  function(dat, arm_name) {
    if (!all(c("dose_before_grouping", "dose_grouped") %in% names(dat))) {
      return(tibble::tibble())
    }

    dat |>
      dplyr::filter(visit > 0) |>
      dplyr::count(dose_before_grouping, dose_grouped, name = "n") |>
      dplyr::arrange(dose_before_grouping, dose_grouped)
  }
)

arm_summary
names(analysis_arms)
