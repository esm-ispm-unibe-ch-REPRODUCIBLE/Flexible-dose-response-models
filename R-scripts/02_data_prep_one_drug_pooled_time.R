# Data preparation for one selected active drug pooled across eligible studies.
#######################################################################################################
###################################### 02_data_prep.R ##################################################
#######################################################################################################

source("01_functions.R")
load_required_packages()

###################################### User settings #######################################

DATA_DIR <- "/Users/kchalkou/Desktop/Projects_ongoing/Dose-response/Saved Data"
CSV_PATTERN <- "\\.csv$"

# Analyse this active drug only. A study is kept only if it also contains placebo.
TARGET_TREATMENT <- "PAROXETINE"
PLACEBO_LABELS <- c("PLACEBO")

# Keep the reported visit time instead of creating integer weeks.
# Use "days" if visit is reported in study days; use "weeks" if it is already in weeks.
VISIT_UNIT <- "days"
ANALYSIS_MAX_DAY <- 42
VISIT_MAX <- ANALYSIS_MAX_DAY / 7

MIN_N_PER_DOSE <- 10

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

if (!VISIT_UNIT %in% c("days", "weeks")) {
  stop("VISIT_UNIT must be 'days' or 'weeks'.")
}

trial_names <- tools::file_path_sans_ext(basename(csv_files))

raw_trials <- purrr::map2(
  csv_files,
  trial_names,
  function(path, trial_name) {
    read.csv(path, stringsAsFactors = FALSE) |>
      standardise_trial_columns(trial_name = trial_name) |>
      dplyr::filter(!is.na(visit), visit >= 0) |>
      dplyr::mutate(
        # Preserve the actual reported time; no weekly binning is performed.
        visit_day = if (VISIT_UNIT == "days") visit else 7 * visit
      ) |>
      dplyr::arrange(pid, visit_day) |>
      dplyr::group_by(pid, visit_day) |>
      dplyr::slice_tail(n = 1) |>
      dplyr::ungroup()
  }
)

names(raw_trials) <- trial_names
raw_all <- dplyr::bind_rows(raw_trials, .id = "trial_file")

###################################### Select eligible studies #######################################

target_key <- toupper(trimws(TARGET_TREATMENT))
placebo_keys <- toupper(trimws(PLACEBO_LABELS))

# Matching is case-insensitive. Keep only studies containing both the selected drug and placebo.
raw_all <- raw_all |>
  dplyr::mutate(treat_key = toupper(trimws(as.character(treat))))

study_treatment_check <- raw_all |>
  dplyr::group_by(studyid) |>
  dplyr::summarise(
    has_target = any(treat_key == target_key, na.rm = TRUE),
    has_placebo = any(treat_key %in% placebo_keys, na.rm = TRUE),
    .groups = "drop"
  )

eligible_studies <- study_treatment_check |>
  dplyr::filter(has_target, has_placebo) |>
  dplyr::pull(studyid)

if (length(eligible_studies) == 0) {
  stop("No study contains both ", TARGET_TREATMENT, " and placebo.")
}

selected_data_full <- raw_all |>
  dplyr::filter(
    studyid %in% eligible_studies,
    treat_key == target_key | treat_key %in% placebo_keys
  ) |>
  dplyr::mutate(
    # Patient IDs can repeat across studies, so make the analysis ID study-specific.
    pid_original = pid,
    pid = paste(studyid, pid_original, sep = "__")
  ) |>
  dplyr::arrange(pid, visit_day) |>
  dplyr::group_by(pid) |>
  dplyr::mutate(
    baseline_day = min(visit_day, na.rm = TRUE),
    analysis_day = visit_day - baseline_day,
    time_weeks = analysis_day / 7,
    # Existing functions expect a variable called visit; it is now continuous time in weeks.
    visit = time_weeks
  ) |>
  dplyr::ungroup()

# Keep the full follow-up object for the later censoring definition.
# The outcome analysis itself uses observations through day 42 only.
selected_data_analysis <- selected_data_full |>
  dplyr::filter(analysis_day <= ANALYSIS_MAX_DAY)

# Useful checks before modelling.
study_treatment_check
eligible_studies

selected_data_analysis |>
  dplyr::count(studyid, treat, name = "n_rows") |>
  dplyr::arrange(studyid, treat)

selected_data_analysis |>
  dplyr::group_by(studyid) |>
  dplyr::summarise(
    n_patients = dplyr::n_distinct(pid),
    min_day = min(analysis_day, na.rm = TRUE),
    max_day = max(analysis_day, na.rm = TRUE),
    .groups = "drop"
  )

###################################### Create pooled active and placebo data #######################################

active_data_full <- selected_data_full |>
  dplyr::filter(treat_key == target_key) |>
  dplyr::select(-treat_key)

placebo_data_full <- selected_data_full |>
  dplyr::filter(treat_key %in% placebo_keys) |>
  dplyr::select(-treat_key)

active_data <- selected_data_analysis |>
  dplyr::filter(treat_key == target_key) |>
  dplyr::select(-treat_key)

placebo_data <- selected_data_analysis |>
  dplyr::filter(treat_key %in% placebo_keys) |>
  dplyr::select(-treat_key)

###################################### Full follow-up status #######################################

# Use all available follow-up to identify each patient's final observed record.
active_followup_status <- active_data_full |>
  dplyr::group_by(studyid, pid) |>
  dplyr::summarise(
    last_observed_day_full = max(analysis_day, na.rm = TRUE),
    .groups = "drop"
  )

# Attach full-follow-up information to the day-42 active dataset.
active_data <- active_data |>
  dplyr::left_join(
    active_followup_status,
    by = c("studyid", "pid")
  ) |>
  dplyr::mutate(
    
    # TRUE if any later record exists anywhere in the full follow-up data.
    has_later_observation_full =
      last_observed_day_full > analysis_day,
    
    # Terminal loss of follow-up before the day-42 analysis horizon.
    no_later_observation_before_horizon =
      analysis_day < ANALYSIS_MAX_DAY &
      !has_later_observation_full
  )
################################################################################
######################## OBSERVED PLACEBO IMPROVEMENT ###########################
################################################################################

placebo_patient_visits <- placebo_data |>
  dplyr::arrange(studyid, pid, visit) |>
  dplyr::group_by(studyid, pid) |>
  dplyr::mutate(
    placebo_outcome_0 = outcome[which.min(visit)],
    placebo_delta_outcome = placebo_outcome_0 - outcome
  ) |>
  dplyr::ungroup()

# These summaries now use the exact observed time in weeks.
placebo_counts_by_trial <- placebo_patient_visits |>
  dplyr::filter(
    visit > 0,
    visit <= VISIT_MAX,
    !is.na(placebo_delta_outcome)
  ) |>
  dplyr::group_by(studyid, visit) |>
  dplyr::summarise(
    placebo_n_patients = dplyr::n_distinct(pid),
    placebo_n_visit_rows = dplyr::n(),
    .groups = "drop"
  )

placebo_observed_by_trial <- placebo_patient_visits |>
  dplyr::filter(
    visit > 0,
    visit <= VISIT_MAX,
    !is.na(placebo_delta_outcome)
  ) |>
  dplyr::group_by(studyid, visit) |>
  dplyr::summarise(
    n_placebo = dplyr::n_distinct(pid),
    placebo_mean_delta = mean(placebo_delta_outcome, na.rm = TRUE),
    placebo_sd_delta = stats::sd(placebo_delta_outcome, na.rm = TRUE),
    placebo_se_delta = placebo_sd_delta / sqrt(n_placebo),
    placebo_lower_95 = placebo_mean_delta - stats::qnorm(0.975) * placebo_se_delta,
    placebo_upper_95 = placebo_mean_delta + stats::qnorm(0.975) * placebo_se_delta,
    .groups = "drop"
  )

placebo_counts_by_trial
placebo_observed_by_trial
###################################### Prepare the two zero-dose analyses #######################################

# Add the side-effect variable before splitting the two analyses.
# When real side effects exist, this function leaves them unchanged.


active_base <- active_data |>
  dplyr::mutate(
    
    # Keep the dose exactly as it appeared in the source data.
    dose_original = dose,
    
    # PAROXETINE baseline dose is 20 mg.
    dose = dplyr::if_else(
      visit == 0,
      20,
      dose
    ),
    
    # Analysis-dose categories.
    dose = dplyr::case_when(
      dose == 10 ~ 20,
      dose == 60 ~ 50,
      TRUE ~ dose
    )
  )


if (SIMULATE_SIDE_EFFECTS_IF_MISSING) {
  
  active_base <- active_base |>
    # Ignore any side-effect values currently present in the source data.
    dplyr::mutate(
      side.effects = NA_real_
    ) |>
    simulate_side_effects_if_needed(
      seed = SIDE_EFFECT_SEED
    )
}
################################################################################
######################## ANALYSIS 1: ZERO IS A DOSE #############################
################################################################################

# Zero remains an observed treatment state.
# Sparse positive doses may still be grouped; zero is never grouped.
active_zero_as_dose <- active_base |>
  group_sparse_doses(min_n_per_dose = MIN_N_PER_DOSE)


################################################################################
################### DATASET 2: ZERO IS DISCONTINUATION ##########################
################################################################################

# A first post-baseline dose = 0 is interpreted as evidence that treatment
# was discontinued sometime after the previous visit and before this visit.
# Therefore the zero-dose row itself is not retained as an on-treatment outcome.
discontinuation_events <- active_base |>
  dplyr::filter(
    visit > 0,
    !is.na(dose),
    dose == 0
  ) |>
  dplyr::arrange(pid, visit) |>
  dplyr::group_by(pid) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup() |>
  dplyr::transmute(
    studyid,
    pid,
    discontinuation_day = analysis_day,
    discontinuation_visit = visit
  )

active_zero_as_discontinuation <- active_base |>
  dplyr::arrange(pid, visit) |>
  dplyr::group_by(pid) |>
  dplyr::mutate(
    
    # First observed post-baseline zero dose.
    zero_postbaseline =
      !is.na(dose) &
      dose == 0 &
      visit > 0,
    
    first_zero_visit = if (any(zero_postbaseline)) {
      min(visit[zero_postbaseline])
    } else {
      Inf
    },
    
    # The row immediately before the first zero is the last on-treatment row.
    next_visit = dplyr::lead(visit),
    
    censor_for_discontinuation =
      is.finite(first_zero_visit) &
      next_visit == first_zero_visit,
    
    # Keep only observations strictly before the first zero-dose visit.
    before_discontinuation =
      visit < first_zero_visit
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(before_discontinuation) |>
  group_sparse_doses(min_n_per_dose = MIN_N_PER_DOSE)


dose_levels_zero_as_dose <- sort(unique(
  active_zero_as_dose$dose[
    active_zero_as_dose$visit > 0 &
      !is.na(active_zero_as_dose$dose)
  ]
))

analysis_zero_as_dose <- prepare_dose_response_arm_data(
  data = active_zero_as_dose,
  visit_min = 0,
  visit_max = VISIT_MAX,
  max_followup_visit = VISIT_MAX,
  dose_levels = dose_levels_zero_as_dose,
  dose_history_levels = sort(unique(c(0, dose_levels_zero_as_dose))),
  include_zero_as_dose = TRUE,
  censoring_scenario = "zero_as_dose"
)

dose_levels_zero_as_discontinuation <- sort(unique(
  active_zero_as_discontinuation$dose[
    active_zero_as_discontinuation$visit > 0 &
      active_zero_as_discontinuation$dose > 0 &
      !is.na(active_zero_as_discontinuation$dose)
  ]
))

analysis_zero_as_discontinuation <- prepare_dose_response_arm_data(
  data = active_zero_as_discontinuation,
  visit_min = 0,
  visit_max = VISIT_MAX,
  max_followup_visit = VISIT_MAX,
  dose_levels = dose_levels_zero_as_discontinuation,
  dose_history_levels = sort(unique(c(
    0,
    dose_levels_zero_as_discontinuation
  ))),
  include_zero_as_dose = FALSE,
  censoring_scenario = "zero_as_discontinuation"
)
