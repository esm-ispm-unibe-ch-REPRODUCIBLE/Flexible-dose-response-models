#step-by-step per-arm master; dose tables and plots are lists by arm.
#######################################################################################################
############################ DR_03_master_stepwise_by_arm.R ###########################################
#######################################################################################################

# Run this script section-by-section in RStudio.
# Results are kept separately per active arm.

source("01_functions.R")
load_required_packages()

# If you have not already sourced the prep script, this will create analysis_arms.
if (!exists("analysis_arms")) {
  source("02_data_prep.R")
}

###################################### User settings #######################################
RESULTS_DIR <- "/Users/kchalkou/Desktop/Projects_ongoing/Dose-response/Results"
VISIT_DF <- 3
TRUNCATION <- c(0.01, 0.99)
INCLUDE_DOSE_TIME_INTERACTION <- TRUE
USE_IPCW <- TRUE

# NULL = use each arm's available grouped dose levels.
# Example: STRATEGY_DOSES <- c(20, 30, 40, 50)
STRATEGY_DOSES <- NULL

# NULL = each target-dose strategy starts at its own target dose.
# If all arms must start at 20 mg, use BASELINE_DOSE_FOR_PREDICTIONS <- 20
BASELINE_DOSE_FOR_PREDICTIONS <- NULL

RUN_MULTINOMIAL_SENSITIVITY <- TRUE
MULTINOMIAL_REF_DOSE <- NULL   # NULL chooses first available dose level as reference; use e.g. 20 if appropriate.

SAVE_PLOTS <- FALSE
PLOT_DIR <- "plots_by_arm"
PREDICTION_Y_LIMITS <- c(-10, 40)
PREDICTION_Y_BREAKS <- seq(-10, 40, by = 10)
################################################################################
######################## STEP 0: DESCRIPTIVE TABLES BY ARM ######################
################################################################################

# Separate tables by arm. Do not row-bind dose factors across arms.
dose_by_visit_tables <- purrr::imap(
  analysis_arms,
  function(dat, arm_name) {
    dat |>
      dplyr::filter(use_treatment_weight) |>
      dplyr::count(visit, dose_f, name = "n") |>
      tidyr::pivot_wider(
        names_from = dose_f,
        values_from = n,
        values_fill = 0
      ) |>
      dplyr::arrange(visit)
  }
)
dose_by_visit_tables

write_xlsx(
  dose_by_visit_tables,
  path = file.path(RESULTS_DIR, "dose_by_visit_tables.xlsx"))

# Separate dose-transition tables by arm.
dose_transition_tables <- purrr::imap(
  analysis_arms,
  function(dat, arm_name) {
    dat |>
      dplyr::filter(use_treatment_weight) |>
      dplyr::count(visit, dose_lag1_f, dose_f, name = "n") |>
      tidyr::pivot_wider(
        names_from = dose_f,
        values_from = n,
        values_fill = 0
      ) |>
      dplyr::arrange(visit, dose_lag1_f)
  }
)
dose_transition_tables


write_xlsx(
  dose_transition_tables,
  path = file.path(RESULTS_DIR,"dose_transition_tables.xlsx"))

# Separate sparse-dose grouping tables by arm.
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
dose_grouping_tables

write_xlsx(
  dose_grouping_tables,
  path = file.path(RESULTS_DIR,"dose_grouping_tables.xlsx"))

################################################################################
######################## STEP 0B: DESCRIPTIVE PLOTS BY ARM ######################
################################################################################

dose_by_visit_plots <- purrr::imap(
  analysis_arms,
  function(dat, arm_name) {
    plot_dat <- dat |>
      dplyr::filter(use_treatment_weight) |>
      dplyr::count(visit, dose_f, name = "n")

    ggplot2::ggplot(plot_dat, ggplot2::aes(x = visit, y = n, fill = dose_f)) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::scale_x_continuous(breaks = sort(unique(plot_dat$visit))) +
      ggplot2::labs(
        x = "Visit",
        y = "Number of patient-visits",
        fill = "Dose",
        title = paste0(arm_name)
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }
)
dose_by_visit_plots
dose_transition_plots <- purrr::imap(
  analysis_arms,
  function(dat, arm_name) {
    plot_dat <- dat |>
      dplyr::filter(use_treatment_weight) |>
      dplyr::count(visit, dose_lag1_f, dose_f, name = "n")

    ggplot2::ggplot(plot_dat, ggplot2::aes(x = dose_f, y = dose_lag1_f, fill = n)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = n), size = 3) +
      ggplot2::facet_wrap(~ visit) +
      ggplot2::labs(
        x = "Current dose",
        y = "Previous dose",
        fill = "n",
        title = paste0( arm_name)
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
  }
)
dose_transition_plots
dose_grouping_plots <- purrr::imap(
  analysis_arms,
  function(dat, arm_name) {
    if (!all(c("dose_before_grouping", "dose_grouped") %in% names(dat))) {
      return(NULL)
    }

    plot_dat <- dat |>
      dplyr::filter(visit > 0) |>
      dplyr::count(dose_before_grouping, dose_grouped, name = "n")

    ggplot2::ggplot(
      plot_dat,
      ggplot2::aes(x = factor(dose_before_grouping), y = n, fill = factor(dose_grouped))
    ) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::labs(
        x = "Original dose before grouping",
        y = "Number of patient-visits",
        fill = "Grouped dose",
        title = paste0(arm_name)
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }
)
dose_grouping_plots

tmp_grouping_plots <- purrr::compact(dose_grouping_plots)
if (length(tmp_grouping_plots) > 0) tmp_grouping_plots[[1]]
tmp_grouping_plots
dose_by_visit_grid <- save_plot_grid_jpg(
  plot_list = dose_by_visit_plots,
  file_path = file.path(RESULTS_DIR, "dose_by_visit_plots_grid.jpg"),
  ncol = 4,
  base_width = 5,
  base_height = 4
)

dose_transition_grid <- save_plot_grid_jpg(
  plot_list = dose_transition_plots,
  file_path = file.path(RESULTS_DIR, "dose_transition_plots_grid.jpg"),
  ncol = 4,
  base_width = 5.5,
  base_height = 5
)

dose_grouping_grid <- save_plot_grid_jpg(
  plot_list = dose_grouping_plots,
  file_path = file.path(RESULTS_DIR, "dose_grouping_plots_grid.jpg"),
  ncol = 4,
  base_width = 5,
  base_height = 4
)
################################################################################
######################## STEP 0C: CHECK WHICH ARMS ARE MODELLABLE ###############
################################################################################

MIN_PATIENTS_FOR_MODEL <- 30
MIN_TREATMENT_ROWS_FOR_MODEL <- 100
MIN_DOSE_LEVELS_FOR_MODEL <- 2

arm_modelability <- purrr::imap_dfr(
  analysis_arms,
  function(dat, arm_name) {
    
    model_dat <- dat %>%
      dplyr::filter(use_treatment_weight)
    
    dose_counts <- model_dat %>%
      dplyr::count(dose_f, name = "n")
    
    tibble::tibble(
      arm_name = arm_name,
      n_patients = dplyr::n_distinct(model_dat$pid),
      n_treatment_rows = nrow(model_dat),
      n_dose_levels = nrow(dose_counts),
      min_n_per_dose = min(dose_counts$n, na.rm = TRUE),
      dose_levels = paste(as.character(dose_counts$dose_f), collapse = ", "),
      modelable =
        dplyr::n_distinct(model_dat$pid) >= MIN_PATIENTS_FOR_MODEL &&
        nrow(model_dat) >= MIN_TREATMENT_ROWS_FOR_MODEL &&
        nrow(dose_counts) >= MIN_DOSE_LEVELS_FOR_MODEL
    )
  }
)

arm_modelability

write_xlsx(
  arm_modelability,
  path = file.path(RESULTS_DIR,"arm_modelability.xlsx"))

modelable_arm_names <- arm_modelability %>%
  dplyr::filter(modelable) %>%
  dplyr::pull(arm_name)
modelable_arm_names
descriptive_only_arm_names <- arm_modelability %>%
  dplyr::filter(!modelable) %>%
  dplyr::pull(arm_name)
descriptive_only_arm_names
analysis_arms_model <- analysis_arms[modelable_arm_names]
analysis_arms_descriptive_only <- analysis_arms[descriptive_only_arm_names]

################################################################################
######################## STEP 1A: IPTW DENOMINATOR BY ARM ######################
################################################################################

iptw_denominator_runs <- safe_imap(
  analysis_arms_model,
  function(dat, arm_name) {
    fit_iptw_denominator_model(data = dat, visit_df = VISIT_DF)
  }
)

iptw_denominator_status <- get_status_table(iptw_denominator_runs)
iptw_denominator_status

iptw_denominator_models <- get_success_results(iptw_denominator_runs)

iptw_denominator_coef_tabs <- purrr::map(
  iptw_denominator_models,
  make_polr_coef_table
)

names(iptw_denominator_coef_tabs)
if (length(iptw_denominator_coef_tabs) > 0) iptw_denominator_coef_tabs[[1]]
iptw_denominator_coef_tabs

coef_tabs_to_excel <- function(coef_tabs) {
  purrr::imap(
    coef_tabs,
    function(tab, arm_name) {
      
      tab_df <- as.data.frame(tab)
      
      if (!is.null(rownames(tab_df))) {
        tab_df <- tab_df |>
          tibble::rownames_to_column("term")
      }
      
      tab_df |>
        dplyr::mutate(
          arm_name = arm_name,
          .before = 1
        )
    }
  )
}

iptw_denominator_coef_tabs_excel <- coef_tabs_to_excel(
  iptw_denominator_coef_tabs
)

writexl::write_xlsx(
  iptw_denominator_coef_tabs_excel,
  path = file.path(RESULTS_DIR, "iptw_denominator_coef_tabs.xlsx")
)

iptw_denominator_coef_tabs_all <- purrr::imap_dfr(
  iptw_denominator_coef_tabs,
  function(tab, arm_name) {
    
    tab_df <- as.data.frame(tab)
    
    if (!is.null(rownames(tab_df))) {
      tab_df <- tab_df |>
        tibble::rownames_to_column("term")
    }
    
    tab_df |>
      dplyr::mutate(
        arm_name = arm_name,
        .before = 1
      )
  }
)

writexl::write_xlsx(
  list(
    "All_arms" = iptw_denominator_coef_tabs_all
  ),
  path = file.path(RESULTS_DIR, "iptw_denominator_coef_tabs_all.xlsx")
)
################################################################################
######################## STEP 1A-PLOT: IPTW DOSE MODEL PLOTS ###################
################################################################################

make_iptw_dose_model_plot_data <- function(model, data) {
  
  model_dat <- data |>
    dplyr::filter(use_treatment_weight) |>
    dplyr::mutate(.row_id = dplyr::row_number())
  
  prob_mat <- predict(
    model,
    newdata = model_dat,
    type = "probs"
  )
  
  if (is.null(dim(prob_mat))) {
    prob_names <- names(prob_mat)
    prob_mat <- matrix(prob_mat, nrow = 1)
    colnames(prob_mat) <- prob_names
  }
  
  dose_levels <- colnames(prob_mat)
  dose_levels_sorted <- as.character(
    sort(as.numeric(as.character(dose_levels)))
  )
  
  prob_df <- as.data.frame(prob_mat)
  prob_df$.row_id <- seq_len(nrow(prob_df))
  
  pred_long <- model_dat |>
    dplyr::select(.row_id, visit, observed_dose = dose_f) |>
    dplyr::left_join(prob_df, by = ".row_id") |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(dose_levels),
      names_to = "dose",
      values_to = "predicted_probability"
    ) |>
    dplyr::mutate(
      dose = factor(dose, levels = dose_levels_sorted)
    )
  
  predicted_summary <- pred_long |>
    dplyr::group_by(visit, dose) |>
    dplyr::summarise(
      mean_predicted_probability = mean(predicted_probability, na.rm = TRUE),
      .groups = "drop"
    )
  
  observed_counts <- model_dat |>
    dplyr::count(
      visit,
      dose = as.character(dose_f),
      name = "observed_n"
    )
  
  visit_totals <- model_dat |>
    dplyr::count(visit, name = "observed_total")
  
  plot_dat <- tidyr::expand_grid(
    visit = sort(unique(model_dat$visit)),
    dose = dose_levels_sorted
  ) |>
    dplyr::mutate(
      dose = factor(dose, levels = dose_levels_sorted)
    ) |>
    dplyr::left_join(
      predicted_summary,
      by = c("visit", "dose")
    ) |>
    dplyr::left_join(
      observed_counts |>
        dplyr::mutate(
          dose = factor(dose, levels = dose_levels_sorted)
        ),
      by = c("visit", "dose")
    ) |>
    dplyr::left_join(
      visit_totals,
      by = "visit"
    ) |>
    dplyr::mutate(
      observed_n = tidyr::replace_na(observed_n, 0L),
      observed_proportion = observed_n / observed_total
    )
  
  plot_dat
}


make_iptw_dose_model_plot <- function(plot_dat, arm_name) {
  
  ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(
      x = visit,
      group = dose,
      color = dose
    )
  ) +
    ggplot2::geom_line(
      ggplot2::aes(y = mean_predicted_probability),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = observed_proportion, shape = dose),
      size = 2,
      alpha = 0.80
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(plot_dat$visit))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    ggplot2::labs(
      x = "Visit",
      y = "Dose probability",
      color = "Dose",
      shape = "Observed dose",
      title = paste0(arm_name),
      subtitle = "Lines = mean model-predicted dose probabilities; points = observed dose proportions"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "right"
    )
}

iptw_dose_model_plot_data <- purrr::imap(
  iptw_denominator_models,
  function(model, arm_name) {
    make_iptw_dose_model_plot_data(
      model = model,
      data = analysis_arms_model[[arm_name]]
    )
  }
)

iptw_dose_model_plots <- purrr::imap(
  iptw_dose_model_plot_data,
  function(plot_dat, arm_name) {
    make_iptw_dose_model_plot(
      plot_dat = plot_dat,
      arm_name = arm_name
    )
  }
)

iptw_dose_model_plots

iptw_dose_model_grid <- save_plot_grid_jpg(
  plot_list = iptw_dose_model_plots,
  file_path = file.path(RESULTS_DIR, "iptw_dose_model_plots_grid.jpg"),
  ncol = 4,
  base_width = 5.5,
  base_height = 4.5
)

iptw_dose_model_grid
################################################################################
################ STEP 1A-PLOT: DOSE-ASSIGNMENT MODEL PROFILES ##################
################################################################################
DOSE_MODEL_PROFILE_PREVIOUS_DOSE_BY_ARM <- list(
  "29060_002_PAROXETINE" = 20,
  "29060_003_PAROXETINE" = 20,
  "29060_003_IMIPRAMINE" = 65
)

DOSE_MODEL_PROFILE_SIDE_EFFECTS <- 4
DOSE_MODEL_PROFILE_VISITS <- 1:8


get_dose_model_profile_previous_dose <- function(data, arm_name) {
  
  model_dat <- data |>
    dplyr::filter(use_treatment_weight)
  
  available_previous_doses <- levels(model_dat$dose_lag1_f) |>
    as.numeric()
  
  available_previous_doses <- sort(
    unique(
      available_previous_doses[
        !is.na(available_previous_doses) &
          available_previous_doses > 0
      ]
    )
  )
  
  if (
    exists("DOSE_MODEL_PROFILE_PREVIOUS_DOSE_BY_ARM") &&
    !is.null(DOSE_MODEL_PROFILE_PREVIOUS_DOSE_BY_ARM) &&
    arm_name %in% names(DOSE_MODEL_PROFILE_PREVIOUS_DOSE_BY_ARM)
  ) {
    previous_dose <- DOSE_MODEL_PROFILE_PREVIOUS_DOSE_BY_ARM[[arm_name]]
  } else {
    previous_dose <- min(available_previous_doses, na.rm = TRUE)
  }
  
  if (!(previous_dose %in% available_previous_doses)) {
    stop(
      "Requested previous dose ", previous_dose,
      " is not available for arm ", arm_name,
      ". Available previous doses are: ",
      paste(available_previous_doses, collapse = ", ")
    )
  }
  
  previous_dose
}
make_dose_model_profile_data <- function(
    model,
    data,
    previous_dose = 20,
    side_effect_value = 4,
    visits = 1:8
) {
  
  model_dat <- data |>
    dplyr::filter(use_treatment_weight)
  
  # If requested previous dose is not available in this arm, use the lowest
  # positive previous-dose level available in that arm.
  dose_lag1_levels <- levels(model_dat$dose_lag1_f)
  dose_lag1_numeric <- as.numeric(as.character(dose_lag1_levels))
  
  available_previous_doses <- dose_lag1_numeric[
    !is.na(dose_lag1_numeric) & dose_lag1_numeric > 0
  ]
  
  if (as.character(previous_dose) %in% dose_lag1_levels) {
    previous_dose_used <- previous_dose
  } else {
    previous_dose_used <- min(available_previous_doses, na.rm = TRUE)
  }
  
  mean_outcome_0 <- mean(model_dat$outcome_0, na.rm = TRUE)
  
  # Define clinical profiles using observed visit-specific quantiles.
  # Improving profile = higher observed HAMD improvement at that visit.
  # Worsening profile = lower observed HAMD improvement at that visit.
  profile_values <- model_dat |>
    dplyr::filter(visit %in% visits) |>
    dplyr::group_by(visit) |>
    dplyr::summarise(
      improving_delta_outcome = stats::quantile(
        delta_outcome,
        probs = 0.75,
        na.rm = TRUE,
        names = FALSE
      ),
      worsening_delta_outcome = stats::quantile(
        delta_outcome,
        probs = 0.25,
        na.rm = TRUE,
        names = FALSE
      ),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(
      cols = c(improving_delta_outcome, worsening_delta_outcome),
      names_to = "profile",
      values_to = "delta_outcome"
    ) |>
    dplyr::mutate(
      profile = dplyr::case_when(
        profile == "improving_delta_outcome" ~ "Improving profile",
        profile == "worsening_delta_outcome" ~ "Worsening profile",
        TRUE ~ profile
      )
    )
  
  newdat <- profile_values |>
    dplyr::mutate(
      outcome_0 = mean_outcome_0,
      side.effects = side_effect_value,
      dose_lag1_f = factor(
        as.character(previous_dose_used),
        levels = levels(model_dat$dose_lag1_f)
      )
    )
  
  prob_mat <- predict(
    model,
    newdata = newdat,
    type = "probs"
  )
  
  if (is.null(dim(prob_mat))) {
    prob_names <- names(prob_mat)
    prob_mat <- matrix(prob_mat, nrow = 1)
    colnames(prob_mat) <- prob_names
  }
  
  dose_levels <- colnames(prob_mat)
  dose_levels_sorted <- as.character(
    sort(as.numeric(as.character(dose_levels)))
  )
  
  prob_df <- as.data.frame(prob_mat)
  names(prob_df) <- dose_levels
  prob_df$.row_id <- seq_len(nrow(prob_df))
  
  plot_dat <- newdat |>
    dplyr::mutate(.row_id = dplyr::row_number()) |>
    dplyr::left_join(prob_df, by = ".row_id") |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(dose_levels),
      names_to = "current_dose",
      values_to = "predicted_probability"
    ) |>
    dplyr::mutate(
      current_dose = factor(
        current_dose,
        levels = dose_levels_sorted
      ),
      profile = factor(
        profile,
        levels = c("Improving profile", "Worsening profile")
      ),
      previous_dose_used = previous_dose_used,
      side_effect_value = side_effect_value,
      mean_outcome_0 = mean_outcome_0
    ) |>
    dplyr::arrange(profile, current_dose, visit)
  
  plot_dat
}


make_dose_model_profile_plot <- function(plot_dat, arm_name) {
  
  previous_dose_used <- unique(plot_dat$previous_dose_used)[1]
  side_effect_value <- unique(plot_dat$side_effect_value)[1]
  mean_outcome_0 <- unique(plot_dat$mean_outcome_0)[1]
  
  ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(
      x = visit,
      y = predicted_probability,
      color = current_dose,
      group = current_dose
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(
      ~ profile,
      ncol = 1
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(plot_dat$visit))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      x = "Visit",
      y = "Predicted probability of assigned dose",
      color = "Current dose",
      ,
      subtitle = paste0(
        arm_name
      )
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_line(color = "grey92", linewidth = 0.25),
      panel.grid.major = ggplot2::element_line(color = "grey88", linewidth = 0.35),
      legend.position = "right"
    )
}

iptw_dose_model_profile_data <- purrr::imap(
  iptw_denominator_models,
  function(model, arm_name) {
    
    previous_dose_this_arm <- get_dose_model_profile_previous_dose(
      data = analysis_arms_model[[arm_name]],
      arm_name = arm_name
    )
    
    make_dose_model_profile_data(
      model = model,
      data = analysis_arms_model[[arm_name]],
      previous_dose = previous_dose_this_arm,
      side_effect_value = DOSE_MODEL_PROFILE_SIDE_EFFECTS,
      visits = DOSE_MODEL_PROFILE_VISITS
    )
  }
)

iptw_dose_model_profile_plots <- purrr::imap(
  iptw_dose_model_profile_data,
  function(plot_dat, arm_name) {
    make_dose_model_profile_plot(
      plot_dat = plot_dat,
      arm_name = arm_name
    )
  }
)


iptw_dose_model_profile_plots

iptw_dose_model_profile_grid <- save_plot_grid_jpg(
  plot_list = iptw_dose_model_profile_plots,
  file_path = file.path(RESULTS_DIR, "iptw_dose_model_profile_plots_grid.jpg"),
  ncol = 4,
  base_width = 5.5,
  base_height = 6
)

iptw_dose_model_profile_grid
################################################################################
######################## STEP 1B: IPTW NUMERATOR BY ARM ########################
################################################################################

iptw_numerator_runs <- safe_imap(
  analysis_arms_model,
  function(dat, arm_name) {
    fit_iptw_numerator_model(data = dat, visit_df = VISIT_DF)
  }
)

iptw_numerator_status <- get_status_table(iptw_numerator_runs)
iptw_numerator_status

iptw_numerator_models <- get_success_results(iptw_numerator_runs)

iptw_numerator_coef_tabs <- purrr::map(
  iptw_numerator_models,
  make_polr_coef_table
)

names(iptw_numerator_coef_tabs)
if (length(iptw_numerator_coef_tabs) > 0) iptw_numerator_coef_tabs[[1]]
iptw_numerator_coef_tabs


# If not already loaded
if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

# Convert each coefficient table to a data frame
iptw_numerator_coef_tabs_excel <- purrr::imap(
  iptw_numerator_coef_tabs,
  function(tab, arm_name) {
    
    tab_df <- as.data.frame(tab)
    
    if (!is.null(rownames(tab_df))) {
      tab_df <- tab_df |>
        tibble::rownames_to_column("term")
    }
    
    tab_df |>
      dplyr::mutate(
        arm_name = arm_name,
        .before = 1
      )
  }
)

# Save one sheet per arm
writexl::write_xlsx(
  iptw_numerator_coef_tabs_excel,
  path = file.path(RESULTS_DIR, "iptw_numerator_coef_tabs.xlsx")
)

iptw_numerator_coef_tabs_all <- purrr::imap_dfr(
  iptw_numerator_coef_tabs,
  function(tab, arm_name) {
    
    tab_df <- as.data.frame(tab)
    
    if (!is.null(rownames(tab_df))) {
      tab_df <- tab_df |>
        tibble::rownames_to_column("term")
    }
    
    tab_df |>
      dplyr::mutate(
        arm_name = arm_name,
        .before = 1
      )
  }
)

writexl::write_xlsx(
  list(
    "All_arms" = iptw_numerator_coef_tabs_all
  ),
  path = file.path(RESULTS_DIR, "iptw_numerator_coef_tabs_all.xlsx")
)
################################################################################
######################## STEP 1C: ADD IPTW WEIGHTS BY ARM ######################
################################################################################
iptw_common_arms <- Reduce(
  intersect,
  list(
    names(analysis_arms_model),
    names(iptw_denominator_models),
    names(iptw_numerator_models)
  )
)

iptw_common_arms

iptw_weight_runs <- safe_imap(
  analysis_arms_model[iptw_common_arms],
  function(dat, arm_name) {
    add_iptw_treatment_weights(
      data = dat,
      denominator_model = iptw_denominator_models[[arm_name]],
      numerator_model = iptw_numerator_models[[arm_name]]
    )
  }
)

iptw_weight_status <- get_status_table(iptw_weight_runs)
iptw_weight_status

analysis_arms_iptw <- get_success_results(iptw_weight_runs)

iptw_weight_summaries <- purrr::map(
  analysis_arms_iptw,
  summarise_iptw_treatment_weights
)

names(iptw_weight_summaries)
iptw_weight_summaries
analysis_arms_iptw <- get_success_results(iptw_weight_runs)

iptw_weight_summaries <- purrr::map(
  analysis_arms_iptw,
  summarise_iptw_treatment_weights
)

names(iptw_weight_summaries)
if (length(iptw_weight_summaries) > 0) iptw_weight_summaries[[1]]
iptw_weight_summaries


write_xlsx(
  iptw_weight_summaries,
  path = file.path(RESULTS_DIR, "iptw_weight_summaries.xlsx"))

################################################################################
######################## STEP 2: IPCW BY ARM ###################################
################################################################################

if (isTRUE(USE_IPCW)) {

  ipcw_denominator_runs <- safe_imap(
    analysis_arms_iptw,
    function(dat, arm_name) {
      fit_ipcw_denominator_model(data = dat)
    }
  )

  ipcw_denominator_status <- get_status_table(ipcw_denominator_runs)
  ipcw_denominator_status

  ipcw_denominator_models <- get_success_results(ipcw_denominator_runs)

  ipcw_denominator_coef_tabs <- purrr::map(
    ipcw_denominator_models,
    make_glm_coef_table
  )
ipcw_denominator_coef_tabs

################################################################################
######################## SAVE IPCW DENOMINATOR COEFFICIENT TABLES ##############
################################################################################

if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

ipcw_denominator_coef_tabs_excel <- purrr::imap(
  ipcw_denominator_coef_tabs,
  function(tab, arm_name) {
    
    tab_df <- as.data.frame(tab)
    
    if (!is.null(rownames(tab_df))) {
      tab_df <- tab_df |>
        tibble::rownames_to_column("term")
    }
    
    tab_df |>
      dplyr::mutate(
        arm_name = arm_name,
        .before = 1
      )
  }
)

writexl::write_xlsx(
  ipcw_denominator_coef_tabs_excel,
  path = file.path(RESULTS_DIR, "ipcw_denominator_coef_tabs.xlsx")
)

ipcw_denominator_coef_tabs_all <- purrr::imap_dfr(
  ipcw_denominator_coef_tabs,
  function(tab, arm_name) {
    
    tab_df <- as.data.frame(tab)
    
    if (!is.null(rownames(tab_df))) {
      tab_df <- tab_df |>
        tibble::rownames_to_column("term")
    }
    
    tab_df |>
      dplyr::mutate(
        arm_name = arm_name,
        .before = 1
      )
  }
)

writexl::write_xlsx(
  list(
    "All_arms" = ipcw_denominator_coef_tabs_all
  ),
  path = file.path(RESULTS_DIR, "ipcw_denominator_coef_tabs_all.xlsx")
)
  
  ipcw_numerator_runs <- safe_imap(
    analysis_arms_iptw,
    function(dat, arm_name) {
      fit_ipcw_numerator_model(data = dat)
    }
  )

  ipcw_numerator_status <- get_status_table(ipcw_numerator_runs)
  ipcw_numerator_status

  ipcw_numerator_models <- get_success_results(ipcw_numerator_runs)

  ipcw_numerator_coef_tabs <- purrr::map(
    ipcw_numerator_models,
    make_glm_coef_table
  )
ipcw_numerator_coef_tabs

################################################################################
######################## SAVE IPCW NUMERATOR COEFFICIENT TABLES ################
################################################################################

if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

ipcw_numerator_coef_tabs_excel <- purrr::imap(
  ipcw_numerator_coef_tabs,
  function(tab, arm_name) {
    
    tab_df <- as.data.frame(tab)
    
    if (!is.null(rownames(tab_df))) {
      tab_df <- tab_df |>
        tibble::rownames_to_column("term")
    }
    
    tab_df |>
      dplyr::mutate(
        arm_name = arm_name,
        .before = 1
      )
  }
)

writexl::write_xlsx(
  ipcw_numerator_coef_tabs_excel,
  path = file.path(RESULTS_DIR, "ipcw_numerator_coef_tabs.xlsx")
)

ipcw_numerator_coef_tabs_all <- purrr::imap_dfr(
  ipcw_numerator_coef_tabs,
  function(tab, arm_name) {
    
    tab_df <- as.data.frame(tab)
    
    if (!is.null(rownames(tab_df))) {
      tab_df <- tab_df |>
        tibble::rownames_to_column("term")
    }
    
    tab_df |>
      dplyr::mutate(
        arm_name = arm_name,
        .before = 1
      )
  }
)

writexl::write_xlsx(
  list(
    "All_arms" = ipcw_numerator_coef_tabs_all
  ),
  path = file.path(RESULTS_DIR, "ipcw_numerator_coef_tabs_all.xlsx")
)
  ipcw_common_arms <- Reduce(intersect, list(
    names(analysis_arms_iptw),
    names(ipcw_denominator_models),
    names(ipcw_numerator_models)
  ))

  ipcw_weight_runs <- safe_imap(
    analysis_arms_iptw[ipcw_common_arms],
    function(dat, arm_name) {
      add_ipcw_censoring_weights(
        data = dat,
        denominator_model = ipcw_denominator_models[[arm_name]],
        numerator_model = ipcw_numerator_models[[arm_name]]
      )
    }
  )

  ipcw_weight_status <- get_status_table(ipcw_weight_runs)
  ipcw_weight_status

  analysis_arms_ipdcw <- get_success_results(ipcw_weight_runs)

} else {
  ipcw_denominator_status <- tibble::tibble()
  ipcw_numerator_status <- tibble::tibble()
  ipcw_weight_status <- tibble::tibble(note = "IPCW not used; censoring weights set to 1.")

  analysis_arms_ipdcw <- purrr::map(
    analysis_arms_iptw,
    add_no_censoring_weights
  )
}

ipcw_weight_summaries <- purrr::map(
  analysis_arms_ipdcw,
  summarise_ipcw_censoring_weights
)

names(ipcw_weight_summaries)
if (length(ipcw_weight_summaries) > 0) ipcw_weight_summaries[[1]]
ipcw_weight_summaries

write_xlsx(
  ipcw_weight_summaries,
  path = file.path(RESULTS_DIR, "ipcw_weight_summaries.xlsx"))

################################################################################
######################## STEP 3: TOTAL WEIGHTS BY ARM ##########################
################################################################################

total_weight_runs <- safe_imap(
  analysis_arms_ipdcw,
  function(dat, arm_name) {
    dat |>
      add_total_weights() |>
      truncate_total_weights(lower = TRUNCATION[1], upper = TRUNCATION[2])
  }
)

total_weight_status <- get_status_table(total_weight_runs)
total_weight_status

analysis_arms_weighted <- get_success_results(total_weight_runs)

total_weight_summaries <- purrr::map(
  analysis_arms_weighted,
  summarise_total_weights
)
total_weight_summaries

write_xlsx(
  total_weight_summaries,
  path = file.path(RESULTS_DIR, "total_weight_summaries.xlsx"))

total_truncated_weight_summaries <- purrr::map(
  analysis_arms_weighted,
  summarise_total_truncated_weights
)

total_weight_truncation_checks <- purrr::map(
  analysis_arms_weighted,
  check_truncation
)
total_truncated_weight_summaries
names(total_truncated_weight_summaries)
if (length(total_truncated_weight_summaries) > 0) total_truncated_weight_summaries[[1]]
if (length(total_weight_truncation_checks) > 0) total_weight_truncation_checks[[1]]

write_xlsx(
  total_truncated_weight_summaries,
  path = file.path(RESULTS_DIR, "total_truncated_weight_summaries.xlsx"))

################################################################################
######################## STEP 4: MSM BY ARM ####################################
################################################################################

msm_runs <- safe_imap(
  analysis_arms_weighted,
  function(dat, arm_name) {
    fit_weighted_msm(
      data = dat,
      weight_var = "SW_total_trunc",
      visit_df = VISIT_DF,
      corstr = "independence",
      include_dose_time_interaction = INCLUDE_DOSE_TIME_INTERACTION
    )
  }
)

msm_status <- get_status_table(msm_runs)
msm_status

msm_models <- get_success_results(msm_runs)

msm_coef_tabs <- purrr::map(msm_models, make_gee_coef_table)
msm_robust_naive_tabs <- purrr::map(msm_models, make_gee_coef_table_robust_naive)

names(msm_coef_tabs)
if (length(msm_coef_tabs) > 0) msm_coef_tabs[[1]]
if (length(msm_robust_naive_tabs) > 0) msm_robust_naive_tabs[[1]]
msm_coef_tabs

################################################################################
######################## SAVE MSM COEFFICIENT TABLES ############################
################################################################################

if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

msm_coef_tabs_excel <- purrr::imap(
  msm_coef_tabs,
  function(tab, arm_name) {
    
    tab_df <- as.data.frame(tab)
    
    if (!is.null(rownames(tab_df))) {
      tab_df <- tab_df |>
        tibble::rownames_to_column("term")
    }
    
    tab_df |>
      dplyr::mutate(
        arm_name = arm_name,
        .before = 1
      )
  }
)

writexl::write_xlsx(
  msm_coef_tabs_excel,
  path = file.path(RESULTS_DIR, "msm_coef_tabs.xlsx")
)

msm_coef_tabs_all <- purrr::imap_dfr(
  msm_coef_tabs,
  function(tab, arm_name) {
    
    tab_df <- as.data.frame(tab)
    
    if (!is.null(rownames(tab_df))) {
      tab_df <- tab_df |>
        tibble::rownames_to_column("term")
    }
    
    tab_df |>
      dplyr::mutate(
        arm_name = arm_name,
        .before = 1
      )
  }
)

writexl::write_xlsx(
  list(
    "All_arms" = msm_coef_tabs_all
  ),
  path = file.path(RESULTS_DIR, "msm_coef_tabs_all.xlsx")
)
################################################################################
######################## STEP 5: PREDICTIONS AND PLOTS BY ARM ##################
################################################################################

MSM_MAX_VISIT <- 9

strip_ansi <- function(x) {
  gsub("\\033\\[[0-9;]*m", "", x)
}

get_prediction_dose_history_levels <- function(model) {
  model_dat <- attr(model, "model_data")
  
  factor_vars <- c("dose_lag1_f", "dose_lag2_f", "dose_lag3_f")
  
  dose_history_levels <- unique(unlist(
    lapply(factor_vars, function(v) {
      if (v %in% names(model_dat)) {
        as.numeric(as.character(levels(model_dat[[v]])))
      } else {
        numeric(0)
      }
    })
  ))
  
  dose_history_levels <- dose_history_levels[!is.na(dose_history_levels)]
  sort(dose_history_levels)
}

get_strategy_doses_for_arm <- function(model) {
  model_dat <- attr(model, "model_data")
  
  if ("dose_f" %in% names(model_dat)) {
    available_doses <- model_dat$dose_f |>
      droplevels() |>
      levels() |>
      as.numeric() |>
      sort()
  } else {
    available_doses <- model_dat$dose_lag1_f |>
      droplevels() |>
      levels() |>
      as.numeric() |>
      sort()
  }
  
  available_doses <- available_doses[!is.na(available_doses) & available_doses > 0]
  
  if (!is.null(STRATEGY_DOSES)) {
    available_doses <- intersect(available_doses, STRATEGY_DOSES)
  }
  
  available_doses
}

predict_gee_ci <- function(model, newdata) {
  beta <- coef(model)
  
  V <- model$geese$vbeta
  dimnames(V) <- list(names(beta), names(beta))
  
  Terms <- stats::delete.response(stats::terms(model))
  
  needed_vars <- intersect(all.vars(Terms), names(newdata))
  
  bad_rows <- !stats::complete.cases(newdata[, needed_vars, drop = FALSE])
  
  if (any(bad_rows)) {
    bad_preview <- newdata[bad_rows, needed_vars, drop = FALSE]
    
    stop(
      "Prediction data contain missing values in model variables. ",
      "This usually means a prediction dose is not supported by the fitted arm-specific factor levels. ",
      "Rows with missing values: ",
      paste(which(bad_rows), collapse = ", "),
      "\nPreview:\n",
      paste(capture.output(print(utils::head(bad_preview, 10))), collapse = "\n")
    )
  }
  
  X <- stats::model.matrix(
    Terms,
    newdata,
    contrasts.arg = model$contrasts,
    xlev = model$xlevels
  )
  
  missing_cols <- setdiff(names(beta), colnames(X))
  
  if (length(missing_cols) > 0) {
    X <- cbind(
      X,
      matrix(
        0,
        nrow = nrow(X),
        ncol = length(missing_cols),
        dimnames = list(NULL, missing_cols)
      )
    )
  }
  
  X <- X[, names(beta), drop = FALSE]
  
  if (nrow(X) != nrow(newdata)) {
    stop(
      "Model matrix has ", nrow(X), " rows but prediction data have ",
      nrow(newdata), " rows. This usually means model.matrix dropped rows."
    )
  }
  
  pred <- as.numeric(X %*% beta)
  se <- sqrt(diag(X %*% V %*% t(X)))
  
  newdata |>
    dplyr::mutate(
      pred_delta_outcome = pred,
      se = se,
      lower_95 = pred_delta_outcome - qnorm(0.975) * se,
      upper_95 = pred_delta_outcome + qnorm(0.975) * se
    )
}

make_prediction_results_one_arm <- function(model, arm_name) {
  model_dat <- attr(model, "model_data")
  trial_this_arm <- unique(as.character(model_dat$trial_name))[1]
  visits <- 1:MSM_MAX_VISIT
  mean_outcome_0 <- mean(model_dat$outcome_0, na.rm = TRUE)
  
  strategy_doses <- get_strategy_doses_for_arm(model)
  
  if (length(strategy_doses) < 2) {
    stop("Fewer than 2 strategy doses available for ", arm_name)
  }
  
  dose_history_levels <- get_prediction_dose_history_levels(model)
  
  strategy_dat <- dplyr::bind_rows(
    lapply(strategy_doses, function(dose_value) {
      
      baseline_dose_this_arm <- if (is.null(BASELINE_DOSE_FOR_PREDICTIONS)) {
        dose_value
      } else if (BASELINE_DOSE_FOR_PREDICTIONS %in% dose_history_levels) {
        BASELINE_DOSE_FOR_PREDICTIONS
      } else {
        dose_value
      }
      
      make_strategy_data(
        strategy_dose = dose_value,
        visits = visits,
        outcome_0_value = mean_outcome_0,
        baseline_dose = baseline_dose_this_arm,
        dose_history_levels = dose_history_levels
      )
    })
  )
  
  factor_vars <- c("dose_lag1_f", "dose_lag2_f", "dose_lag3_f")
  
  for (v in factor_vars) {
    if (v %in% names(strategy_dat) && v %in% names(model_dat)) {
      strategy_dat[[v]] <- factor(
        as.character(strategy_dat[[v]]),
        levels = levels(model_dat[[v]])
      )
    }
  }
  
  # Explicit check before prediction
  factor_na_check <- strategy_dat |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(factor_vars),
        ~ sum(is.na(.x))
      )
    )
  
  if (any(factor_na_check > 0)) {
    stop(
      "Some prediction dose-history factors became NA for arm ",
      arm_name,
      ". This means the requested prediction strategy is not supported by the fitted model levels.\n",
      paste(capture.output(print(factor_na_check)), collapse = "\n")
    )
  }
  
  strategy_dat <- predict_gee_ci(
    model = model,
    newdata = strategy_dat
  )
  
  strategy_dat <- strategy_dat |>
    dplyr::mutate(
      arm_name = arm_name,
      trial_name = trial_this_arm
    )
  
  placebo_dat <- placebo_observed_by_trial |>
    dplyr::filter(
      trial_name == trial_this_arm,
      visit %in% visits
    )
  strategy_dat <- strategy_dat %>%
    dplyr::mutate(
      strategy_dose_num = as.numeric(gsub("[^0-9.]", "", as.character(strategy)))
    ) %>%
    dplyr::arrange(strategy_dose_num, visit) %>%
    dplyr::mutate(
      strategy = factor(
        strategy,
        levels = unique(strategy[order(strategy_dose_num)])
      )
    )
  n_dat <- model_dat |>
    dplyr::filter(visit %in% visits) |>
    dplyr::group_by(visit) |>
    dplyr::summarise(
      active_n_patients = dplyr::n_distinct(pid),
      active_n_visit_rows = dplyr::n(),
      .groups = "drop"
    )
  
  y_range <- range(
    c(
      strategy_dat$lower_95,
      strategy_dat$upper_95,
      placebo_dat$placebo_lower_95,
      placebo_dat$placebo_upper_95
    ),
    na.rm = TRUE
  )
  
  y_range <- PREDICTION_Y_LIMITS
  
  # Put the active-arm n just below the visible y-axis.
  # It will still show because coord_cartesian uses clip = "off".
  y_n <- -8.5
  
  strategy_plot <- ggplot2::ggplot(
    strategy_dat,
    ggplot2::aes(
      x = visit,
      y = pred_delta_outcome,
      group = strategy,
      color = strategy,
      linetype = strategy,
      shape = strategy
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4,
      color = "grey60"
    ) +
    
    # Observed placebo uncertainty band
    ggplot2::geom_ribbon(
      data = placebo_dat,
      ggplot2::aes(
        x = visit,
        ymin = placebo_lower_95,
        ymax = placebo_upper_95
      ),
      inherit.aes = FALSE,
      fill = "grey70",
      alpha = 0.30
    ) +
    
    # Active-dose model-based prediction uncertainty bands
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = lower_95,
        ymax = upper_95,
        fill = strategy
      ),
      alpha = 0.12,
      color = NA,
      show.legend = FALSE
    ) +
    
    # Active-dose model-based prediction lines
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    
    # Observed placebo line
    ggplot2::geom_line(
      data = placebo_dat,
      ggplot2::aes(
        x = visit,
        y = placebo_mean_delta
      ),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 1,
      linetype = "longdash"
    ) +
    ggplot2::geom_point(
      data = placebo_dat,
      ggplot2::aes(
        x = visit,
        y = placebo_mean_delta
      ),
      inherit.aes = FALSE,
      color = "black",
      shape = 17,
      size = 2.5
    ) +
   ggplot2::geom_text(
      data = n_dat,
      ggplot2::aes(x = visit, y = y_n, label = active_n_patients),
      inherit.aes = FALSE,
      size = 3
    ) +
    ggplot2::annotate(
      "text",
      x = min(n_dat$visit) - 0.35,
      y = y_n,
      label = "n=",
      hjust = 1,
      size = 3
    ) +
    ggplot2::scale_x_continuous(breaks = visits) +
    ggplot2::scale_y_continuous(
      breaks = PREDICTION_Y_BREAKS
    ) +
    ggplot2::labs(
      x = "Weeks from randomization",
      y = "Predicted HAMD improvement from baseline",
      color = "Dose strategy",
      linetype = "Dose strategy",
      shape = "Dose strategy",
      title = paste0(arm_name)
    ) +
    ggplot2::coord_cartesian(
      ylim = PREDICTION_Y_LIMITS,
      clip = "off"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(10, 10, 45, 45),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = c(0.03, 0.97),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(fill = "white", colour = "black"),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA)
    )
  
  list(
    arm_name = arm_name,
    strategy_doses = strategy_doses,
    dose_history_levels = dose_history_levels,
    strategy_dat = strategy_dat,
    n_dat = n_dat,
    plot = strategy_plot
  )
}

prediction_runs <- safe_imap(
  msm_models,
  function(model, arm_name) {
    make_prediction_results_one_arm(
      model = model,
      arm_name = arm_name
    )
  }
)

prediction_status <- get_status_table(prediction_runs)
prediction_status

# If any prediction fails, print full readable error
if (any(!prediction_status$ok)) {
  failed_prediction_arms <- prediction_status |>
    dplyr::filter(!ok) |>
    dplyr::pull(arm_name)
  
  for (arm_name in failed_prediction_arms) {
    cat("\n\n==============================\n")
    cat("Prediction error:", arm_name, "\n")
    cat("==============================\n")
    cat(strip_ansi(prediction_runs[[arm_name]]$error), "\n")
  }
}

prediction_results <- get_success_results(prediction_runs)

strategy_prediction_data <- purrr::map(
  prediction_results,
  "strategy_dat"
)

strategy_prediction_plots <- purrr::map(
  prediction_results,
  "plot"
)

strategy_prediction_plots

strategy_prediction_grid <- save_plot_grid_jpg(
  plot_list = strategy_prediction_plots,
  file_path = file.path(RESULTS_DIR, "strategy_prediction_plots_grid.jpg"),
  ncol = 4,
  base_width = 6,
  base_height = 5
)

strategy_prediction_grid

################################################################################
######################## STEP 6: PAIRWISE CONTRASTS BY ARM #####################
################################################################################

pairwise_runs <- safe_imap(
  msm_models,
  function(model, arm_name) {
    strategy_doses <- get_strategy_doses_for_arm(model)

    if (length(strategy_doses) < 2) {
      stop("Fewer than 2 strategy doses available for ", arm_name)
    }

    estimate_all_pairwise_strategy_contrasts(
      model = model,
      target_visit = max(attr(model, "model_data")$visit, na.rm = TRUE),
      strategy_doses = strategy_doses,
      baseline_dose = BASELINE_DOSE_FOR_PREDICTIONS
    )
  }
)

pairwise_status <- get_status_table(pairwise_runs)
pairwise_status

pairwise_effects <- get_success_results(pairwise_runs)

pairwise_effects_all <- purrr::imap_dfr(
  pairwise_effects,
  function(tab, arm_name) {
    tab |> dplyr::mutate(arm_name = arm_name, .before = 1)
  }
)

pairwise_effects_all

write_xlsx(
  pairwise_effects_all,
  path = file.path(RESULTS_DIR, "pairwise_effects_all.xlsx"))

################################################################################
################ FINAL TABLE: DOSE STRATEGIES VERSUS OBSERVED PLACEBO ##########
################################################################################
make_dose_vs_placebo_table <- function(
    prediction_results,
    placebo_observed_by_trial,
    placebo_counts_by_trial = NULL,
    max_visit = 8,
    weight_model = "Ordinal IPTW"
) {
  
  if (is.null(placebo_counts_by_trial)) {
    placebo_counts_by_trial <- placebo_observed_by_trial |>
      dplyr::transmute(
        trial_name,
        visit,
        placebo_n_patients = n_placebo,
        placebo_n_visit_rows = n_placebo
      )
  }
  
  active_counts_by_arm <- purrr::imap_dfr(
    prediction_results,
    function(res, arm_name) {
      
      n_dat <- res$n_dat
      
      if ("n" %in% names(n_dat) && !"active_n_patients" %in% names(n_dat)) {
        n_dat <- n_dat |>
          dplyr::mutate(active_n_patients = n)
      }
      
      if (!"active_n_visit_rows" %in% names(n_dat)) {
        n_dat <- n_dat |>
          dplyr::mutate(active_n_visit_rows = active_n_patients)
      }
      
      n_dat |>
        dplyr::transmute(
          arm_name = arm_name,
          visit,
          active_n_patients,
          active_n_visit_rows
        )
    }
  )
  
  placebo_for_join <- placebo_observed_by_trial |>
    dplyr::select(
      trial_name,
      visit,
      n_placebo,
      placebo_mean_delta,
      placebo_se_delta,
      placebo_lower_95,
      placebo_upper_95
    ) |>
    dplyr::left_join(
      placebo_counts_by_trial,
      by = c("trial_name", "visit")
    ) |>
    dplyr::mutate(
      placebo_n_patients = dplyr::coalesce(placebo_n_patients, n_placebo),
      placebo_n_visit_rows = dplyr::coalesce(placebo_n_visit_rows, n_placebo)
    )
  
  out <- purrr::imap_dfr(
    prediction_results,
    function(res, arm_name) {
      
      pred_dat <- res$strategy_dat |>
        dplyr::filter(
          visit >= 1,
          visit <= max_visit
        )
      
      if (!"trial_name" %in% names(pred_dat)) {
        pred_dat <- pred_dat |>
          dplyr::mutate(
            trial_name = sub("_[^_]+$", "", arm_name)
          )
      }
      
      if (!"strategy_dose" %in% names(pred_dat)) {
        if ("strategy_dose_num" %in% names(pred_dat)) {
          pred_dat <- pred_dat |>
            dplyr::mutate(strategy_dose = strategy_dose_num)
        } else {
          pred_dat <- pred_dat |>
            dplyr::mutate(
              strategy_dose = as.numeric(
                gsub("[^0-9.]", "", as.character(strategy))
              )
            )
        }
      }
      
      pred_dat |>
        dplyr::left_join(
          placebo_for_join,
          by = c("trial_name", "visit")
        ) |>
        dplyr::left_join(
          active_counts_by_arm,
          by = c("arm_name", "visit")
        ) |>
        dplyr::mutate(
          dose_vs_placebo_difference =
            pred_delta_outcome - placebo_mean_delta,
          
          SE =
            sqrt(se^2 + placebo_se_delta^2),
          
          lower_95 =
            dose_vs_placebo_difference -
            stats::qnorm(0.975) * SE,
          
          upper_95 =
            dose_vs_placebo_difference +
            stats::qnorm(0.975) * SE
        ) |>
        dplyr::transmute(
          weight_model = weight_model,
          arm_name = arm_name,
          trial_name = trial_name,
          week = visit,
          dose_strategy = as.character(strategy),
          strategy_dose = strategy_dose,
          
          active_n_patients = active_n_patients,
          active_n_visit_rows = active_n_visit_rows,
          
          placebo_n_patients = placebo_n_patients,
          placebo_n_visit_rows = placebo_n_visit_rows,
          
          predicted_active_improvement = pred_delta_outcome,
          SE_active_prediction = se,
          
          observed_placebo_improvement = placebo_mean_delta,
          SE_observed_placebo = placebo_se_delta,
          
          dose_vs_placebo_difference = dose_vs_placebo_difference,
          SE = SE,
          lower_95 = lower_95,
          upper_95 = upper_95
        )
    }
  )
  
  out |>
    dplyr::arrange(
      trial_name,
      arm_name,
      strategy_dose,
      week
    )
}

dose_vs_placebo_predictions_all <- make_dose_vs_placebo_table(
  prediction_results = prediction_results,
  placebo_observed_by_trial = placebo_observed_by_trial,
  max_visit = 8,
  weight_model = "Ordinal IPTW"
)

dose_vs_placebo_predictions_all

write_xlsx(
  dose_vs_placebo_predictions_all,
  path = file.path(RESULTS_DIR, "dose_vs_placebo_predictions_all.xlsx"))


gee_and_placebo_means_all <- make_gee_and_placebo_mean_table(
  prediction_results = prediction_results,
  placebo_observed_by_trial = placebo_observed_by_trial,
  placebo_counts_by_trial = placebo_counts_by_trial,
  max_visit = 8,
  weight_model = "Ordinal IPTW",
  repeat_placebo_by_active_arm = TRUE
)

gee_and_placebo_means_all

write_xlsx(
  gee_and_placebo_means_all,
  path = file.path(RESULTS_DIR, "gee_and_placebo_means_all.xlsx"))


################################################################################
######################## MULTINOMIAL IPTW SENSITIVITY BY ARM ###################
################################################################################

if (isTRUE(RUN_MULTINOMIAL_SENSITIVITY)) {
  
  ################################################################################
  ######################## SETTINGS USED BY MULTINOMIAL SECTION ##################
  ################################################################################
  
  TABLE_MAX_VISIT <- 8
  
  if (!exists("DOSE_MODEL_PROFILE_PREVIOUS_DOSE")) {
    DOSE_MODEL_PROFILE_PREVIOUS_DOSE <- 20
  }
  
  if (!exists("DOSE_MODEL_PROFILE_SIDE_EFFECTS")) {
    DOSE_MODEL_PROFILE_SIDE_EFFECTS <- 4
  }
  
  if (!exists("DOSE_MODEL_PROFILE_VISITS")) {
    DOSE_MODEL_PROFILE_VISITS <- 1:8
  }
  
  placebo_counts_for_table <- if (exists("placebo_counts_by_trial")) {
    placebo_counts_by_trial
  } else {
    NULL
  }
  
  ################################################################################
  ######################## MULTINOMIAL STEP 1A: DENOMINATOR #####################
  ################################################################################
  
  multinom_denominator_runs <- safe_imap(
    analysis_arms_model,
    function(dat, arm_name) {
      fit_iptw_multinom_denominator_model(
        data = dat,
        visit_df = VISIT_DF,
        ref_dose = MULTINOMIAL_REF_DOSE
      )
    }
  )
  
  multinom_denominator_status <- get_status_table(multinom_denominator_runs)
  multinom_denominator_status
  
  multinom_denominator_models <- get_success_results(multinom_denominator_runs)
  
  ################################################################################
  ################ MULTINOMIAL STEP 1A-PLOT: DOSE MODEL PROFILE PLOTS ###########
  ################################################################################
  
  multinom_dose_model_profile_runs <- safe_imap(
    multinom_denominator_models,
    function(model, arm_name) {
      
      plot_dat <- make_dose_model_profile_data(
        model = model,
        data = analysis_arms_model[[arm_name]],
        previous_dose = DOSE_MODEL_PROFILE_PREVIOUS_DOSE,
        side_effect_value = DOSE_MODEL_PROFILE_SIDE_EFFECTS,
        visits = DOSE_MODEL_PROFILE_VISITS
      )
      
      make_dose_model_profile_plot(
        plot_dat = plot_dat,
        arm_name = paste0(arm_name, " - multinomial IPTW")
      )
    }
  )
  
  multinom_dose_model_profile_status <- get_status_table(
    multinom_dose_model_profile_runs
  )
  multinom_dose_model_profile_status
  
  multinom_dose_model_profile_plots <- get_success_results(
    multinom_dose_model_profile_runs
  )
  
  names(multinom_dose_model_profile_plots)
  
  if (length(multinom_dose_model_profile_plots) > 0) {
    multinom_dose_model_profile_plots[[1]]
  }
  
  ################################################################################
  ######################## MULTINOMIAL STEP 1B: NUMERATOR #######################
  ################################################################################
  
  multinom_numerator_runs <- safe_imap(
    analysis_arms_model,
    function(dat, arm_name) {
      fit_iptw_multinom_numerator_model(
        data = dat,
        visit_df = VISIT_DF,
        ref_dose = MULTINOMIAL_REF_DOSE
      )
    }
  )
  
  multinom_numerator_status <- get_status_table(multinom_numerator_runs)
  multinom_numerator_status
  
  multinom_numerator_models <- get_success_results(multinom_numerator_runs)
  
  ################################################################################
  ######################## MULTINOMIAL STEP 1C: ADD WEIGHTS #####################
  ################################################################################
  
  multinom_common_arms <- Reduce(
    intersect,
    list(
      names(analysis_arms_model),
      names(multinom_denominator_models),
      names(multinom_numerator_models),
      names(analysis_arms_ipdcw)
    )
  )
  
  multinom_common_arms
  
  multinom_weight_runs <- safe_imap(
    analysis_arms_ipdcw[multinom_common_arms],
    function(dat, arm_name) {
      add_iptw_treatment_weights_multinom(
        data = dat,
        denominator_model = multinom_denominator_models[[arm_name]],
        numerator_model = multinom_numerator_models[[arm_name]]
      ) |>
        add_total_weights_multinom() |>
        truncate_total_weights_multinom(
          lower = TRUNCATION[1],
          upper = TRUNCATION[2]
        )
    }
  )
  
  multinom_weight_status <- get_status_table(multinom_weight_runs)
  multinom_weight_status
  
  analysis_arms_multinom_weighted <- get_success_results(multinom_weight_runs)
  
  multinom_weight_summaries <- purrr::map(
    analysis_arms_multinom_weighted,
    summarise_iptw_multinom_treatment_weights
  )
  
  multinom_total_weight_summaries <- purrr::map(
    analysis_arms_multinom_weighted,
    summarise_total_multinom_weights
  )
  
  names(multinom_weight_summaries)
  
  if (length(multinom_weight_summaries) > 0) {
    multinom_weight_summaries[[1]]
  }
  
  if (length(multinom_total_weight_summaries) > 0) {
    multinom_total_weight_summaries[[1]]
  }
  
  ################################################################################
  ######################## MULTINOMIAL STEP 2: WEIGHTED MSM #####################
  ################################################################################
  
  multinom_msm_runs <- safe_imap(
    analysis_arms_multinom_weighted,
    function(dat, arm_name) {
      fit_weighted_msm(
        data = dat,
        weight_var = "SW_total_multinom_trunc",
        visit_df = VISIT_DF,
        corstr = "independence",
        include_dose_time_interaction = INCLUDE_DOSE_TIME_INTERACTION
      )
    }
  )
  
  multinom_msm_status <- get_status_table(multinom_msm_runs)
  multinom_msm_status
  
  multinom_msm_models <- get_success_results(multinom_msm_runs)
  
  multinom_msm_coef_tabs <- purrr::map(
    multinom_msm_models,
    make_gee_coef_table
  )
  
  multinom_msm_robust_naive_tabs <- purrr::map(
    multinom_msm_models,
    make_gee_coef_table_robust_naive
  )
  
  names(multinom_msm_coef_tabs)
  
  if (length(multinom_msm_coef_tabs) > 0) {
    multinom_msm_coef_tabs[[1]]
  }
  
  ################################################################################
  ######################## MULTINOMIAL STEP 3: PREDICTIONS ######################
  ################################################################################
  
  multinom_prediction_runs <- safe_imap(
    multinom_msm_models,
    function(model, arm_name) {
      make_prediction_results_one_arm(
        model = model,
        arm_name = arm_name
      )
    }
  )
  
  multinom_prediction_status <- get_status_table(multinom_prediction_runs)
  multinom_prediction_status
  
  if (any(!multinom_prediction_status$ok)) {
    
    failed_multinom_prediction_arms <- multinom_prediction_status |>
      dplyr::filter(!ok) |>
      dplyr::pull(arm_name)
    
    for (arm_name in failed_multinom_prediction_arms) {
      cat("\n\n==============================\n")
      cat("Multinomial prediction error:", arm_name, "\n")
      cat("==============================\n")
      cat(strip_ansi(multinom_prediction_runs[[arm_name]]$error), "\n")
    }
  }
  
  multinom_prediction_results <- get_success_results(multinom_prediction_runs)
  
  multinom_strategy_prediction_data <- purrr::map(
    multinom_prediction_results,
    "strategy_dat"
  )
  
  multinom_strategy_prediction_plots <- purrr::imap(
    multinom_prediction_results,
    function(res, arm_name) {
      res$plot +
        ggplot2::labs(
          title = paste0(
            arm_name,
            " - multinomial IPTW"
          )
        )
    }
  )
  
  names(multinom_strategy_prediction_plots)
  
  if (length(multinom_strategy_prediction_plots) > 0) {
    multinom_strategy_prediction_plots[[1]]
  }
  
  ################################################################################
  ################ MULTINOMIAL STEP 4: DOSE VERSUS OBSERVED PLACEBO #############
  ################################################################################
  
  multinom_dose_vs_placebo_predictions_all <- make_dose_vs_placebo_table(
    prediction_results = multinom_prediction_results,
    placebo_observed_by_trial = placebo_observed_by_trial,
    placebo_counts_by_trial = placebo_counts_for_table,
    max_visit = TABLE_MAX_VISIT,
    weight_model = "Multinomial IPTW"
  )
  
  multinom_dose_vs_placebo_predictions_all
  
  ################################################################################
  ################ COMPARE ORDINAL VERSUS MULTINOMIAL AGAINST PLACEBO ###########
  ################################################################################
  
  if (exists("dose_vs_placebo_predictions_all")) {
    
    compare_ordinal_vs_multinom_placebo <- dplyr::bind_rows(
      dose_vs_placebo_predictions_all,
      multinom_dose_vs_placebo_predictions_all
    )
    
  } else {
    
    compare_ordinal_vs_multinom_placebo <- multinom_dose_vs_placebo_predictions_all
  }
  
  compare_ordinal_vs_multinom_placebo
}

multinom_gee_and_placebo_means_all <- make_gee_and_placebo_mean_table(
  prediction_results = multinom_prediction_results,
  placebo_observed_by_trial = placebo_observed_by_trial,
  placebo_counts_by_trial = placebo_counts_by_trial,
  max_visit = 8,
  weight_model = "Multinomial IPTW",
  repeat_placebo_by_active_arm = TRUE
)

multinom_gee_and_placebo_means_all


write_xlsx(
  multinom_gee_and_placebo_means_all,
  path = file.path(RESULTS_DIR, "multinom_gee_and_placebo_means_all.xlsx"))

gee_and_placebo_means_compare <- dplyr::bind_rows(
  gee_and_placebo_means_all,
  multinom_gee_and_placebo_means_all
)

gee_and_placebo_means_compare

write_xlsx(
  gee_and_placebo_means_compare,
  path = file.path(RESULTS_DIR, "gee_and_placebo_means_compare.xlsx"))

################################################################################
######################## PRINT ALL MULTINOMIAL PLOTS ############################
################################################################################

if (exists("multinom_strategy_prediction_plots")) {
  purrr::iwalk(
    multinom_strategy_prediction_plots,
    function(plot_obj, arm_name) {
      print(plot_obj)
    }
  )
}

if (exists("multinom_dose_model_profile_plots")) {
  purrr::iwalk(
    multinom_dose_model_profile_plots,
    function(plot_obj, arm_name) {
      print(plot_obj)
    }
  )
}

################################################################################
################ SAVE MULTINOMIAL PLOTS AS TWO JPG GRIDS ########################
################################################################################

# 1) Multinomial final strategy prediction plots
if (exists("multinom_strategy_prediction_plots")) {
  
  multinom_strategy_prediction_grid <- save_plot_grid_jpg(
    plot_list = multinom_strategy_prediction_plots,
    file_path = file.path(RESULTS_DIR, "multinom_strategy_prediction_plots_grid.jpg"),
    ncol = 4,
    base_width = 6,
    base_height = 5
  )
  
  multinom_strategy_prediction_grid
}

# 2) Multinomial dose-model profile plots
if (exists("multinom_dose_model_profile_plots")) {
  
  multinom_dose_model_profile_grid <- save_plot_grid_jpg(
    plot_list = multinom_dose_model_profile_plots,
    file_path = file.path(RESULTS_DIR, "multinom_dose_model_profile_plots_grid.jpg"),
    ncol = 4,
    base_width = 5.5,
    base_height = 6
  )
  
  multinom_dose_model_profile_grid
}
