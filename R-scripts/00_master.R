#step-by-step per-arm master; dose tables and plots are lists by arm.
#######################################################################################################
############################ DR_03_master_stepwise_by_arm.R ###########################################
#######################################################################################################

# Run this script section-by-section in RStudio.
# Results are kept separately per active arm.

source("01_functions.R")
load_required_packages()

source("02_data_prep_one_drug_pooled_time.R")

################################################################################
############################ USER SETTINGS #####################################
################################################################################

RESULTS_DIR <- "/Users/kchalkou/Desktop/Projects_ongoing/Dose-response/Results"
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Restricted cubic spline specification for continuous time.
VISIT_DF <- 3

# Weight truncation.
TRUNCATION <- c(0.01, 0.99)

# MSM specification.
INCLUDE_DOSE_TIME_INTERACTION <- TRUE

# Use censoring weights.
USE_IPCW <- TRUE

# NULL = use all supported grouped doses.
STRATEGY_DOSES_ZERO_AS_DOSE <- c(20, 30, 40, 50)
PLOT_STRATEGY_DOSES <- STRATEGY_DOSES_ZERO_AS_DOSE

STRATEGY_DOSES_ZERO_AS_DISCONTINUATION <- c(20, 30, 40, 50)


# Real PAROXETINE data start at 20 mg.
BASELINE_DOSE_FOR_PREDICTIONS <- 20

# Multinomial sensitivity will be turned back on after its
# treatment-weight timing is updated to match the ordinal analysis.
RUN_MULTINOMIAL_SENSITIVITY <- TRUE
MULTINOMIAL_REF_DOSE <- 20

SAVE_PLOTS <- FALSE

PREDICTION_Y_LIMITS <- c(-10, 20)
PREDICTION_Y_BREAKS <- seq(-10, 40, by = 10)
DOSE_COLORS <- c(
  "0"  = "#7F7F7F",
  "20" = "#F8766D",
  "30" = "#7CAE00",
  "40" = "#00BFC4",
  "50" = "#C77CFF"
)



################################################################################
######################## CREATE MODEL-READY DATASETS ############################
################################################################################

# Zero as a valid dose.
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
  dose_history_levels = sort(unique(c(
    0,
    dose_levels_zero_as_dose
  ))),
  include_zero_as_dose = TRUE,
  censoring_scenario = "zero_as_dose"
)


# Zero as treatment discontinuation.
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


# Put the two analyses in the list expected by the master.
analysis_arms <- list(
  zero_as_dose = analysis_zero_as_dose
)

analysis_arms



################################################################################
######################## STEP 0: DESCRIPTIVE TABLES #############################
################################################################################

# Broad windows are used only for descriptive summaries.
# Continuous visit time is retained for modelling.
add_analysis_time_window <- function(dat) {
  
  dat |>
    dplyr::mutate(
      time_window = dplyr::case_when(
        visit == 0 ~ "Baseline",
        visit > 0 & visit <= 2 ~ ">0-2 weeks",
        visit > 2 & visit <= 4 ~ ">2-4 weeks",
        visit > 4 & visit <= VISIT_MAX ~ ">4-6 weeks",
        TRUE ~ NA_character_
      ),
      time_window = factor(
        time_window,
        levels = c(
          "Baseline",
          ">0-2 weeks",
          ">2-4 weeks",
          ">4-6 weeks"
        )
      )
    )
}


######################## DOSE DISTRIBUTION #####################################

dose_by_window_tables <- purrr::imap(
  analysis_arms,
  function(dat, analysis_name) {
    
    dat |>
      add_analysis_time_window() |>
      dplyr::filter(
        visit > 0,
        visit <= VISIT_MAX,
        !is.na(dose_f)
      ) |>
      dplyr::count(
        time_window,
        dose_f,
        name = "n"
      ) |>
      tidyr::pivot_wider(
        names_from = dose_f,
        values_from = n,
        values_fill = 0
      ) |>
      dplyr::arrange(time_window)
  }
)

dose_by_window_tables

writexl::write_xlsx(
  dose_by_window_tables,
  path = file.path(
    RESULTS_DIR,
    "dose_by_time_window_tables.xlsx"
  )
)


######################## DOSE TRANSITIONS ######################################

dose_transition_tables <- purrr::imap(
  analysis_arms,
  function(dat, analysis_name) {
    
    dat |>
      add_analysis_time_window() |>
      dplyr::filter(
        visit > 0,
        visit <= VISIT_MAX,
        !is.na(dose_f),
        !is.na(dose_lag1_f)
      ) |>
      dplyr::count(
        time_window,
        dose_lag1_f,
        dose_f,
        name = "n"
      ) |>
      tidyr::pivot_wider(
        names_from = dose_f,
        values_from = n,
        values_fill = 0
      ) |>
      dplyr::arrange(
        time_window,
        dose_lag1_f
      )
  }
)

dose_transition_tables

writexl::write_xlsx(
  dose_transition_tables,
  path = file.path(
    RESULTS_DIR,
    "dose_transition_tables.xlsx"
  )
)


######################## SPARSE-DOSE GROUPING ##################################

dose_grouping_tables <- purrr::imap(
  analysis_arms,
  function(dat, analysis_name) {
    
    if (!all(
      c(
        "dose_before_grouping",
        "dose_grouped"
      ) %in% names(dat)
    )) {
      return(tibble::tibble())
    }
    
    dat |>
      dplyr::filter(
        visit > 0,
        visit <= VISIT_MAX
      ) |>
      dplyr::count(
        dose_before_grouping,
        dose_grouped,
        name = "n"
      ) |>
      dplyr::arrange(
        dose_before_grouping,
        dose_grouped
      )
  }
)

dose_grouping_tables

writexl::write_xlsx(
  dose_grouping_tables,
  path = file.path(
    RESULTS_DIR,
    "dose_grouping_tables.xlsx"
  )
)


################################################################################
######################## STEP 0B: DESCRIPTIVE PLOTS #############################
################################################################################

dose_by_window_plots <- purrr::imap(
  analysis_arms,
  function(dat, analysis_name) {
    
    plot_dat <- dat |>
      add_analysis_time_window() |>
      dplyr::filter(
        visit > 0,
        visit <= VISIT_MAX,
        !is.na(dose_f)
      ) |>
      dplyr::count(
        time_window,
        dose_f,
        name = "n"
      )
    
    ggplot2::ggplot(
      plot_dat,
      ggplot2::aes(
        x = time_window,
        y = n,
        fill = dose_f
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::labs(
        x = "Time from baseline",
        y = "Number of patient-visits",
        fill = "Dose",
        title = analysis_name
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title =
          ggplot2::element_text(face = "bold")
      )
  }
)

dose_by_window_plots


dose_transition_plots <- purrr::imap(
  analysis_arms,
  function(dat, analysis_name) {
    
    plot_dat <- dat |>
      add_analysis_time_window() |>
      dplyr::filter(
        visit > 0,
        visit <= VISIT_MAX,
        !is.na(dose_f),
        !is.na(dose_lag1_f)
      ) |>
      dplyr::count(
        time_window,
        dose_lag1_f,
        dose_f,
        name = "n"
      )
    
    ggplot2::ggplot(
      plot_dat,
      ggplot2::aes(
        x = dose_f,
        y = dose_lag1_f,
        fill = n
      )
    ) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(
        ggplot2::aes(label = n),
        size = 3
      ) +
      ggplot2::facet_wrap(
        ~ time_window
      ) +
      ggplot2::labs(
        x = "Current dose",
        y = "Previous observed dose",
        fill = "n",
        title = analysis_name
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title =
          ggplot2::element_text(face = "bold"),
        axis.text.x =
          ggplot2::element_text(
            angle = 45,
            hjust = 1
          )
      )
  }
)

dose_transition_plots


dose_grouping_plots <- purrr::imap(
  analysis_arms,
  function(dat, analysis_name) {
    
    if (!all(
      c(
        "dose_before_grouping",
        "dose_grouped"
      ) %in% names(dat)
    )) {
      return(NULL)
    }
    
    plot_dat <- dat |>
      dplyr::filter(
        visit > 0,
        visit <= VISIT_MAX
      ) |>
      dplyr::count(
        dose_before_grouping,
        dose_grouped,
        name = "n"
      )
    
    ggplot2::ggplot(
      plot_dat,
      ggplot2::aes(
        x = factor(dose_before_grouping),
        y = n,
        fill = factor(dose_grouped)
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::labs(
        x = "Original dose",
        y = "Number of patient-visits",
        fill = "Grouped dose",
        title = analysis_name
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title =
          ggplot2::element_text(face = "bold")
      )
  }
)

dose_grouping_plots


# Save only when requested.
if (isTRUE(SAVE_PLOTS)) {
  
  save_plot_grid_jpg(
    plot_list = dose_by_window_plots,
    file_path = file.path(
      RESULTS_DIR,
      "dose_by_time_window_plots.jpg"
    ),
    ncol = 2
  )
  
  save_plot_grid_jpg(
    plot_list = dose_transition_plots,
    file_path = file.path(
      RESULTS_DIR,
      "dose_transition_plots.jpg"
    ),
    ncol = 2
  )
  
  save_plot_grid_jpg(
    plot_list = dose_grouping_plots,
    file_path = file.path(
      RESULTS_DIR,
      "dose_grouping_plots.jpg"
    ),
    ncol = 2
  )
}

positivity_by_study_window <- purrr::imap(
  analysis_arms,
  function(dat, analysis_name) {
    
    dat |>
      add_analysis_time_window() |>
      dplyr::filter(
        use_treatment_weight,
        !is.na(dose_f)
      ) |>
      dplyr::count(
        studyid,
        time_window,
        dose_f,
        name = "n"
      ) |>
      tidyr::complete(
        studyid,
        time_window,
        dose_f,
        fill = list(n = 0)
      ) |>
      dplyr::arrange(
        studyid,
        time_window,
        dose_f
      )
  }
)

positivity_by_study_window

positivity_by_study_window_wide <- purrr::map(
  positivity_by_study_window,
  ~ .x |>
    tidyr::pivot_wider(
      names_from = dose_f,
      values_from = n,
      values_fill = 0
    )
)

positivity_by_study_window_wide
################################################################################
######################## STEP 0C: MODELABILITY #################################
################################################################################

MIN_PATIENTS_FOR_MODEL <- 30
MIN_TREATMENT_ROWS_FOR_MODEL <- 100
MIN_DOSE_LEVELS_FOR_MODEL <- 2

analysis_modelability <- purrr::imap_dfr(
  analysis_arms,
  function(dat, analysis_name) {
    
    # These are exactly the rows entering the treatment-assignment model.
    model_dat <- dat |>
      dplyr::filter(use_treatment_weight)
    
    dose_counts <- model_dat |>
      dplyr::count(
        dose_f,
        name = "n"
      )
    
    tibble::tibble(
      analysis_name = analysis_name,
      
      n_patients =
        dplyr::n_distinct(model_dat$pid),
      
      n_treatment_rows =
        nrow(model_dat),
      
      n_dose_levels =
        nrow(dose_counts),
      
      min_n_per_dose =
        if (nrow(dose_counts) > 0) {
          min(dose_counts$n)
        } else {
          0
        },
      
      dose_levels =
        paste(
          as.character(dose_counts$dose_f),
          collapse = ", "
        ),
      
      modelable =
        dplyr::n_distinct(model_dat$pid) >=
        MIN_PATIENTS_FOR_MODEL &&
        nrow(model_dat) >=
        MIN_TREATMENT_ROWS_FOR_MODEL &&
        nrow(dose_counts) >=
        MIN_DOSE_LEVELS_FOR_MODEL
    )
  }
)

analysis_modelability

writexl::write_xlsx(
  analysis_modelability,
  path = file.path(
    RESULTS_DIR,
    "analysis_modelability.xlsx"
  )
)

modelable_analysis_names <- analysis_modelability |>
  dplyr::filter(modelable) |>
  dplyr::pull(analysis_name)

modelable_analysis_names

descriptive_only_analysis_names <- analysis_modelability |>
  dplyr::filter(!modelable) |>
  dplyr::pull(analysis_name)

descriptive_only_analysis_names

analysis_arms_model <-
  analysis_arms[modelable_analysis_names]

analysis_arms_descriptive_only <-
  analysis_arms[descriptive_only_analysis_names]


################################################################################
######################## STEP 1A: IPTW DENOMINATOR BY ARM ######################
################################################################################

# Fit the ordinal treatment-assignment denominator model
# separately for the two zero-dose analyses.
iptw_denominator_runs <- safe_imap(
  analysis_arms_model,
  function(dat, analysis_name) {
    
    fit_iptw_denominator_model(
      data = dat,
      visit_df = VISIT_DF
    )
  }
)


################################################################################
######################## MODEL STATUS ##########################################
################################################################################

iptw_denominator_status <- get_status_table(
  iptw_denominator_runs
)

iptw_denominator_status


# Keep successfully fitted models.
iptw_denominator_models <- get_success_results(
  iptw_denominator_runs
)

names(iptw_denominator_models)


################################################################################
######################## COEFFICIENT TABLES ####################################
################################################################################

iptw_denominator_coef_tabs <- purrr::map(
  iptw_denominator_models,
  make_polr_coef_table
)

iptw_denominator_coef_tabs


# Add the analysis name to each coefficient table.
iptw_denominator_coef_tabs_excel <- purrr::imap(
  iptw_denominator_coef_tabs,
  function(tab, analysis_name) {
    
    tab_df <- as.data.frame(tab)
    
    if (!is.null(rownames(tab_df))) {
      tab_df <- tab_df |>
        tibble::rownames_to_column("term")
    }
    
    tab_df |>
      dplyr::mutate(
        analysis_name = analysis_name,
        .before = 1
      )
  }
)


################################################################################
######################## SAVE COEFFICIENT TABLES ###############################
################################################################################

writexl::write_xlsx(
  iptw_denominator_coef_tabs_excel,
  path = file.path(
    RESULTS_DIR,
    "iptw_denominator_coef_tabs.xlsx"
  )
)


################################################################################
######## STEP 1A-PLOT: OBSERVED VS PREDICTED DOSE PROBABILITIES ################
################################################################################

make_iptw_dose_model_plot_data <- function(model, data) {
  
  model_dat <- data |>
    dplyr::filter(use_treatment_weight) |>
    add_analysis_time_window() |>
    dplyr::filter(!is.na(time_window)) |>
    dplyr::mutate(
      .row_id = dplyr::row_number(),
      studyid = as.character(studyid),
      observed_dose = as.character(dose_f)
    )
  
  # Predicted probability of every possible dose
  # for every treatment-decision row.
  prob_mat <- predict(
    model,
    newdata = model_dat,
    type = "probs"
  )
  
  if (is.null(dim(prob_mat))) {
    
    prob_names <- names(prob_mat)
    
    prob_mat <- matrix(
      prob_mat,
      nrow = 1
    )
    
    colnames(prob_mat) <- prob_names
  }
  
  dose_levels <- colnames(prob_mat)
  
  
  ##############################################################################
  # MODEL-PREDICTED PROBABILITIES
  ##############################################################################
  
  prob_long <- as.data.frame(prob_mat) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      .row_id = dplyr::row_number()
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(dose_levels),
      names_to = "dose",
      values_to = "predicted_probability"
    ) |>
    dplyr::left_join(
      model_dat |>
        dplyr::select(
          .row_id,
          studyid,
          time_window
        ),
      by = ".row_id"
    )
  
  
  predicted_summary <- prob_long |>
    dplyr::group_by(
      studyid,
      time_window,
      dose
    ) |>
    dplyr::summarise(
      predicted_probability =
        mean(
          predicted_probability,
          na.rm = TRUE
        ),
      .groups = "drop"
    )
  
  
  ##############################################################################
  # OBSERVED DOSE PROPORTIONS
  ##############################################################################
  
  observed_counts <- model_dat |>
    dplyr::count(
      studyid,
      time_window,
      observed_dose,
      name = "observed_n"
    )
  
  
  observed_totals <- model_dat |>
    dplyr::count(
      studyid,
      time_window,
      name = "observed_total"
    )
  
  
  ##############################################################################
  # COMBINE OBSERVED AND PREDICTED
  ##############################################################################
  
  tidyr::expand_grid(
    studyid =
      sort(unique(model_dat$studyid)),
    
    time_window =
      levels(model_dat$time_window),
    
    dose =
      dose_levels
  ) |>
    
    dplyr::filter(
      time_window != "Baseline"
    ) |>
    
    dplyr::left_join(
      predicted_summary,
      by = c(
        "studyid",
        "time_window",
        "dose"
      )
    ) |>
    
    dplyr::left_join(
      observed_counts |>
        dplyr::rename(
          dose = observed_dose
        ),
      by = c(
        "studyid",
        "time_window",
        "dose"
      )
    ) |>
    
    dplyr::left_join(
      observed_totals,
      by = c(
        "studyid",
        "time_window"
      )
    ) |>
    
    dplyr::mutate(
      
      observed_n =
        tidyr::replace_na(
          observed_n,
          0L
        ),
      
      observed_probability =
        observed_n / observed_total,
      
      time_window = factor(
        time_window,
        levels = c(
          ">0-2 weeks",
          ">2-4 weeks",
          ">4-6 weeks"
        )
      ),
      
      dose = factor(
        dose,
        levels = dose_levels
      )
    )
}

make_iptw_dose_model_plot <- function(
    plot_dat,
    analysis_name
) {
  
  plot_long <- plot_dat |>
    dplyr::select(
      studyid,
      time_window,
      dose,
      observed_probability,
      predicted_probability
    ) |>
    
    tidyr::pivot_longer(
      cols = c(
        observed_probability,
        predicted_probability
      ),
      names_to = "source",
      values_to = "probability"
    ) |>
    
    dplyr::mutate(
      source = dplyr::recode(
        source,
        observed_probability =
          "Observed",
        predicted_probability =
          "Model predicted"
      )
    )
  
  
  ggplot2::ggplot(
    plot_long,
    ggplot2::aes(
      x = dose,
      y = probability,
      color = dose,
      shape = source
    )
  ) +
    
    ggplot2::geom_point(
      position =
        ggplot2::position_dodge(
          width = 0.25
        ),
      size = 3
    ) +
    
    # Same dose colors as your other figures.
    ggplot2::scale_color_manual(
      values = DOSE_COLORS
    ) +
    
    ggplot2::facet_grid(
      studyid ~ time_window
    ) +
    
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        by = 0.2
      ),
      labels =
        scales::percent_format(
          accuracy = 1
        )
    ) +
    
    ggplot2::labs(
      x = "Prescribed dose (mg)",
      y = "Probability",
      color = "Dose (mg)",
      shape = NULL,
      
      title = paste0(
        "Observed vs model-predicted dose probabilities: ",
        analysis_name
      )
    ) +
    
    ggplot2::theme_minimal(
      base_size = 11
    ) +
    
    ggplot2::theme(
      plot.title =
        ggplot2::element_text(
          face = "bold"
        ),
      legend.position =
        "bottom"
    )
}

iptw_dose_model_plot_data <- purrr::imap(
  iptw_denominator_models,
  function(model, analysis_name) {
    
    make_iptw_dose_model_plot_data(
      model = model,
      data =
        analysis_arms_model[[analysis_name]]
    )
  }
)


iptw_dose_model_plots <- purrr::imap(
  iptw_dose_model_plot_data,
  function(plot_dat, analysis_name) {
    
    make_iptw_dose_model_plot(
      plot_dat = plot_dat,
      analysis_name = analysis_name
    )
  }
)


iptw_dose_model_plots
################################################################################
######## STEP 1A-PLOT: RESPONSE-PROFILE DOSE PROBABILITIES WITH FIXED SIDE EFFECTS ############
################################################################################
#Among patient-visits that were actually improving versus those 
# that were actually not improving/worsening, what dose probabilities 
# does the fitted model predict if side effects are standardized to 4?

PROFILE_SIDE_EFFECT_VALUE <- 4
make_response_profile_probability_data <- function(
    model,
    data,
    side_effect_value = 4
) {
  
  model_dat <- data |>
    dplyr::filter(
      use_treatment_weight,
      !is.na(delta_outcome_locf)
    ) |>
    
    dplyr::mutate(
      
      response_profile = dplyr::case_when(
        
        delta_outcome_locf > 2 ~
          "Improving",
        
        delta_outcome_locf < -2 ~
          "Worsening",
        
        TRUE ~
          NA_character_
      ),
      
      response_profile = factor(
        response_profile,
        levels = c(
          "Improving",
          "Worsening"
        )
      )
    ) |>
    
    # Exclude the middle range: -2 to +2
    dplyr::filter(
      !is.na(response_profile)
    ) |>
    
    dplyr::mutate(
      
      # Standardize side effects to the same value
      # in both response groups.
      side.effects_model_locf =
        side_effect_value,
      
      .row_id =
        dplyr::row_number()
    )
  
  ##############################################################################
  # Predicted probability of each possible prescribed dose
  ##############################################################################
  
  prob_mat <- predict(
    model,
    newdata = model_dat,
    type = "probs"
  )
  
  if (is.null(dim(prob_mat))) {
    
    prob_names <- names(prob_mat)
    
    prob_mat <- matrix(
      prob_mat,
      nrow = 1
    )
    
    colnames(prob_mat) <- prob_names
  }
  
  dose_levels <- colnames(prob_mat)
  
  
  ##############################################################################
  # Long-format prediction dataset
  ##############################################################################
  
  probability_dat <- as.data.frame(prob_mat) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      .row_id = dplyr::row_number()
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(dose_levels),
      names_to = "dose",
      values_to = "predicted_probability"
    ) |>
    dplyr::left_join(
      model_dat |>
        dplyr::select(
          .row_id,
          pid,
          studyid,
          visit,
          response_profile
        ),
      by = ".row_id"
    ) |>
    dplyr::mutate(
      dose = factor(
        dose,
        levels = dose_levels
      )
    )
  
  probability_dat
}

make_response_profile_probability_plot <- function(
    plot_dat,
    analysis_name,
    side_effect_value = 4
) {
  
  ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(
      x = visit,
      y = predicted_probability,
      color = dose,
      group = dose
    )
  ) +
    
    ggplot2::geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      span = 0.8,
      linewidth = 1
    ) +
    
    ggplot2::facet_wrap(
      ~ response_profile,
      nrow = 1
    ) +
    
    # SAME DOSE COLORS AS EVERY OTHER PLOT
    ggplot2::scale_color_manual(
      values = DOSE_COLORS
    ) +
    
    ggplot2::scale_x_continuous(
      breaks = 0:VISIT_MAX,
      limits = c(0, VISIT_MAX)
    ) +
    
    ggplot2::scale_y_continuous(
      breaks = seq(
        0,
        1,
        by = 0.25
      ),
      labels =
        scales::percent_format(
          accuracy = 1
        )
    ) +
    
    ggplot2::coord_cartesian(
      ylim = c(0, 1)
    ) +
    
    ggplot2::labs(
      x = "Time from baseline (weeks)",
      y = "Predicted probability of prescribed dose",
      color = "Dose (mg)",
      title = paste0(
        "Dose assignment by response profile: ",
        analysis_name
      ),
      subtitle = paste0(
        "Side-effect score fixed at ",
        side_effect_value,
        " in both response groups"
      )
    ) +
    
    ggplot2::theme_minimal(
      base_size = 12
    ) +
    
    ggplot2::theme(
      plot.title =
        ggplot2::element_text(
          face = "bold"
        ),
      legend.position = "right"
    )
}

iptw_response_profile_data <- purrr::imap(
  iptw_denominator_models,
  function(model, analysis_name) {
    
    make_response_profile_probability_data(
      model = model,
      data =
        analysis_arms_model[[analysis_name]],
      side_effect_value =
        PROFILE_SIDE_EFFECT_VALUE
    )
  }
)

iptw_response_profile_plots <- purrr::imap(
  iptw_response_profile_data,
  function(plot_dat, analysis_name) {
    
    make_response_profile_probability_plot(
      plot_dat = plot_dat,
      analysis_name = analysis_name,
      side_effect_value =
        PROFILE_SIDE_EFFECT_VALUE
    )
  }
)

iptw_response_profile_plots

################################################################################
######################## STEP 1B: IPTW NUMERATOR ###############################
################################################################################

iptw_numerator_runs <- safe_imap(
  analysis_arms_model,
  function(dat, analysis_name) {
    
    fit_iptw_numerator_model(
      data = dat,
      visit_df = VISIT_DF
    )
  }
)


################################################################################
######################## MODEL STATUS ##########################################
################################################################################

iptw_numerator_status <- get_status_table(
  iptw_numerator_runs
)

iptw_numerator_status


# Keep successfully fitted models.
iptw_numerator_models <- get_success_results(
  iptw_numerator_runs
)

names(iptw_numerator_models)


################################################################################
######################## COEFFICIENT TABLES ####################################
################################################################################

iptw_numerator_coef_tabs <- purrr::map(
  iptw_numerator_models,
  make_polr_coef_table
)

iptw_numerator_coef_tabs


# Prepare one Excel sheet per analysis scenario.
iptw_numerator_coef_tabs_excel <- purrr::imap(
  iptw_numerator_coef_tabs,
  function(tab, analysis_name) {
    
    tab_df <- as.data.frame(tab)
    
    if (!is.null(rownames(tab_df))) {
      tab_df <- tab_df |>
        tibble::rownames_to_column("term")
    }
    
    tab_df |>
      dplyr::mutate(
        analysis_name = analysis_name,
        .before = 1
      )
  }
)


writexl::write_xlsx(
  iptw_numerator_coef_tabs_excel,
  path = file.path(
    RESULTS_DIR,
    "iptw_numerator_coef_tabs.xlsx"
  )
)


################################################################################
######################## STEP 1C: ADD IPTW WEIGHTS #############################
################################################################################

# Analyses for which both numerator and denominator models fitted successfully.
iptw_common_analysis_names <- Reduce(
  intersect,
  list(
    names(analysis_arms_model),
    names(iptw_denominator_models),
    names(iptw_numerator_models)
  )
)

iptw_common_analysis_names


################################################################################
######################## CALCULATE IPTW ########################################
################################################################################

iptw_weight_runs <- safe_imap(
  analysis_arms_model[iptw_common_analysis_names],
  function(dat, analysis_name) {
    
    add_iptw_treatment_weights(
      data = dat,
      denominator_model =
        iptw_denominator_models[[analysis_name]],
      numerator_model =
        iptw_numerator_models[[analysis_name]]
    )
  }
)


iptw_weight_status <- get_status_table(
  iptw_weight_runs
)

iptw_weight_status


analysis_arms_iptw <- get_success_results(
  iptw_weight_runs
)

names(analysis_arms_iptw)


################################################################################
######################## IPTW SUMMARIES ########################################
################################################################################

iptw_weight_summaries <- purrr::map(
  analysis_arms_iptw,
  summarise_iptw_treatment_weights
)

iptw_weight_summaries

print(
  iptw_weight_summaries$zero_as_dose,
  width = Inf
)


writexl::write_xlsx(
  iptw_weight_summaries,
  path = file.path(
    RESULTS_DIR,
    "iptw_weight_summaries.xlsx"
  )
)

################################################################################
################### IPTW PROBABILITY / POSITIVITY DIAGNOSTICS ##################
################################################################################

iptw_probability_diagnostics <- purrr::imap_dfr(
  analysis_arms_iptw,
  function(dat, analysis_name) {
    
    dat |>
      dplyr::filter(
        use_treatment_weight,
        !is.na(p_dose_denominator)
      ) |>
      dplyr::summarise(
        analysis_name = analysis_name,
        
        n = dplyr::n(),
        n_patients = dplyr::n_distinct(pid),
        
        min_p_denominator =
          min(p_dose_denominator),
        
        p1_p_denominator =
          as.numeric(
            stats::quantile(
              p_dose_denominator,
              0.01
            )
          ),
        
        p5_p_denominator =
          as.numeric(
            stats::quantile(
              p_dose_denominator,
              0.05
            )
          ),
        
        median_p_denominator =
          median(p_dose_denominator),
        
        mean_p_denominator =
          mean(p_dose_denominator),
        
        n_p_below_0_01 =
          sum(p_dose_denominator < 0.01),
        
        percent_p_below_0_01 =
          100 * mean(
            p_dose_denominator < 0.01
          ),
        
        n_p_below_0_05 =
          sum(p_dose_denominator < 0.05),
        
        percent_p_below_0_05 =
          100 * mean(
            p_dose_denominator < 0.05
          )
      )
  }
)

iptw_probability_diagnostics

iptw_probability_by_study_dose <- purrr::imap_dfr(
  analysis_arms_iptw,
  function(dat, analysis_name) {
    
    dat |>
      dplyr::filter(
        use_treatment_weight,
        !is.na(p_dose_denominator)
      ) |>
      dplyr::mutate(
        studyid = as.character(studyid),
        dose = as.character(dose_f)
      ) |>
      dplyr::group_by(
        studyid,
        dose
      ) |>
      dplyr::summarise(
        n = dplyr::n(),
        
        min_p =
          min(p_dose_denominator),
        
        p5_p =
          as.numeric(
            stats::quantile(
              p_dose_denominator,
              0.05
            )
          ),
        
        median_p =
          median(p_dose_denominator),
        
        mean_p =
          mean(p_dose_denominator),
        
        max_p =
          max(p_dose_denominator),
        
        .groups = "drop"
      ) |>
      dplyr::mutate(
        analysis_name = analysis_name,
        .before = 1
      )
  }
)

iptw_probability_by_study_dose

print(
  iptw_probability_diagnostics,
  width = Inf
)

print(
  iptw_probability_by_study_dose,
  width = Inf
)

print(
  iptw_probability_by_study_dose,
  n = Inf,
  width = Inf
)

iptw_low_probability_counts <- purrr::imap_dfr(
  analysis_arms_iptw,
  function(dat, analysis_name) {
    
    dat |>
      dplyr::filter(
        use_treatment_weight,
        !is.na(p_dose_denominator)
      ) |>
      dplyr::mutate(
        studyid = as.character(studyid),
        dose = as.character(dose_f)
      ) |>
      dplyr::group_by(
        studyid,
        dose
      ) |>
      dplyr::summarise(
        n = dplyr::n(),
        
        n_below_0_01 =
          sum(p_dose_denominator < 0.01),
        
        pct_below_0_01 =
          100 * mean(p_dose_denominator < 0.01),
        
        n_below_0_05 =
          sum(p_dose_denominator < 0.05),
        
        pct_below_0_05 =
          100 * mean(p_dose_denominator < 0.05),
        
        .groups = "drop"
      ) |>
      dplyr::mutate(
        analysis_name = analysis_name,
        .before = 1
      )
  }
)

print(
  iptw_low_probability_counts,
  n = Inf,
  width = Inf
)

################################################################################
######################## STEP 2: IPCW ##########################################
################################################################################


  ##############################################################################
  ################ STEP 2A: IPCW DENOMINATOR ###################################
  ##############################################################################
  
  ipcw_denominator_runs <- safe_imap(
    analysis_arms_iptw,
    function(dat, analysis_name) {
      
      fit_ipcw_denominator_model(
        data = dat,
        visit_df = VISIT_DF
      )
    }
  )
  
  
  ipcw_denominator_status <- get_status_table(
    ipcw_denominator_runs
  )
  
  ipcw_denominator_status
  
  
  # Keep successfully fitted models.
  ipcw_denominator_models <- get_success_results(
    ipcw_denominator_runs
  )
  
  names(ipcw_denominator_models)
  
  
  ##############################################################################
  ################ COEFFICIENT TABLES ##########################################
  ##############################################################################
  
  ipcw_denominator_coef_tabs <- purrr::map(
    ipcw_denominator_models,
    make_glm_coef_table
  )
  
  ipcw_denominator_coef_tabs
  
  ################################################################################
  ######## STEP 2A-PLOT: OBSERVED VS PREDICTED CENSORING #########################
  ################################################################################
  make_ipcw_denominator_plot_data <- function(
    model,
    data
  ) {
    
    ##############################################################################
    # Rows used in the censoring denominator model
    ##############################################################################
    
    model_dat <- data |>
      dplyr::filter(
        use_censoring_weight,
        !is.na(R_next)
      ) |>
      add_analysis_time_window() |>
      dplyr::filter(
        !is.na(time_window)
      )
    
    
    ##############################################################################
    # MODEL-PREDICTED PROBABILITY OF REMAINING UNCENSORED
    #
    # R_next = 1 -> remains uncensored
    # R_next = 0 -> censoring event
    ##############################################################################
    
    model_dat$predicted_remain_uncensored <-
      stats::predict(
        model,
        newdata = model_dat,
        type = "response"
      )
    
    
    ##############################################################################
    # Convert to probability of censoring
    ##############################################################################
    
    model_dat <- model_dat |>
      dplyr::mutate(
        
        observed_censor =
          1 - R_next,
        
        predicted_censor =
          1 - predicted_remain_uncensored
      )
    
    
    ##############################################################################
    # Observed vs predicted censoring by study and broad time window
    ##############################################################################
    
    model_dat |>
      dplyr::group_by(
        studyid,
        time_window
      ) |>
      dplyr::summarise(
        
        n_at_risk =
          dplyr::n(),
        
        n_censored =
          sum(
            observed_censor,
            na.rm = TRUE
          ),
        
        observed_censor_probability =
          mean(
            observed_censor,
            na.rm = TRUE
          ),
        
        predicted_censor_probability =
          mean(
            predicted_censor,
            na.rm = TRUE
          ),
        
        .groups = "drop"
      )
  }
  
  
  make_ipcw_denominator_plot <- function(
    plot_dat,
    analysis_name
  ) {
    
    plot_long <- plot_dat |>
      dplyr::select(
        studyid,
        time_window,
        n_at_risk,
        n_censored,
        observed_censor_probability,
        predicted_censor_probability
      ) |>
      tidyr::pivot_longer(
        cols = c(
          observed_censor_probability,
          predicted_censor_probability
        ),
        names_to = "source",
        values_to = "censor_probability"
      ) |>
      dplyr::mutate(
        source = dplyr::recode(
          source,
          
          observed_censor_probability =
            "Observed",
          
          predicted_censor_probability =
            "Model predicted"
        )
      )
    
    
    ggplot2::ggplot(
      plot_long,
      ggplot2::aes(
        x = time_window,
        y = censor_probability,
        group = source,
        shape = source
      )
    ) +
      
      ggplot2::geom_point(
        size = 3,
        position =
          ggplot2::position_dodge(
            width = 0.20
          )
      ) +
      
      ggplot2::geom_line(
        position =
          ggplot2::position_dodge(
            width = 0.20
          )
      ) +
      
      ggplot2::facet_wrap(
        ~ studyid
      ) +
      
      ggplot2::scale_y_continuous(
        labels =
          scales::percent_format(
            accuracy = 1
          )
      ) +
      
      ggplot2::labs(
        x = "Time from baseline",
        y = "Probability of censoring",
        shape = NULL,
        
        title = paste0(
          "Observed vs model-predicted censoring: ",
          analysis_name
        )
      ) +
      
      ggplot2::theme_minimal(
        base_size = 12
      ) +
      
      ggplot2::theme(
        plot.title =
          ggplot2::element_text(
            face = "bold"
          ),
        legend.position = "bottom"
      )
  }
  
  ipcw_denominator_plot_data <- purrr::imap(
    ipcw_denominator_models,
    function(model, analysis_name) {
      
      make_ipcw_denominator_plot_data(
        model = model,
        data =
          analysis_arms_iptw[[analysis_name]]
      )
    }
  )
  
  ipcw_denominator_plots <- purrr::imap(
    ipcw_denominator_plot_data,
    function(plot_dat, analysis_name) {
      
      make_ipcw_denominator_plot(
        plot_dat = plot_dat,
        analysis_name = analysis_name
      )
    }
  )
  
  ipcw_denominator_plots
  ##############################################################################
  ################ SAVE COEFFICIENT TABLES #####################################
  ##############################################################################
  
  ipcw_denominator_coef_tabs_excel <- purrr::imap(
    ipcw_denominator_coef_tabs,
    function(tab, analysis_name) {
      
      tab_df <- as.data.frame(tab)
      
      if (!is.null(rownames(tab_df))) {
        tab_df <- tab_df |>
          tibble::rownames_to_column("term")
      }
      
      tab_df |>
        dplyr::mutate(
          analysis_name = analysis_name,
          .before = 1
        )
    }
  )
  
  
  writexl::write_xlsx(
    ipcw_denominator_coef_tabs_excel,
    path = file.path(
      RESULTS_DIR,
      "ipcw_denominator_coef_tabs.xlsx"
    )
  )
  
  ipcw_censoring_by_dose <- purrr::imap(
    analysis_arms_iptw,
    function(dat, analysis_name) {
      
      dat |>
        dplyr::filter(
          use_censoring_weight
        ) |>
        dplyr::count(
          studyid,
          dose_censor_f,
          R_next,
          name = "n"
        ) |>
        dplyr::arrange(
          studyid,
          dose_censor_f,
          R_next
        )
    }
  )
  
  ipcw_censoring_by_dose

  ################################################################################
  ######################## STEP 2B: IPCW NUMERATOR ###############################
  ################################################################################
  
  if (isTRUE(USE_IPCW)) {
    
    ipcw_numerator_runs <- safe_imap(
      analysis_arms_iptw,
      function(dat, analysis_name) {
        
        fit_ipcw_numerator_model(
          data = dat,
          visit_df = VISIT_DF
        )
      }
    )
    
    
    ##############################################################################
    ################ MODEL STATUS ################################################
    ##############################################################################
    
    ipcw_numerator_status <- get_status_table(
      ipcw_numerator_runs
    )
    
    ipcw_numerator_status
    
    
    # Keep successfully fitted models.
    ipcw_numerator_models <- get_success_results(
      ipcw_numerator_runs
    )
    
    names(ipcw_numerator_models)
    
    
    ##############################################################################
    ################ COEFFICIENT TABLES ##########################################
    ##############################################################################
    
    ipcw_numerator_coef_tabs <- purrr::map(
      ipcw_numerator_models,
      make_glm_coef_table
    )
    
    ipcw_numerator_coef_tabs
    
    
    ##############################################################################
    ################ SAVE COEFFICIENT TABLES #####################################
    ##############################################################################
    
    ipcw_numerator_coef_tabs_excel <- purrr::imap(
      ipcw_numerator_coef_tabs,
      function(tab, analysis_name) {
        
        tab_df <- as.data.frame(tab)
        
        if (!is.null(rownames(tab_df))) {
          tab_df <- tab_df |>
            tibble::rownames_to_column("term")
        }
        
        tab_df |>
          dplyr::mutate(
            analysis_name = analysis_name,
            .before = 1
          )
      }
    )
    
    
    writexl::write_xlsx(
      ipcw_numerator_coef_tabs_excel,
      path = file.path(
        RESULTS_DIR,
        "ipcw_numerator_coef_tabs.xlsx"
      )
    )
  }
  
  ipcw_numerator_status
  
  ipcw_numerator_coef_tabs
  
  ################################################################################
  ######################## STEP 2C: ADD IPCW WEIGHTS #############################
  ################################################################################
  
  if (isTRUE(USE_IPCW)) {
    
    # Keep only analyses where both censoring models fitted successfully.
    ipcw_common_analysis_names <- Reduce(
      intersect,
      list(
        names(analysis_arms_iptw),
        names(ipcw_denominator_models),
        names(ipcw_numerator_models)
      )
    )
    
    ipcw_common_analysis_names
    
    
    ##############################################################################
    ################ CALCULATE IPCW ##############################################
    ##############################################################################
    
    ipcw_weight_runs <- safe_imap(
      analysis_arms_iptw[ipcw_common_analysis_names],
      function(dat, analysis_name) {
        
        add_ipcw_censoring_weights(
          data = dat,
          denominator_model =
            ipcw_denominator_models[[analysis_name]],
          numerator_model =
            ipcw_numerator_models[[analysis_name]]
        )
      }
    )
    
    
    ipcw_weight_status <- get_status_table(
      ipcw_weight_runs
    )
    
    ipcw_weight_status
    
    
    analysis_arms_ipcw <- get_success_results(
      ipcw_weight_runs
    )
    
    names(analysis_arms_ipcw)
    
    
  } else {
    
    ipcw_weight_status <- tibble::tibble(
      note = "IPCW not used; censoring weights set to 1."
    )
    
    analysis_arms_ipcw <- purrr::map(
      analysis_arms_iptw,
      add_no_censoring_weights
    )
  }
  
  ################################################################################
  ######################## IPCW SUMMARIES ########################################
  ################################################################################
  
  ipcw_weight_summaries <- purrr::map(
    analysis_arms_ipcw,
    summarise_ipcw_censoring_weights
  )
  
  print(
    ipcw_weight_summaries$zero_as_dose,
    width = Inf
  )
  
 
  
  
  ################################################################################
  ######################## STEP 3: TOTAL WEIGHTS #################################
  ################################################################################
  
  total_weight_runs <- safe_imap(
    analysis_arms_ipcw,
    function(dat, analysis_name) {
      
      dat |>
        add_total_weights() |>
        truncate_total_weights(
          lower = TRUNCATION[1],
          upper = TRUNCATION[2]
        )
    }
  )
  
  total_weight_status <- get_status_table(
    total_weight_runs
  )
  
  total_weight_status
  
  
  analysis_arms_weighted <- get_success_results(
    total_weight_runs
  )
  
  names(analysis_arms_weighted)
  
  ################################################################################
  ######################## TOTAL WEIGHT SUMMARIES ################################
  ################################################################################
  
  total_weight_summaries <- purrr::map(
    analysis_arms_weighted,
    summarise_total_weights
  )
  
  print(
    total_weight_summaries$zero_as_dose,
    width = Inf
  )
  
 
  
  truncated_weight_summaries <- purrr::map(
    analysis_arms_weighted,
    summarise_total_truncated_weights
  )
  
  print(
    truncated_weight_summaries$zero_as_dose,
    width = Inf
  )
  

  
  truncation_checks <- purrr::map(
    analysis_arms_weighted,
    check_truncation
  )
  
  truncation_checks

  ################################################################################
  ######################## STEP 4: WEIGHTED MSM ##################################
  ################################################################################
  
  msm_runs <- safe_imap(
    analysis_arms_weighted,
    function(dat, analysis_name) {
      
      fit_weighted_msm(
        data = dat,
        weight_var = "SW_total_trunc",
        visit_df = VISIT_DF,
        corstr = "independence",
        include_dose_time_interaction =
          INCLUDE_DOSE_TIME_INTERACTION
      )
    }
  )
  
  
  ################################################################################
  ######################## MSM STATUS ############################################
  ################################################################################
  
  msm_status <- get_status_table(
    msm_runs
  )
  
  msm_status
  
  
  # Keep successfully fitted models.
  msm_models <- get_success_results(
    msm_runs
  )
  
  names(msm_models)
  
  ################################################################################
  ######################## MSM COEFFICIENT TABLES ################################
  ################################################################################
  
  msm_coef_tabs <- purrr::map(
    msm_models,
    make_gee_coef_table_robust_naive
  )
  
  
  print(
    msm_coef_tabs$zero_as_dose
  )
 

  
  #### UNTIL HERE
  ################################################################################
  ################ STEP 5A: COUNTERFACTUAL DOSE STRATEGIES #######################
  ################################################################################
  
  make_counterfactual_strategy_data <- function(
    model,
    analysis_data,
    strategy_doses = c(0, 20, 30, 40, 50),
    baseline_dose = 20
  ) {
    
    ##############################################################################
    # MSM data
    ##############################################################################
    
    model_dat <- attr(
      model,
      "model_data"
    )
    
    if (is.null(model_dat)) {
      stop(
        "The MSM does not contain the model_data attribute."
      )
    }
    
    
    # Rows that actually contributed to the MSM.
    msm_rows <- model_dat |>
      dplyr::distinct(
        pid,
        visit
      )
    
    
    ##############################################################################
    # Factor levels used in fitted MSM
    ##############################################################################
    
    dose_lag1_levels <-
      levels(
        model_dat$dose_lag1_f
      )
    
    dose_lag2_levels <-
      levels(
        model_dat$dose_lag2_f
      )
    
    dose_lag3_levels <-
      levels(
        model_dat$dose_lag3_f
      )
    
    study_levels <-
      levels(
        model_dat$studyid
      )
    
    
    ##############################################################################
    # Actual patient visit structure
    ##############################################################################
    
    base_dat <- analysis_data |>
      dplyr::filter(
        pid %in% model_dat$pid
      ) |>
      dplyr::arrange(
        pid,
        visit
      )
    
    
    ##############################################################################
    # Construct one counterfactual longitudinal dataset per sustained dose
    ##############################################################################
    
    strategy_dat <- purrr::map_dfr(
      strategy_doses,
      function(target_dose) {
        
        tmp <- base_dat |>
          dplyr::group_by(pid) |>
          dplyr::arrange(
            visit,
            .by_group = TRUE
          ) |>
          dplyr::mutate(
            
            visit_order =
              dplyr::row_number(),
            
            
            ######################################################################
            # Counterfactual current prescribed dose
            #
            # Everyone starts 20 mg at baseline.
            # Thereafter the assigned strategy is maintained.
            ######################################################################
            
            strategy_current_dose =
              dplyr::case_when(
                
                is_baseline ~
                  as.numeric(
                    baseline_dose
                  ),
                
                TRUE ~
                  as.numeric(
                    target_dose
                  )
              ),
            
            
            ######################################################################
            # Counterfactual dose history
            #
            # Before baseline patients were untreated = 0 mg.
            ######################################################################
            
            strategy_dose_lag1 =
              dplyr::case_when(
                
                visit_order == 1 ~
                  0,
                
                TRUE ~
                  dplyr::lag(
                    strategy_current_dose,
                    1
                  )
              ),
            
            
            strategy_dose_lag2 =
              dplyr::case_when(
                
                visit_order <= 2 ~
                  0,
                
                TRUE ~
                  dplyr::lag(
                    strategy_current_dose,
                    2
                  )
              ),
            
            
            strategy_dose_lag3 =
              dplyr::case_when(
                
                visit_order <= 3 ~
                  0,
                
                TRUE ~
                  dplyr::lag(
                    strategy_current_dose,
                    3
                  )
              ),
            
            
            ######################################################################
            # Older treatment history
            #
            # EXACTLY the same definition used in the fitted MSM.
            ######################################################################
            
            strategy_avg_dose_before_lag3 =
              vapply(
                seq_along(
                  strategy_current_dose
                ),
                function(j) {
                  
                  if (j <= 4) {
                    return(0)
                  }
                  
                  early_doses <-
                    strategy_current_dose[
                      seq_len(
                        j - 4
                      )
                    ]
                  
                  if (
                    all(
                      is.na(
                        early_doses
                      )
                    )
                  ) {
                    
                    return(0)
                    
                  } else {
                    
                    return(
                      mean(
                        early_doses,
                        na.rm = TRUE
                      )
                    )
                  }
                },
                numeric(1)
              )
          ) |>
          
          dplyr::ungroup() |>
          
          
          ########################################################################
        # Put counterfactual history into variables expected by MSM
        ########################################################################
        
        dplyr::mutate(
          
          strategy_dose =
            target_dose,
          
          strategy =
            paste0(
              target_dose,
              " mg"
            ),
          
          dose_lag1_f =
            factor(
              as.character(
                strategy_dose_lag1
              ),
              levels =
                dose_lag1_levels
            ),
          
          dose_lag2_f =
            factor(
              as.character(
                strategy_dose_lag2
              ),
              levels =
                dose_lag2_levels
            ),
          
          dose_lag3_f =
            factor(
              as.character(
                strategy_dose_lag3
              ),
              levels =
                dose_lag3_levels
            ),
          
          avg_dose_before_lag3 =
            strategy_avg_dose_before_lag3,
          
          studyid =
            factor(
              as.character(
                studyid
              ),
              levels =
                study_levels
            )
        ) |>
          
          
          ########################################################################
        # Keep only rows corresponding to MSM outcome rows
        ########################################################################
        
        dplyr::semi_join(
          msm_rows,
          by = c(
            "pid",
            "visit"
          )
        )
        
        tmp
      }
    )
    
    
    ##############################################################################
    # Safety check
    ##############################################################################
    
    if (
      any(is.na(strategy_dat$dose_lag1_f)) ||
      any(is.na(strategy_dat$dose_lag2_f)) ||
      any(is.na(strategy_dat$dose_lag3_f))
    ) {
      
      stop(
        paste0(
          "At least one counterfactual dose-history value ",
          "is absent from the fitted MSM factor levels."
        )
      )
    }
    
    
    strategy_dat
  }
  
  
  ################################################################################
  # Create strategies
  ################################################################################
  
  strategy_data_zero_as_dose <-
    make_counterfactual_strategy_data(
      
      model =
        msm_models$zero_as_dose,
      
      analysis_data =
        analysis_zero_as_dose,
      
      strategy_doses =
        STRATEGY_DOSES_ZERO_AS_DOSE,
      
      baseline_dose =
        BASELINE_DOSE_FOR_PREDICTIONS
    )
  
  
  ################################################################################
  # Check strategy datasets
  ################################################################################
  
  strategy_data_zero_as_dose |>
    dplyr::count(
      strategy
    )
  
  
  strategy_data_zero_as_dose |>
    dplyr::group_by(
      strategy
    ) |>
    dplyr::summarise(
      
      n_rows =
        dplyr::n(),
      
      n_patients =
        dplyr::n_distinct(
          pid
        ),
      
      min_visit =
        min(visit),
      
      max_visit =
        max(visit),
      
      .groups = "drop"
    )
  
  
  ################################################################################
  ######## STEP 5B: MSM PREDICTIONS AT ACTUAL CONTINUOUS VISIT TIMES #############
  ################################################################################
  
  
  ################################################################################
  # Build MSM design matrix and obtain robust covariance
  ################################################################################
  
  make_msm_design_matrix <- function(
    model,
    newdata
  ) {
    
    beta <-
      stats::coef(
        model
      )
    
    # Robust sandwich covariance matrix.
    V <-
      model$geese$vbeta
    
    dimnames(V) <- list(
      names(beta),
      names(beta)
    )
    
    
    ##############################################################################
    # Same design matrix as fitted MSM
    ##############################################################################
    
    Terms <-
      stats::delete.response(
        stats::terms(
          model
        )
      )
    
    X <-
      stats::model.matrix(
        Terms,
        data = newdata,
        contrasts.arg = model$contrasts,
        xlev = model$xlevels
      )
    
    
    ##############################################################################
    # Align matrix columns with fitted coefficients
    ##############################################################################
    
    missing_cols <-
      setdiff(
        names(beta),
        colnames(X)
      )
    
    if (
      length(
        missing_cols
      ) > 0
    ) {
      
      X <- cbind(
        X,
        matrix(
          0,
          nrow = nrow(X),
          ncol = length(
            missing_cols
          ),
          dimnames = list(
            NULL,
            missing_cols
          )
        )
      )
    }
    
    
    X <- X[
      ,
      names(beta),
      drop = FALSE
    ]
    
    
    ##############################################################################
    # Verify that matrix prediction reproduces predict.geeglm()
    ##############################################################################
    
    pred_from_X <-
      as.numeric(
        X %*% beta
      )
    
    pred_from_model <-
      as.numeric(
        stats::predict(
          model,
          newdata = newdata,
          type = "response"
        )
      )
    
    max_prediction_difference <-
      max(
        abs(
          pred_from_X -
            pred_from_model
        ),
        na.rm = TRUE
      )
    
    
    if (
      is.finite(
        max_prediction_difference
      ) &&
      max_prediction_difference > 1e-6
    ) {
      
      warning(
        paste0(
          "Model-matrix predictions differ from predict(). ",
          "Maximum absolute difference = ",
          signif(
            max_prediction_difference,
            4
          )
        )
      )
    }
    
    
    list(
      
      X =
        X,
      
      beta =
        beta,
      
      V =
        V,
      
      max_prediction_difference =
        max_prediction_difference
    )
  }
  
  
  ################################################################################
  # Standardized strategy predictions at each observed continuous visit time
  ################################################################################
  
  estimate_strategy_means_by_visit <- function(
    model,
    strategy_data
  ) {
    
    strategy_data <- strategy_data |>
      dplyr::mutate(
        .row_id =
          dplyr::row_number()
      )
    
    
    design <-
      make_msm_design_matrix(
        model = model,
        newdata = strategy_data
      )
    
    X <-
      design$X
    
    beta <-
      design$beta
    
    V <-
      design$V
    
    
    ##############################################################################
    # Strategy × exact observed visit-time groups
    ##############################################################################
    
    groups <- strategy_data |>
      dplyr::group_by(
        strategy_dose,
        strategy,
        visit
      ) |>
      dplyr::summarise(
        
        row_ids =
          list(
            .row_id
          ),
        
        n_rows =
          dplyr::n(),
        
        n_patients =
          dplyr::n_distinct(
            pid
          ),
        
        .groups = "drop"
      )
    
    
    ##############################################################################
    # Marginal mean and robust 95% CI
    ##############################################################################
    
    estimates <- purrr::pmap_dfr(
      groups,
      function(
    strategy_dose,
    strategy,
    visit,
    row_ids,
    n_rows,
    n_patients
      ) {
        
        X_this <-
          X[
            row_ids,
            ,
            drop = FALSE
          ]
        
        X_bar <-
          colMeans(
            X_this
          )
        
        estimate <-
          as.numeric(
            X_bar %*%
              beta
          )
        
        variance <-
          as.numeric(
            t(X_bar) %*%
              V %*%
              X_bar
          )
        
        se <-
          sqrt(
            pmax(
              variance,
              0
            )
          )
        
        
        tibble::tibble(
          
          strategy_dose =
            strategy_dose,
          
          strategy =
            strategy,
          
          visit =
            visit,
          
          n_rows =
            n_rows,
          
          n_patients =
            n_patients,
          
          predicted_improvement =
            estimate,
          
          robust_se =
            se,
          
          lower_95 =
            estimate -
            stats::qnorm(0.975) * se,
          
          upper_95 =
            estimate +
            stats::qnorm(0.975) * se
        )
      }
    )
    
    
    attr(
      estimates,
      "max_prediction_difference"
    ) <-
      design$max_prediction_difference
    
    
    estimates |>
      dplyr::arrange(
        strategy_dose,
        visit
      )
  }
  
  
  ################################################################################
  # Run trajectory predictions
  ################################################################################
  
  strategy_predictions_by_visit <-
    estimate_strategy_means_by_visit(
      
      model =
        msm_models$zero_as_dose,
      
      strategy_data =
        strategy_data_zero_as_dose
    )
  
  
  ################################################################################
  # Validation
  ################################################################################
  
  attr(
    strategy_predictions_by_visit,
    "max_prediction_difference"
  )
  
  
  ################################################################################
  # Inspect predictions
  ################################################################################
  
  print(
    strategy_predictions_by_visit,
    n = 100,
    width = Inf
  )

  ################################################################################
  ################ STEP 5C: FINAL PREDICTIONS + PLACEBO ###########################
  ################################################################################
  
  # Same role as Step 3 in the Flexible-dose-response-models GitHub analysis:
  #
  # 1. Predicted HAMD improvement under sustained dose strategies
  # 2. Observed placebo HAMD improvement
  # 3. Dose-strategy vs placebo table
  # 4. Active MSM + placebo means table
  # 5. Final trajectory graph with placebo
  
  
  ################################################################################
  ######################## PREDICTION WEEKS ######################################
  ################################################################################
  
  PREDICTION_WEEKS <- 1:VISIT_MAX
  
  
  ################################################################################
  ######## CREATE STANDARDIZED STRATEGY DATA AT A TARGET WEEK ####################
  ################################################################################
  
  make_target_strategy_data <- function(
    model,
    analysis_data,
    strategy_doses = c(0, 20, 30, 40, 50),
    baseline_dose = 20,
    target_visit
  ) {
    
    ##############################################################################
    # MSM population
    ##############################################################################
    
    model_dat <- attr(
      model,
      "model_data"
    )
    
    if (is.null(model_dat)) {
      stop(
        "The MSM does not contain the model_data attribute."
      )
    }
    
    
    ##############################################################################
    # One row per MSM patient
    #
    # Keep actual:
    # - patient
    # - study
    # - baseline HAMD
    ##############################################################################
    
    patient_dat <- model_dat |>
      dplyr::group_by(pid) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::select(
        pid,
        studyid,
        outcome_0
      )
    
    
    ##############################################################################
    # Actual treatment-decision schedule before target time
    #
    # We retain actual continuous visit times.
    ##############################################################################
    
    schedule_dat <- analysis_data |>
      dplyr::filter(
        pid %in% patient_dat$pid,
        visit < target_visit
      ) |>
      dplyr::select(
        pid,
        visit,
        is_baseline
      ) |>
      dplyr::distinct(
        pid,
        visit,
        .keep_all = TRUE
      ) |>
      dplyr::left_join(
        patient_dat,
        by = "pid"
      ) |>
      dplyr::mutate(
        target_row = FALSE
      )
    
    
    ##############################################################################
    # Add one synthetic outcome row at the target week for every patient
    ##############################################################################
    
    target_rows <- patient_dat |>
      dplyr::mutate(
        
        visit =
          as.numeric(target_visit),
        
        is_baseline =
          FALSE,
        
        target_row =
          TRUE
      )
    
    
    ##############################################################################
    # Factor levels required by the fitted MSM
    ##############################################################################
    
    dose_lag1_levels <-
      levels(
        model_dat$dose_lag1_f
      )
    
    dose_lag2_levels <-
      levels(
        model_dat$dose_lag2_f
      )
    
    dose_lag3_levels <-
      levels(
        model_dat$dose_lag3_f
      )
    
    study_levels <-
      levels(
        model_dat$studyid
      )
    
    
    ##############################################################################
    # Construct sustained treatment strategies
    ##############################################################################
    
    target_dat <- purrr::map_dfr(
      strategy_doses,
      function(target_dose) {
        
        tmp <- dplyr::bind_rows(
          schedule_dat,
          target_rows
        ) |>
          
          dplyr::group_by(pid) |>
          
          dplyr::arrange(
            visit,
            .by_group = TRUE
          ) |>
          
          dplyr::mutate(
            
            visit_order =
              dplyr::row_number(),
            
            
            ######################################################################
            # Treatment strategy
            #
            # Before study = 0 mg
            # Baseline = 20 mg
            # Postbaseline = sustained target dose
            ######################################################################
            
            strategy_current_dose =
              dplyr::case_when(
                
                is_baseline ~
                  as.numeric(
                    baseline_dose
                  ),
                
                TRUE ~
                  as.numeric(
                    target_dose
                  )
              ),
            
            
            ######################################################################
            # Dose lags
            ######################################################################
            
            strategy_dose_lag1 =
              dplyr::case_when(
                
                visit_order == 1 ~
                  0,
                
                TRUE ~
                  dplyr::lag(
                    strategy_current_dose,
                    1
                  )
              ),
            
            
            strategy_dose_lag2 =
              dplyr::case_when(
                
                visit_order <= 2 ~
                  0,
                
                TRUE ~
                  dplyr::lag(
                    strategy_current_dose,
                    2
                  )
              ),
            
            
            strategy_dose_lag3 =
              dplyr::case_when(
                
                visit_order <= 3 ~
                  0,
                
                TRUE ~
                  dplyr::lag(
                    strategy_current_dose,
                    3
                  )
              ),
            
            
            ######################################################################
            # Older dose history
            #
            # EXACTLY the same definition as in the fitted MSM.
            ######################################################################
            
            strategy_avg_dose_before_lag3 =
              vapply(
                seq_along(
                  strategy_current_dose
                ),
                function(j) {
                  
                  if (j <= 4) {
                    return(0)
                  }
                  
                  early_doses <-
                    strategy_current_dose[
                      seq_len(
                        j - 4
                      )
                    ]
                  
                  if (
                    all(
                      is.na(
                        early_doses
                      )
                    )
                  ) {
                    
                    return(0)
                    
                  }
                  
                  mean(
                    early_doses,
                    na.rm = TRUE
                  )
                },
                numeric(1)
              )
          ) |>
          
          dplyr::ungroup() |>
          
          
          ########################################################################
        # Keep only synthetic target-time outcome row
        ########################################################################
        
        dplyr::filter(
          target_row
        ) |>
          
          
          ########################################################################
        # Put histories into MSM variable names
        ########################################################################
        
        dplyr::mutate(
          
          strategy_dose =
            target_dose,
          
          strategy =
            paste0(
              target_dose,
              " mg"
            ),
          
          dose_lag1_f =
            factor(
              as.character(
                strategy_dose_lag1
              ),
              levels =
                dose_lag1_levels
            ),
          
          dose_lag2_f =
            factor(
              as.character(
                strategy_dose_lag2
              ),
              levels =
                dose_lag2_levels
            ),
          
          dose_lag3_f =
            factor(
              as.character(
                strategy_dose_lag3
              ),
              levels =
                dose_lag3_levels
            ),
          
          avg_dose_before_lag3 =
            strategy_avg_dose_before_lag3,
          
          studyid =
            factor(
              as.character(
                studyid
              ),
              levels =
                study_levels
            )
        )
        
        
        tmp
      }
    )
    
    
    ##############################################################################
    # Safety check
    ##############################################################################
    
    if (
      any(is.na(target_dat$dose_lag1_f)) ||
      any(is.na(target_dat$dose_lag2_f)) ||
      any(is.na(target_dat$dose_lag3_f))
    ) {
      
      stop(
        "Counterfactual history contains a dose level absent from the fitted MSM."
      )
    }
    
    
    target_dat
  }
  
  
  ################################################################################
  ######## CREATE STANDARDIZED DATA FOR WEEKS 1-6 ################################
  ################################################################################
  
  weekly_strategy_data <- purrr::map_dfr(
    PREDICTION_WEEKS,
    function(target_week) {
      
      make_target_strategy_data(
        
        model =
          msm_models$zero_as_dose,
        
        analysis_data =
          analysis_zero_as_dose,
        
        strategy_doses =
          STRATEGY_DOSES_ZERO_AS_DOSE,
        
        baseline_dose =
          BASELINE_DOSE_FOR_PREDICTIONS,
        
        target_visit =
          target_week
      )
    }
  )
  
  
  ################################################################################
  # Check
  #
  # Each strategy should have:
  # 416 patients x 6 weeks = 2496 rows
  ################################################################################
  
  weekly_strategy_data |>
    dplyr::count(
      strategy
    )
  
  
  ################################################################################
  ######## STANDARDIZED ACTIVE MSM PREDICTIONS ###################################
  ################################################################################
  
  # We already created this function in Step 5B:
  #
  # estimate_strategy_means_by_visit()
  #
  # It averages the MSM design matrix over the empirical patients/studies
  # and uses the robust GEE covariance matrix.
  
  weekly_strategy_predictions <-
    estimate_strategy_means_by_visit(
      
      model =
        msm_models$zero_as_dose,
      
      strategy_data =
        weekly_strategy_data
    )
  
  
  ################################################################################
  # Check prediction calculation
  ################################################################################
  
  attr(
    weekly_strategy_predictions,
    "max_prediction_difference"
  )
  
  
  ################################################################################
  # Active strategy predictions
  ################################################################################
  
  weekly_strategy_predictions

  ################################################################################
  ############ STEP 5D: STUDY-STANDARDIZED PLACEBO GEE ###########################
  ################################################################################
  
  # The active MSM curves are standardized over the 416 patients
  # represented in the active MSM.
  #
  # We fit the placebo trajectory using ALL observed placebo data
  # at their actual continuous visit times.
  #
  # Then placebo is predicted at weeks 1-6 for the SAME empirical
  # distribution of:
  #
  #   - study
  #   - baseline HAMD
  #
  # as the active MSM population.
  #
  # This avoids sparse exact-week placebo cells.
  
  
  ################################################################################
  ######################## PREPARE PLACEBO DATA ##################################
  ################################################################################
  
  placebo_model_dat <-
    placebo_data |>
    
    dplyr::mutate(
      
      pid =
        as.character(pid),
      
      studyid =
        as.factor(studyid),
      
      visit =
        as.numeric(visit),
      
      outcome =
        as.numeric(outcome)
    ) |>
    
    dplyr::arrange(
      studyid,
      pid,
      visit
    ) |>
    
    dplyr::group_by(
      studyid,
      pid
    ) |>
    
    dplyr::mutate(
      
      ##########################################################################
      # Baseline HAMD
      ##########################################################################
      
      outcome_0 =
        outcome[
          which.min(visit)
        ],
      
      ##########################################################################
      # Positive value = HAMD improvement
      ##########################################################################
      
      delta_outcome =
        outcome_0 -
        outcome
    ) |>
    
    dplyr::ungroup() |>
    
    dplyr::mutate(
      
      ##########################################################################
      # Unique placebo patient ID across studies
      ##########################################################################
      
      placebo_pid =
        factor(
          paste(
            studyid,
            pid,
            sep = "__"
          )
        )
    ) |>
    
    dplyr::filter(
      
      visit > 0,
      
      visit <= VISIT_MAX,
      
      !is.na(
        delta_outcome
      ),
      
      !is.na(
        outcome_0
      ),
      
      !is.na(
        studyid
      )
    )
  
  
  ################################################################################
  ######################## FIT PLACEBO GEE #######################################
  ################################################################################
  
  # Mirrors the non-treatment part of the active MSM:
  #
  #   nonlinear continuous time
  #   baseline HAMD
  #   baseline HAMD x time
  #   study fixed effect
  #
  # No treatment variables are needed because everyone here is placebo.
  
  placebo_gee <-
    geepack::geeglm(
      
      delta_outcome ~
        rms::rcs(
          visit,
          VISIT_DF
        ) * outcome_0 +
        studyid,
      
      id =
        placebo_pid,
      
      waves =
        visit,
      
      data =
        placebo_model_dat,
      
      family =
        gaussian(
          link = "identity"
        ),
      
      corstr =
        "independence",
      
      std.err =
        "san.se"
    )
  
  
  ################################################################################
  ######################## PLACEBO MODEL TABLE ###################################
  ################################################################################
  
  placebo_gee_coef_tab <-
    make_gee_coef_table_robust_naive(
      placebo_gee
    )
  
  placebo_gee_coef_tab
  
  
  ################################################################################
  ################ ACTIVE STANDARDIZATION POPULATION #############################
  ################################################################################
  
  # Exactly the same 416 patients used to standardize the active MSM predictions.
  
  active_standardization_population <-
    attr(
      msm_models$zero_as_dose,
      "model_data"
    ) |>
    
    dplyr::group_by(pid) |>
    
    dplyr::slice_head(
      n = 1
    ) |>
    
    dplyr::ungroup() |>
    
    dplyr::select(
      pid,
      studyid,
      outcome_0
    )
  
  
  ################################################################################
  # Check
  ################################################################################
  
  active_standardization_population |>
    dplyr::summarise(
      
      n_patients =
        dplyr::n_distinct(
          pid
        )
    )
  
  
  ################################################################################
  ################ PLACEBO PREDICTION DATA: WEEKS 1-6 ############################
  ################################################################################
  
  placebo_prediction_data <-
    tidyr::crossing(
      
      active_standardization_population,
      
      week =
        PREDICTION_WEEKS
    ) |>
    
    dplyr::mutate(
      
      visit =
        as.numeric(
          week
        ),
      
      ##########################################################################
      # Match factor levels used in placebo model
      ##########################################################################
      
      studyid =
        factor(
          as.character(
            studyid
          ),
          levels =
            levels(
              placebo_model_dat$studyid
            )
        )
    )
  
  
  ################################################################################
  # Safety check:
  # every active study must also have placebo data
  ################################################################################
  
  if (
    any(
      is.na(
        placebo_prediction_data$studyid
      )
    )
  ) {
    
    stop(
      paste0(
        "At least one study represented in the active MSM ",
        "has no corresponding placebo data."
      )
    )
  }
  
  
  ################################################################################
  ################ BUILD PLACEBO DESIGN MATRIX ###################################
  ################################################################################
  
  # The function created in Step 5B is generic enough to use for this GEE too.
  
  placebo_design <-
    make_msm_design_matrix(
      
      model =
        placebo_gee,
      
      newdata =
        placebo_prediction_data
    )
  
  
  X_placebo <-
    placebo_design$X
  
  beta_placebo <-
    placebo_design$beta
  
  V_placebo <-
    placebo_design$V
  
  
  ################################################################################
  ################ STANDARDIZED PLACEBO MEAN BY WEEK #############################
  ################################################################################
  
  placebo_prediction_data <-
    placebo_prediction_data |>
    
    dplyr::mutate(
      .row_id =
        dplyr::row_number()
    )
  
  
  placebo_week_groups <-
    placebo_prediction_data |>
    
    dplyr::group_by(
      week
    ) |>
    
    dplyr::summarise(
      
      row_ids =
        list(
          .row_id
        ),
      
      n_standardized_patients =
        dplyr::n_distinct(
          pid
        ),
      
      .groups =
        "drop"
    )
  
  
  placebo_standardized_weekly <-
    purrr::map_dfr(
      
      seq_len(
        nrow(
          placebo_week_groups
        )
      ),
      
      function(i) {
        
        ids <-
          placebo_week_groups$row_ids[[i]]
        
        
        X_this <-
          X_placebo[
            ids,
            ,
            drop = FALSE
          ]
        
        
        ##########################################################################
        # Average design vector over SAME 416 active MSM patients
        ##########################################################################
        
        X_bar <-
          colMeans(
            X_this
          )
        
        
        placebo_estimate <-
          as.numeric(
            X_bar %*%
              beta_placebo
          )
        
        
        placebo_variance <-
          as.numeric(
            t(
              X_bar
            ) %*%
              V_placebo %*%
              X_bar
          )
        
        
        placebo_se <-
          sqrt(
            pmax(
              placebo_variance,
              0
            )
          )
        
        
        tibble::tibble(
          
          week =
            placebo_week_groups$week[i],
          
          n_standardized_patients =
            placebo_week_groups$n_standardized_patients[i],
          
          placebo_mean_delta =
            placebo_estimate,
          
          placebo_se_delta =
            placebo_se,
          
          placebo_lower_95 =
            placebo_estimate -
            stats::qnorm(
              0.975
            ) *
            placebo_se,
          
          placebo_upper_95 =
            placebo_estimate +
            stats::qnorm(
              0.975
            ) *
            placebo_se
        )
      }
    )
  
  
  ################################################################################
  ################ ADD DESCRIPTIVE PLACEBO MODEL SAMPLE SIZE #####################
  ################################################################################
  
  n_placebo_model_patients <-
    dplyr::n_distinct(
      placebo_model_dat$placebo_pid
    )
  
  
  n_placebo_model_rows <-
    nrow(
      placebo_model_dat
    )
  
  
  placebo_standardized_weekly <-
    placebo_standardized_weekly |>
    
    dplyr::mutate(
      
      placebo_model_n_patients =
        n_placebo_model_patients,
      
      placebo_model_n_rows =
        n_placebo_model_rows
    )
  
  
  placebo_standardized_weekly
  
  
  ################################################################################
  # Validation:
  # model-matrix prediction must reproduce predict.geeglm()
  ################################################################################
  
  placebo_design$max_prediction_difference
  
  ################################################################################
  ######## TABLE 1: DOSE STRATEGIES VS STANDARDIZED PLACEBO ######################
  ################################################################################
  
  dose_vs_placebo_predictions_all <-
    weekly_strategy_predictions |>
    
    dplyr::mutate(
      
      week =
        as.integer(
          round(
            visit
          )
        )
    ) |>
    
    dplyr::left_join(
      placebo_standardized_weekly,
      by = "week"
    ) |>
    
    dplyr::mutate(
      
      ##########################################################################
      # Positive difference = more improvement under active strategy
      ##########################################################################
      
      dose_vs_placebo_difference =
        predicted_improvement -
        placebo_mean_delta,
      
      
      ##########################################################################
      # Approximate SE of active minus placebo
      ##########################################################################
      
      SE_difference =
        sqrt(
          robust_se^2 +
            placebo_se_delta^2
        ),
      
      
      difference_lower_95 =
        dose_vs_placebo_difference -
        stats::qnorm(
          0.975
        ) *
        SE_difference,
      
      
      difference_upper_95 =
        dose_vs_placebo_difference +
        stats::qnorm(
          0.975
        ) *
        SE_difference
    ) |>
    
    dplyr::transmute(
      
      week =
        week,
      
      dose_strategy =
        strategy,
      
      strategy_dose =
        strategy_dose,
      
      
      ##########################################################################
      # Standardization population
      ##########################################################################
      
      active_standardization_n =
        n_patients,
      
      placebo_standardization_n =
        n_standardized_patients,
      
      
      ##########################################################################
      # Active MSM
      ##########################################################################
      
      predicted_active_improvement =
        predicted_improvement,
      
      SE_active_prediction =
        robust_se,
      
      active_lower_95 =
        lower_95,
      
      active_upper_95 =
        upper_95,
      
      
      ##########################################################################
      # Placebo GEE
      ##########################################################################
      
      standardized_placebo_improvement =
        placebo_mean_delta,
      
      SE_standardized_placebo =
        placebo_se_delta,
      
      placebo_lower_95 =
        placebo_lower_95,
      
      placebo_upper_95 =
        placebo_upper_95,
      
      
      ##########################################################################
      # Difference
      ##########################################################################
      
      dose_vs_placebo_difference =
        dose_vs_placebo_difference,
      
      SE_difference =
        SE_difference,
      
      difference_lower_95 =
        difference_lower_95,
      
      difference_upper_95 =
        difference_upper_95
    ) |>
    
    dplyr::arrange(
      strategy_dose,
      week
    )
  
  
  dose_vs_placebo_predictions_all
  ################################################################################
  ######## TABLE 2: ACTIVE MSM + STANDARDIZED PLACEBO MEANS ######################
  ################################################################################
  
  active_mean_table <-
    weekly_strategy_predictions |>
    
    dplyr::transmute(
      
      week =
        as.integer(
          round(
            visit
          )
        ),
      
      source =
        "Active GEE/MSM prediction",
      
      dose_strategy =
        strategy,
      
      strategy_dose =
        strategy_dose,
      
      mean_improvement =
        predicted_improvement,
      
      SE =
        robust_se,
      
      lower_95 =
        lower_95,
      
      upper_95 =
        upper_95,
      
      n_standardized =
        n_patients
    )
  
  
  placebo_mean_table <-
    placebo_standardized_weekly |>
    
    dplyr::transmute(
      
      week =
        week,
      
      source =
        "Study-standardized placebo GEE",
      
      dose_strategy =
        "Placebo",
      
      strategy_dose =
        NA_real_,
      
      mean_improvement =
        placebo_mean_delta,
      
      SE =
        placebo_se_delta,
      
      lower_95 =
        placebo_lower_95,
      
      upper_95 =
        placebo_upper_95,
      
      n_standardized =
        n_standardized_patients
    )
  
  
  gee_and_placebo_means_all <-
    dplyr::bind_rows(
      
      active_mean_table,
      
      placebo_mean_table
    ) |>
    
    dplyr::mutate(
      
      source_order =
        dplyr::if_else(
          source ==
            "Study-standardized placebo GEE",
          0L,
          1L
        )
    ) |>
    
    dplyr::arrange(
      week,
      source_order,
      strategy_dose
    ) |>
    
    dplyr::select(
      -source_order
    )
  
  
  gee_and_placebo_means_all
  
  ################################################################################
  ######## FINAL GRAPH: DOSE STRATEGIES + STANDARDIZED PLACEBO ###################
  ################################################################################
  
  strategy_plot_dat <-
    weekly_strategy_predictions |>
    dplyr::filter(
      strategy_dose %in% PLOT_STRATEGY_DOSES
    ) |>
    dplyr::mutate(
      dose = factor(
        as.character(
          strategy_dose
        ),
        levels = as.character(
          PLOT_STRATEGY_DOSES
        )
      )
    )
  
  strategy_prediction_plot <-
    ggplot2::ggplot(
      
      strategy_plot_dat,
      
      ggplot2::aes(
        
        x =
          visit,
        
        y =
          predicted_improvement,
        
        color =
          dose,
        
        group =
          dose
      )
    ) +
    
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4,
      color = "grey60"
    ) +
    
    
    ##############################################################################
  # Standardized placebo CI
  ##############################################################################
  
  ggplot2::geom_ribbon(
    
    data =
      placebo_standardized_weekly,
    
    ggplot2::aes(
      
      x =
        week,
      
      ymin =
        placebo_lower_95,
      
      ymax =
        placebo_upper_95
    ),
    
    inherit.aes =
      FALSE,
    
    fill =
      "grey70",
    
    alpha =
      0.30
  ) +
    
    
    ##############################################################################
  # Active MSM CIs
  ##############################################################################
  
  ggplot2::geom_ribbon(
    
    ggplot2::aes(
      
      ymin =
        lower_95,
      
      ymax =
        upper_95,
      
      fill =
        dose
    ),
    
    alpha =
      0.10,
    
    color =
      NA,
    
    show.legend =
      FALSE
  ) +
    
    
    ##############################################################################
  # Active strategies
  ##############################################################################
  
  ggplot2::geom_line(
    linewidth = 1
  ) +
    
    ggplot2::geom_point(
      size = 2
    ) +
    
    
    ##############################################################################
  # Placebo
  ##############################################################################
  
  ggplot2::geom_line(
    
    data =
      placebo_standardized_weekly,
    
    ggplot2::aes(
      
      x =
        week,
      
      y =
        placebo_mean_delta,
      
      linetype =
        "Study-standardized placebo"
    ),
    
    inherit.aes =
      FALSE,
    
    color =
      "black",
    
    linewidth =
      1
  ) +
    
    
    ggplot2::geom_point(
      
      data =
        placebo_standardized_weekly,
      
      ggplot2::aes(
        
        x =
          week,
        
        y =
          placebo_mean_delta,
        
        shape =
          "Study-standardized placebo"
      ),
      
      inherit.aes =
        FALSE,
      
      color =
        "black",
      
      size =
        2.8
    ) +
    
    
    ggplot2::scale_color_manual(
      values =
        DOSE_COLORS
    ) +
    
    ggplot2::scale_fill_manual(
      values =
        DOSE_COLORS
    ) +
    
    ggplot2::scale_linetype_manual(
      values =
        c(
          "Study-standardized placebo" =
            "longdash"
        )
    ) +
    
    ggplot2::scale_shape_manual(
      values =
        c(
          "Study-standardized placebo" =
            17
        )
    ) +
    
    ggplot2::scale_x_continuous(
      breaks =
        PREDICTION_WEEKS
    ) +
    
    ggplot2::scale_y_continuous(
      breaks =
        PREDICTION_Y_BREAKS
    ) +
    
    ggplot2::labs(
      
      x =
        "Weeks from baseline",
      
      y =
        "Predicted HAMD improvement from baseline",
      
      color =
        "Dose strategy (mg)",
      
      fill =
        "Dose strategy (mg)",
      
      linetype =
        NULL,
      
      shape =
        NULL,
      
      title =
        "Paroxetine dose strategies and standardized placebo",
      
      subtitle =
        paste0(
          "Active and placebo predictions standardized over the same ",
          dplyr::n_distinct(
            active_standardization_population$pid
          ),
          " patients' baseline HAMD and study distribution"
        )
    ) +
    
    ggplot2::coord_cartesian(
      ylim =
        PREDICTION_Y_LIMITS
    ) +
    
    ggplot2::theme_classic(
      base_size = 12
    ) +
    
    ggplot2::theme(
      
      plot.title =
        ggplot2::element_text(
          face = "bold"
        ),
      
      panel.grid.major.y =
        ggplot2::element_line(
          color = "grey90",
          linewidth = 0.3
        ),
      
      legend.position =
        "right"
    )
  
  
  strategy_prediction_plot
  
  
### UNTIL HERE

  ################################################################################
  ################ MULTINOMIAL IPTW SENSITIVITY ANALYSIS #########################
  ################################################################################
  
  if (isTRUE(RUN_MULTINOMIAL_SENSITIVITY)) {
    
    
    ##############################################################################
    ######## MULTINOMIAL STEP 1A: TREATMENT DENOMINATOR ##########################
    ##############################################################################
    
    multinom_denominator_runs <-
      safe_imap(
        
        analysis_arms_model,
        
        function(
    dat,
    analysis_name
        ) {
          
          fit_iptw_multinom_denominator_model(
            
            data =
              dat,
            
            visit_df =
              VISIT_DF,
            
            ref_dose =
              MULTINOMIAL_REF_DOSE
          )
        }
      )
    
    
    multinom_denominator_status <-
      get_status_table(
        multinom_denominator_runs
      )
    
    multinom_denominator_status
    
    
    multinom_denominator_models <-
      get_success_results(
        multinom_denominator_runs
      )
    
    
    names(
      multinom_denominator_models
    )
    
    
    ##############################################################################
    # Convergence check
    ##############################################################################
    
    multinom_denominator_convergence <-
      purrr::imap_dfr(
        
        multinom_denominator_models,
        
        function(
    model,
    analysis_name
        ) {
          
          tibble::tibble(
            
            analysis_name =
              analysis_name,
            
            convergence_code =
              model$convergence,
            
            converged =
              model$convergence == 0
          )
        }
      )
    
    
    multinom_denominator_convergence
    
    
    ##############################################################################
    # Denominator coefficient tables
    ##############################################################################
    
    multinom_denominator_coef_tabs <-
      purrr::map(
        
        multinom_denominator_models,
        
        make_multinom_coef_table
      )
    
    
    multinom_denominator_coef_tabs
    
    
    writexl::write_xlsx(
      
      multinom_denominator_coef_tabs,
      
      path =
        file.path(
          RESULTS_DIR,
          "multinom_denominator_coef_tabs.xlsx"
        )
    )
    
    
    ##############################################################################
    ######## MULTINOMIAL OBSERVED VS PREDICTED DOSE DIAGNOSTIC ##################
    ##############################################################################
    
    # Reuse the exact same diagnostic function as the ordinal analysis.
    # predict(..., type = "probs") also works for multinom.
    
    multinom_dose_model_plot_data <-
      purrr::imap(
        
        multinom_denominator_models,
        
        function(
    model,
    analysis_name
        ) {
          
          make_iptw_dose_model_plot_data(
            
            model =
              model,
            
            data =
              analysis_arms_model[[analysis_name]]
          )
        }
      )
    
    
    multinom_dose_model_plots <-
      purrr::imap(
        
        multinom_dose_model_plot_data,
        
        function(
    plot_dat,
    analysis_name
        ) {
          
          make_iptw_dose_model_plot(
            
            plot_dat =
              plot_dat,
            
            analysis_name =
              paste0(
                analysis_name,
                " - multinomial IPTW"
              )
          )
        }
      )
    
    
    multinom_dose_model_plots
    
    
    ##############################################################################
    ######## MULTINOMIAL STEP 1B: TREATMENT NUMERATOR ############################
    ##############################################################################
    
    multinom_numerator_runs <-
      safe_imap(
        
        analysis_arms_model,
        
        function(
    dat,
    analysis_name
        ) {
          
          fit_iptw_multinom_numerator_model(
            
            data =
              dat,
            
            visit_df =
              VISIT_DF,
            
            ref_dose =
              MULTINOMIAL_REF_DOSE
          )
        }
      )
    
    
    multinom_numerator_status <-
      get_status_table(
        multinom_numerator_runs
      )
    
    multinom_numerator_status
    
    
    multinom_numerator_models <-
      get_success_results(
        multinom_numerator_runs
      )
    
    
    multinom_numerator_convergence <-
      purrr::imap_dfr(
        
        multinom_numerator_models,
        
        function(
    model,
    analysis_name
        ) {
          
          tibble::tibble(
            
            analysis_name =
              analysis_name,
            
            convergence_code =
              model$convergence,
            
            converged =
              model$convergence == 0
          )
        }
      )
    
    
    multinom_numerator_convergence
    
    
    multinom_numerator_coef_tabs <-
      purrr::map(
        
        multinom_numerator_models,
        
        make_multinom_coef_table
      )
    
    
    multinom_numerator_coef_tabs
    
    
    writexl::write_xlsx(
      
      multinom_numerator_coef_tabs,
      
      path =
        file.path(
          RESULTS_DIR,
          "multinom_numerator_coef_tabs.xlsx"
        )
    )
    
    
    ##############################################################################
    ######## MULTINOMIAL STEP 1C: TREATMENT WEIGHTS ##############################
    ##############################################################################
    
    # IMPORTANT:
    #
    # Use analysis_arms_ipcw because the censoring weights have already
    # been correctly estimated.
    #
    # The SAME IPCW is used in ordinal and multinomial analyses.
    
    multinom_common_analysis_names <-
      Reduce(
        
        intersect,
        
        list(
          
          names(
            analysis_arms_model
          ),
          
          names(
            multinom_denominator_models
          ),
          
          names(
            multinom_numerator_models
          ),
          
          names(
            analysis_arms_ipcw
          )
        )
      )
    
    
    multinom_common_analysis_names
    
    
    multinom_weight_runs <-
      safe_imap(
        
        analysis_arms_ipcw[
          multinom_common_analysis_names
        ],
        
        function(
    dat,
    analysis_name
        ) {
          
          dat |>
            
            add_iptw_treatment_weights_multinom(
              
              denominator_model =
                multinom_denominator_models[[analysis_name]],
              
              numerator_model =
                multinom_numerator_models[[analysis_name]]
            ) |>
            
            add_total_weights_multinom() |>
            
            truncate_total_weights_multinom(
              
              lower =
                TRUNCATION[1],
              
              upper =
                TRUNCATION[2]
            )
        }
      )
    
    
    multinom_weight_status <-
      get_status_table(
        multinom_weight_runs
      )
    
    multinom_weight_status
    
    
    analysis_arms_multinom_weighted <-
      get_success_results(
        multinom_weight_runs
      )
    
    
    names(
      analysis_arms_multinom_weighted
    )
    
    
    ##############################################################################
    ################ MULTINOMIAL IPTW SUMMARY ####################################
    ##############################################################################
    
    multinom_weight_summaries <-
      purrr::map(
        
        analysis_arms_multinom_weighted,
        
        summarise_iptw_multinom_treatment_weights
      )
    
    
    print(
      multinom_weight_summaries$zero_as_dose,
      width = Inf
    )
    
    
    ##############################################################################
    ################ MULTINOMIAL POSITIVITY DIAGNOSTIC ###########################
    ##############################################################################
    
    multinom_probability_diagnostics <-
      purrr::imap_dfr(
        
        analysis_arms_multinom_weighted,
        
        function(
    dat,
    analysis_name
        ) {
          
          dat |>
            
            dplyr::filter(
              
              use_treatment_weight,
              
              !is.na(
                p_dose_denominator_multinom
              )
            ) |>
            
            dplyr::summarise(
              
              analysis_name =
                analysis_name,
              
              n =
                dplyr::n(),
              
              n_patients =
                dplyr::n_distinct(
                  pid
                ),
              
              min_p_denominator =
                min(
                  p_dose_denominator_multinom
                ),
              
              p1_p_denominator =
                as.numeric(
                  stats::quantile(
                    p_dose_denominator_multinom,
                    0.01
                  )
                ),
              
              p5_p_denominator =
                as.numeric(
                  stats::quantile(
                    p_dose_denominator_multinom,
                    0.05
                  )
                ),
              
              median_p_denominator =
                stats::median(
                  p_dose_denominator_multinom
                ),
              
              mean_p_denominator =
                mean(
                  p_dose_denominator_multinom
                ),
              
              n_p_below_0_01 =
                sum(
                  p_dose_denominator_multinom <
                    0.01
                ),
              
              percent_p_below_0_01 =
                100 *
                mean(
                  p_dose_denominator_multinom <
                    0.01
                ),
              
              n_p_below_0_05 =
                sum(
                  p_dose_denominator_multinom <
                    0.05
                ),
              
              percent_p_below_0_05 =
                100 *
                mean(
                  p_dose_denominator_multinom <
                    0.05
                )
            )
        }
      )
    
    
    print(
      multinom_probability_diagnostics,
      width = Inf
    )
    
    
    ##############################################################################
    ################ MULTINOMIAL TOTAL WEIGHTS ###################################
    ##############################################################################
    
    multinom_total_weight_summaries <-
      purrr::map(
        
        analysis_arms_multinom_weighted,
        
        summarise_total_multinom_weights
      )
    
    
    multinom_total_truncated_weight_summaries <-
      purrr::map(
        
        analysis_arms_multinom_weighted,
        
        summarise_total_multinom_truncated_weights
      )
    
    
    multinom_truncation_checks <-
      purrr::map(
        
        analysis_arms_multinom_weighted,
        
        check_multinom_truncation
      )
    
    
    print(
      multinom_total_weight_summaries$zero_as_dose,
      width = Inf
    )
    
    
    print(
      multinom_total_truncated_weight_summaries$zero_as_dose,
      width = Inf
    )
    
    
    multinom_truncation_checks
    
    
    ##############################################################################
    ######## MULTINOMIAL STEP 2: SAME WEIGHTED MSM ###############################
    ##############################################################################
    
    # The OUTCOME MODEL is exactly the same.
    #
    # Only the treatment weights differ.
    
    multinom_msm_runs <-
      safe_imap(
        
        analysis_arms_multinom_weighted,
        
        function(
    dat,
    analysis_name
        ) {
          
          fit_weighted_msm(
            
            data =
              dat,
            
            weight_var =
              "SW_total_multinom_trunc",
            
            visit_df =
              VISIT_DF,
            
            corstr =
              "independence",
            
            include_dose_time_interaction =
              INCLUDE_DOSE_TIME_INTERACTION
          )
        }
      )
    
    
    multinom_msm_status <-
      get_status_table(
        multinom_msm_runs
      )
    
    multinom_msm_status
    
    
    multinom_msm_models <-
      get_success_results(
        multinom_msm_runs
      )
    
    
    names(
      multinom_msm_models
    )
    
    
    multinom_msm_coef_tabs <-
      purrr::map(
        
        multinom_msm_models,
        
        make_gee_coef_table_robust_naive
      )
    
    
    print(
      multinom_msm_coef_tabs$zero_as_dose
    )
    
    
    ##############################################################################
    ######## CHECK THAT ORDINAL AND MULTINOMIAL USE SAME MSM POPULATION ##########
    ##############################################################################
    
    ordinal_model_dat <-
      attr(
        msm_models$zero_as_dose,
        "model_data"
      )
    
    
    multinom_model_dat <-
      attr(
        multinom_msm_models$zero_as_dose,
        "model_data"
      )
    
    
    ordinal_pid_set <-
      sort(
        unique(
          as.character(
            ordinal_model_dat$pid
          )
        )
      )
    
    
    multinom_pid_set <-
      sort(
        unique(
          as.character(
            multinom_model_dat$pid
          )
        )
      )
    
    
    ordinal_row_set <-
      paste(
        ordinal_model_dat$pid,
        ordinal_model_dat$visit,
        sep = "__"
      )
    
    
    multinom_row_set <-
      paste(
        multinom_model_dat$pid,
        multinom_model_dat$visit,
        sep = "__"
      )
    
    
    multinom_population_check <-
      tibble::tibble(
        
        ordinal_n_patients =
          length(
            ordinal_pid_set
          ),
        
        multinom_n_patients =
          length(
            multinom_pid_set
          ),
        
        same_patient_set =
          setequal(
            ordinal_pid_set,
            multinom_pid_set
          ),
        
        ordinal_n_rows =
          nrow(
            ordinal_model_dat
          ),
        
        multinom_n_rows =
          nrow(
            multinom_model_dat
          ),
        
        same_msm_rows =
          setequal(
            ordinal_row_set,
            multinom_row_set
          )
      )
    
    
    multinom_population_check
    
    
    if (
      !multinom_population_check$same_patient_set ||
      !multinom_population_check$same_msm_rows
    ) {
      
      stop(
        paste0(
          "Ordinal and multinomial MSMs are not using the same population/rows. ",
          "Do not compare them until this is resolved."
        )
      )
    }
    
    
    ##############################################################################
    ######## MULTINOMIAL STEP 3: SAME STANDARDIZED STRATEGY PREDICTIONS ##########
    ##############################################################################
    
    # Reuse the FINAL prediction functions from the ordinal analysis:
    #
    # make_target_strategy_data()
    # estimate_strategy_means_by_visit()
    #
    # Same:
    # - 416 patients
    # - baseline HAMD
    # - study distribution
    # - empirical treatment-decision structure
    # - 20 mg at baseline
    # - weeks 1-6
    
    multinom_weekly_strategy_data <-
      purrr::map_dfr(
        
        PREDICTION_WEEKS,
        
        function(
    target_week
        ) {
          
          make_target_strategy_data(
            
            model =
              multinom_msm_models$zero_as_dose,
            
            analysis_data =
              analysis_zero_as_dose,
            
            strategy_doses =
              STRATEGY_DOSES_ZERO_AS_DOSE,
            
            baseline_dose =
              BASELINE_DOSE_FOR_PREDICTIONS,
            
            target_visit =
              target_week
          )
        }
      )
    
    
    ##############################################################################
    # Expected:
    # 416 patients x 6 weeks = 2496 rows per strategy.
    ##############################################################################
    
    multinom_weekly_strategy_data |>
      dplyr::count(
        strategy
      )
    
    
    multinom_weekly_strategy_predictions <-
      estimate_strategy_means_by_visit(
        
        model =
          multinom_msm_models$zero_as_dose,
        
        strategy_data =
          multinom_weekly_strategy_data
      )
    
    
    multinom_weekly_strategy_predictions
    
    
    attr(
      multinom_weekly_strategy_predictions,
      "max_prediction_difference"
    )
    
    
    ##############################################################################
    ######## MULTINOMIAL TABLE 1: DOSE STRATEGY VS SAME PLACEBO ##################
    ##############################################################################
    
    # IMPORTANT:
    #
    # Do NOT refit placebo.
    #
    # placebo_standardized_weekly is identical for the ordinal and
    # multinomial sensitivity analyses.
    #
    # That way the only difference is the IPTW treatment model.
    
    multinom_dose_vs_placebo_predictions_all <-
      multinom_weekly_strategy_predictions |>
      
      dplyr::mutate(
        
        week =
          as.integer(
            round(
              visit
            )
          )
      ) |>
      
      dplyr::left_join(
        
        placebo_standardized_weekly,
        
        by =
          "week"
      ) |>
      
      dplyr::mutate(
        
        dose_vs_placebo_difference =
          predicted_improvement -
          placebo_mean_delta,
        
        SE_difference =
          sqrt(
            robust_se^2 +
              placebo_se_delta^2
          ),
        
        difference_lower_95 =
          dose_vs_placebo_difference -
          stats::qnorm(
            0.975
          ) *
          SE_difference,
        
        difference_upper_95 =
          dose_vs_placebo_difference +
          stats::qnorm(
            0.975
          ) *
          SE_difference
      ) |>
      
      dplyr::transmute(
        
        week =
          week,
        
        dose_strategy =
          strategy,
        
        strategy_dose =
          strategy_dose,
        
        active_standardization_n =
          n_patients,
        
        placebo_standardization_n =
          n_standardized_patients,
        
        predicted_active_improvement =
          predicted_improvement,
        
        SE_active_prediction =
          robust_se,
        
        active_lower_95 =
          lower_95,
        
        active_upper_95 =
          upper_95,
        
        standardized_placebo_improvement =
          placebo_mean_delta,
        
        SE_standardized_placebo =
          placebo_se_delta,
        
        placebo_lower_95 =
          placebo_lower_95,
        
        placebo_upper_95 =
          placebo_upper_95,
        
        dose_vs_placebo_difference =
          dose_vs_placebo_difference,
        
        SE_difference =
          SE_difference,
        
        difference_lower_95 =
          difference_lower_95,
        
        difference_upper_95 =
          difference_upper_95
      ) |>
      
      dplyr::arrange(
        strategy_dose,
        week
      )
    
    
    multinom_dose_vs_placebo_predictions_all
    
    
    ##############################################################################
    ######## MULTINOMIAL TABLE 2: MSM + SAME STANDARDIZED PLACEBO ################
    ##############################################################################
    
    multinom_active_mean_table <-
      multinom_weekly_strategy_predictions |>
      
      dplyr::transmute(
        
        week =
          as.integer(
            round(
              visit
            )
          ),
        
        source =
          "Active multinomial-IPTW MSM prediction",
        
        dose_strategy =
          strategy,
        
        strategy_dose =
          strategy_dose,
        
        mean_improvement =
          predicted_improvement,
        
        SE =
          robust_se,
        
        lower_95 =
          lower_95,
        
        upper_95 =
          upper_95,
        
        n_standardized =
          n_patients
      )
    
    
    multinom_placebo_mean_table <-
      placebo_standardized_weekly |>
      
      dplyr::transmute(
        
        week =
          week,
        
        source =
          "Study-standardized placebo GEE",
        
        dose_strategy =
          "Placebo",
        
        strategy_dose =
          NA_real_,
        
        mean_improvement =
          placebo_mean_delta,
        
        SE =
          placebo_se_delta,
        
        lower_95 =
          placebo_lower_95,
        
        upper_95 =
          placebo_upper_95,
        
        n_standardized =
          n_standardized_patients
      )
    
    
    multinom_gee_and_placebo_means_all <-
      dplyr::bind_rows(
        
        multinom_active_mean_table,
        
        multinom_placebo_mean_table
      ) |>
      
      dplyr::mutate(
        
        source_order =
          dplyr::if_else(
            
            source ==
              "Study-standardized placebo GEE",
            
            0L,
            
            1L
          )
      ) |>
      
      dplyr::arrange(
        
        week,
        
        source_order,
        
        strategy_dose
      ) |>
      
      dplyr::select(
        -source_order
      )
    
    
    multinom_gee_and_placebo_means_all
    
    
    ##############################################################################
    ######## DIRECT ORDINAL VS MULTINOMIAL RESULTS TABLE #########################
    ##############################################################################
    
    compare_ordinal_vs_multinom_placebo <-
      dplyr::bind_rows(
        
        dose_vs_placebo_predictions_all |>
          dplyr::mutate(
            weight_model =
              "Ordinal IPTW + IPCW",
            .before = 1
          ),
        
        multinom_dose_vs_placebo_predictions_all |>
          dplyr::mutate(
            weight_model =
              "Multinomial IPTW + IPCW",
            .before = 1
          )
      ) |>
      
      dplyr::arrange(
        
        strategy_dose,
        
        week,
        
        weight_model
      )
    
    
    compare_ordinal_vs_multinom_placebo
    
    
    ##############################################################################
    # Compare full active/placebo mean tables
    ##############################################################################
    
    gee_and_placebo_means_compare <-
      dplyr::bind_rows(
        
        gee_and_placebo_means_all |>
          dplyr::mutate(
            weight_model =
              "Ordinal IPTW + IPCW",
            .before = 1
          ),
        
        multinom_gee_and_placebo_means_all |>
          dplyr::mutate(
            weight_model =
              "Multinomial IPTW + IPCW",
            .before = 1
          )
      )
    
    
    gee_and_placebo_means_compare
    
    
    ##############################################################################
    ######## MULTINOMIAL FINAL STRATEGY + PLACEBO GRAPH ##########################
    ##############################################################################
    
    multinom_strategy_plot_dat <-
      multinom_weekly_strategy_predictions |>
      dplyr::filter(
        strategy_dose %in% PLOT_STRATEGY_DOSES
      ) |>
      dplyr::mutate(
        dose = factor(
          as.character(strategy_dose),
          levels = as.character(PLOT_STRATEGY_DOSES)
        )
      )
    
    multinom_strategy_prediction_plot <-
      ggplot2::ggplot(
        
        multinom_strategy_plot_dat,
        
        ggplot2::aes(
          
          x =
            visit,
          
          y =
            predicted_improvement,
          
          color =
            dose,
          
          group =
            dose
        )
      ) +
      
      ggplot2::geom_hline(
        
        yintercept =
          0,
        
        linetype =
          "dashed",
        
        linewidth =
          0.4,
        
        color =
          "grey60"
      ) +
      
      
      ##########################################################################
    # SAME placebo curve
    ##########################################################################
    
    ggplot2::geom_ribbon(
      
      data =
        placebo_standardized_weekly,
      
      ggplot2::aes(
        
        x =
          week,
        
        ymin =
          placebo_lower_95,
        
        ymax =
          placebo_upper_95
      ),
      
      inherit.aes =
        FALSE,
      
      fill =
        "grey70",
      
      alpha =
        0.30
    ) +
      
      
      ##########################################################################
    # Multinomial MSM uncertainty
    ##########################################################################
    
    ggplot2::geom_ribbon(
      
      ggplot2::aes(
        
        ymin =
          lower_95,
        
        ymax =
          upper_95,
        
        fill =
          dose
      ),
      
      alpha =
        0.10,
      
      color =
        NA,
      
      show.legend =
        FALSE
    ) +
      
      
      ggplot2::geom_line(
        linewidth = 1
      ) +
      
      ggplot2::geom_point(
        size = 2
      ) +
      
      
      ##########################################################################
    # Placebo curve
    ##########################################################################
    
    ggplot2::geom_line(
      
      data =
        placebo_standardized_weekly,
      
      ggplot2::aes(
        
        x =
          week,
        
        y =
          placebo_mean_delta,
        
        linetype =
          "Study-standardized placebo"
      ),
      
      inherit.aes =
        FALSE,
      
      color =
        "black",
      
      linewidth =
        1
    ) +
      
      
      ggplot2::geom_point(
        
        data =
          placebo_standardized_weekly,
        
        ggplot2::aes(
          
          x =
            week,
          
          y =
            placebo_mean_delta,
          
          shape =
            "Study-standardized placebo"
        ),
        
        inherit.aes =
          FALSE,
        
        color =
          "black",
        
        size =
          2.8
      ) +
      
      
      ggplot2::scale_color_manual(
        values =
          DOSE_COLORS
      ) +
      
      ggplot2::scale_fill_manual(
        values =
          DOSE_COLORS
      ) +
      
      ggplot2::scale_linetype_manual(
        
        values =
          c(
            "Study-standardized placebo" =
              "longdash"
          )
      ) +
      
      ggplot2::scale_shape_manual(
        
        values =
          c(
            "Study-standardized placebo" =
              17
          )
      ) +
      
      ggplot2::scale_x_continuous(
        breaks =
          PREDICTION_WEEKS
      ) +
      
      ggplot2::scale_y_continuous(
        breaks =
          PREDICTION_Y_BREAKS
      ) +
      
      ggplot2::labs(
        
        x =
          "Weeks from baseline",
        
        y =
          "Predicted HAMD improvement from baseline",
        
        color =
          "Dose strategy (mg)",
        
        fill =
          "Dose strategy (mg)",
        
        linetype =
          NULL,
        
        shape =
          NULL,
        
        title =
          "Paroxetine dose strategies: multinomial IPTW sensitivity",
        
        subtitle =
          "Same MSM, IPCW, target population and standardized placebo as primary analysis"
      ) +
      
      ggplot2::coord_cartesian(
        ylim =
          PREDICTION_Y_LIMITS
      ) +
      
      ggplot2::theme_classic(
        base_size = 12
      ) +
      
      ggplot2::theme(
        
        plot.title =
          ggplot2::element_text(
            face = "bold"
          ),
        
        panel.grid.major.y =
          ggplot2::element_line(
            color = "grey90",
            linewidth = 0.3
          ),
        
        legend.position =
          "right"
      )
    
    
    multinom_strategy_prediction_plot
    
    
    ##############################################################################
    ######## DIRECT ORDINAL VS MULTINOMIAL ACTIVE-CURVE GRAPH ####################
    ##############################################################################
    
    ordinal_vs_multinom_plot_dat <-
      dplyr::bind_rows(
        weekly_strategy_predictions |>
          dplyr::mutate(
            weight_model = "Ordinal IPTW"
          ),
        multinom_weekly_strategy_predictions |>
          dplyr::mutate(
            weight_model = "Multinomial IPTW"
          )
      ) |>
      dplyr::filter(
        strategy_dose %in% PLOT_STRATEGY_DOSES
      ) |>
      dplyr::mutate(
        dose = factor(
          as.character(
            strategy_dose
          ),
          levels = as.character(
            PLOT_STRATEGY_DOSES
          )
        )
      )
    
    
    ordinal_vs_multinom_strategy_plot <-
      ggplot2::ggplot(
        
        ordinal_vs_multinom_plot_dat,
        
        ggplot2::aes(
          
          x =
            visit,
          
          y =
            predicted_improvement,
          
          color =
            dose,
          
          linetype =
            weight_model,
          
          group =
            interaction(
              dose,
              weight_model
            )
        )
      ) +
      
      ggplot2::geom_line(
        linewidth = 1
      ) +
      
      ggplot2::geom_point(
        size = 1.8
      ) +
      
      ggplot2::geom_line(
        
        data =
          placebo_standardized_weekly,
        
        ggplot2::aes(
          
          x =
            week,
          
          y =
            placebo_mean_delta
        ),
        
        inherit.aes =
          FALSE,
        
        color =
          "black",
        
        linewidth =
          1,
        
        linetype =
          "longdash"
      ) +
      
      ggplot2::scale_color_manual(
        values =
          DOSE_COLORS
      ) +
      
      ggplot2::scale_x_continuous(
        breaks =
          PREDICTION_WEEKS
      ) +
      
      ggplot2::scale_y_continuous(
        breaks =
          PREDICTION_Y_BREAKS
      ) +
      
      ggplot2::labs(
        
        x =
          "Weeks from baseline",
        
        y =
          "Predicted HAMD improvement from baseline",
        
        color =
          "Dose strategy (mg)",
        
        linetype =
          "Treatment-weight model",
        
        title =
          "Ordinal versus multinomial IPTW sensitivity",
        
        subtitle =
          "Black dashed line = study-standardized placebo"
      ) +
      
      ggplot2::coord_cartesian(
        ylim =
          PREDICTION_Y_LIMITS
      ) +
      
      ggplot2::theme_classic(
        base_size = 12
      ) +
      
      ggplot2::theme(
        
        plot.title =
          ggplot2::element_text(
            face = "bold"
          ),
        
        panel.grid.major.y =
          ggplot2::element_line(
            color = "grey90",
            linewidth = 0.3
          ),
        
        legend.position =
          "right"
      )
    
    
    ordinal_vs_multinom_strategy_plot
    
    
    ##############################################################################
    ################ SAVE MULTINOMIAL SENSITIVITY RESULTS ########################
    ##############################################################################
    
    writexl::write_xlsx(
      
      list(
        
        "denom_coef" =
          multinom_denominator_coef_tabs$zero_as_dose,
        
        "num_coef" =
          multinom_numerator_coef_tabs$zero_as_dose,
        
        "treat_weights" =
          multinom_weight_summaries$zero_as_dose,
        
        "positivity" =
          multinom_probability_diagnostics,
        
        "total_weights" =
          multinom_total_weight_summaries$zero_as_dose,
        
        "trunc_weights" =
          multinom_total_truncated_weight_summaries$zero_as_dose,
        
        "trunc_check" =
          multinom_truncation_checks$zero_as_dose,
        
        "msm_coef" =
          multinom_msm_coef_tabs$zero_as_dose,
        
        "dose_vs_placebo" =
          multinom_dose_vs_placebo_predictions_all,
        
        "means_placebo" =
          multinom_gee_and_placebo_means_all,
        
        "compare_dose_placebo" =
          compare_ordinal_vs_multinom_placebo,
        
        "compare_means" =
          gee_and_placebo_means_compare
        
      ),
      
      path =
        file.path(
          RESULTS_DIR,
          "multinomial_IPTW_sensitivity.xlsx"
        )
    )
    
    
    ##############################################################################
    ################ SAVE MULTINOMIAL PLOTS ######################################
    ##############################################################################
    
    if (isTRUE(SAVE_PLOTS)) {
      
      ggplot2::ggsave(
        
        filename =
          file.path(
            RESULTS_DIR,
            "multinomial_IPTW_strategy_placebo.jpg"
          ),
        
        plot =
          multinom_strategy_prediction_plot,
        
        width =
          10,
        
        height =
          7,
        
        dpi =
          300
      )
      
      
      ggplot2::ggsave(
        
        filename =
          file.path(
            RESULTS_DIR,
            "ordinal_vs_multinomial_IPTW.jpg"
          ),
        
        plot =
          ordinal_vs_multinom_strategy_plot,
        
        width =
          10,
        
        height =
          7,
        
        dpi =
          300
      )
      
      
      if (
        length(
          multinom_dose_model_plots
        ) > 0
      ) {
        
        save_plot_grid_jpg(
          
          plot_list =
            multinom_dose_model_plots,
          
          file_path =
            file.path(
              RESULTS_DIR,
              "multinom_observed_vs_predicted_dose.jpg"
            ),
          
          ncol =
            1
        )
      }
    }
    
  }
  
  ################################################################################
  ######## MULTINOMIAL FINAL GRAPH: DOSE STRATEGIES + PLACEBO ####################
  ################################################################################
  
  multinom_strategy_plot_dat <-
    multinom_weekly_strategy_predictions |>
    
    dplyr::mutate(
      
      dose =
        factor(
          as.character(
            strategy_dose
          ),
          levels =
            as.character(
              STRATEGY_DOSES_ZERO_AS_DOSE
            )
        )
    )
  
  
  multinom_strategy_prediction_plot <-
    ggplot2::ggplot(
      
      multinom_strategy_plot_dat,
      
      ggplot2::aes(
        
        x =
          visit,
        
        y =
          predicted_improvement,
        
        color =
          dose,
        
        group =
          dose
      )
    ) +
    
    
    ##############################################################################
  # Zero reference
  ##############################################################################
  
  ggplot2::geom_hline(
    
    yintercept =
      0,
    
    linetype =
      "dashed",
    
    linewidth =
      0.4,
    
    color =
      "grey60"
  ) +
    
    
    ##############################################################################
  # Standardized placebo 95% CI
  ##############################################################################
  
  ggplot2::geom_ribbon(
    
    data =
      placebo_standardized_weekly,
    
    ggplot2::aes(
      
      x =
        week,
      
      ymin =
        placebo_lower_95,
      
      ymax =
        placebo_upper_95
    ),
    
    inherit.aes =
      FALSE,
    
    fill =
      "grey70",
    
    alpha =
      0.30
  ) +
    
    
    ##############################################################################
  # Multinomial MSM 95% CIs
  ##############################################################################
  
  ggplot2::geom_ribbon(
    
    ggplot2::aes(
      
      ymin =
        lower_95,
      
      ymax =
        upper_95,
      
      fill =
        dose
    ),
    
    alpha =
      0.10,
    
    color =
      NA,
    
    show.legend =
      FALSE
  ) +
    
    
    ##############################################################################
  # Multinomial dose-strategy trajectories
  ##############################################################################
  
  ggplot2::geom_line(
    linewidth =
      1
  ) +
    
    ggplot2::geom_point(
      size =
        2
    ) +
    
    
    ##############################################################################
  # Same standardized placebo curve
  ##############################################################################
  
  ggplot2::geom_line(
    
    data =
      placebo_standardized_weekly,
    
    ggplot2::aes(
      
      x =
        week,
      
      y =
        placebo_mean_delta,
      
      linetype =
        "Study-standardized placebo"
    ),
    
    inherit.aes =
      FALSE,
    
    color =
      "black",
    
    linewidth =
      1
  ) +
    
    
    ggplot2::geom_point(
      
      data =
        placebo_standardized_weekly,
      
      ggplot2::aes(
        
        x =
          week,
        
        y =
          placebo_mean_delta,
        
        shape =
          "Study-standardized placebo"
      ),
      
      inherit.aes =
        FALSE,
      
      color =
        "black",
      
      size =
        2.8
    ) +
    
    
    ##############################################################################
  # Same dose colors as primary analysis
  ##############################################################################
  
  ggplot2::scale_color_manual(
    
    values =
      DOSE_COLORS
  ) +
    
    ggplot2::scale_fill_manual(
      
      values =
        DOSE_COLORS
    ) +
    
    
    ggplot2::scale_linetype_manual(
      
      values =
        c(
          "Study-standardized placebo" =
            "longdash"
        )
    ) +
    
    
    ggplot2::scale_shape_manual(
      
      values =
        c(
          "Study-standardized placebo" =
            17
        )
    ) +
    
    
    ggplot2::scale_x_continuous(
      
      breaks =
        PREDICTION_WEEKS
    ) +
    
    
    ggplot2::scale_y_continuous(
      
      breaks =
        PREDICTION_Y_BREAKS
    ) +
    
    
    ##############################################################################
  # Labels
  ##############################################################################
  
  ggplot2::labs(
    
    x =
      "Weeks from baseline",
    
    y =
      "Predicted HAMD improvement from baseline",
    
    color =
      "Dose strategy (mg)",
    
    fill =
      "Dose strategy (mg)",
    
    linetype =
      NULL,
    
    shape =
      NULL,
    
    title =
      "Paroxetine dose strategies and standardized placebo",
    
    subtitle =
      paste0(
        "Multinomial IPTW sensitivity analysis; ",
        "predictions standardized over ",
        dplyr::n_distinct(
          attr(
            multinom_msm_models$zero_as_dose,
            "model_data"
          )$pid
        ),
        " active MSM patients"
      )
  ) +
    
    
    ggplot2::coord_cartesian(
      
      ylim =
        PREDICTION_Y_LIMITS
    ) +
    
    
    ggplot2::theme_classic(
      base_size =
        12
    ) +
    
    
    ggplot2::theme(
      
      plot.title =
        ggplot2::element_text(
          face =
            "bold"
        ),
      
      panel.grid.major.y =
        ggplot2::element_line(
          color =
            "grey90",
          linewidth =
            0.3
        ),
      
      legend.position =
        "right"
    )
  
  
  ################################################################################
  # Print graph
  ################################################################################
  
  multinom_strategy_prediction_plot
  
  