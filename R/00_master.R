
#######################################################################################################
#################################    MASTER R-script for flexible dose-response analysis ###########################
####################################################################################################################################3

rm(list = ls())
####################################### Load needed variables ###############################################
library(dplyr)
library(forcats)
library(splines)
library(MASS)
library(rms)
library(purrr)
library(ggplot2)
# source the needed functions
source("Functions.R")

############################### Data preparation ####################################################

### 1st in 01_data_prep.R the real data need to be loaded and to include these variables named as:
#"pid", "studyid","age",
#"sex","dose", "outcome" , "visit","treat" and "side.effects" (number of reposrted: 0 to 10)
source("01_data_prep.R")

### Create the needed vars for the analyses
trials <- list(
  trial_1 = trial_1
)

analysis_trials <- lapply(
  trials,
  prepare_dose_response_data,
  treatment_name = "PAROXETINE",
  visit_min = 0,
  visit_max = 8,
  dose_levels = c(0,10,20, 30, 40, 50),
  dose_history_levels = c(0, 10, 20, 30, 40, 50)
)
##make a summary of included rows in the analysis - some may be dropped due to NAs
lapply(analysis_trials, summarise_dose_response_data)
#test that your IPTW needed vars are correct
test<-lapply(analysis_trials, make_vars_check)
#test[["trial_1"]]

############################### Some checks for visits and dosages ####################################################
#### test doses per visit
analysis_trials$trial_1 %>%
  dplyr::filter(use_treatment_weight) %>%
  dplyr::count(visit, dose_f) %>%
  tidyr::pivot_wider(
    names_from = dose_f,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::arrange(visit)

## test doses titration per visit

check<-analysis_trials$trial_1 %>%
  dplyr::filter(use_treatment_weight) %>%
  dplyr::count(visit, dose_lag1_f, dose_f) %>%
  tidyr::pivot_wider(
    names_from = dose_f,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::arrange(visit, dose_lag1_f)

print(check, n=Inf)
### check how the visits are distributes across rows
visit_day_by_patient_row <- trial_1 %>%
  arrange(pid, visit_day) %>%
  group_by(pid) %>%
  mutate(patient_row = row_number()) %>%
  ungroup() %>%
  group_by(patient_row) %>%
  summarise(
    n = n(),
    min_day = min(visit_day, na.rm = TRUE),
    p25_day = quantile(visit_day, 0.25, na.rm = TRUE),
    median_day = median(visit_day, na.rm = TRUE),
    mean_day = mean(visit_day, na.rm = TRUE),
    p75_day = quantile(visit_day, 0.75, na.rm = TRUE),
    max_day = max(visit_day, na.rm = TRUE),
    unique_days = paste(sort(unique(visit_day)), collapse = ", "),
    .groups = "drop"
  )

visit_day_by_patient_row

###################################### STEP 1: IPTW #######################################


#### Step 1a; IPTW for dose


#denominator model
iptw_denominator_models <- lapply(
  analysis_trials,
  fit_iptw_denominator_model
)
   #coefficients table
iptw_denominator_coef_tabs <- lapply(
  iptw_denominator_models,
  make_polr_coef_table
)
   #results
#lapply(iptw_denominator_models, summary)

#numerator model
iptw_numerator_models <- lapply(
  analysis_trials,
  fit_iptw_numerator_model
)
   #coefficients
iptw_numerator_coef_tabs <- lapply(
  iptw_numerator_models,
  make_polr_coef_table
)
   #results
#lapply(iptw_numerator_models, summary)

analysis_trials$trial_1 <- add_iptw_treatment_weights(
  data = analysis_trials$trial_1,
  denominator_model = iptw_denominator_models$trial_1,
  numerator_model = iptw_numerator_models$trial_1
)
# 
# ### test
# test_trial_1 <- analysis_trials$trial_1 %>%
#   dplyr::select(
#     pid,
#     studyid,
#     age,
#     sex,
#     dose,
#     outcome,
#     visit,
#     treat,
#     side.effects,
#     dose_f,
#     outcome_0,
#     delta_outcome,
#     dose_lag1,
#     p_dose_denominator,
#     p_dose_numerator,
#     SW_treatment,
#     cSW_treatment
#   )
# 
# analysis_trials$trial_1 %>%
#   dplyr::filter(use_treatment_weight) %>%
#   group_by(dose_f) %>%
#   summarise(
#     n = n(),
#     min_p_den = min(p_dose_denominator, na.rm = TRUE),
#     p1_p_den = quantile(p_dose_denominator, 0.01, na.rm = TRUE),
#     median_p_den = median(p_dose_denominator, na.rm = TRUE),
#     max_SW = max(SW_treatment, na.rm = TRUE),
#     p99_SW = quantile(SW_treatment, 0.99, na.rm = TRUE),
#     .groups = "drop"
#   )



#### Step 1b; IPTW for censoring

#denominator
ipcw_denominator_models <- lapply(
  analysis_trials,
  fit_ipcw_denominator_model
)
#numerator
ipcw_numerator_models <- lapply(
  analysis_trials,
  fit_ipcw_numerator_model
)
#results
ipcw_denominator_coef_tabs <- lapply(
  ipcw_denominator_models,
  function(model) {
    coef_tab <- coef(summary(model))
    cbind(
      coef_tab,
      p_value = 2 * pnorm(abs(coef_tab[, "z value"]), lower.tail = FALSE)
    )
  }
)

ipcw_numerator_coef_tabs <- lapply(
  ipcw_numerator_models,
  function(model) {
    coef_tab <- coef(summary(model))
    cbind(
      coef_tab,
      p_value = 2 * pnorm(abs(coef_tab[, "z value"]), lower.tail = FALSE)
    )
  }
)

# add the weigths in the data
analysis_trials$trial_1 <- add_ipcw_censoring_weights(
  data = analysis_trials$trial_1,
  denominator_model = ipcw_denominator_models$trial_1,
  numerator_model = ipcw_numerator_models$trial_1
)


# ## test
 # test_trial_1 <- analysis_trials$trial_1 %>%
 #   dplyr::select(
 #     pid,
 #     studyid,
 #     age,
 #     sex,
 #     dose,
 #     outcome,
 #     visit,
 #     side.effects,
 #     side.effects_lag1,
 #     dose_f,
 #     outcome_0,
 #     delta_outcome,
 #     delta_outcome_lag1,
 #     dose_lag1,
 #     p_censor_denominator,
 #     p_censor_numerator,
 #     SW_censoring,
 #     cSW_censoring,
 #     p_dose_denominator,
 #     p_dose_numerator,
 #     SW_treatment,
 #     cSW_treatment
 #   )
#

##test that all the lags are correct
# test_ipcw_denominator_trial_1 <- analysis_trials$trial_1 %>%
#    dplyr::select(
#      pid,
#      studyid,
#      age,
#      sex,
#      dose,
#      outcome,
#      visit,
#      side.effects,
#      R_next,
#      outcome_0,
#      delta_outcome,
#      delta_outcome_lag1,
#      side.effects_lag1,
#      p_censor_denominator,
#      p_censor_numerator,
#      SW_censoring,
#      cSW_censoring
#    )


#### Step 1c: IPTW combined for dose and censoring

analysis_trials <- lapply(analysis_trials, add_total_weights)

#check missingness
total_weight_missingness <- lapply(
  analysis_trials,
  check_total_weight_missingness
)

total_weight_missingness

##summary of weights

total_weight_summaries <- lapply(
  analysis_trials,
  summarise_total_weights
)

total_weight_summaries

### truncate at 1st and 99th
analysis_trials <- lapply(
  analysis_trials,
  truncate_total_weights
)
##summarise truncated values

total_truncated_weight_summaries <- lapply(
  analysis_trials,
  summarise_total_truncated_weights
)

total_truncated_weight_summaries
##check how many lines were truncated
total_weight_truncation_checks <- lapply(
  analysis_trials,
  check_truncation
)

total_weight_truncation_checks

##### Check for IPTW alone

iptw_weight_summaries <- lapply(
  analysis_trials,
  summarise_iptw_treatment_weights
)


iptw_weight_summaries

## check for IPCW alone

ipcw_weight_summaries <- lapply(
  analysis_trials,
  summarise_ipcw_censoring_weights
)

ipcw_weight_summaries


#### Plots

## plot 1  Censoring per week (to see if we will use it with rcs or with linear model)

observed_censor_week <- analysis_trials$trial_1 %>%
  filter(use_censoring_weight, !is.na(R_next)) %>%
  mutate(
    observed_censor = 1 - R_next
  ) %>%
  group_by(visit) %>%
  summarise(
    n_at_risk = n(),
    n_censored = sum(observed_censor == 1),
    observed_p_censor = mean(observed_censor),
    .groups = "drop"
  )

observed_censor_week

censoring_plot<-ggplot(observed_censor_week, aes(x = visit, y = observed_p_censor)) +
  geom_line() +
  geom_point() +
  geom_text(
    aes(label = paste0(n_censored, "/", n_at_risk)),
    vjust = -0.7,
    size = 3
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Week",
    y = "Censoring proportion",
    title = "Censoring by week"
  )
censoring_plot

### Plot 2

model_dat <- analysis_trials$trial_1 %>%
  dplyr::filter(use_treatment_weight)

mean_outcome_0 <- mean(model_dat$outcome_0, na.rm = TRUE)

profile_dat <- tidyr::expand_grid(
  visit = 1:8,
  profile = c("Improving profile", "Worsening profile"),
  dose_lag1 = c(20, 30, 40, 50)
) %>%
  group_by(profile, dose_lag1) %>%
  mutate(
    delta_outcome = if_else(
      profile == "Improving profile",
      seq(0, 30, length.out = n()),
      seq(0, -15, length.out = n())
    ),
    side.effects = 4,
    outcome_0 = mean_outcome_0
  ) %>%
  ungroup()

pred_probs <- predict(
  iptw_denominator_models$trial_1,
  newdata = profile_dat,
  type = "probs"
)

plot_dat <- cbind(profile_dat, as.data.frame(pred_probs)) %>%
  tidyr::pivot_longer(
    cols = c("20", "30", "40", "50"),
    names_to = "dose",
    values_to = "probability"
  ) %>%
  mutate(
    dose = factor(dose, levels = c("20", "30", "40", "50")),
    previous_dose = paste0("Previous dose: ", dose_lag1, " mg")
  )

dose_plot <- ggplot(plot_dat, aes(x = visit, y = probability, color = dose, group = dose)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(profile ~ previous_dose) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_x_continuous(breaks = 1:8) +
  labs(
    x = "Week",
    y = "Predicted probability of prescribed dose",
    color = "Current dose",
    title = "Predicted dose-assignment probabilities by clinical profile and previous dose",
    subtitle = "Side-effect burden fixed at 4; baseline HAMD fixed at the sample mean"
  ) +
  theme_minimal()

dose_plot





























##check for all trials
weight_summary_all_trials <- bind_rows(
  lapply(names(analysis_trials), function(trial_name) {
    analysis_trials[[trial_name]] %>%
      dplyr::filter(use_msm, !is.na(SW_total_trunc)) %>%
      summarise(
        trial = trial_name,
        n = n(),
        n_patients = n_distinct(pid),
        mean_SW_total = mean(SW_total),
        sd_SW_total = sd(SW_total),
        p1_SW_total = quantile(SW_total, 0.01),
        p50_SW_total = quantile(SW_total, 0.50),
        p99_SW_total = quantile(SW_total, 0.99),
        max_SW_total = max(SW_total),
        ESS_SW_total = sum(SW_total)^2 / sum(SW_total^2),
        
        mean_SW_total_trunc = mean(SW_total_trunc),
        sd_SW_total_trunc = sd(SW_total_trunc),
        p1_SW_total_trunc = quantile(SW_total_trunc, 0.01),
        p50_SW_total_trunc = quantile(SW_total_trunc, 0.50),
        p99_SW_total_trunc = quantile(SW_total_trunc, 0.99),
        max_SW_total_trunc = max(SW_total_trunc),
        ESS_SW_total_trunc = sum(SW_total_trunc)^2 / sum(SW_total_trunc^2)
      )
  })
)

weight_summary_all_trials







############## Not ready yet!!
##check balance
analysis_trials$trial_1 <- analysis_trials$trial_1 %>%
  mutate(
    visit_window = factor(visit)
  )

balance_unweighted <- assess_balance_by_window(
  analysis_trials$trial_1,
  weight_var = NULL
)

balance_treatment <- assess_balance_by_window(
  analysis_trials$trial_1,
  weight_var = "cSW_treatment"
)

balance_total <- assess_balance_by_window(
  analysis_trials$trial_1,
  weight_var = "SW_total_trunc"
)

balance_all <- bind_rows(
  balance_unweighted,
  balance_treatment,
  balance_total
)
balance_summary <- balance_all %>%
  group_by(weight_type, covariate) %>%
  summarise(
    n_cells = n(),
    median_abs_smd = median(abs_smd, na.rm = TRUE),
    p75_abs_smd = quantile(abs_smd, 0.75, na.rm = TRUE),
    p90_abs_smd = quantile(abs_smd, 0.90, na.rm = TRUE),
    max_abs_smd = max(abs_smd, na.rm = TRUE),
    percent_below_0_1 = 100 * mean(abs_smd < 0.1, na.rm = TRUE),
    percent_below_0_2 = 100 * mean(abs_smd < 0.2, na.rm = TRUE),
    .groups = "drop"
  )

balance_summary

balance_all %>%
  dplyr::filter(weight_type %in% c("unweighted", "cSW_treatment")) %>%
  group_by(weight_type, covariate) %>%
  summarise(
    median_abs_smd = median(abs_smd, na.rm = TRUE),
    p75_abs_smd = quantile(abs_smd, 0.75, na.rm = TRUE),
    p90_abs_smd = quantile(abs_smd, 0.90, na.rm = TRUE),
    max_abs_smd = max(abs_smd, na.rm = TRUE),
    .groups = "drop"
  )

rm(list = ls())
