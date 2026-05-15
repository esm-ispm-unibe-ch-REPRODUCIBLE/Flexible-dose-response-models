######################################### R-script with all functions needed for ################################
###########################          Flexible dose-response models                 #####################


######################### Functions for data preperation#################################
########################################################################################

### Function that prepares the needed variable for the IPTW/MSM analysis 

max_followup_visit <- 8
prepare_dose_response_data <- function(
  data,
  treatment_name = "PAROXETINE",
  visit_min = 0,
  visit_max = 8,
  dose_levels = c(10, 20, 30, 40, 50),
  dose_history_levels = c(0, 10, 20, 30, 40, 50)
) {
  data %>%
    filter(
      treat == treatment_name,
      visit >= visit_min,
      visit <= visit_max
    ) %>%
    mutate(
      pid = as.factor(pid),
      studyid = as.factor(studyid),
      sex = as.factor(sex),
      
      age = as.numeric(age),
      visit = as.numeric(visit),
      dose = as.numeric(dose),
      outcome = as.numeric(outcome),
      side.effects = as.numeric(side.effects),
      side.effects_model = case_when(
        visit == 0 ~ 0,
        TRUE ~ side.effects
      ),
      side.effects = case_when(
        visit == 0 ~ 0,
        TRUE ~ side.effects
      ),
      dose_f = factor(dose, levels = dose_levels, ordered = TRUE)
    ) %>%
    arrange(pid, visit) %>%
    group_by(pid) %>%
    mutate(
      baseline_visit = min(visit, na.rm = TRUE),
      is_baseline = visit == baseline_visit,
      
      outcome_0 = outcome[which.min(visit)],
      delta_outcome = outcome_0 - outcome,
      
      dose_lag1 = coalesce(lag(dose, 1), 0),
      dose_lag2 = coalesce(lag(dose, 2), 0),
      dose_lag3 = coalesce(lag(dose, 3), 0),
      
      delta_outcome_locf = na.locf(delta_outcome, na.rm = FALSE),
      side.effects_model_locf = na.locf(side.effects_model, na.rm = FALSE),
      
      delta_outcome_lag1 = coalesce(lag(delta_outcome_locf), 0),
      side.effects_lag1 = coalesce(lag(side.effects_model_locf), 0),
      
      
      next_visit = lead(visit),
      R_next = case_when(
        visit >= max_followup_visit ~ NA_integer_,
        !is.na(next_visit) ~ 1L,
        is.na(next_visit) ~ 0L
      ),
      
      avg_dose_before_lag3 = sapply(seq_along(dose), function(j) {
        if (j <= 3) {
          0
        } else {
          mean(dose[seq_len(j - 3)], na.rm = TRUE)
        }
      })
    ) %>%
    ungroup() %>%
    mutate(
      dose_f = factor(dose, levels = dose_levels, ordered = TRUE),
      
      dose_lag1_f = factor(dose_lag1, levels = dose_history_levels),
      dose_lag2_f = factor(dose_lag2, levels = dose_history_levels),
      dose_lag3_f = factor(dose_lag3, levels = dose_history_levels),
      
      avg_dose_before_lag3_f = factor(
        avg_dose_before_lag3,
        levels = dose_history_levels
      ),
      
      use_treatment_weight = !is_baseline &
        !is.na(dose_f) &
        !is.na(delta_outcome) &
        !is.na(side.effects) &
        !is.na(outcome_0),
      
      use_censoring_weight = !is_baseline &
        !is.na(R_next) &
        !is.na(dose) &
        !is.na(delta_outcome) &
        !is.na(side.effects) &
        !is.na(outcome_0) &
        !is.na(age) &
        !is.na(sex),
      
      use_msm = !is_baseline &
        !is.na(delta_outcome) &
        !is.na(outcome_0) &
        !is.na(age) &
        !is.na(sex)
    )
}

## check how many visits are used in each model
summarise_dose_response_data <- function(data) {
  data %>%
    summarise(
      n_rows = n(),
      n_patients = n_distinct(pid),
      n_treatment_weight_rows = sum(use_treatment_weight),
      n_censoring_weight_rows = sum(use_censoring_weight),
      n_msm_rows = sum(use_msm)
    )
}

## Inspect also visually what you want
make_vars_check <- function(data) {
  data %>%
    dplyr::select(
      pid,
      studyid,
      age,
      sex,
      dose,
      outcome,
      visit,
      treat,
      side.effects,
      dose_f,
      outcome_0,
      delta_outcome,
      dose_lag1_f
    )
}



###################################### STEP 1: IPTW #######################################


#### Step 1a; IPTW for dose
#denominator model
fit_iptw_denominator_model <- function(data, visit_df = 3) {
  model_dat <- data %>%
    dplyr::filter(use_treatment_weight) %>%
    mutate(
      dose_f = droplevels(dose_f),
      dose_lag1<-as.numeric(dose_lag1)
    )
  
  polr(
    dose_f ~ rcs(visit, visit_df) +
      delta_outcome +
      side.effects +
      dose_lag1 +
      dose_lag1:delta_outcome +
      dose_lag1:side.effects +
      outcome_0,
    data = model_dat,
    Hess = TRUE,
    method = "logistic"
  )
}

#export table for results
make_polr_coef_table <- function(model) {
  coef_tab <- coef(summary(model))
  
  cbind(
    coef_tab,
    p_value = 2 * pnorm(abs(coef_tab[, "t value"]), lower.tail = FALSE)
  )
}

## numerator model

fit_iptw_numerator_model <- function(data, visit_df = 3) {
  model_dat <- data %>%
    filter(use_treatment_weight) %>%
    mutate(
      dose_f = droplevels(dose_f),
      dose_lag1 = as.numeric(dose_lag1),
      dose_lag1<-as.numeric(dose_lag1)
    )
  
  polr(
    dose_f ~ rcs(visit, visit_df) +
      dose_lag1 +
      outcome_0,
    data = model_dat,
    Hess = TRUE,
    method = "logistic"
  )
}


### Function to extract the observed-dose probability
get_observed_dose_prob <- function(model, data, outcome_var = "dose_f") {
  prob_mat <- predict(model, newdata = data, type = "probs")
  
  observed_dose <- as.character(data[[outcome_var]])
  
  prob_mat[
    cbind(
      seq_len(nrow(prob_mat)),
      match(observed_dose, colnames(prob_mat))
    )
  ]
}
### Function to add IPTW treatment weights
add_iptw_treatment_weights <- function(
  data,
  denominator_model,
  numerator_model,
  min_prob = 1e-6
) {
  weight_dat <- data %>%
    filter(use_treatment_weight) %>%
    mutate(
      dose_f = droplevels(dose_f),
      dose_lag1 = as.numeric(dose_lag1)
    )
  
  weight_dat$p_dose_denominator <- get_observed_dose_prob(
    model = denominator_model,
    data = weight_dat,
    outcome_var = "dose_f"
  )
  
  weight_dat$p_dose_numerator <- get_observed_dose_prob(
    model = numerator_model,
    data = weight_dat,
    outcome_var = "dose_f"
  )
  
  weight_dat <- weight_dat %>%
    mutate(
      p_dose_denominator = pmax(p_dose_denominator, min_prob),
      p_dose_numerator = pmax(p_dose_numerator, min_prob),
      SW_treatment = p_dose_numerator / p_dose_denominator
    ) %>%
    arrange(pid, visit) %>%
    group_by(pid) %>%
    mutate(
      cSW_treatment = cumprod(SW_treatment)
    ) %>%
    ungroup()
  
  data %>%
    left_join(
      weight_dat %>%
        dplyr::select(
          pid,
          visit,
          p_dose_denominator,
          p_dose_numerator,
          SW_treatment,
          cSW_treatment
        ),
      by = c("pid", "visit")
    )
}

#### Step 1b; IPTW for censoring
#denominator model

fit_ipcw_denominator_model <- function(data) {
  model_dat <- data %>%
    dplyr::filter(use_censoring_weight) %>%
    mutate(
      sex = droplevels(sex),
      dose_lag1<-as.numeric(dose_lag1)
    )
  
  glm(
    R_next ~ visit +
      dose +
      delta_outcome +
      delta_outcome_lag1 +
      side.effects +
      side.effects_lag1 +
      outcome_0 +
      age +
      sex,
    data = model_dat,
    family = binomial()
  )
}

#numerator model
fit_ipcw_numerator_model <- function(data) {
  model_dat <- data %>%
    dplyr::filter(use_censoring_weight) %>%
    mutate(
      sex = droplevels(sex),
      dose_lag1<-as.numeric(dose_lag1)
    )
  
  glm(
    R_next ~ visit +
      dose +
      outcome_0 +
      age +
      sex,
    data = model_dat,
    family = binomial()
  )
}
# add the weigths in the data
add_ipcw_censoring_weights <- function(
  data,
  denominator_model,
  numerator_model,
  min_prob = 1e-6
) {
  data <- data %>%
    dplyr::select(
      -any_of(c(
        "p_censor_denominator_raw",
        "p_censor_numerator_raw",
        "p_censor_denominator",
        "p_censor_numerator",
        "SW_censoring",
        "cSW_censoring"
      ))
    )
  
  weight_dat <- data %>%
    dplyr::filter(use_censoring_weight) %>%
    mutate(
      sex = droplevels(sex),
      dose_f = droplevels(dose_f)
    )
  
  weight_dat$p_censor_denominator_raw <- predict(
    denominator_model,
    newdata = weight_dat,
    type = "response"
  )
  
  weight_dat$p_censor_numerator_raw <- predict(
    numerator_model,
    newdata = weight_dat,
    type = "response"
  )
  
  weight_dat <- weight_dat %>%
    mutate(
      p_censor_denominator = pmax(p_censor_denominator_raw, min_prob),
      p_censor_numerator = pmax(p_censor_numerator_raw, min_prob),
      SW_censoring = p_censor_numerator / p_censor_denominator
    ) %>%
    arrange(pid, visit) %>%
    group_by(pid) %>%
    mutate(
      cSW_censoring = cumprod(SW_censoring)
    ) %>%
    ungroup()
  
  data %>%
    left_join(
      weight_dat %>%
        dplyr::select(
          pid,
          visit,
          p_censor_denominator_raw,
          p_censor_numerator_raw,
          p_censor_denominator,
          p_censor_numerator,
          SW_censoring,
          cSW_censoring
        ),
      by = c("pid", "visit")
    )
}

#### Step 1c: IPTW combined for dose and censoring

## multiply the weigths
add_total_weights <- function(data) {
  data %>%
    mutate(
      SW_total = case_when(
        use_msm &
          !is.na(cSW_treatment) &
          !is.na(cSW_censoring) ~ cSW_treatment * cSW_censoring,
        TRUE ~ NA_real_
      )
    )
}

#check missingness
check_total_weight_missingness <- function(data) {
  data %>%
    dplyr::filter(use_msm) %>%
    summarise(
      n_msm_rows = n(),
      n_patients = n_distinct(pid),
      missing_cSW_treatment = sum(is.na(cSW_treatment)),
      missing_cSW_censoring = sum(is.na(cSW_censoring)),
      missing_SW_total = sum(is.na(SW_total)),
      percent_missing_SW_total = 100 * mean(is.na(SW_total))
    )
}
##summary of weights
summarise_total_weights <- function(data) {
  data %>%
    dplyr::filter(use_msm, !is.na(SW_total)) %>%
    summarise(
      n = n(),
      n_patients = n_distinct(pid),
      
      mean_SW_total = mean(SW_total),
      sd_SW_total = sd(SW_total),
      min_SW_total = min(SW_total),
      p1_SW_total = quantile(SW_total, 0.01),
      p5_SW_total = quantile(SW_total, 0.05),
      p25_SW_total = quantile(SW_total, 0.25),
      p50_SW_total = quantile(SW_total, 0.50),
      p75_SW_total = quantile(SW_total, 0.75),
      p95_SW_total = quantile(SW_total, 0.95),
      p99_SW_total = quantile(SW_total, 0.99),
      max_SW_total = max(SW_total),
      
      ESS_SW_total = sum(SW_total)^2 / sum(SW_total^2)
    )
}
## truncate to 1st and 99th
truncate_total_weights <- function(data, lower = 0.01, upper = 0.99) {
  cutoffs <- quantile(
    data$SW_total,
    probs = c(lower, upper),
    na.rm = TRUE
  )
  
  data %>%
    mutate(
      SW_total_trunc = case_when(
        is.na(SW_total) ~ NA_real_,
        SW_total < cutoffs[[1]] ~ as.numeric(cutoffs[[1]]),
        SW_total > cutoffs[[2]] ~ as.numeric(cutoffs[[2]]),
        TRUE ~ SW_total
      )
    )
}

### summarise truncated weights
summarise_total_truncated_weights <- function(data) {
  data %>%
    dplyr::filter(use_msm, !is.na(SW_total_trunc)) %>%
    summarise(
      n = n(),
      n_patients = n_distinct(pid),
      
      mean_SW_total_trunc = mean(SW_total_trunc),
      sd_SW_total_trunc = sd(SW_total_trunc),
      min_SW_total_trunc = min(SW_total_trunc),
      p1_SW_total_trunc = quantile(SW_total_trunc, 0.01),
      p5_SW_total_trunc = quantile(SW_total_trunc, 0.05),
      p25_SW_total_trunc = quantile(SW_total_trunc, 0.25),
      p50_SW_total_trunc = quantile(SW_total_trunc, 0.50),
      p75_SW_total_trunc = quantile(SW_total_trunc, 0.75),
      p95_SW_total_trunc = quantile(SW_total_trunc, 0.95),
      p99_SW_total_trunc = quantile(SW_total_trunc, 0.99),
      max_SW_total_trunc = max(SW_total_trunc),
      
      ESS_SW_total_trunc = sum(SW_total_trunc)^2 / sum(SW_total_trunc^2)
    )
}

##check how many lines were truncated

check_truncation <- function(data) {
  data %>%
    dplyr::filter(use_msm, !is.na(SW_total)) %>%
    summarise(
      n = n(),
      n_lower_truncated = sum(SW_total != SW_total_trunc & SW_total < SW_total_trunc),
      n_upper_truncated = sum(SW_total != SW_total_trunc & SW_total > SW_total_trunc),
      n_any_truncated = sum(SW_total != SW_total_trunc),
      percent_any_truncated = 100 * mean(SW_total != SW_total_trunc)
    )
}

#### check of the balance before and after
weighted_var <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  x <- x[ok]
  w <- w[ok]
  
  if (length(x) <= 1 || sum(w) == 0) return(NA_real_)
  
  mu <- weighted.mean(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

assess_balance_by_window <- function(
  data,
  covariates = c("delta_outcome", "side.effects"),
  weight_var = NULL
) {
  dat <- data %>%
    dplyr::filter(
      use_msm,
      !is.na(visit_window),
      !is.na(dose_f)
    )
  
  purrr::map_dfr(covariates, function(covariate) {
    dat %>%
      dplyr::group_by(visit_window) %>%
      dplyr::group_modify(function(window_dat, key) {
        x_all <- window_dat[[covariate]]
        
        if (is.null(weight_var)) {
          mean_all <- mean(x_all, na.rm = TRUE)
          var_all <- var(x_all, na.rm = TRUE)
          
          window_dat %>%
            dplyr::group_by(dose_f) %>%
            dplyr::summarise(
              n = dplyr::n(),
              mean_cov = mean(.data[[covariate]], na.rm = TRUE),
              smd = (mean_cov - mean_all) / sqrt(var_all),
              .groups = "drop"
            )
        } else {
          w_all <- window_dat[[weight_var]]
          mean_all <- weighted.mean(x_all, w_all, na.rm = TRUE)
          var_all <- weighted_var(x_all, w_all)
          
          window_dat %>%
            dplyr::group_by(dose_f) %>%
            dplyr::summarise(
              n = dplyr::n(),
              mean_cov = weighted.mean(
                .data[[covariate]],
                .data[[weight_var]],
                na.rm = TRUE
              ),
              smd = (mean_cov - mean_all) / sqrt(var_all),
              .groups = "drop"
            )
        }
      }) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        covariate = covariate,
        abs_smd = abs(smd),
        weight_type = ifelse(is.null(weight_var), "unweighted", weight_var)
      ) %>%
      dplyr::select(
        weight_type,
        covariate,
        visit_window,
        dose = dose_f,
        n,
        mean_cov,
        smd,
        abs_smd
      )
  })
}

###only for IPTW
summarise_iptw_treatment_weights <- function(data) {
  data %>%
    dplyr::filter(use_treatment_weight, !is.na(cSW_treatment)) %>%
    summarise(
      n = n(),
      n_patients = n_distinct(pid),
      
      mean_SW_treatment = mean(SW_treatment),
      sd_SW_treatment = sd(SW_treatment),
      min_SW_treatment = min(SW_treatment),
      p1_SW_treatment = quantile(SW_treatment, 0.01),
      p5_SW_treatment = quantile(SW_treatment, 0.05),
      p25_SW_treatment = quantile(SW_treatment, 0.25),
      p50_SW_treatment = quantile(SW_treatment, 0.50),
      p75_SW_treatment = quantile(SW_treatment, 0.75),
      p95_SW_treatment = quantile(SW_treatment, 0.95),
      p99_SW_treatment = quantile(SW_treatment, 0.99),
      max_SW_treatment = max(SW_treatment),
      
      mean_cSW_treatment = mean(cSW_treatment),
      sd_cSW_treatment = sd(cSW_treatment),
      min_cSW_treatment = min(cSW_treatment),
      p1_cSW_treatment = quantile(cSW_treatment, 0.01),
      p5_cSW_treatment = quantile(cSW_treatment, 0.05),
      p25_cSW_treatment = quantile(cSW_treatment, 0.25),
      p50_cSW_treatment = quantile(cSW_treatment, 0.50),
      p75_cSW_treatment = quantile(cSW_treatment, 0.75),
      p95_cSW_treatment = quantile(cSW_treatment, 0.95),
      p99_cSW_treatment = quantile(cSW_treatment, 0.99),
      max_cSW_treatment = max(cSW_treatment),
      
      ESS_cSW_treatment = sum(cSW_treatment)^2 / sum(cSW_treatment^2)
    )
}

###Only for IPCW

summarise_ipcw_censoring_weights <- function(data) {
  data %>%
    dplyr::filter(use_censoring_weight, !is.na(cSW_censoring)) %>%
    summarise(
      n = n(),
      n_patients = n_distinct(pid),
      
      mean_SW_censoring = mean(SW_censoring),
      sd_SW_censoring = sd(SW_censoring),
      min_SW_censoring = min(SW_censoring),
      p1_SW_censoring = quantile(SW_censoring, 0.01),
      p5_SW_censoring = quantile(SW_censoring, 0.05),
      p25_SW_censoring = quantile(SW_censoring, 0.25),
      p50_SW_censoring = quantile(SW_censoring, 0.50),
      p75_SW_censoring = quantile(SW_censoring, 0.75),
      p95_SW_censoring = quantile(SW_censoring, 0.95),
      p99_SW_censoring = quantile(SW_censoring, 0.99),
      max_SW_censoring = max(SW_censoring),
      
      mean_cSW_censoring = mean(cSW_censoring),
      sd_cSW_censoring = sd(cSW_censoring),
      min_cSW_censoring = min(cSW_censoring),
      p1_cSW_censoring = quantile(cSW_censoring, 0.01),
      p5_cSW_censoring = quantile(cSW_censoring, 0.05),
      p25_cSW_censoring = quantile(cSW_censoring, 0.25),
      p50_cSW_censoring = quantile(cSW_censoring, 0.50),
      p75_cSW_censoring = quantile(cSW_censoring, 0.75),
      p95_cSW_censoring = quantile(cSW_censoring, 0.95),
      p99_cSW_censoring = quantile(cSW_censoring, 0.99),
      max_cSW_censoring = max(cSW_censoring),
      
      ESS_cSW_censoring = sum(cSW_censoring)^2 / sum(cSW_censoring^2)
    )
}
