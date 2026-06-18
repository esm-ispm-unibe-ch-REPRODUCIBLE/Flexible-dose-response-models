max_followup_visit <- 8

locf_simple <- function(x) {
  out <- x
  last <- NA
  
  for (i in seq_along(out)) {
    if (!is.na(out[i])) {
      last <- out[i]
    } else if (!is.na(last)) {
      out[i] <- last
    }
  }
  
  out
}

categorise_to_dose_level <- function(x, levels) {
  levels <- sort(as.numeric(levels))
  
  if (length(levels) < 2) {
    return(factor(x, levels = levels))
  }
  
  cut_points <- (levels[-1] + levels[-length(levels)]) / 2
  breaks <- c(-Inf, cut_points, Inf)
  
  factor(
    as.numeric(as.character(cut(
      x,
      breaks = breaks,
      labels = levels,
      include.lowest = TRUE,
      right = TRUE
    ))),
    levels = levels
  )
}

add_dose_history_factors <- function(data, dose_history_levels = c(0, 10, 20, 30, 40, 50)) {
  data %>%
    dplyr::mutate(
      dose_lag1_f = factor(dose_lag1, levels = dose_history_levels),
      dose_lag2_f = factor(dose_lag2, levels = dose_history_levels),
      dose_lag3_f = factor(dose_lag3, levels = dose_history_levels),
      avg_dose_before_lag3_f = categorise_to_dose_level(
        avg_dose_before_lag3,
        levels = dose_history_levels
      )
    )
}

prepare_dose_response_data <- function(
    data,
    treatment_name = "PAROXETINE",
    visit_min = 0,
    visit_max = 8,
    dose_levels = c(10, 20, 30, 40, 50),
    dose_history_levels = c(0, 10, 20, 30, 40, 50)
) {
  data %>%
    dplyr::filter(
      treat == treatment_name,
      visit >= visit_min,
      visit <= visit_max
    ) %>%
    dplyr::mutate(
      pid = as.factor(pid),
      studyid = as.factor(studyid),
      sex = as.factor(sex),
      
      age = as.numeric(age),
      visit = as.numeric(visit),
      dose = as.numeric(dose),
      outcome = as.numeric(outcome),
      side.effects = as.numeric(side.effects),
      
      side.effects_model = dplyr::case_when(
        visit == 0 ~ 0,
        TRUE ~ side.effects
      ),
      side.effects = dplyr::case_when(
        visit == 0 ~ 0,
        TRUE ~ side.effects
      ),
      
      dose_f = factor(dose, levels = dose_levels, ordered = TRUE),
      dose_current_f = factor(dose, levels = dose_levels)
    ) %>%
    dplyr::arrange(pid, visit) %>%
    dplyr::group_by(pid) %>%
    dplyr::mutate(
      baseline_visit = min(visit, na.rm = TRUE),
      is_baseline = visit == baseline_visit,
      
      outcome_0 = outcome[which.min(visit)],
      delta_outcome = outcome_0 - outcome,
      
      dose_lag1 = dplyr::coalesce(dplyr::lag(dose, 1), 0),
      dose_lag2 = dplyr::coalesce(dplyr::lag(dose, 2), 0),
      dose_lag3 = dplyr::coalesce(dplyr::lag(dose, 3), 0),
      
      delta_outcome_locf = locf_simple(delta_outcome),
      side.effects_model_locf = locf_simple(side.effects_model),
      
      delta_outcome_lag1 = dplyr::coalesce(dplyr::lag(delta_outcome_locf), 0),
      side.effects_lag1 = dplyr::coalesce(dplyr::lag(side.effects_model_locf), 0),
      
      next_visit = dplyr::lead(visit),
      R_next = dplyr::case_when(
        visit >= max_followup_visit ~ NA_integer_,
        !is.na(next_visit) ~ 1L,
        is.na(next_visit) ~ 0L
      ),
      
      avg_dose_before_lag3 = vapply(seq_along(dose), function(j) {
        if (j <= 4) {
          return(0)
        }
        
        early_doses <- dose[seq_len(j - 4)]
        
        if (all(is.na(early_doses))) {
          0
        } else {
          mean(early_doses, na.rm = TRUE)
        }
      }, numeric(1))
    ) %>%
    dplyr::ungroup() %>%
    add_dose_history_factors(dose_history_levels = dose_history_levels) %>%
    dplyr::mutate(
      use_treatment_weight = !is_baseline &
        !is.na(dose_f) &
        !is.na(delta_outcome) &
        !is.na(side.effects) &
        !is.na(outcome_0),
      
      use_censoring_weight = !is_baseline &
        !is.na(R_next) &
        !is.na(dose_current_f) &
        !is.na(delta_outcome) &
        !is.na(side.effects) &
        !is.na(outcome_0) &
        !is.na(age) &
        !is.na(sex),
      
      use_msm = !is_baseline &
        !is.na(delta_outcome) &
        !is.na(outcome_0) &
        !is.na(age) &
        !is.na(sex) &
        !is.na(dose_lag1_f) &
        !is.na(dose_lag2_f) &
        !is.na(dose_lag3_f) &
        !is.na(avg_dose_before_lag3_f)
    )
}

summarise_dose_response_data <- function(data) {
  data %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_patients = dplyr::n_distinct(pid),
      n_treatment_weight_rows = sum(use_treatment_weight),
      n_censoring_weight_rows = sum(use_censoring_weight),
      n_msm_rows = sum(use_msm)
    )
}

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
      dose_current_f,
      outcome_0,
      delta_outcome,
      dose_lag1_f,
      dose_lag2_f,
      dose_lag3_f,
      avg_dose_before_lag3,
      avg_dose_before_lag3_f
    )
}

fit_iptw_denominator_model <- function(data, visit_df = 3) {
  model_dat <- data %>%
    dplyr::filter(use_treatment_weight) %>%
    dplyr::mutate(
      dose_f = droplevels(dose_f),
      dose_lag1_f = droplevels(dose_lag1_f)
    )
  
  MASS::polr(
    dose_f ~
      rms::rcs(visit, visit_df) +
      delta_outcome +
      side.effects +
      dose_lag1_f +
      dose_lag1_f:delta_outcome +
      dose_lag1_f:side.effects +
      outcome_0,
    data = model_dat,
    Hess = TRUE,
    method = "logistic"
  )
}

fit_iptw_numerator_model <- function(data, visit_df = 3) {
  model_dat <- data %>%
    dplyr::filter(use_treatment_weight) %>%
    dplyr::mutate(
      dose_f = droplevels(dose_f),
      dose_lag1_f = droplevels(dose_lag1_f)
    )
  
  MASS::polr(
    dose_f ~
      rms::rcs(visit, visit_df) +
      dose_lag1_f +
      outcome_0,
    data = model_dat,
    Hess = TRUE,
    method = "logistic"
  )
}

make_polr_coef_table <- function(model) {
  coef_tab <- coef(summary(model))
  
  cbind(
    coef_tab,
    p_value = 2 * pnorm(abs(coef_tab[, "t value"]), lower.tail = FALSE)
  )
}

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

add_iptw_treatment_weights <- function(
    data,
    denominator_model,
    numerator_model,
    min_prob = 1e-6
) {
  weight_dat <- data %>%
    dplyr::filter(use_treatment_weight) %>%
    dplyr::mutate(
      dose_f = droplevels(dose_f),
      dose_lag1_f = droplevels(dose_lag1_f)
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
    dplyr::mutate(
      p_dose_denominator = pmax(p_dose_denominator, min_prob),
      p_dose_numerator = pmax(p_dose_numerator, min_prob),
      SW_treatment = p_dose_numerator / p_dose_denominator
    ) %>%
    dplyr::arrange(pid, visit) %>%
    dplyr::group_by(pid) %>%
    dplyr::mutate(
      cSW_treatment = cumprod(SW_treatment)
    ) %>%
    dplyr::ungroup()
  
  data %>%
    dplyr::left_join(
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

fit_ipcw_denominator_model <- function(data) {
  model_dat <- data %>%
    dplyr::filter(use_censoring_weight) %>%
    dplyr::mutate(
      sex = droplevels(sex),
      dose_current_f = droplevels(dose_current_f)
    )
  
  glm(
    R_next ~
      visit +
      dose_current_f +
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

fit_ipcw_numerator_model <- function(data) {
  model_dat <- data %>%
    dplyr::filter(use_censoring_weight) %>%
    dplyr::mutate(
      sex = droplevels(sex),
      dose_current_f = droplevels(dose_current_f)
    )
  
  glm(
    R_next ~
      visit +
      dose_current_f +
      outcome_0 +
      age +
      sex,
    data = model_dat,
    family = binomial()
  )
}

add_ipcw_censoring_weights <- function(
    data,
    denominator_model,
    numerator_model,
    min_prob = 1e-6
) {
  data <- data %>%
    dplyr::select(
      -dplyr::any_of(c(
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
    dplyr::mutate(
      sex = droplevels(sex),
      dose_current_f = droplevels(dose_current_f)
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
    dplyr::mutate(
      p_censor_denominator = pmax(p_censor_denominator_raw, min_prob),
      p_censor_numerator = pmax(p_censor_numerator_raw, min_prob),
      SW_censoring = p_censor_numerator / p_censor_denominator
    ) %>%
    dplyr::arrange(pid, visit) %>%
    dplyr::group_by(pid) %>%
    dplyr::mutate(
      cSW_censoring = cumprod(SW_censoring)
    ) %>%
    dplyr::ungroup()
  
  data %>%
    dplyr::left_join(
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

add_total_weights <- function(data) {
  data %>%
    dplyr::mutate(
      SW_total = dplyr::case_when(
        use_msm &
          !is.na(cSW_treatment) &
          !is.na(cSW_censoring) ~ cSW_treatment * cSW_censoring,
        TRUE ~ NA_real_
      )
    )
}

check_total_weight_missingness <- function(data) {
  data %>%
    dplyr::filter(use_msm) %>%
    dplyr::summarise(
      n_msm_rows = dplyr::n(),
      n_patients = dplyr::n_distinct(pid),
      missing_cSW_treatment = sum(is.na(cSW_treatment)),
      missing_cSW_censoring = sum(is.na(cSW_censoring)),
      missing_SW_total = sum(is.na(SW_total)),
      percent_missing_SW_total = 100 * mean(is.na(SW_total))
    )
}

summarise_total_weights <- function(data) {
  data %>%
    dplyr::filter(use_msm, !is.na(SW_total)) %>%
    dplyr::summarise(
      n = dplyr::n(),
      n_patients = dplyr::n_distinct(pid),
      
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

truncate_total_weights <- function(data, lower = 0.01, upper = 0.99) {
  cutoffs <- quantile(
    data$SW_total,
    probs = c(lower, upper),
    na.rm = TRUE
  )
  
  data %>%
    dplyr::mutate(
      SW_total_trunc = dplyr::case_when(
        is.na(SW_total) ~ NA_real_,
        SW_total < cutoffs[[1]] ~ as.numeric(cutoffs[[1]]),
        SW_total > cutoffs[[2]] ~ as.numeric(cutoffs[[2]]),
        TRUE ~ SW_total
      )
    )
}

summarise_total_truncated_weights <- function(data) {
  data %>%
    dplyr::filter(use_msm, !is.na(SW_total_trunc)) %>%
    dplyr::summarise(
      n = dplyr::n(),
      n_patients = dplyr::n_distinct(pid),
      
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

check_truncation <- function(data) {
  data %>%
    dplyr::filter(use_msm, !is.na(SW_total)) %>%
    dplyr::summarise(
      n = dplyr::n(),
      n_lower_truncated = sum(SW_total != SW_total_trunc & SW_total < SW_total_trunc),
      n_upper_truncated = sum(SW_total != SW_total_trunc & SW_total > SW_total_trunc),
      n_any_truncated = sum(SW_total != SW_total_trunc),
      percent_any_truncated = 100 * mean(SW_total != SW_total_trunc)
    )
}

weighted_var <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  x <- x[ok]
  w <- w[ok]
  
  if (length(x) <= 1 || sum(w) == 0) {
    return(NA_real_)
  }
  
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

summarise_iptw_treatment_weights <- function(data) {
  data %>%
    dplyr::filter(use_treatment_weight, !is.na(cSW_treatment)) %>%
    dplyr::summarise(
      n = dplyr::n(),
      n_patients = dplyr::n_distinct(pid),
      
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

summarise_ipcw_censoring_weights <- function(data) {
  data %>%
    dplyr::filter(use_censoring_weight, !is.na(cSW_censoring)) %>%
    dplyr::summarise(
      n = dplyr::n(),
      n_patients = dplyr::n_distinct(pid),
      
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
fit_weighted_msm <- function(
    data,
    weight_var = "SW_total_trunc",
    visit_df = 3,
    corstr = "independence"
) {
  needed_vars <- c(
    "pid",
    "visit",
    "delta_outcome",
    "outcome_0",
    "dose_lag1_f",
    "dose_lag2_f",
    "dose_lag3_f",
    "avg_dose_before_lag3",
    "use_msm",
    weight_var
  )
  
  missing_vars <- setdiff(needed_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop("Missing variables: ", paste(missing_vars, collapse = ", "))
  }
  
  model_dat <- data %>%
    dplyr::filter(
      use_msm,
      !is.na(.data[[weight_var]]),
      .data[[weight_var]] > 0,
      !is.na(delta_outcome),
      !is.na(outcome_0),
      !is.na(visit),
      !is.na(dose_lag1_f),
      !is.na(dose_lag2_f),
      !is.na(dose_lag3_f),
      !is.na(avg_dose_before_lag3)
    ) %>%
    dplyr::mutate(
      final_weight = .data[[weight_var]],
      dose_lag1_f = droplevels(dose_lag1_f),
      dose_lag2_f = droplevels(dose_lag2_f),
      dose_lag3_f = droplevels(dose_lag3_f)
    ) %>%
    dplyr::arrange(pid, visit)
  
  msm_formula <- stats::as.formula(
    paste0(
      "delta_outcome ~ ",
      "rms::rcs(visit, ", visit_df, ") * (outcome_0 + dose_lag1_f) + ",
      "dose_lag2_f + ",
      "dose_lag3_f + ",
      "avg_dose_before_lag3"
    )
  )
  
  fit <- geepack::geeglm(
    msm_formula,
    id = pid,
    waves = visit,
    data = model_dat,
    weights = final_weight,
    family = gaussian(link = "identity"),
    corstr = corstr,
    std.err = "san.se"
  )
  
  attr(fit, "model_data") <- model_dat
  attr(fit, "dose_history_levels") <- levels(model_dat$dose_lag1_f)
  
  fit
}

make_gee_coef_table <- function(model) {
  coef_tab <- as.data.frame(coef(summary(model)))
  
  data.frame(
    Term = rownames(coef_tab),
    coef_tab,
    row.names = NULL,
    check.names = FALSE
  )
}

make_gee_coef_table_robust_naive <- function(model) {
  beta <- coef(model)
  
  V_robust <- model$geese$vbeta
  V_naive <- model$geese$vbeta.naiv
  
  dimnames(V_robust) <- list(names(beta), names(beta))
  dimnames(V_naive) <- list(names(beta), names(beta))
  
  robust_se <- sqrt(diag(V_robust))
  naive_se <- sqrt(diag(V_naive))
  z_robust <- beta / robust_se
  
  data.frame(
    Term = names(beta),
    Estimate = beta,
    Naive_SE = naive_se,
    Robust_SE = robust_se,
    Wald_robust = z_robust^2,
    p_value_robust = 2 * pnorm(abs(z_robust), lower.tail = FALSE),
    row.names = NULL,
    check.names = FALSE
  )
}

make_strategy_data <- function(
    strategy_dose,
    visits,
    outcome_0_value,
    baseline_dose = strategy_dose,
    dose_history_levels = c(0, 10, 20, 30, 40, 50)
) {
  visit_order <- seq_along(visits)
  
  dat <- tibble::tibble(
    visit = visits,
    visit_order = visit_order,
    strategy_dose = strategy_dose,
    strategy = paste0("Always ", strategy_dose, " mg"),
    outcome_0 = outcome_0_value
  ) %>%
    dplyr::mutate(
      dose_lag1 = dplyr::if_else(
        visit_order == 1,
        baseline_dose,
        strategy_dose
      ),
      dose_lag2 = dplyr::case_when(
        visit_order == 1 ~ 0,
        visit_order == 2 ~ baseline_dose,
        TRUE ~ strategy_dose
      ),
      dose_lag3 = dplyr::case_when(
        visit_order <= 2 ~ 0,
        visit_order == 3 ~ baseline_dose,
        TRUE ~ strategy_dose
      ),
      avg_dose_before_lag3 = dplyr::case_when(
        visit_order <= 2 ~ 0,
        TRUE ~ (baseline_dose + pmax(visit_order - 3, 0) * strategy_dose) /
          (visit_order - 2)
      )
    ) %>%
    add_dose_history_factors(dose_history_levels = dose_history_levels)
  
  dat
}

make_strategy_history <- function(
    strategy_dose,
    target_visit = 8,
    baseline_dose = strategy_dose,
    outcome_0_value = NA_real_,
    dose_history_levels = c(0, 10, 20, 30, 40, 50)
) {
  make_strategy_data(
    strategy_dose = strategy_dose,
    visits = seq_len(target_visit),
    outcome_0_value = outcome_0_value,
    baseline_dose = baseline_dose,
    dose_history_levels = dose_history_levels
  ) %>%
    dplyr::filter(visit == target_visit) %>%
    dplyr::slice_tail(n = 1)
}

estimate_all_pairwise_strategy_contrasts <- function(
    model,
    target_visit = 8,
    strategy_doses = c(20, 30, 40, 50),
    baseline_dose = NULL
) {
  model_dat <- attr(model, "model_data")
  
  if (is.null(model_dat)) {
    stop("Model is missing model_data attribute.")
  }
  
  outcome_0_value <- mean(model_dat$outcome_0, na.rm = TRUE)
  dose_history_levels <- attr(model, "dose_history_levels")
  
  if (is.null(dose_history_levels)) {
    dose_history_levels <- c(0, 10, 20, 30, 40, 50)
  }
  
  strategy_histories <- dplyr::bind_rows(
    lapply(strategy_doses, function(dose_value) {
      make_strategy_history(
        strategy_dose = dose_value,
        target_visit = target_visit,
        baseline_dose = if (is.null(baseline_dose)) dose_value else baseline_dose,
        outcome_0_value = outcome_0_value,
        dose_history_levels = dose_history_levels
      )
    })
  )
  
  beta <- coef(model)
  
  V <- tryCatch(
    vcov(model),
    error = function(e) model$geese$vbeta
  )
  
  dimnames(V) <- list(names(beta), names(beta))
  
  Terms <- stats::delete.response(stats::terms(model))
  X <- stats::model.matrix(Terms, strategy_histories)
  
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
  
  pairs <- t(combn(seq_along(strategy_doses), 2))
  
  pairwise_contrasts <- lapply(seq_len(nrow(pairs)), function(i) {
    lower_id <- pairs[i, 1]
    higher_id <- pairs[i, 2]
    
    lower_dose <- strategy_doses[lower_id]
    higher_dose <- strategy_doses[higher_id]
    
    contrast_vector <- X[higher_id, ] - X[lower_id, ]
    
    estimate <- as.numeric(sum(contrast_vector * beta))
    se <- sqrt(as.numeric(t(contrast_vector) %*% V %*% contrast_vector))
    
    tibble::tibble(
      contrast = paste0(
        "Always ", higher_dose,
        " mg vs always ", lower_dose,
        " mg"
      ),
      target_visit = target_visit,
      estimate = estimate,
      se = se,
      lower_95 = estimate - qnorm(0.975) * se,
      upper_95 = estimate + qnorm(0.975) * se,
      p_value = 2 * pnorm(abs(estimate / se), lower.tail = FALSE)
    )
  })
  
  dplyr::bind_rows(pairwise_contrasts)
}

fit_weighted_lmm_msm <- function(
    data,
    weight_var = "SW_total_trunc",
    visit_df = 3
) {
  needed_vars <- c(
    "pid",
    "visit",
    "delta_outcome",
    "outcome_0",
    "dose_lag1_f",
    "dose_lag2_f",
    "dose_lag3_f",
    "avg_dose_before_lag3_f",
    "use_msm",
    weight_var
  )
  
  missing_vars <- setdiff(needed_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop("Missing variables: ", paste(missing_vars, collapse = ", "))
  }
  
  model_dat <- data %>%
    dplyr::filter(
      use_msm,
      !is.na(.data[[weight_var]]),
      .data[[weight_var]] > 0,
      !is.na(delta_outcome),
      !is.na(outcome_0),
      !is.na(visit),
      !is.na(dose_lag1_f),
      !is.na(dose_lag2_f),
      !is.na(dose_lag3_f),
      !is.na(avg_dose_before_lag3_f)
    ) %>%
    dplyr::mutate(
      final_weight = .data[[weight_var]],
      dose_lag1_f = droplevels(dose_lag1_f),
      dose_lag2_f = droplevels(dose_lag2_f),
      dose_lag3_f = droplevels(dose_lag3_f),
      avg_dose_before_lag3_f = droplevels(avg_dose_before_lag3_f)
    ) %>%
    dplyr::arrange(pid, visit)
  
  lmm_formula <- stats::as.formula(
    paste0(
      "delta_outcome ~ rms::rcs(visit, ", visit_df, ") * (",
      "outcome_0 + ",
      "dose_lag1_f + ",
      "dose_lag2_f + ",
      "dose_lag3_f + ",
      "avg_dose_before_lag3_f",
      ") + (1 | pid)"
    )
  )
  
  fit <- lmerTest::lmer(
    lmm_formula,
    data = model_dat,
    weights = final_weight,
    REML = FALSE,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 100000)
    )
  )
  
  attr(fit, "model_data") <- model_dat
  attr(fit, "dose_history_levels") <- levels(model_dat$dose_lag1_f)
  
  fit
}