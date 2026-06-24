# step-by-step per-arm scripts; no row-binding of arm-specific dose-factor tables; fixed visit conversion scalar if_else issue.
#######################################################################################################
############################ DR_01_functions_stepwise_by_arm.R ########################################
#######################################################################################################

# General helper and modelling functions for flexible-dose IPTW/MSM analyses.
# This file contains functions only. It does not run the analysis.

safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

load_required_packages <- function() {
  pkgs <- c(
    "dplyr", "tidyr", "purrr", "tibble", "ggplot2", "zoo",
    "MASS", "rms", "geepack", "nnet", "lme4", "lmerTest", "writexl"
  )
  invisible(lapply(pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package not installed: ", pkg)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

strip_ansi <- function(x) {
  gsub("\\033\\[[0-9;]*m", "", x)
}

safe_imap <- function(.x, .f) {
  purrr::imap(
    .x,
    function(dat, arm_name) {
      tryCatch(
        {
          out <- .f(dat, arm_name)
          list(ok = TRUE, result = out, error = NA_character_)
        },
        error = function(e) {
          list(ok = FALSE, result = NULL, error = strip_ansi(conditionMessage(e)))
        }
      )
    }
  )
}

get_success_results <- function(runs) {
  purrr::map(
    purrr::keep(runs, ~ isTRUE(.x$ok)),
    "result"
  )
}

get_status_table <- function(runs) {
  purrr::imap_dfr(
    runs,
    function(x, arm_name) {
      tibble::tibble(
        arm_name = arm_name,
        ok = isTRUE(x$ok),
        error = ifelse(isTRUE(x$ok), NA_character_, x$error)
      )
    }
  )
}

make_numeric_dose_factor <- function(x, include_zero = FALSE, ordered = FALSE) {
  
  x_num <- as.numeric(as.character(x))
  
  lev <- sort(unique(x_num[!is.na(x_num)]))
  
  if (!include_zero) {
    lev <- lev[lev > 0]
  }
  
  factor(
    as.character(x_num),
    levels = as.character(lev),
    ordered = ordered
  )
}
# ---------------------------- data standardisation ----------------------------

standardise_trial_columns <- function(data, trial_name = NA_character_) {
  names(data) <- trimws(names(data))

  rename_if_present <- function(dat, standard_name, aliases) {
    aliases <- unique(c(standard_name, aliases))
    hit <- aliases[aliases %in% names(dat)]
    if (length(hit) > 0 && standard_name != hit[1]) {
      names(dat)[names(dat) == hit[1]] <- standard_name
    }
    dat
  }

  data <- data |>
    rename_if_present("pid", c("patient", "patient_id", "subject", "subject_id", "USUBJID", "id")) |>
    rename_if_present("studyid", c("study", "trial", "trial_id", "study_id")) |>
    rename_if_present("age", c("AGE", "Age")) |>
    rename_if_present("sex", c("SEX", "Sex", "gender")) |>
    rename_if_present("dose", c("DOSE", "Dose", "dose_mg", "daily_dose")) |>
    rename_if_present("outcome", c("HAMD", "hamd", "score", "total_score", "Y")) |>
    rename_if_present("visit", c("VISIT", "Visit", "day", "visit_day", "time")) |>
    rename_if_present("treat", c("TRT", "treatment", "arm", "drug", "treatment_name")) |>
    rename_if_present("side.effects", c("side_effects", "side_effect", "ae_score", "safety"))

  required <- c("pid", "age", "sex", "dose", "outcome", "visit", "treat")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns in ", trial_name, ": ", paste(missing, collapse = ", "))
  }

  if (!"studyid" %in% names(data)) {
    data$studyid <- trial_name
  }

  if (!"side.effects" %in% names(data)) {
    data$side.effects <- NA_real_
  }

  data |>
    dplyr::mutate(
      trial_name = as.character(trial_name),
      pid = as.character(pid),
      studyid = as.character(studyid),
      treat = as.character(treat),
      sex = as.character(sex),
      age = suppressWarnings(as.numeric(age)),
      dose = suppressWarnings(as.numeric(dose)),
      outcome = suppressWarnings(as.numeric(outcome)),
      visit = suppressWarnings(as.numeric(visit)),
      side.effects = suppressWarnings(as.numeric(side.effects))
    )
}

convert_visit_to_week <- function(data, visit_max = 9, assume_visit_is_days = NULL) {
  dat <- data |> dplyr::filter(!is.na(visit), visit >= 0)

  if (is.null(assume_visit_is_days)) {
    positive_visits <- dat$visit[dat$visit > 0]
    assume_visit_is_days <- length(positive_visits) > 0 && max(positive_visits, na.rm = TRUE) > visit_max + 2
  }

  if (isTRUE(assume_visit_is_days)) {
    dat |>
      dplyr::mutate(
        visit_day = visit,
        visit = dplyr::case_when(
          visit_day == 0 ~ 0,
          visit_day > 0 ~ ceiling(visit_day / 7),
          TRUE ~ NA_real_
        ),
        visit = dplyr::if_else(visit > visit_max, as.numeric(visit_max), as.numeric(visit))
      )
  } else {
    dat |>
      dplyr::mutate(
        visit_day = visit,
        visit = as.numeric(visit)
      )
  }
}

keep_last_record_per_patient_visit <- function(data) {
  data |>
    dplyr::arrange(pid, visit, visit_day) |>
    dplyr::group_by(pid, visit) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::ungroup()
}

impute_active_zero_dose <- function(data, method = c("nearest", "previous", "next")) {
  method <- match.arg(method)

  dat <- data |>
    dplyr::arrange(pid, visit) |>
    dplyr::group_by(pid) |>
    dplyr::mutate(
      dose_before_zero_impute = dose,
      dose_nonzero = dplyr::if_else(dose > 0, dose, NA_real_),
      prev_dose = zoo::na.locf(dose_nonzero, na.rm = FALSE),
      next_dose = rev(zoo::na.locf(rev(dose_nonzero), na.rm = FALSE))
    ) |>
    dplyr::ungroup()

  dat <- dat |>
    dplyr::mutate(
      dose = dplyr::case_when(
        dose != 0 | is.na(dose) ~ dose,
        method == "previous" & !is.na(prev_dose) ~ prev_dose,
        method == "next" & !is.na(next_dose) ~ next_dose,
        method == "nearest" & !is.na(prev_dose) ~ prev_dose,
        method == "nearest" & is.na(prev_dose) & !is.na(next_dose) ~ next_dose,
        TRUE ~ dose
      )
    ) |>
    dplyr::select(-dose_nonzero, -prev_dose, -next_dose)

  dat
}

# Group sparse doses into the nearest dose with adequate support.
# This is arm-specific and uses post-baseline active-arm rows.
group_sparse_doses <- function(data, min_n_per_dose = 10) {
  dat <- data |>
    dplyr::mutate(
      dose_before_grouping = dose
    )

  counts <- dat |>
    dplyr::filter(visit > 0, !is.na(dose), dose > 0) |>
    dplyr::count(dose, name = "n") |>
    dplyr::arrange(dose)

  if (nrow(counts) == 0) {
    dat$dose_grouped <- dat$dose
    return(dat)
  }

  good <- counts |> dplyr::filter(n >= min_n_per_dose)

  # If no dose has enough support, keep all observed doses. The model may fail later,
  # but the data-prep step should not invent a category.
  if (nrow(good) == 0) {
    dat$dose_grouped <- dat$dose
    return(dat)
  }

  good_doses <- good$dose

  map_one <- function(x) {
    if (is.na(x) || x <= 0) return(x)
    if (x %in% good_doses) return(x)
    good_doses[which.min(abs(good_doses - x))]
  }

  dat |>
    dplyr::mutate(
      dose_grouped = vapply(dose, map_one, numeric(1)),
      dose = dose_grouped
    )
}

simulate_side_effects_if_needed <- function(data, seed = 2025) {
  if ("side.effects" %in% names(data) && any(!is.na(data$side.effects[data$visit > 0]))) {
    return(data)
  }

  set.seed(seed)
  dose_range <- range(data$dose[data$dose >= 0], na.rm = TRUE)
  if (!all(is.finite(dose_range)) || diff(dose_range) == 0) {
    dose_range <- c(0, 1)
  }

  data |>
    dplyr::mutate(
      dose_scaled_tmp = dplyr::if_else(
        dose >= 0,
        (dose - dose_range[1]) / (dose_range[2] - dose_range[1]),
        NA_real_
      ),
      inv_effect_tmp = dplyr::if_else(dose > 0, 1 - dose_scaled_tmp, NA_real_),
      side_raw_tmp = inv_effect_tmp + rnorm(dplyr::n(), mean = 0, sd = 0.35),
      side_scaled_tmp = pmin(pmax(side_raw_tmp, 0), 1),
      side.effects = dplyr::if_else(visit > 0, as.numeric(round(10 * side_scaled_tmp)), NA_real_)
    ) |>
    dplyr::select(-dose_scaled_tmp, -inv_effect_tmp, -side_raw_tmp, -side_scaled_tmp)
}

# ---------------------------- longitudinal variables ----------------------------

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
  if (length(levels) < 2) return(factor(x, levels = levels))
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

add_dose_history_factors <- function(data, dose_history_levels = NULL) {
  if (is.null(dose_history_levels)) {
    dose_history_levels <- sort(unique(c(0, data$dose[data$dose > 0])))
  }

  data |>
    dplyr::mutate(
      dose_lag1_f = make_numeric_dose_factor(
        dose_lag1,
        include_zero = TRUE,
        ordered = FALSE
      ),
      dose_lag2_f = make_numeric_dose_factor(
        dose_lag2,
        include_zero = TRUE,
        ordered = FALSE
      ),
      dose_lag3_f = make_numeric_dose_factor(
        dose_lag3,
        include_zero = TRUE,
        ordered = FALSE
      ),
      avg_dose_before_lag3_f = categorise_to_dose_level(avg_dose_before_lag3, levels = dose_history_levels)
    )
}

prepare_dose_response_arm_data <- function(
    data,
    visit_min = 0,
    visit_max = 9,
    max_followup_visit = visit_max,
    dose_levels = NULL,
    dose_history_levels = NULL
) {
  dat <- data |>
    dplyr::filter(visit >= visit_min, visit <= visit_max)

  if (is.null(dose_levels)) {
    dose_levels <- sort(unique(dat$dose[dat$visit > 0 & dat$dose > 0 & !is.na(dat$dose)]))
  }
  if (is.null(dose_history_levels)) {
    dose_history_levels <- sort(unique(c(0, dose_levels)))
  }

  dat |>
    dplyr::mutate(
      pid = as.factor(pid),
      studyid = as.factor(studyid),
      sex = as.factor(sex),
      treat = as.factor(treat),
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
      dose_f = make_numeric_dose_factor(
        dose,
        include_zero = FALSE,
        ordered = TRUE
      ),
      dose_current_f = make_numeric_dose_factor(
        dose,
        include_zero = FALSE,
        ordered = FALSE
      )
    ) |>
    dplyr::arrange(pid, visit) |>
    dplyr::group_by(pid) |>
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
        if (j <= 4) return(0)
        early_doses <- dose[seq_len(j - 4)]
        if (all(is.na(early_doses))) 0 else mean(early_doses, na.rm = TRUE)
      }, numeric(1))
    ) |>
    dplyr::ungroup() |>
    add_dose_history_factors(dose_history_levels = dose_history_levels) |>
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
        !is.na(avg_dose_before_lag3)
    )
}

summarise_dose_response_data <- function(data) {
  data |>
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_patients = dplyr::n_distinct(pid),
      n_treatment_weight_rows = sum(use_treatment_weight, na.rm = TRUE),
      n_censoring_weight_rows = sum(use_censoring_weight, na.rm = TRUE),
      n_msm_rows = sum(use_msm, na.rm = TRUE),
      dose_levels = paste(levels(droplevels(dose_f)), collapse = ", ")
    )
}

# ---------------------------- IPTW models ----------------------------

make_formula_from_terms <- function(outcome, terms) {
  stats::as.formula(paste(outcome, "~", paste(terms, collapse = " + ")))
}

fit_iptw_denominator_model <- function(data, visit_df = 3) {
  model_dat <- data |>
    dplyr::filter(use_treatment_weight) |>
    dplyr::mutate(
      dose_f = droplevels(dose_f),
      dose_lag1_f = droplevels(dose_lag1_f)
    )

  if (nlevels(model_dat$dose_f) < 2) stop("IPTW denominator: fewer than 2 dose levels.")

  terms <- c(
    paste0("rms::rcs(visit, ", visit_df, ")"),
    "delta_outcome",
    "side.effects",
    "dose_lag1_f",
    "dose_lag1_f:delta_outcome",
    "dose_lag1_f:side.effects",
    "outcome_0"
  )

  MASS::polr(
    make_formula_from_terms("dose_f", terms),
    data = model_dat,
    Hess = TRUE,
    method = "logistic"
  )
}

fit_iptw_numerator_model <- function(data, visit_df = 3) {
  model_dat <- data |>
    dplyr::filter(use_treatment_weight) |>
    dplyr::mutate(
      dose_f = droplevels(dose_f),
      dose_lag1_f = droplevels(dose_lag1_f)
    )

  if (nlevels(model_dat$dose_f) < 2) stop("IPTW numerator: fewer than 2 dose levels.")

  terms <- c(
    paste0("rms::rcs(visit, ", visit_df, ")"),
    "dose_lag1_f",
    "outcome_0"
  )

  MASS::polr(
    make_formula_from_terms("dose_f", terms),
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
  
  if (is.null(dim(prob_mat))) {
    prob_mat <- matrix(
      prob_mat,
      nrow = 1,
      dimnames = list(NULL, model$lev)
    )
  }
  
  observed_dose <- as.character(data[[outcome_var]])
  
  matched_col <- match(observed_dose, colnames(prob_mat))
  
  if (any(is.na(matched_col))) {
    missing_doses <- unique(observed_dose[is.na(matched_col)])
    stop(
      "Some observed doses are not present in the model-predicted probability columns: ",
      paste(missing_doses, collapse = ", ")
    )
  }
  
  prob_mat[
    cbind(
      seq_len(nrow(prob_mat)),
      matched_col
    )
  ]
}
add_iptw_treatment_weights <- function(
    data,
    denominator_model,
    numerator_model,
    min_prob = 1e-6
) {
  data <- data |>
    dplyr::select(
      -dplyr::any_of(c(
        "p_dose_denominator",
        "p_dose_numerator",
        "SW_treatment",
        "cSW_treatment"
      ))
    )
  
  dose_levels_model <- denominator_model$lev
  
  weight_dat <- data |>
    dplyr::filter(use_treatment_weight) |>
    dplyr::mutate(
      dose_f = factor(
        as.character(dose_f),
        levels = dose_levels_model,
        ordered = TRUE
      )
    )
  
  if (any(is.na(weight_dat$dose_f))) {
    missing_doses <- data |>
      dplyr::filter(use_treatment_weight) |>
      dplyr::mutate(dose_chr = as.character(dose_f)) |>
      dplyr::filter(!(dose_chr %in% dose_levels_model)) |>
      dplyr::distinct(dose_chr) |>
      dplyr::pull(dose_chr)
    
    stop(
      "Some observed dose levels are not present in the IPTW model: ",
      paste(missing_doses, collapse = ", ")
    )
  }
  
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
  
  weight_dat <- weight_dat |>
    dplyr::mutate(
      p_dose_denominator = pmax(p_dose_denominator, min_prob),
      p_dose_numerator = pmax(p_dose_numerator, min_prob),
      SW_treatment = p_dose_numerator / p_dose_denominator
    ) |>
    dplyr::arrange(pid, visit) |>
    dplyr::group_by(pid) |>
    dplyr::mutate(
      cSW_treatment = cumprod(SW_treatment)
    ) |>
    dplyr::ungroup()
  
  data |>
    dplyr::left_join(
      weight_dat |>
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

# ---------------------------- IPCW models ----------------------------

fit_ipcw_denominator_model <- function(data) {
  model_dat <- data |>
    dplyr::filter(use_censoring_weight) |>
    dplyr::mutate(
      sex = droplevels(sex),
      dose_current_f = droplevels(dose_current_f)
    )

  if (length(unique(model_dat$R_next)) < 2) stop("IPCW denominator: R_next has fewer than 2 values.")

  glm(
    R_next ~ visit + dose_current_f + delta_outcome + delta_outcome_lag1 +
      side.effects + side.effects_lag1 + outcome_0 + age + sex,
    data = model_dat,
    family = binomial()
  )
}

fit_ipcw_numerator_model <- function(data) {
  model_dat <- data |>
    dplyr::filter(use_censoring_weight) |>
    dplyr::mutate(
      sex = droplevels(sex),
      dose_current_f = droplevels(dose_current_f)
    )

  if (length(unique(model_dat$R_next)) < 2) stop("IPCW numerator: R_next has fewer than 2 values.")

  glm(
    R_next ~ visit + dose_current_f + outcome_0 + age + sex,
    data = model_dat,
    family = binomial()
  )
}

make_glm_coef_table <- function(model) {
  coef_tab <- coef(summary(model))
  cbind(
    coef_tab,
    p_value = 2 * pnorm(abs(coef_tab[, "z value"]), lower.tail = FALSE)
  )
}

add_ipcw_censoring_weights <- function(data, denominator_model, numerator_model, min_prob = 1e-6) {
  data <- data |>
    dplyr::select(-dplyr::any_of(c(
      "p_censor_denominator_raw", "p_censor_numerator_raw",
      "p_censor_denominator", "p_censor_numerator",
      "SW_censoring", "cSW_censoring"
    )))

  weight_dat <- data |>
    dplyr::filter(use_censoring_weight) |>
    dplyr::mutate(
      sex = droplevels(sex),
      dose_current_f = droplevels(dose_current_f)
    )

  weight_dat$p_censor_denominator_raw <- predict(denominator_model, newdata = weight_dat, type = "response")
  weight_dat$p_censor_numerator_raw <- predict(numerator_model, newdata = weight_dat, type = "response")

  weight_dat <- weight_dat |>
    dplyr::mutate(
      p_censor_denominator = pmax(p_censor_denominator_raw, min_prob),
      p_censor_numerator = pmax(p_censor_numerator_raw, min_prob),
      SW_censoring = p_censor_numerator / p_censor_denominator
    ) |>
    dplyr::arrange(pid, visit) |>
    dplyr::group_by(pid) |>
    dplyr::mutate(cSW_censoring = cumprod(SW_censoring)) |>
    dplyr::ungroup()

  data |>
    dplyr::left_join(
      weight_dat |>
        dplyr::select(
          pid, visit,
          p_censor_denominator_raw, p_censor_numerator_raw,
          p_censor_denominator, p_censor_numerator,
          SW_censoring, cSW_censoring
        ),
      by = c("pid", "visit")
    )
}

add_no_censoring_weights <- function(data) {
  data |>
    dplyr::mutate(
      p_censor_denominator_raw = dplyr::if_else(use_msm, 1, NA_real_),
      p_censor_numerator_raw = dplyr::if_else(use_msm, 1, NA_real_),
      p_censor_denominator = dplyr::if_else(use_msm, 1, NA_real_),
      p_censor_numerator = dplyr::if_else(use_msm, 1, NA_real_),
      SW_censoring = dplyr::if_else(use_msm, 1, NA_real_),
      cSW_censoring = dplyr::if_else(use_msm, 1, NA_real_)
    )
}

# ---------------------------- weights summaries ----------------------------

add_total_weights <- function(data) {
  data |>
    dplyr::mutate(
      SW_total = dplyr::case_when(
        use_msm & !is.na(cSW_treatment) & !is.na(cSW_censoring) ~ cSW_treatment * cSW_censoring,
        TRUE ~ NA_real_
      )
    )
}

truncate_total_weights <- function(data, lower = 0.01, upper = 0.99) {
  cutoffs <- quantile(data$SW_total, probs = c(lower, upper), na.rm = TRUE)
  data |>
    dplyr::mutate(
      SW_total_trunc = dplyr::case_when(
        is.na(SW_total) ~ NA_real_,
        SW_total < cutoffs[[1]] ~ as.numeric(cutoffs[[1]]),
        SW_total > cutoffs[[2]] ~ as.numeric(cutoffs[[2]]),
        TRUE ~ SW_total
      )
    )
}

summarise_iptw_treatment_weights <- function(data) {
  data |>
    dplyr::filter(use_treatment_weight, !is.na(cSW_treatment)) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_patients = dplyr::n_distinct(pid),
      mean_SW_treatment = mean(SW_treatment),
      sd_SW_treatment = sd(SW_treatment),
      min_SW_treatment = min(SW_treatment),
      p1_SW_treatment = quantile(SW_treatment, 0.01),
      p50_SW_treatment = quantile(SW_treatment, 0.50),
      p99_SW_treatment = quantile(SW_treatment, 0.99),
      max_SW_treatment = max(SW_treatment),
      mean_cSW_treatment = mean(cSW_treatment),
      sd_cSW_treatment = sd(cSW_treatment),
      min_cSW_treatment = min(cSW_treatment),
      p1_cSW_treatment = quantile(cSW_treatment, 0.01),
      p50_cSW_treatment = quantile(cSW_treatment, 0.50),
      p99_cSW_treatment = quantile(cSW_treatment, 0.99),
      max_cSW_treatment = max(cSW_treatment),
      ESS_cSW_treatment = sum(cSW_treatment)^2 / sum(cSW_treatment^2)
    )
}

summarise_ipcw_censoring_weights <- function(data) {
  data |>
    dplyr::filter(use_censoring_weight, !is.na(cSW_censoring)) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_patients = dplyr::n_distinct(pid),
      mean_SW_censoring = mean(SW_censoring),
      sd_SW_censoring = sd(SW_censoring),
      min_SW_censoring = min(SW_censoring),
      p1_SW_censoring = quantile(SW_censoring, 0.01),
      p50_SW_censoring = quantile(SW_censoring, 0.50),
      p99_SW_censoring = quantile(SW_censoring, 0.99),
      max_SW_censoring = max(SW_censoring),
      mean_cSW_censoring = mean(cSW_censoring),
      sd_cSW_censoring = sd(cSW_censoring),
      min_cSW_censoring = min(cSW_censoring),
      p1_cSW_censoring = quantile(cSW_censoring, 0.01),
      p50_cSW_censoring = quantile(cSW_censoring, 0.50),
      p99_cSW_censoring = quantile(cSW_censoring, 0.99),
      max_cSW_censoring = max(cSW_censoring),
      ESS_cSW_censoring = sum(cSW_censoring)^2 / sum(cSW_censoring^2)
    )
}

summarise_total_weights <- function(data) {
  data |>
    dplyr::filter(use_msm, !is.na(SW_total)) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_patients = dplyr::n_distinct(pid),
      mean_SW_total = mean(SW_total),
      sd_SW_total = sd(SW_total),
      min_SW_total = min(SW_total),
      p1_SW_total = quantile(SW_total, 0.01),
      p50_SW_total = quantile(SW_total, 0.50),
      p99_SW_total = quantile(SW_total, 0.99),
      max_SW_total = max(SW_total),
      ESS_SW_total = sum(SW_total)^2 / sum(SW_total^2)
    )
}

summarise_total_truncated_weights <- function(data) {
  data |>
    dplyr::filter(use_msm, !is.na(SW_total_trunc)) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_patients = dplyr::n_distinct(pid),
      mean_SW_total_trunc = mean(SW_total_trunc),
      sd_SW_total_trunc = sd(SW_total_trunc),
      min_SW_total_trunc = min(SW_total_trunc),
      p1_SW_total_trunc = quantile(SW_total_trunc, 0.01),
      p50_SW_total_trunc = quantile(SW_total_trunc, 0.50),
      p99_SW_total_trunc = quantile(SW_total_trunc, 0.99),
      max_SW_total_trunc = max(SW_total_trunc),
      ESS_SW_total_trunc = sum(SW_total_trunc)^2 / sum(SW_total_trunc^2)
    )
}

check_truncation <- function(data) {
  data |>
    dplyr::filter(use_msm, !is.na(SW_total), !is.na(SW_total_trunc)) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_lower_truncated = sum(SW_total != SW_total_trunc & SW_total < SW_total_trunc),
      n_upper_truncated = sum(SW_total != SW_total_trunc & SW_total > SW_total_trunc),
      n_any_truncated = sum(SW_total != SW_total_trunc),
      percent_any_truncated = 100 * mean(SW_total != SW_total_trunc)
    )
}

# ---------------------------- MSM and predictions ----------------------------

fit_weighted_msm <- function(
    data,
    weight_var = "SW_total_trunc",
    visit_df = 3,
    corstr = "independence",
    include_dose_time_interaction = TRUE
) {
  needed_vars <- c(
    "pid", "visit", "delta_outcome", "outcome_0",
    "dose_lag1_f", "dose_lag2_f", "dose_lag3_f",
    "avg_dose_before_lag3", "use_msm", weight_var
  )
  missing_vars <- setdiff(needed_vars, names(data))
  if (length(missing_vars) > 0) stop("Missing variables: ", paste(missing_vars, collapse = ", "))

  model_dat <- data |>
    dplyr::filter(
      use_msm,
      !is.na(.data[[weight_var]]),
      .data[[weight_var]] > 0,
      !is.na(delta_outcome), !is.na(outcome_0), !is.na(visit),
      !is.na(dose_lag1_f), !is.na(dose_lag2_f), !is.na(dose_lag3_f),
      !is.na(avg_dose_before_lag3)
    ) |>
    dplyr::mutate(
      final_weight = .data[[weight_var]],
      dose_lag1_f = droplevels(dose_lag1_f),
      dose_lag2_f = droplevels(dose_lag2_f),
      dose_lag3_f = droplevels(dose_lag3_f)
    ) |>
    dplyr::arrange(pid, visit)

  if (include_dose_time_interaction) {
    rhs <- paste0(
      "rms::rcs(visit, ", visit_df, ") * (outcome_0 + dose_lag1_f) + ",
      "dose_lag2_f + dose_lag3_f + avg_dose_before_lag3"
    )
  } else {
    rhs <- paste0(
      "rms::rcs(visit, ", visit_df, ") * outcome_0 + ",
      "dose_lag1_f + dose_lag2_f + dose_lag3_f + avg_dose_before_lag3"
    )
  }

  msm_formula <- stats::as.formula(paste("delta_outcome ~", rhs))

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
  data.frame(Term = rownames(coef_tab), coef_tab, row.names = NULL, check.names = FALSE)
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
    dose_history_levels = NULL
) {
  if (is.null(dose_history_levels)) dose_history_levels <- sort(unique(c(0, strategy_dose, baseline_dose)))
  visit_order <- seq_along(visits)

  tibble::tibble(
    visit = visits,
    visit_order = visit_order,
    strategy_dose = strategy_dose,
    strategy = paste0("Target ", strategy_dose, " mg"),
    outcome_0 = outcome_0_value
  ) |>
    dplyr::mutate(
      dose_lag1 = dplyr::if_else(visit_order == 1, baseline_dose, strategy_dose),
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
        TRUE ~ (baseline_dose + pmax(visit_order - 3, 0) * strategy_dose) / (visit_order - 2)
      )
    ) |>
    add_dose_history_factors(dose_history_levels = dose_history_levels)
}

predict_gee_ci <- function(model, newdata) {
  beta <- coef(model)
  V <- model$geese$vbeta
  dimnames(V) <- list(names(beta), names(beta))
  Terms <- stats::delete.response(stats::terms(model))

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
      matrix(0, nrow = nrow(X), ncol = length(missing_cols), dimnames = list(NULL, missing_cols))
    )
  }
  X <- X[, names(beta), drop = FALSE]
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

estimate_all_pairwise_strategy_contrasts <- function(
    model,
    target_visit = NULL,
    strategy_doses,
    baseline_dose = NULL
) {
  model_dat <- attr(model, "model_data")
  if (is.null(model_dat)) stop("Model is missing model_data attribute.")
  if (is.null(target_visit)) target_visit <- max(model_dat$visit, na.rm = TRUE)

  outcome_0_value <- mean(model_dat$outcome_0, na.rm = TRUE)
  dose_history_levels <- attr(model, "dose_history_levels")

  strategy_histories <- dplyr::bind_rows(
    lapply(strategy_doses, function(dose_value) {
      make_strategy_data(
        strategy_dose = dose_value,
        visits = seq_len(target_visit),
        outcome_0_value = outcome_0_value,
        baseline_dose = if (is.null(baseline_dose)) dose_value else baseline_dose,
        dose_history_levels = dose_history_levels
      ) |>
        dplyr::filter(visit == target_visit) |>
        dplyr::slice_tail(n = 1)
    })
  )

  beta <- coef(model)
  V <- tryCatch(vcov(model), error = function(e) model$geese$vbeta)
  dimnames(V) <- list(names(beta), names(beta))

  Terms <- stats::delete.response(stats::terms(model))
  X <- stats::model.matrix(Terms, strategy_histories)
  missing_cols <- setdiff(names(beta), colnames(X))
  if (length(missing_cols) > 0) {
    X <- cbind(X, matrix(0, nrow = nrow(X), ncol = length(missing_cols), dimnames = list(NULL, missing_cols)))
  }
  X <- X[, names(beta), drop = FALSE]

  pairs <- t(combn(seq_along(strategy_doses), 2))
  dplyr::bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
    lower_id <- pairs[i, 1]
    higher_id <- pairs[i, 2]
    lower_dose <- strategy_doses[lower_id]
    higher_dose <- strategy_doses[higher_id]
    contrast_vector <- X[higher_id, ] - X[lower_id, ]
    estimate <- as.numeric(sum(contrast_vector * beta))
    se <- sqrt(as.numeric(t(contrast_vector) %*% V %*% contrast_vector))
    tibble::tibble(
      contrast = paste0("Target ", higher_dose, " mg vs target ", lower_dose, " mg"),
      target_visit = target_visit,
      estimate = estimate,
      se = se,
      lower_95 = estimate - qnorm(0.975) * se,
      upper_95 = estimate + qnorm(0.975) * se,
      p_value = 2 * pnorm(abs(estimate / se), lower.tail = FALSE)
    )
  }))
}

################################################################################
######## TABLE: ACTIVE GEE MEANS AND OBSERVED PLACEBO MEANS BY WEEK #############
################################################################################

make_gee_and_placebo_mean_table <- function(
    prediction_results,
    placebo_observed_by_trial,
    placebo_counts_by_trial = NULL,
    max_visit = 8,
    weight_model = "Ordinal IPTW",
    repeat_placebo_by_active_arm = TRUE
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
  
  ##############################################################################
  # Active-arm counts from fitted MSM data
  ##############################################################################
  
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
  
  ##############################################################################
  # Active GEE/MSM predicted means
  ##############################################################################
  
  active_mean_table <- purrr::imap_dfr(
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
          active_counts_by_arm,
          by = c("arm_name", "visit")
        ) |>
        dplyr::transmute(
          weight_model = weight_model,
          trial_name = trial_name,
          arm_name = arm_name,
          week = visit,
          source = "Active GEE/MSM prediction",
          dose_strategy = as.character(strategy),
          strategy_dose = strategy_dose,
          
          mean_improvement = pred_delta_outcome,
          SE = se,
          lower_95 = lower_95,
          upper_95 = upper_95,
          
          n_patients = active_n_patients,
          n_visit_rows = active_n_visit_rows
        )
    }
  )
  
  ##############################################################################
  # Observed placebo means
  ##############################################################################
  
  placebo_for_table <- placebo_observed_by_trial |>
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
    ) |>
    dplyr::filter(
      visit >= 1,
      visit <= max_visit
    )
  
  if (isTRUE(repeat_placebo_by_active_arm)) {
    
    active_arm_trial_map <- active_mean_table |>
      dplyr::distinct(
        weight_model,
        trial_name,
        arm_name
      )
    
    placebo_mean_table <- active_arm_trial_map |>
      dplyr::left_join(
        placebo_for_table,
        by = "trial_name"
      ) |>
      dplyr::transmute(
        weight_model = weight_model,
        trial_name = trial_name,
        arm_name = arm_name,
        week = visit,
        source = "Observed placebo",
        dose_strategy = "Observed placebo",
        strategy_dose = NA_real_,
        
        mean_improvement = placebo_mean_delta,
        SE = placebo_se_delta,
        lower_95 = placebo_lower_95,
        upper_95 = placebo_upper_95,
        
        n_patients = placebo_n_patients,
        n_visit_rows = placebo_n_visit_rows
      )
    
  } else {
    
    placebo_mean_table <- placebo_for_table |>
      dplyr::transmute(
        weight_model = weight_model,
        trial_name = trial_name,
        arm_name = paste0(trial_name, "_PLACEBO"),
        week = visit,
        source = "Observed placebo",
        dose_strategy = "Observed placebo",
        strategy_dose = NA_real_,
        
        mean_improvement = placebo_mean_delta,
        SE = placebo_se_delta,
        lower_95 = placebo_lower_95,
        upper_95 = placebo_upper_95,
        
        n_patients = placebo_n_patients,
        n_visit_rows = placebo_n_visit_rows
      )
  }
  
  ##############################################################################
  # Combine active and placebo rows
  ##############################################################################
  
  dplyr::bind_rows(
    active_mean_table,
    placebo_mean_table
  ) |>
    dplyr::mutate(
      source_order = dplyr::case_when(
        source == "Observed placebo" ~ 0L,
        TRUE ~ 1L
      )
    ) |>
    dplyr::arrange(
      trial_name,
      arm_name,
      week,
      source_order,
      strategy_dose
    ) |>
    dplyr::select(-source_order)
}
# ---------------------------- multinomial sensitivity ----------------------------

make_multinom_dose_factor <- function(dose_f, ref_dose = NULL) {
  
  dose_num <- as.numeric(as.character(dose_f))
  dose_levels <- sort(unique(dose_num[!is.na(dose_num)]))
  
  out <- factor(
    as.character(dose_num),
    levels = as.character(dose_levels)
  )
  
  if (!is.null(ref_dose) && as.character(ref_dose) %in% levels(out)) {
    out <- stats::relevel(out, ref = as.character(ref_dose))
  }
  
  out
}

fit_iptw_multinom_denominator_model <- function(data, visit_df = 3, ref_dose = NULL) {
  model_dat <- data |>
    dplyr::filter(use_treatment_weight) |>
    dplyr::mutate(
      dose_multinom_f = make_multinom_dose_factor(dose_f, ref_dose = ref_dose),
      dose_lag1_f = droplevels(dose_lag1_f)
    )

  if (nlevels(model_dat$dose_multinom_f) < 2) stop("Multinomial denominator: fewer than 2 dose levels.")

  nnet::multinom(
    dose_multinom_f ~ rms::rcs(visit, visit_df) + delta_outcome + side.effects +
      dose_lag1_f + dose_lag1_f:delta_outcome + dose_lag1_f:side.effects + outcome_0,
    data = model_dat,
    trace = FALSE,
    maxit = 1000,
    MaxNWts = 10000,
    Hess = TRUE
  )
}

fit_iptw_multinom_numerator_model <- function(data, visit_df = 3, ref_dose = NULL) {
  model_dat <- data |>
    dplyr::filter(use_treatment_weight) |>
    dplyr::mutate(
      dose_multinom_f = make_multinom_dose_factor(dose_f, ref_dose = ref_dose),
      dose_lag1_f = droplevels(dose_lag1_f)
    )

  if (nlevels(model_dat$dose_multinom_f) < 2) stop("Multinomial numerator: fewer than 2 dose levels.")

  nnet::multinom(
    dose_multinom_f ~ rms::rcs(visit, visit_df) + dose_lag1_f + outcome_0,
    data = model_dat,
    trace = FALSE,
    maxit = 1000,
    MaxNWts = 10000,
    Hess = TRUE
  )
}

get_observed_dose_prob_multinom <- function(model, data, outcome_var = "dose_multinom_f") {
  prob_mat <- predict(model, newdata = data, type = "probs")
  if (is.null(dim(prob_mat))) {
    lev <- model$lev
    prob_mat <- cbind(1 - prob_mat, prob_mat)
    colnames(prob_mat) <- lev
  }
  observed_dose <- as.character(data[[outcome_var]])
  prob_mat[cbind(seq_len(nrow(prob_mat)), match(observed_dose, colnames(prob_mat)))]
}

add_iptw_treatment_weights_multinom <- function(data, denominator_model, numerator_model, min_prob = 1e-6) {
  data <- data |>
    dplyr::select(-dplyr::any_of(c(
      "p_dose_denominator_multinom", "p_dose_numerator_multinom",
      "SW_treatment_multinom", "cSW_treatment_multinom"
    )))

  weight_dat <- data |>
    dplyr::filter(use_treatment_weight) |>
    dplyr::mutate(
      dose_multinom_f = factor(as.character(dose_f), levels = denominator_model$lev),
      dose_lag1_f = droplevels(dose_lag1_f)
    )

  weight_dat$p_dose_denominator_multinom <- get_observed_dose_prob_multinom(denominator_model, weight_dat)
  weight_dat$p_dose_numerator_multinom <- get_observed_dose_prob_multinom(numerator_model, weight_dat)

  weight_dat <- weight_dat |>
    dplyr::mutate(
      p_dose_denominator_multinom = pmax(p_dose_denominator_multinom, min_prob),
      p_dose_numerator_multinom = pmax(p_dose_numerator_multinom, min_prob),
      SW_treatment_multinom = p_dose_numerator_multinom / p_dose_denominator_multinom
    ) |>
    dplyr::arrange(pid, visit) |>
    dplyr::group_by(pid) |>
    dplyr::mutate(cSW_treatment_multinom = cumprod(SW_treatment_multinom)) |>
    dplyr::ungroup()

  data |>
    dplyr::left_join(
      weight_dat |>
        dplyr::select(pid, visit, p_dose_denominator_multinom, p_dose_numerator_multinom, SW_treatment_multinom, cSW_treatment_multinom),
      by = c("pid", "visit")
    )
}

add_total_weights_multinom <- function(data) {
  data |>
    dplyr::mutate(
      SW_total_multinom = dplyr::case_when(
        use_msm & !is.na(cSW_treatment_multinom) & !is.na(cSW_censoring) ~ cSW_treatment_multinom * cSW_censoring,
        TRUE ~ NA_real_
      )
    )
}

truncate_total_weights_multinom <- function(data, lower = 0.01, upper = 0.99) {
  cutoffs <- quantile(data$SW_total_multinom, probs = c(lower, upper), na.rm = TRUE)
  data |>
    dplyr::mutate(
      SW_total_multinom_trunc = dplyr::case_when(
        is.na(SW_total_multinom) ~ NA_real_,
        SW_total_multinom < cutoffs[[1]] ~ as.numeric(cutoffs[[1]]),
        SW_total_multinom > cutoffs[[2]] ~ as.numeric(cutoffs[[2]]),
        TRUE ~ SW_total_multinom
      )
    )
}

summarise_iptw_multinom_treatment_weights <- function(data) {
  data |>
    dplyr::filter(use_treatment_weight, !is.na(cSW_treatment_multinom)) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_patients = dplyr::n_distinct(pid),
      mean_SW_treatment_multinom = mean(SW_treatment_multinom),
      sd_SW_treatment_multinom = sd(SW_treatment_multinom),
      min_SW_treatment_multinom = min(SW_treatment_multinom),
      p1_SW_treatment_multinom = quantile(SW_treatment_multinom, 0.01),
      p50_SW_treatment_multinom = quantile(SW_treatment_multinom, 0.50),
      p99_SW_treatment_multinom = quantile(SW_treatment_multinom, 0.99),
      max_SW_treatment_multinom = max(SW_treatment_multinom),
      mean_cSW_treatment_multinom = mean(cSW_treatment_multinom),
      sd_cSW_treatment_multinom = sd(cSW_treatment_multinom),
      min_cSW_treatment_multinom = min(cSW_treatment_multinom),
      p1_cSW_treatment_multinom = quantile(cSW_treatment_multinom, 0.01),
      p50_cSW_treatment_multinom = quantile(cSW_treatment_multinom, 0.50),
      p99_cSW_treatment_multinom = quantile(cSW_treatment_multinom, 0.99),
      max_cSW_treatment_multinom = max(cSW_treatment_multinom),
      ESS_cSW_treatment_multinom = sum(cSW_treatment_multinom)^2 / sum(cSW_treatment_multinom^2)
    )
}

summarise_total_multinom_weights <- function(data) {
  data |>
    dplyr::filter(use_msm, !is.na(SW_total_multinom_trunc)) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_patients = dplyr::n_distinct(pid),
      mean_SW_total_multinom_trunc = mean(SW_total_multinom_trunc),
      sd_SW_total_multinom_trunc = sd(SW_total_multinom_trunc),
      min_SW_total_multinom_trunc = min(SW_total_multinom_trunc),
      p1_SW_total_multinom_trunc = quantile(SW_total_multinom_trunc, 0.01),
      p50_SW_total_multinom_trunc = quantile(SW_total_multinom_trunc, 0.50),
      p99_SW_total_multinom_trunc = quantile(SW_total_multinom_trunc, 0.99),
      max_SW_total_multinom_trunc = max(SW_total_multinom_trunc),
      ESS_SW_total_multinom_trunc = sum(SW_total_multinom_trunc)^2 / sum(SW_total_multinom_trunc^2)
    )
}

make_multinom_coef_table <- function(model) {
  s <- summary(model)
  coef_mat <- s$coefficients
  se_mat <- s$standard.errors
  if (is.null(dim(coef_mat))) {
    coef_mat <- matrix(coef_mat, nrow = 1, dimnames = list(model$lev[2], names(coef_mat)))
    se_mat <- matrix(se_mat, nrow = 1, dimnames = list(model$lev[2], names(se_mat)))
  }
  dplyr::bind_rows(lapply(rownames(coef_mat), function(dose_level) {
    tibble::tibble(
      Dose_level_vs_reference = dose_level,
      Term = colnames(coef_mat),
      Estimate = as.numeric(coef_mat[dose_level, ]),
      SE = as.numeric(se_mat[dose_level, ]),
      z_value = Estimate / SE,
      p_value = 2 * pnorm(abs(z_value), lower.tail = FALSE)
    )
  }))
}

################################################################################
######################## SAVE PLOT LISTS AS JPG GRIDS ###########################
################################################################################

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}

dir.create("plot_grids", showWarnings = FALSE)
################################################################################
######################## SAVE PLOT LISTS AS JPG GRIDS ###########################
################################################################################

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}

save_plot_grid_jpg <- function(
    plot_list,
    file_path,
    ncol = 4,
    base_width = 5,
    base_height = 4,
    dpi = 300
) {
  
  plot_list <- plot_list[!vapply(plot_list, is.null, logical(1))]
  
  if (length(plot_list) == 0) {
    message("No plots to save for: ", file_path)
    return(NULL)
  }
  
  out_dir <- dirname(file_path)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  nrow <- ceiling(length(plot_list) / ncol)
  
  combined_plot <- patchwork::wrap_plots(
    plot_list,
    ncol = ncol
  )
  
  ggplot2::ggsave(
    filename = file_path,
    plot = combined_plot,
    width = ncol * base_width,
    height = nrow * base_height,
    units = "in",
    dpi = dpi,
    device = "jpeg",
    bg = "white",
    limitsize = FALSE
  )
  
  combined_plot
}
