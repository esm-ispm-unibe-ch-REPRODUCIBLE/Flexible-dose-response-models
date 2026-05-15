library(dplyr)
library(zoo)
### STUDIES with long format 
# variable visit indicates the visit time of the patient

### three studies comparing PAROXETINE vs PLACEBO
# trial 1
trial_1_before<- read.csv("O:/PETRUSHA - dose response/Saved Data/29060_002.csv")
table(trial_1_before$treat)
nrow(trial_1_before[trial_1_before$visit>0 & trial_1_before$treat=="PAROXETINE" & trial_1_before$dose==0,])

## Fix the dose, it should be 20 at visit=0 for PAROXETINE
trial_1_before$dose[trial_1_before$visit==0 & trial_1_before$treat=="PAROXETINE"]<-20
table(trial_1_before$dose[trial_1_before$visit==0])
table(trial_1_before$dose[trial_1_before$visit==0 & trial_1_before$treat=="PAROXETINE"])

## delete prebaseline visits
trial_1_before <- trial_1_before %>%
  dplyr::filter(!is.na(visit), visit >= 0)
#make days in weeks
trial_1_before <- trial_1_before %>%
  mutate(
    visit_day = visit,
    visit = case_when(
      visit_day == 0 ~ 0,
      visit_day >= 1  & visit_day <= 7  ~ 1,
      visit_day >= 8  & visit_day <= 14 ~ 2,
      visit_day >= 15 & visit_day <= 21 ~ 3,
      visit_day >= 22 & visit_day <=28 ~ 4,
      visit_day >= 29 & visit_day <= 35 ~ 5,
      visit_day >= 36 & visit_day <= 42 ~ 6,
      visit_day >= 43 & visit_day <= 49 ~ 7,
      visit_day >= 50 & visit_day <= 57 ~ 8,
      visit_day>57~ 9,
      TRUE ~ NA_real_
    ))
    
trial_1_before <- trial_1_before %>%
  filter(!is.na(visit)) %>%
  arrange(pid, visit, visit_day) %>%
  group_by(pid, visit) %>%
  slice_tail(n = 1) %>%
  ungroup()

trial_1_before <- trial_1_before %>%
  mutate(
    dose_original = dose,
    dose = if_else(dose == 10, 20, dose)
  )

### If dose becomes 0 in between visits dplyr::select one between the previous and the next one
set.seed(123)

trial_1_before <- trial_1_before %>%
  arrange(pid, visit) %>%                 # ensure visits are in order within pid
  group_by(pid) %>%
  mutate(
    # optional: enforce baseline dose = 20 for PAROXETINE
    dose = if_else(visit == 0 & treat == "PAROXETINE", 20, dose),
    
    # helper: NA instead of 0 so we can carry-forward / backward non-zero doses
    dose_nozero = if_else(dose == 0, NA_real_, dose),
    
    # previous non-zero dose within pid
    prev_dose = na.locf(dose_nozero, na.rm = FALSE),
    
    # next non-zero dose within pid
    next_dose = rev(na.locf(rev(dose_nozero), na.rm = FALSE))
  ) %>%
  rowwise() %>%
  mutate(
    dose = if (visit > 0 && treat == "PAROXETINE" && dose == 0) {
      candidates <- c(prev_dose, next_dose)
      candidates <- candidates[!is.na(candidates)]
      if (length(candidates) == 0) dose else candidates[sample.int(length(candidates), 1)]
    } else {
      dose
    }
  ) %>%
  ungroup() %>%
  dplyr::select(-dose_nozero, -prev_dose, -next_dose) %>%
  ## === here comes the sorting you want ===
  mutate(
    pid_num = as.numeric(sub(".*_", "", pid)),                 # numeric suffix
    treat   = factor(treat, levels = c("PAROXETINE", "PLACEBO"))
  ) %>%
  arrange(treat, pid_num, visit) %>%
  dplyr::select(-pid_num)

## Create a variable for side effects. It needs to be between 0 and 10 and been associated (weakly) with dose
### lower dose means higher side effects. It needs to be weak and with random noise
set.seed(2025)  # for reproducibility

# dose range for positive doses
d_range <- range(trial_1_before$dose[trial_1_before$dose >= 0], na.rm = TRUE)

trial_1 <- trial_1_before %>%
  mutate(
    # scale dose to [0,1] for positive doses
    dose_scaled = if_else(
      dose >= 0,
      (dose - d_range[1]) / (d_range[2] - d_range[1]),
      NA_real_
    ),
    # inverse: lower dose -> higher "effect"
    inv_effect = if_else(dose > 0, 1 - dose_scaled, NA_real_),
    
    # add random noise (weak association)
    side_raw = inv_effect + rnorm(n(), mean = 0, sd = 0.35),
    
    # clamp to [0,1]
    side_scaled = pmin(pmax(side_raw, 0), 1),
    
    # continuous 0–10
    side_cont = 10 * side_scaled,
    
    # integer 0–10
    side_int = as.integer(round(side_cont)),
    
    # final variable: integer, NA for visit <= 0
    side.effects = if_else(visit > 0, side_int, NA_integer_)
  ) %>%
  dplyr::select(-dose_scaled, -inv_effect, -side_raw, -side_scaled, -side_cont, -side_int)



#check 
ggplot(trial_1, aes(x = side.effects, y = dose)) +
      geom_jitter(width = 0.25, height = 1.5, alpha = 0.5) +
      geom_smooth(method = "lm", se = FALSE) +
      labs(
            x = "Side effects",
            y = "Dose"
       )


# trial 2 - trial 2 is multi-arm for now I excluded the IMIPRAMINE arm and included only PAROXETINE vs PLACEBO
trial_2<- read.csv("O:/PETRUSHA - dose response/Saved Data/29060_003.csv")
table(trial_2$treat)
trial_2<-trial_2[trial_2$treat!="IMIPRAMINE",]
table(trial_2$treat)
# trial 3 - trial 3 is multi-arm for now I excluded the AMITRIPTYLINE arm and included only PAROXETINE vs PLACEBO
trial_3<- read.csv("O:/PETRUSHA - dose response/Saved Data/29060_007.csv")
table(trial_3$treat)
trial_3<-trial_3[trial_3$treat!="AMITRIPTYLINE",]
table(trial_3$treat)
