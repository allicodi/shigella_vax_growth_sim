# ---------------------------------------------------------------------------
# Function to estimate hazard ratio for baseline growth --> Shigella in MAL-ED
# ---------------------------------------------------------------------------

library(survival)

hr_maled <- function(dose_schedule = "6mo"){
  # Load data 
  tac_data <- read.csv(here::here("data/maled/Tac_Sep2018.csv")) %>%
    filter(country_id != "PK")
  growth_data <- readxl::read_xlsx(here::here("data/maled/growth.xlsx")) %>%
    filter(country_id != "PK")
  
  # based on dose schedule get hazard ratio in MAL-ED
  # if dose schedule finish by 6mo:
  #   bl_haz = height closest to 6mo
  #   fit cox model for 6mo growth measurement start time, shigella 6-12 mo (censor at 6mo post 6mo growth measurement) + stratify hazard by site
  #   fit cox model for 6mo growth measurement start time, shigella 12-18 mo only in kids who survived model 1, 6mo growth measurement censor at 12mo past 6mo growth measurement and stratify hazard by site
  
  # if dose schedule finish by 12mo:
  #   bl_haz = height closest to 12mo
  #   fit cox model for 12mo growth measurement start time, shigella 12-18 mo (censor at 6mo post 12mo growth measurement) + stratify hazard by site
  #   fit cox model for 12mo growth measurement start time, shigella 18-24 mo only in kids who survived model 1, 12mo growth measurement censor at 12mo past 12mo growth measurement and stratify hazard by site
  
  # results in two coef for HAZ feed into parameter grid below & find intercept 
  
  if(dose_schedule == "6mo"){
    bl_mo <- 6
    max_mo_1 <- 12
    max_mo_2 <- 18
  } else{
    bl_mo <- 12
    max_mo_1 <- 18
    max_mo_2 <- 24
  }
  
  # Make baseline data
  bl_data <- growth_data %>%
    filter(month_ss == bl_mo) %>%
    arrange(pid, date) %>%  # optional: ensures ordering
    distinct(pid, .keep_all = TRUE) %>%  # keep first per pid
    filter(!is.na(zlen)) %>%
    select(pid, country_id, date, zlen) %>%
    rename('bl_haz' = zlen,
           'bl_date' = date) 
  
  # Get censoring date for all (date of measurement closest to month max_mo_1)
  end_date_1 <- growth_data %>%
    filter(month_ss > bl_mo & month_ss <= max_mo_1) %>%
    filter(!is.na(date)) %>%
    arrange(pid, desc(month_ss), date) %>%  # prioritize higher month_ss first
    distinct(pid, .keep_all = TRUE) %>%     # keep the latest available month ≤ 18
    select(pid, date) %>%
    rename(cens_date_1 = date)
  
  # Get censoring date for all (date of measurement closest to month 18)
  end_date_2 <- growth_data %>%
    filter(month_ss > bl_mo & month_ss <= max_mo_2) %>%
    filter(!is.na(date)) %>%
    arrange(pid, desc(month_ss), date) %>%  # prioritize higher month_ss first
    distinct(pid, .keep_all = TRUE) %>%     # keep the latest available month ≤ 18
    select(pid, date) %>%
    rename(cens_date_2 = date)
  
  # Get diarrheal stool data month bl to month max_mo_2
  stool_data <- tac_data %>%
    filter(stooltype == "D1", month_ss >= bl_mo, month_ss <= max_mo_2) %>%
    select(pid, dob, date, shigella_eiec, shigella_eiec_afe) %>%
    mutate(shigella_pos_overall = if_else(shigella_eiec_afe > 0.5, 1, 0)) %>%
    filter(shigella_pos_overall == 1) %>%               # Keep only rows with positive Shigella
    group_by(pid) %>%
    slice_min(order_by = date, n = 1, with_ties = FALSE) %>%  # Keep first positive episode
    ungroup() %>%
    rename('diar_date' = date) %>%
    mutate(
      dob = as.Date(dob, format = "%m/%d/%Y"),
      diar_date = as.Date(diar_date, , format = "%m/%d/%Y")
    )
  
  # Join data & make dataset for first 6mo
  data_1 <- left_join(bl_data, end_date_1, by = "pid") %>%
    left_join(end_date_2, by = "pid") %>%
    left_join(stool_data, by = "pid") %>%
    mutate(shigella_pos_overall = if_else(is.na(shigella_pos_overall), 0, shigella_pos_overall),
           shigella_pos_1 = if_else(!is.na(diar_date) & diar_date <= cens_date_1, 1, 0), 
           ftime1 = if_else(!is.na(diar_date) & diar_date <= cens_date_1,
                           as.numeric(difftime(diar_date, bl_date, units = "weeks")),
                           as.numeric(difftime(cens_date_1, bl_date, units = "weeks")))) %>%
    filter(!is.na(ftime1)) 
  
  # Make dataset for first 6mo, excluding people who had shigella or dropped out of study in the first six months
  data_2 <- data_1 %>%
    filter(!(shigella_pos_1 == 1)) %>%
    mutate(shigella_pos_2 = if_else(!is.na(diar_date) & diar_date <= cens_date_2, 1, 0),
           ftime2 = if_else(!is.na(diar_date) & diar_date <= cens_date_2,
                            as.numeric(difftime(diar_date, bl_date, units = "weeks")),
                            as.numeric(difftime(cens_date_2, bl_date, units = "weeks")))) %>%
    filter(!(cens_date_1 == cens_date_2))
  
  
  # Cox model for first 6 months, adjust for country_id
  cox_model_1 <- coxph(
    Surv(time = ftime1, event = shigella_pos_1) ~ bl_haz + strata(country_id),
    data = data_1
  )
  
  # Cox model for last 6 months, adjust for country id
  cox_model_2 <- coxph(
    Surv(time = ftime2, event = shigella_pos_2) ~ bl_haz + strata(country_id),
    data = data_2
  )
  
  # ^ note some of these are opposite the expected effect (hazard > 1... but not significant)
  # so not sure what to do there
  
  return(data.frame(coef_0_6 = cox_model_1$coefficients['bl_haz'],
                    hazard_0_6 = exp(cox_model_1$coefficients['bl_haz']),
                    coef_6_12 = cox_model_2$coefficients['bl_haz'],
                    hazard_6_12 = exp(cox_model_2$coefficients['bl_haz'])))
  
}
