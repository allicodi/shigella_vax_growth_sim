
get_baseline_growth <- function(dose_schedule = "6mo",
                                enroll_site = "Overall",
                                enr_haz_var_name = "enr_haz", 
                                age_var_name = "enr_age_months"){
  
  # Estimate from EFGH data in age group min_age_bl to max_age_bl
  data <- readRDS(here::here("data/efgh/efgh_data.Rds"))
  
  if(dose_schedule == "6mo"){
    data <- data[data[[age_var_name]] >= 6 & data[[age_var_name]] <= 8,]
  } else{
    data <- data[data[[age_var_name]] >= 9 & data[[age_var_name]] < 12,]
  }
  
  # If site specific analysis
  if(enroll_site != "Overall"){
    data <- data[data$enroll_site == enroll_site,]
  }
  
  mean_enr_haz <- mean(data[[enr_haz_var_name]], na.rm = TRUE)
  sd_enr_haz <- sd(data[[enr_haz_var_name]], na.rm = TRUE)
  
  return(data.frame(mean_enr_haz = mean_enr_haz,
                    sd_enr_haz = sd_enr_haz))
  
}

get_effect_baseline_growth_sev_inf <- function(dose_schedule = "6mo",
                                               enroll_site = "Overall",
                                               enr_haz_var_name = "enr_haz", 
                                               inf_var_name = "positive_tac_or_culture",
                                               severity_var_name = "gems_msd_bin",
                                               age_var_name = "enr_age_months",
                                               site_var_name = "enroll_site",
                                               site_adj = TRUE){
  
  # Estimate from EFGH data in age group min_age_bl to max_age_bl
  data <- readRDS(here::here("data/efgh/efgh_data.Rds"))
  data$gems_msd_bin <- ifelse(data$gems_msd == "Less-severe", 0, 1)
  
  if(dose_schedule == "6mo"){
    data <- data[data[[age_var_name]] >= 6 & data[[age_var_name]] <= 8,]
  } else{
    data <- data[data[[age_var_name]] >= 9 & data[[age_var_name]] < 12,]
  }
  
  # Fit Severe Infection ~ Baseline Growth
  # MSD ~ BL Growth | Infected
  inf_data <- data[data[[inf_var_name]] == 1, ]
  
  # if dose schedule = 6mo, fit model 6-12 & 12-18
  # if dose schedule = 12mo, fit model 12-18 & 18-24
  
  if(site_adj){
    # depending on the dose schedule needs to subset to 6-18 or 12 - 24 
    fit.sevinf.enr_haz <- glm(as.formula(paste0(severity_var_name, " ~ ", enr_haz_var_name, " + ", site_var_name)), 
                              data = inf_data,
                              family = "binomial")
    
    sev_coef <- coef(fit.sevinf.enr_haz)[[enr_haz_var_name]]
    
    # just use intercept for the reference site as is if overall (assume they're all approx same)
    if(enroll_site == "Overall"){
      sev_int <- coef(fit.sevinf.enr_haz)[1]
    } else{
      sev_int <- coef(fit.sevinf.enr_haz)[[paste0("enroll_site",enroll_site)]]
    }
    
  } else{
    fit.sevinf.enr_haz <- glm(as.formula(paste0(severity_var_name, " ~ ", enr_haz_var_name)), 
                              data = inf_data,
                              family = "binomial")
    
    sev_coef <- coef(fit.sevinf.enr_haz)[[enr_haz_var_name]]
    sev_int <- coef(fit.sevinf.enr_haz)[1]
  }
  
  return(list(sev_coef = sev_coef,
              sev_int = sev_int))
  
}

get_incidence <- function(dose_schedule = "6mo",
                          enroll_site = "Total"){
  
  efgh_mad <- readRDS(here::here("data/efgh/boot_all_agegrps_2025-06-24.RDS"))
  efgh_sev <- readRDS(here::here("data/efgh/boot_all_agegrps_msd_2025-07-30.RDS"))
  
  # All diarrhea
  mad_inc <- efgh_mad[efgh_mad$country == enroll_site &
                        efgh_mad$tac_culture == "tac" &
                        efgh_mad$inc_type == "quad",]
  
  msd_inc <- efgh_sev[efgh_sev$country == enroll_site &
                        efgh_sev$tac_culture == "tac" &
                        efgh_sev$inc_type == "quad",]
  
  if(dose_schedule == "6mo"){
    mad_inc_0_6 <- mean(as.numeric(mad_inc$est[mad_inc$agegrp == "6-8 months"]),
                        as.numeric(mad_inc$est[mad_inc$agegrp == "9-11 months"])) / 100 / 2
    mad_inc_6_12 <- as.numeric(mad_inc$est[mad_inc$agegrp == "12-17 months"]) / 100 / 2
    
    msd_inc_0_6 <- mean(as.numeric(msd_inc$est[mad_inc$agegrp == "6-8 months"]),
                        as.numeric(msd_inc$est[mad_inc$agegrp == "9-11 months"])) / 100 / 2
    msd_inc_6_12 <- as.numeric(msd_inc$est[msd_inc$agegrp == "12-17 months"]) / 100 / 2
    
  } else{
    mad_inc_0_6 <- as.numeric(mad_inc$est[mad_inc$agegrp == "12-17 months"]) / 100 / 2
    mad_inc_6_12 <- as.numeric(mad_inc$est[mad_inc$agegrp == "18-23 months"]) / 100 / 2
    
    msd_inc_0_6 <- as.numeric(msd_inc$est[mad_inc$agegrp == "12-17 months"]) / 100 / 2
    msd_inc_6_12 <- as.numeric(msd_inc$est[msd_inc$agegrp == "18-23 months"]) / 100 / 2
  }

  return(data.frame(mad_inc_0_6 = mad_inc_0_6,
                    mad_inc_6_12 = mad_inc_6_12,
                    msd_inc_0_6 = msd_inc_0_6,
                    msd_inc_6_12 = msd_inc_6_12))
}

