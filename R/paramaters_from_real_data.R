
get_baseline_growth <- function(data, 
                                enr_haz_var_name = "enr_haz", 
                                age_var_name = "enr_age_months",
                                min_age = 9, 
                                max_age = 12,
                                enroll_site = "Total"){
  
  # Subset data to age range if applicable; default age 9-12 months
  if(!is.na(min_age)){
    data <- data[data[[age_var_name]] >= min_age,]
  }
  
  if(!is.na(max_age)){
    data <- data[data[[age_var_name]] <= max_age,]
  }
  
  # If site specific analysis
  if(enroll_site != "Total"){
    data <- data[data$enroll_site == enroll_site,]
  }
  
  mean_enr_haz <- mean(data[[enr_haz_var_name]], na.rm = TRUE)
  sd_enr_haz <- sd(data[[enr_haz_var_name]], na.rm = TRUE)
  
  return(data.frame(mean_enr_haz = mean_enr_haz,
                    sd_enr_haz = sd_enr_haz))
  
}

get_effect_baseline_growth_inf <- function(data,
                                           enr_haz_var_name = "enr_haz", 
                                           inf_var_name = "shigella_attributable",
                                           severity_var_name = "gems_msd",
                                           age_var_name = "enr_age_months",
                                           site_var_name = "enroll_site",
                                           site_adj = TRUE,
                                           min_age = 9, 
                                           max_age = 12,
                                           enroll_site = "Total"){
  
  # ORIGINAL SUBSETTING
  # Subset data to age range if applicable; default age 9-12 months
  if(!is.na(min_age)){
    data <- data[data[[age_var_name]] >= min_age,]
  }

  if(!is.na(max_age)){
    data <- data[data[[age_var_name]] <= max_age,]
  }
  
  # NEW???
  # if(!is.na(max_age)){
  #   data <- data[data[[age_var_name]] >= max_age,]
  # }

  # Fit Infection ~ Baseline Growth
  if(site_adj){
    fit.inf.enr_haz <- glm(as.formula(paste0(inf_var_name, " ~ 0 + ", enr_haz_var_name, " + ", site_var_name)), 
                           data = data,
                           family = "binomial")
    
    inf_coef <- coef(fit.inf.enr_haz)[[enr_haz_var_name]]
    
    # QUESTION taking the intercept in general is site specific, 
    # so center at 0 to have int for each site, then avg by site if not doing site specific analysis?? idk 
    if(enroll_site == "Total"){
      inf_int <- mean(coef(fit.inf.enr_haz)[2:length(unique(data$enroll_site))+1])
    } else{
      inf_int <- coef(fit.inf.enr_haz)[[paste0("enroll_site",enroll_site)]]
    }
    
  } else{
    fit.inf.enr_haz <- glm(as.formula(paste0(inf_var_name, " ~ ", enr_haz_var_name)), 
                           data = data,
                           family = "binomial")
    
    inf_coef <- coef(fit.inf.enr_haz)[[enr_haz_var_name]]
    inf_int <- coef(fit.inf.enr_haz)[1]
  }
  

  # Fit Severe Infection ~ Baseline Growth
  # MSD ~ BL Growth | Infected
  inf_data <- data[data[[inf_var_name]] == 1, ]
  
  if(site_adj){
    fit.sevinf.enr_haz <- glm(as.formula(paste0(severity_var_name, " ~ 0 + ", enr_haz_var_name, " + ", site_var_name)), 
                              data = inf_data,
                              family = "binomial")
    
    sev_coef <- coef(fit.sevinf.enr_haz)[[enr_haz_var_name]]
    
    # QUESTION taking the intercept in general is site specific, 
    # so center at 0 to have int for each site, then avg by site if not doing site specific analysis?? idk 
    if(enroll_site == "Total"){
      sev_int <- mean(coef(fit.sevinf.enr_haz)[2:length(unique(data$enroll_site))+1])
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
  
  return(list(inf_coef = inf_coef,
              inf_int = inf_int,
              sev_coef = sev_coef,
              sev_int = sev_int))

}


# ----------------------------------------------------------------------------
# Explore parameters (will be automated in simulation code)

# NOTE this is the cleaned meta-analysis no etiology data
# n = 9131 (eliminated people without TAC results)
# is it ok to use this? or should i go back to rawest?

# data <- readRDS(here::here("data/efgh/efgh_data.Rds"))
# # 
# data$gems_msd_bin <- ifelse(data$gems_msd == "Less-severe", 0, 1)
# # 
# # get_baseline_growth(data, 
# #                     enr_haz_var_name = "enr_haz", 
# #                     age_var_name = "enr_age_months", 
# #                     min_age = 9, 
# #                     max_age = 12)
# 
# # Mean enr_haz: -0.94
# # SD enr_haz: 1.17
# 
# get_effect_baseline_growth_inf(data,
#                                enr_haz_var_name = "enr_haz",
#                                inf_var_name = "positive_tac_or_culture",
#                                severity_var_name = "gems_msd_bin",
#                                age_var_name = "enr_age_months",
#                                site_var_name = "enroll_site",
#                                site_adj = TRUE,
#                                min_age = 9,
#                                max_age = 12)
# 
# # Adjusting for site:
# # Infection: -0.064
# # Severe infection: -0.222
# 
# # Not adjusting for site:
# # Infection: -0.119
# # Severe infection: -0.145
# 
# # ----------------------------------------------------------------------------
# # Incidence from EFGH final data adjusted estimates
# 
# # Overall, all MAD, quadrivalent (2a, 3a, 6, S. sonnei)
# efgh_age <- readRDS(here::here("data/efgh/boot_all_agegrps_2025-06-24.RDS"))
# 
# mad_inc <- as.numeric(efgh_age$est[efgh_age$country == "Total" &
#                                      efgh_age$tac_culture == "tac" &
#                                      efgh_age$inc_type == "quad" &
#                                      efgh_age$agegrp == "6-35 months"])
# 
# # final est = 6.72
# 
# # Overall, GEMS MSD, quadrivalent (2a, 3a, 6, S. sonnei)
# efgh_sev <- readRDS(here::here("data/efgh/boot_severity_subsets_2025-06-24.RDS"))
# 
# msd_inc <- as.numeric(efgh_sev$est[efgh_sev$country == "Total" &
#                           efgh_sev$tac_culture == "tac" &
#                           efgh_sev$inc_type == "quad" & 
#                           efgh_sev$sevgrp == "gems_msd"]) 
# 
# # final est = 3.18
