# ------------------------------------------------------------------------------
# Function to simulate data based on geometric distribution, other misc changes
# ------------------------------------------------------------------------------

here::i_am("R/simulate_data_short_term.R")

source(here::here("R/paramaters_from_real_data.R"))
source(here::here("R/tune_parameters.R"))

# Data structure: 

# X = covariate (baseline HAZ)
# Z = treatment (vaccine)
# S_inf = infection (Shigella)
# S_sev = severe infection (moderate to severe Shigella)
# Y = outcome (lfazd90)

simulate_data <- function(n = 1e6,
                          seed = 12345,
                          min_age_bl = 9,
                          max_age_bl = 12,
                          VE_mild = 0.4,
                          VE_severe = 0.6,
                          effect_lsd_inf_growth = -0.056,
                          effect_msd_inf_growth = -0.089,
                          enroll_site = "Total"){
  set.seed(seed)
  # ---------------------------------------------------------------------------
  # Z: Vaccination (RCT 1:1) --------------------------------------------------
  # ---------------------------------------------------------------------------
  Z <- rbinom(n, size=1, prob=0.5)
  
  # ---------------------------------------------------------------------------
  # X: Baseline HAZ -----------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  # Estimate from EFGH data in age group min_age_bl to max_age_bl
  efgh_data <- readRDS(here::here("data/efgh/efgh_data.Rds"))
  efgh_data$gems_msd_bin <- ifelse(efgh_data$gems_msd == "Less-severe", 0, 1)
  
  bl_growth_param <- get_baseline_growth(data = efgh_data, 
                                         enr_haz_var_name = "enr_haz", 
                                         age_var_name = "enr_age_months", 
                                         min_age = 9, 
                                         max_age = 12,
                                         enroll_site = enroll_site)
  
  # Simulate growth based on real EFGH data
  X <- rnorm(n, mean = bl_growth_param$mean, sd = bl_growth_param$sd)
  
  # ---------------------------------------------------------------------------
  # S: Shigella infection (S_inf / S_sev for regular and severe; S_inf_time) --
  # ---------------------------------------------------------------------------
  
  # 1) Incidence of Shigella: MAD & GEMS MSD from EFGH data -------------------
  
  # All MAD, TAC attributable, quadrivalent targets (2a, 3a, 6, S. sonnei)
  efgh_age <- readRDS(here::here("data/efgh/boot_all_agegrps_2025-06-24.RDS"))
  
  # Reported in incidence / 100 child years, get in years
  mad_inc <- as.numeric(efgh_age$est[efgh_age$country == enroll_site &
                                       efgh_age$tac_culture == "tac" &
                                       efgh_age$inc_type == "quad" &
                                       efgh_age$agegrp == "6-35 months"]) / 100
  
  # GEMS MSD, TAC attributable, quadrivalent targets (2a, 3a, 6, S. sonnei)
  efgh_sev <- readRDS(here::here("data/efgh/boot_severity_subsets_2025-06-24.RDS"))
  
  msd_inc <- as.numeric(efgh_sev$est[efgh_sev$country == enroll_site &
                                       efgh_sev$tac_culture == "tac" &
                                       efgh_sev$inc_type == "quad" & 
                                       efgh_sev$sevgrp == "gems_msd"]) / 100
  
  # 2) Effect of baseline HAZ on infection from EFGH data ---------------------
  
  # question abt min/max age subsetting
  effect_enr_haz_inf <- get_effect_baseline_growth_inf(data = efgh_data, 
                                                       enr_haz_var_name = "enr_haz",
                                                       inf_var_name = "positive_tac_or_culture",
                                                       severity_var_name = "gems_msd_bin",
                                                       age_var_name = "enr_age_months",
                                                       site_var_name = "enroll_site",
                                                       site_adj = TRUE,
                                                       min_age = min_age_bl,
                                                       max_age = max_age_bl,
                                                       enroll_site = enroll_site)
  
  # 3) Tune intercept & LAZ coefficient for appropriate cumulative incidence, 
  #    and get hazard & risk ratios at those parameters
  
  ### Infection (general medically attended diarrhea) 
  # QUESTION incidence is never remotely close unless drop intercept by a lot...
  # is that because hazard has t0=12 monthly? 
  intercept_S_inf_range <- seq(effect_enr_haz_inf$inf_int - 5, 
                               effect_enr_haz_inf$inf_int + 0.5, 
                               by = 0.01)
  coef_S_inf_range <- seq(effect_enr_haz_inf$inf_coef - 0.1, 
                         effect_enr_haz_inf$inf_coef + 0.1, 
                         by = 0.01)
  
  S_inf_grid <- expand.grid(intercept = intercept_S_inf_range,
                            laz_coef = coef_S_inf_range,
                            mean_X = bl_growth_param$mean)
  
  # Get hazard for LAZ -0.5 & LAZ  + 0.5 --> hazard ratio for each intercept / LAZ coefficient combination 
  hazard_ratio_combos_inf <- cbind(S_inf_grid, do.call(rbind, apply(S_inf_grid, 1, hazard_ratio)))
  # Use hazard results to iterate over again for risk ratio
  risk_ratio_combos_inf <- cbind(hazard_ratio_combos_inf, do.call(rbind, apply(hazard_ratio_combos_inf, 1, risk_ratio)))
  
  # Tune such that cumulative incidence is comparable to mad_inc/year MAD
  rr_params_inf <- risk_ratio_combos_inf %>%
    filter(cum_inc_laz_mean > mad_inc - 0.005 & cum_inc_laz_mean < mad_inc + 0.005) %>% # cumulative incidence within +-0.005 of cumulative incidence from data
    mutate(dist = sqrt((cum_inc_laz_mean - mad_inc)^2 + # get distance from cum incidence & mean X, find row with the smallest distance
                         (laz_coef - effect_enr_haz_inf$inf_coef)^2)) %>%
    filter(dist == min(dist)) # select combination with minimum distance
  
  ### Moderate to severe diarrhea (MSD) infection
  # same issues as above with intercept except more extreme...
  intercept_S_sev_range <- seq(effect_enr_haz_inf$sev_int - 8, 
                               effect_enr_haz_inf$sev_int + 0.5, 
                               by = 0.01)
  coef_S_sev_range <- seq(effect_enr_haz_inf$sev_coef - 0.1, 
                          effect_enr_haz_inf$sev_coef + 0.1, 
                          by = 0.01)
  
  S_sev_grid <- expand.grid(intercept = intercept_S_sev_range,
                            laz_coef = coef_S_sev_range,
                            mean_X = bl_growth_param$mean)
  
  # Get hazard for LAZ -0.5 & LAZ  + 0.5 --> hazard ratio for each intercept / LAZ coefficient combination 
  hazard_ratio_combos_sev <- cbind(S_sev_grid, do.call(rbind, apply(S_sev_grid, 1, hazard_ratio)))
  # Use hazard results to iterate over again for risk ratio
  risk_ratio_combos_sev <- cbind(hazard_ratio_combos_sev, do.call(rbind, apply(hazard_ratio_combos_sev, 1, risk_ratio)))
  
  # Tune such that cumulative incidence is comparable to msd_inc/year MSD
  rr_params_sev <- risk_ratio_combos_sev %>%
    filter(cum_inc_laz_mean > msd_inc - 0.005 & cum_inc_laz_mean < msd_inc + 0.005) %>% # cumulative incidence within +-0.005 of cumulative incidence from data
    mutate(dist = sqrt((cum_inc_laz_mean - msd_inc)^2 + # get distance from cum incidence & mean X, find row with the smallest distance
                         (laz_coef - effect_enr_haz_inf$sev_coef)^2)) %>%
    filter(dist == min(dist)) # select combination with minimum distance
  
  # 4) Use the conditional hazard with intercept + coef found above to model infection / severe infection time??? in absence of vaccine?
  # help idk what im doing
  
  hazard_inf_i <- hazard(intercept = rr_params_inf$intercept, 
                         laz_coef = rr_params_inf$laz_coef,
                         laz = X)
  
  S_inf_time_Z0 <- rgeom(n, prob = hazard_inf_i)
  #S_inf_time_Z0 <- ifelse(S_inf_time_Z0 <= 12, S_inf_time_Z0, 0)
  # changed hazard ratio function t0 to 52? so it's in weeks?
  S_inf_time_Z0 <- ifelse(S_inf_time_Z0 <= 52, S_inf_time_Z0, 0) # if infected after study end, 0
  S_inf_Z0 <- ifelse(S_inf_time_Z0 != 0, 1, 0)
  
  hazard_sev_i <- rep(0,n)
  hazard_sev_i[S_inf_Z0 == 1] <- hazard(intercept = rr_params_sev$intercept,
                                     laz_coef = rr_params_sev$laz_coef,
                                     laz = X[S_inf_Z0 == 1])
  
  S_sev_Z0 <- rbinom(n, 1, prob = hazard_sev_i / hazard_inf_i)
  S_sev_Z0 <- ifelse(is.na(S_sev_Z0), 0, S_sev_Z0)
  
  # 5) Now add in vaccine? Just doing the way we were before, not sure if this should impact hazard, timing, etc 
  #    or if that was the more complex way with additional strata
  prob_S_inf_Z1 <- rep(0,n)
  prob_S_inf_Z1[S_inf_Z0 == 1] <- 1 - VE_mild
  S_inf_Z1 <- rbinom(n, 1, prob_S_inf_Z1)
  
  S_inf_time_Z1 <- ifelse(S_inf_Z0 == 1 & S_inf_Z1 == 0, 0, S_inf_time_Z0)
  
  prob_S_sev_Z1 <- rep(0, n)
  prob_S_sev_Z1[S_sev_Z0 == 1 & S_inf_Z1 == 1] <- (1 - VE_severe) / (1 - VE_mild)
  S_sev_Z1 <- rbinom(n, 1, prob_S_sev_Z1)
  
  # ---------------------------------------------------------------------------
  # X_1 through X_12: Monthly HAZ outcome 
  # ---------------------------------------------------------------------------
  
  maled_data <- readRDS(here::here("data/maled/growth_wide_maled.Rds"))
  
  # QUESTION subset to bangladesh?? that's what we did before but not sure why
  # also get rid of pakistan? 
  
  maled_data <- maled_data[maled_data$country_id == "BG",]
  
  # Model from MAL-ED
  get_sd_monthly <- function(y_idx, x_idx, data){
    data <- data[,c(x_idx, y_idx)] %>% drop_na()
    
    fit <- glm(data[[2]] ~ data[[1]], data = data, family = "gaussian")
    rss_fit <- sum((fit$residuals)^2)
    sd_fit <- sqrt(rss_fit / nrow(data))
    return(data.frame(y = names(data)[2],
                      x = names(data)[1],
                      beta_0 = fit$coefficients[1],
                      beta_1 = fit$coefficients[2],
                      sd = sd_fit))
  }
  
  # get sd for each month to month combo
  cols <- seq(3, 14)
  sd_monthly <- data.frame()
  for(i in cols){
    sd_monthly <- rbind(sd_monthly, get_sd_monthly(i + 1, i, maled_data))
  }
  
  # Get infection adjustments from meta-analysis
  # 'catch up growth' from other version? aka subtract from one timepoint only?
  growth_adj_Z0 <- ifelse(S_sev_Z0, effect_msd_inf_growth,
                          ifelse(S_inf_Z0, effect_lsd_inf_growth, 0))
  growth_adj_Z1 <- ifelse(S_sev_Z1, effect_msd_inf_growth,
                          ifelse(S_inf_Z1, effect_lsd_inf_growth, 0))
  
  # Simulate monthly growth based on sd_monthly
  
  # noise first so same for each Z0 and Z1
  noise_matrix <- sapply(1:12, function(i) rnorm(n, mean = 0, sd = sd_monthly$sd[i]))
  
  X_vec_Z0 <- vector(mode = "list", length = 12) 
  X_vec_Z0[[1]] <- sd_monthly$beta_0[1] + sd_monthly$beta_1[1] * X + as.numeric(S_inf_time_Z0 == 1)*growth_adj_Z0 + noise_matrix[,1]
  for(i in 2:length(X_vec_Z0)){
    X_vec_Z0[[i]] <- sd_monthly$beta_0[i] + sd_monthly$beta_1[i] * X + as.numeric(S_inf_time_Z0 == i)*growth_adj_Z0 + noise_matrix[,i]
  }
  
  X_df_Z0 <- data.frame(X_vec_Z0)
  colnames(X_df_Z0) <- c("X_1_Z0", "X_2_Z0", "X_3_Z0", "X_4_Z0", "X_5_Z0", "X_6_Z0",
                      "X_7_Z0", "X_8_Z0", "X_9_Z0", "X_10_Z0", "X_11_Z0", "X_12_Z0")
  
  X_vec_Z1 <- vector(mode = "list", length = 12) 
  X_vec_Z1[[1]] <- sd_monthly$beta_0[1] + sd_monthly$beta_1[1] * X + as.numeric(S_inf_time_Z1 == 1)*growth_adj_Z1 + noise_matrix[,1]
  for(i in 2:length(X_vec_Z0)){
    X_vec_Z1[[i]] <- sd_monthly$beta_0[i] + sd_monthly$beta_1[i] * X + as.numeric(S_inf_time_Z1 == i)*growth_adj_Z1 + noise_matrix[,i]
  }
  
  X_df_Z1 <- data.frame(X_vec_Z1)
  colnames(X_df_Z1) <- c("X_1_Z1", "X_2_Z1", "X_3_Z1", "X_4_Z1", "X_5_Z1", "X_6_Z1",
                         "X_7_Z1", "X_8_Z1", "X_9_Z1", "X_10_Z1", "X_11_Z1", "X_12_Z1")
  
  # ---------------------------------------------------------------------------
  # Construct final dataset, paramaters, effects 
  # ---------------------------------------------------------------------------
  
  # Counterfactual
  truth_data <- data.frame(id = 1:n,
                         X = X,
                         Z = Z,
                         S_inf_Z0 = S_inf_Z0,
                         S_sev_Z0 = S_sev_Z0,
                         S_inf_time_Z0 = S_inf_time_Z0,
                         S_inf_Z1 = S_inf_Z1,
                         S_sev_Z1 = S_sev_Z1,
                         S_inf_time_Z1 = S_inf_time_Z1)
  
  truth_data <- cbind(truth_data, X_df_Z0, X_df_Z1)
  
  # Observed only
  analysis_data <- data.frame(id = 1:n,
                              X = X,
                              Z = Z)
  
  analysis_data$S_inf <- ifelse(Z == 0, truth_data$S_inf_Z0, truth_data$S_inf_Z1)
  analysis_data$S_sev <- ifelse(Z == 0, truth_data$S_sev_Z0, truth_data$S_sev_Z1)
  analysis_data$S_inf_time <- ifelse(Z == 0, truth_data$S_inf_time_Z0, truth_data$S_inf_time_Z1)
  
  for (i in 1:12) {
    colname <- paste0("X_", i)
    analysis_data[[colname]] <- ifelse(Z == 0, 
                                       truth_data[[paste0(colname, "_Z0")]], 
                                       truth_data[[paste0(colname, "_Z1")]])
  }
  
  # Parameters
  params <- list(mean_X = bl_growth_param$mean, 
                 sd_X = bl_growth_param$sd,
                 efgh_mad_inc = mad_inc,
                 efgh_msd_inc = msd_inc,
                 efgh_enr_haz_inf_coef = effect_enr_haz_inf$inf_coef,
                 efgh_enr_haz_inf_int = effect_enr_haz_inf$inf_int,
                 efgh_enr_haz_sev_coef = effect_enr_haz_inf$sev_coef,
                 efgh_enr_haz_sev_int = effect_enr_haz_inf$sev_int,
                 final_inf_params = rr_params_inf, 
                 final_sev_params = rr_params_sev,
                 monthly_growth_model = sd_monthly)
  
  return(list(truth_data = truth_data,
              analysis_data = analysis_data, 
              params = params))

}

sim <- simulate_data()
