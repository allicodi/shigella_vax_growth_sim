# ------------------------------------------------------------------------------
# Script to find intercept & haz coefficient for desired parameter combinations
# ------------------------------------------------------------------------------

library(tidyverse)

here::i_am("R/tune_parameters.R")

# Functions for finding coefficients  based on incidence, hazard ---------------

cum_inc <- function(hazard, t0){
  sum(hazard * (1 - hazard)^((1:t0) - 1))
}

#expit
hazard <- function(intercept, haz_coef, haz){
  plogis(intercept + haz_coef * haz)
}

cum_inc_by_row <- function(row){

  hazard_haz_mean <- hazard(intercept = row[['intercept']],
                            haz_coef = row[['haz_coef']],
                            haz = row[['mean_X']])
  
  cum_inc_haz_mean <- cum_inc(hazard = hazard_haz_mean, t0 = 26)
   
  return(cum_inc_haz_mean)
  
}

# ------------------------------------------------------------------------------

# Main set of parameters from 'default settings' on poster

# Cumulative incidence of (all) Shigella - 0.0618 / year
# Cumulative incidence of moderate to severe Shigella - 0.0246 / year
# NOTE: these should be updated based on final EFGH data

# effect of baseline haz on (all) Shigella = haz coefficient?? should reflect Liz risk ratio 0.73
# NOTE Liz would rather get this from data as we did for MSD below

# updated from efgh data: -0.064

# effect of baseline haz on moderate to severe Shigella = -0.199
# NOTE this should be updated based on final EFGH data

# updated from efgh data: -0.222

# ------------------------------------------------------------------------------
# 1) All Shigella
# ------------------------------------------------------------------------------
# 
# # Incidence = 0.068 / year Shigella, Risk Ratio = 0.73
# 
# # QUESTION - cumulative incidence should = 0.06 at what value of haz?? mean?? 
# 
# # Grid to search over to best match risk ratio / cumulative incidence
# intercept <- seq(-8, 5, by = 0.05)
# haz_coef <- log(seq(0.01, 1, by = 0.01))
# mean_X <- -0.9395
# 
# combos <- expand.grid(intercept = intercept, 
#                       haz_coef = haz_coef,
#                       mean_X = mean_X)
# 
# # Get hazard for haz = -0.5 & haz = 0.5, hazard ratio for each intercept / haz coefficient combination 
# hazard_ratio_combos <- cbind(combos, do.call(rbind, apply(combos, 1, hazard_ratio)))
# 
# # Use hazard results to iterate over again for risk ratio
# risk_ratio_combos <- cbind(hazard_ratio_combos, do.call(rbind, apply(hazard_ratio_combos, 1, risk_ratio)))
# 
# # Find risk ratio close to 0.73 & cumulative incidence at mean close to 0.0618
# subset_rr <- risk_ratio_combos %>%
#   filter(risk_ratio > 0.72 & risk_ratio < 0.74) %>%
#   filter(cum_inc_haz_mean > 0.055 & cum_inc_haz_mean < 0.065)
# 
# # Manually choose closest: Intercept = -5.55, haz coef = -0.3285041 (log(0.72))
# 
# # intercept   haz_coef  mean_X hazard_haz_mean hazard_haz_minus0.5 hazard_haz_0.5 hazard_ratio cum_inc_haz_mean cum_inc_haz_minus0.5 cum_inc_haz_0.5 risk_ratio
# #  -5.55 -0.3285041    -0.9395      0.00526513         0.006199189    0.004471177    0.7212519       0.06138367           0.07190557      0.05235415  0.7280959
# 
# # ------------------------------------------------------------------------------
# # 2) Moderate-to-severe Shigella
# # ------------------------------------------------------------------------------
# 
# # Incidence = 0.0246 / year moderate to severe Shigella, Risk Ratio = ???
# # define Y_sev using GEMS MSD definition, fit severe diarrhea ~ enrollment HAZ + site only in kids age 9-12 months (the age they'd be when they finish getting vaccinated)
# 
# # From shigella_ve/efgh_data/explore_real_data.R
# # Line 206-207
# # summary(fit.severity.9.12)
# # 
# # Call:
# #   glm(formula = gems_msd ~ enr_haz + enroll_site, family = "binomial", 
# #       data = severe_haz_df_age_9_12)
# # 
# # Coefficients:
# #   Estimate Std. Error z value Pr(>|z|)    
# # (Intercept)   -0.1753     0.2823  -0.621   0.5347    
# # enr_haz       -0.1990     0.1114  -1.787   0.0740 .  
# # enroll_site2   0.3144     0.4854   0.648   0.5172    
# # enroll_site3  -1.3068     0.7090  -1.843   0.0653 .  
# # enroll_site4  -1.2595     0.5008  -2.515   0.0119 *  
# # enroll_site5  -0.4605     0.3488  -1.320   0.1868    
# # enroll_site6   2.0085     0.5078   3.955 7.65e-05 ***
# # enroll_site7  -0.7847     0.4246  -1.848   0.0646 . 
# 
# # Grid to search over to best match risk ratio / cumulative incidence
# intercept <- seq(-8, 2, by = 0.01)       # should be ~ -0.175 ??? except cumulative incidence high around that so drop intercept?
# haz_coef <- seq(-0.25, -0.15, by = 0.01) # should be ~ -0.199
# mean_X <- -0.9395
# 
# combos <- expand.grid(intercept = intercept, 
#                       haz_coef = haz_coef,
#                       mean_X = mean_X)
# 
# # Get hazard for haz = -0.5 & haz = 0.5, hazard ratio for each intercept / haz coefficient combination 
# hazard_ratio_combos <- cbind(combos, do.call(rbind, apply(combos, 1, hazard_ratio)))
# 
# # Use hazard results to iterate over again for risk ratio
# risk_ratio_combos <- cbind(hazard_ratio_combos, do.call(rbind, apply(hazard_ratio_combos, 1, risk_ratio)))
# 
# # Find risk ratio close to 0.73 & cumulative incidence at mean close to 0.0618
# subset_rr <- risk_ratio_combos %>%
#   filter(cum_inc_haz_mean > 0.022 & cum_inc_haz_mean < 0.026)

# Manually choose closest: Intercept = -6.36, haz coef = -0.20

# intercept haz_coef  mean_X hazard_haz_mean hazard_haz_minus0.5 hazard_haz_0.5 hazard_ratio cum_inc_haz_mean cum_inc_haz_minus0.5 cum_inc_haz_0.5 risk_ratio
#   -6.36     -0.2 -0.9395     0.002082503         0.002301018      0.0018847    0.8190724       0.02470579           0.02726543      0.02238343  0.8209454
