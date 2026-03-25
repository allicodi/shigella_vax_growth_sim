# ----------------------------------------------------
# Simple power calculation for primary endpoint
# ----------------------------------------------------

library(gsDesign)

here::i_am("misc/primary_endpoint_power.R")

# 90% power to detect a vaccine efficacy of 60% with 10% dropout, 
# 1:1 group allocation, a two-sided test, and a null hypothesis of 20%

all_settings <- c("default",
                  "base_12mo",
                  "optimistic_6mo",
                  "optimistic_12mo",
                  "peru_6mo",
                  "peru_12mo")

sample_size_df <- data.frame(setting = all_settings,
                             inc_0_6 = rep(NA, length(all_settings)),
                             inc_6_12 = rep(NA, length(all_settings)),
                             inc_0_12 = rep(NA, length(all_settings)),
                             
                             inc_sev_0_6 = rep(NA, length(all_settings)),
                             inc_sev_6_12 = rep(NA, length(all_settings)),
                             inc_sev_0_12 = rep(NA, length(all_settings)),
                             n = rep(NA, length(all_settings)))

for(setting in all_settings){
  parameters <- readRDS(here::here(paste0("parameters/parameters_", setting, ".Rds")))
  
  # always assume 60% efficacy against severe disease, lower bound has to be >= 20%
  VE_sev <- 0.60
  VE_null <- 0.20
  
  inc_0_6 <- parameters$incidence_S_mean_X_0_6   # 6 month incidence of MAD in first half of trial (to pull from df later, not needed for primary power)
  inc_6_12 <- parameters$incidence_S_mean_X_6_12 # 6 month incidence of MAD in second half of trial
  
  sev_inc_0_6 <- parameters$incidence_S_sev_mean_X_0_6   # 6 month incidence of severe disease in first half of trial
  sev_inc_6_12 <- parameters$incidence_S_sev_mean_X_6_12 # 6 month incidence of severe disease in second half of trial
  
  p1 <- sev_inc_0_6 + sev_inc_6_12 
  
  n_total <- nBinomial(
    p1 = p1,
    p2 = p1 * (1 - VE_sev),
    scale = "RR",
    delta0 = log(p1/(p1 * (1 - VE_null))),
    alpha = 0.05,
    beta = 0.10,
    sided = 2     
  )
  
  n_total <- ceiling(n_total / (1 - 0.1))  # adjust for 10% dropout
  
  sample_size_df$inc_0_6[sample_size_df$setting == setting] <- inc_0_6
  sample_size_df$inc_6_12[sample_size_df$setting == setting] <- inc_6_12
  sample_size_df$inc_0_12[sample_size_df$setting == setting] <- inc_0_6 + inc_6_12
  
  sample_size_df$inc_sev_0_6[sample_size_df$setting == setting] <- sev_inc_0_6
  sample_size_df$inc_sev_6_12[sample_size_df$setting == setting] <- sev_inc_6_12
  sample_size_df$inc_sev_0_12[sample_size_df$setting == setting] <- p1
  
  sample_size_df$n[sample_size_df$setting == setting] <- n_total
}

saveRDS(sample_size_df, here::here("results/primary_endpoint_power.Rds"))

#######################################

library(gsDesign)

ref <- c(0.01292, 0.03116, 0.00836, 0.01976, 0.01596,
         0.015844, 0.038212, 0.010252, 0.024232, 0.019572,
         0.01666, 0.04018, 0.01078, 0.02548, 0.02058)

n_total <- nBinomial(
  p1 = ref,
  p2 = ref * (1 - 0.6),
  scale = "RR",
  delta0 = log(ref/(ref * (1 - 0.2))),
  alpha = 0.05,
  beta = 0.10,
  sided = 2     
)

n_total <- ceiling(n_total * 1.1) # matches numbers in table yay

# n_total <- ceiling(n_total / (1 - 0.1))  # adjust for 10% dropout
