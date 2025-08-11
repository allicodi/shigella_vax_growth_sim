# ------------------------------------------------------------------------------
# Function to get truth long term growth effects
# ------------------------------------------------------------------------------

here::i_am("R/get_truth.R")

source(here::here("R/simulate_parameters.R"))
source(here::here("R/simulate_data_short_term.R"))

# type = counterfactual, observed, both
params <- simulate_parameters(dose_schedule = "6mo")
data <- simulate_data(seed = 54321, parameters = params, n = 1e7, type = "counterfactual")

# Long-term growth effect (12 month growth)
true_ge <- mean((data$Y_12_Z1 - data$Y_12_Z0)[data$S_inf_Z0 == 1])

# Population growth effect (12 month growth; naive way)
true_pop_ge <- mean(data$Y_12_Z1 - data$Y_12_Z0)

# Short-term growth effect 
umax <- 52
E_wt_df <- data.frame(u = 1:umax,
                      wt = rep(NA, umax),
                      E_hat = rep(NA, umax))

for(u in 1:umax){
  E_wt_df$wt[u] <- mean(as.numeric(data$S_inf_time_Z0[data$S_inf_Z0 == 1] == u))
  
  idx <- which(data$S_inf_Z0 == 1 & data$S_inf_time_Z0 == u)

  if(u <= 4.3){
    E_wt_df$E_hat[u] <- mean(data$Y_3_Z1[idx] - data$Y_3_Z0[idx])
  } else if(u > 4.3 & u <= 17.3){
    E_wt_df$E_hat[u] <- mean(data$Y_6_Z1[idx] - data$Y_6_Z0[idx])
  } else if(u > 17.3 & u <= 30.4){
    E_wt_df$E_hat[u] <- mean(data$Y_9_Z1[idx] - data$Y_9_Z0[idx])
  } else if(u > 30.4 & u <= umax){
    E_wt_df$E_hat[u] <- mean(data$Y_12_Z1[idx] - data$Y_12_Z0[idx])
  } 
}

E_wt_df$wt_x_E_hat <- E_wt_df$wt * E_wt_df$E_hat
true_short_term_ge <- sum(E_wt_df$wt_x_E_hat)

# true_ge = 0.038
# true pop = 0.0024
# short = 0.026

# --------------------------

# type = counterfactual, observed, both
params <- simulate_parameters(dose_schedule = "6mo", 
                              effect_shigella_growth_formula = y ~ -1 + x + I(pmax(0, x - 3)) + I(pmax(0, x - 6)) + I(pmax(0, x - 9)))
data <- simulate_data(parameters = params, n = 1e7, type = "counterfactual")

# Long-term growth effect (12 month growth)
true_ge <- mean((data$Y_12_Z1 - data$Y_12_Z0)[data$S_inf_Z0 == 1])

# Population growth effect (12 month growth; naive way)
true_pop_ge <- mean(data$Y_12_Z1 - data$Y_12_Z0)

# Short-term growth effect 
umax <- 52
E_wt_df <- data.frame(u = 1:umax,
                      wt = rep(NA, umax),
                      E_hat = rep(NA, umax))

for(u in 1:umax){
  E_wt_df$wt[u] <- mean(as.numeric(data$S_inf_time_Z0[data$S_inf_Z0 == 1] == u))
  
  idx <- which(data$S_inf_Z0 == 1 & data$S_inf_time_Z0 == u)
  
  if(u <= 4.3){
    E_wt_df$E_hat[u] <- mean(data$Y_3_Z1[idx] - data$Y_3_Z0[idx])
  } else if(u > 4.3 & u <= 17.3){
    E_wt_df$E_hat[u] <- mean(data$Y_6_Z1[idx] - data$Y_6_Z0[idx])
  } else if(u > 17.3 & u <= 30.4){
    E_wt_df$E_hat[u] <- mean(data$Y_9_Z1[idx] - data$Y_9_Z0[idx])
  } else if(u > 30.4 & u <= umax){
    E_wt_df$E_hat[u] <- mean(data$Y_12_Z1[idx] - data$Y_12_Z0[idx])
  } 
}

E_wt_df$wt_x_E_hat <- E_wt_df$wt * E_wt_df$E_hat
true_short_term_ge <- sum(E_wt_df$wt_x_E_hat)

# true ge = 0.039
# true pop ge = 0.00248
# true short term ge = 0.031

# --------------------------

# type = counterfactual, observed, both
params <- simulate_parameters(dose_schedule = "12mo", 
                              incidence_shigella_0_6 = 0.03731, 
                              incidence_shigella_6_12 = 0.05796, 
                              incidence_severe_shigella_0_6 = 0.02025,
                              incidence_severe_shigella_6_12 = 0.02548)
data <- simulate_data(parameters = params, n = 1e7, type = "counterfactual")

# Long-term growth effect (12 month growth)
true_ge <- mean((data$Y_12_Z1 - data$Y_12_Z0)[data$S_inf_Z0 == 1])

# Population growth effect (12 month growth; naive way)
true_pop_ge <- mean(data$Y_12_Z1 - data$Y_12_Z0)

# Short-term growth effect 
umax <- 52
E_wt_df <- data.frame(u = 1:umax,
                      wt = rep(NA, umax),
                      E_hat = rep(NA, umax))

for(u in 1:umax){
  E_wt_df$wt[u] <- mean(as.numeric(data$S_inf_time_Z0[data$S_inf_Z0 == 1] == u))
  
  idx <- which(data$S_inf_Z0 == 1 & data$S_inf_time_Z0 == u)
  
  if(u <= 4.3){
    E_wt_df$E_hat[u] <- mean(data$Y_3_Z1[idx] - data$Y_3_Z0[idx])
  } else if(u > 4.3 & u <= 17.3){
    E_wt_df$E_hat[u] <- mean(data$Y_6_Z1[idx] - data$Y_6_Z0[idx])
  } else if(u > 17.3 & u <= 30.4){
    E_wt_df$E_hat[u] <- mean(data$Y_9_Z1[idx] - data$Y_9_Z0[idx])
  } else if(u > 30.4 & u <= umax){
    E_wt_df$E_hat[u] <- mean(data$Y_12_Z1[idx] - data$Y_12_Z0[idx])
  } 
}

E_wt_df$wt_x_E_hat <- E_wt_df$wt * E_wt_df$E_hat
true_short_term_ge <- sum(E_wt_df$wt_x_E_hat)

# true ge = 0.036
# true pop ge = 0.0035
# true short term = 0.025

# type = counterfactual, observed, both
params <- simulate_parameters(dose_schedule = "12mo", 
                              incidence_shigella_0_6 = 0.03731, 
                              incidence_shigella_6_12 = 0.05796, 
                              incidence_severe_shigella_0_6 = 0.02025,
                              incidence_severe_shigella_6_12 = 0.02548,
                              effect_shigella_growth_formula = y ~ -1 + x + I(pmax(0, x - 3)) + I(pmax(0, x - 6)) + I(pmax(0, x - 9)))
data <- simulate_data(parameters = params, n = 1e7, type = "counterfactual")

# Long-term growth effect (12 month growth)
true_ge <- mean((data$Y_12_Z1 - data$Y_12_Z0)[data$S_inf_Z0 == 1])

# Population growth effect (12 month growth; naive way)
true_pop_ge <- mean(data$Y_12_Z1 - data$Y_12_Z0)

# Short-term growth effect 
umax <- 52
E_wt_df <- data.frame(u = 1:umax,
                      wt = rep(NA, umax),
                      E_hat = rep(NA, umax))

for(u in 1:umax){
  E_wt_df$wt[u] <- mean(as.numeric(data$S_inf_time_Z0[data$S_inf_Z0 == 1] == u))
  
  idx <- which(data$S_inf_Z0 == 1 & data$S_inf_time_Z0 == u)
  
  if(u <= 4.3){
    E_wt_df$E_hat[u] <- mean(data$Y_3_Z1[idx] - data$Y_3_Z0[idx])
  } else if(u > 4.3 & u <= 17.3){
    E_wt_df$E_hat[u] <- mean(data$Y_6_Z1[idx] - data$Y_6_Z0[idx])
  } else if(u > 17.3 & u <= 30.4){
    E_wt_df$E_hat[u] <- mean(data$Y_9_Z1[idx] - data$Y_9_Z0[idx])
  } else if(u > 30.4 & u <= umax){
    E_wt_df$E_hat[u] <- mean(data$Y_12_Z1[idx] - data$Y_12_Z0[idx])
  } 
}

E_wt_df$wt_x_E_hat <- E_wt_df$wt * E_wt_df$E_hat
true_short_term_ge <- sum(E_wt_df$wt_x_E_hat)

# true ge = 0.037
# true pop ge = 0.0036
# true short term = 0.029
