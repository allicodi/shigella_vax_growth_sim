# ------------------------------------------------------------------------------
# Function to get truth long term growth effects
# ------------------------------------------------------------------------------

here::i_am("R/estimation_long_term.R")

source(here::here("R/simulate_parameters.R"))

# type = counterfactual, observed, both
params <- simulate_parameters(dose_schedule = "6mo")
data <- simulate_data(parameters = params, n = 1e7, type = "counterfactual")

true_ge <- mean((data$Y_12_Z1 - data$Y_12_Z0)[data$S_inf_Z0 == 1])
true_pop_ge <- mean(data$Y_12_Z1 - data$Y_12_Z0)

# > true_ge
# [1] 0.03812903
# 
# > true_pop_ge
# [1] 0.00241075

params <- simulate_parameters(dose_schedule = "12mo", 
                              incidence_shigella_0_6 = 0.03731, 
                              incidence_shigella_6_12 = 0.05796, 
                              incidence_severe_shigella_0_6 = 0.02025,
                              incidence_severe_shigella_6_12 = 0.02548)
data <- simulate_data(parameters = params, n = 1e7, type = "counterfactual")

true_ge <- mean((data$Y_12_Z1 - data$Y_12_Z0)[data$S_inf_Z0 == 1])
true_pop_ge <- mean(data$Y_12_Z1 - data$Y_12_Z0)

# > true_ge
# [1] 0.0359065

# > true_pop_ge
# [1] 0.003459764
