# ------------------------------------------------------------------------------
# Function to get truth long term growth effects
# ------------------------------------------------------------------------------

here::i_am("R/get_truth.R")

source(here::here("R/simulate_parameters.R"))
source(here::here("R/simulate_data_short_term.R"))

#' Function to get truth in naturally infected, long-term (12mo) estimand
#' 
#' @param data dataset of counterfactual data from simulate_data_short_term
#' 
#' @returns long-term effect in naturally infected
truth_nat_inf_12mo <- function(data){
  mean((data$Y_12_Z1 - data$Y_12_Z0)[data$S_inf_Z0 == 1])
}

#' Function to get truth in population, long-term (12mo)
#' 
#' @param data dataset of counterfactual data from simulate_data_short_term
#' 
#' @returns population-level effect
truth_pop_12mo <- function(data){
  mean(data$Y_12_Z1 - data$Y_12_Z0)
}


#' Function to get truth in naturally infected, short-term (V_u_months intervals)
#' 
#' @param data dataset of counterfactual data from simulate_data_short_term
#' @param V_u_months months to measure growth 
#' @param V_u_week_interval intervals in weeks corresponding to the infection times that fall into each V_u_month measurement. Length should be length(V_u_months) + 1
#' 
#' @returns population-level effect
truth_short_term <- function(data,
                             V_u_months = c(3, 6, 9, 12),
                             V_u_week_interval = c(0, 4.3, 17.3, 30.4, 52)){
  
  umax <- V_u_week_interval[length(V_u_week_interval)]
  
  E_wt_df <- data.frame(u = 1:umax,
                        wt = rep(NA, umax),
                        E_hat = rep(NA, umax))
  
  for(u in 1:umax){
    E_wt_df$wt[u] <- mean(as.numeric(data$S_inf_time_Z0[data$S_inf_Z0 == 1] == u))
    
    idx <- which(data$S_inf_Z0 == 1 & data$S_inf_time_Z0 == u)
    
    Y_V_u_name_Z0 <- paste0("Y_", V_u_months[findInterval(u, V_u_week_interval, rightmost.closed = FALSE, left.open = TRUE)], "_Z0")
    Y_V_u_name_Z1 <- paste0("Y_", V_u_months[findInterval(u, V_u_week_interval, rightmost.closed = FALSE, left.open = TRUE)], "_Z1")
    
    E_wt_df$E_hat[u] <- mean(data[[Y_V_u_name_Z1]][idx] - data[[Y_V_u_name_Z0]][idx]) 
    
  }
  
  E_wt_df$wt_x_E_hat <- E_wt_df$wt * E_wt_df$E_hat
  return(sum(E_wt_df$wt_x_E_hat))
  
}

# ------------------------------------------------------------------------------
# 6 month schedule
# ------------------------------------------------------------------------------

params <- simulate_parameters(dose_schedule = "6mo")
data <- simulate_data(parameters = params, n = 1e7, type = "counterfactual")

long_term <- truth_nat_inf_12mo(data = data)
pop <- truth_pop_12mo(data = data)
short_term <- truth_short_term(data = data)

# seed = 54321
# > long_term
# [1] 0.03826897
# > pop
# [1] 0.002413157
# > short_term
# [1] 0.02617366

# seed = 12345
# > long_term
# [1] 0.03812903
# > pop
# [1] 0.00241075
# > short_term
# [1] 0.02608022

short_term <- truth_short_term(data = data, 
                               V_u_months = c(3,6,9,12),
                               V_u_week_interval = c(0, 4.3, 17.3, 30.4, 43.5))
# > short_term
# [1] 0.02453296

short_term <- truth_short_term(data = data, 
                               V_u_months = c(6,12),
                               V_u_week_interval = c(0, 17.3, 52))
# short_term
# [1] 0.03305071

short_term <- truth_short_term(data = data, 
                               V_u_months = c(6,12),
                               V_u_week_interval = c(0, 26.5, 52))

# short_term
# [1] 0.02517426

short_term <- truth_short_term(data = data, 
                               V_u_months = c(6,12),
                               V_u_week_interval = c(0, 17.3, 43.5))

# short_term
# [1] 0.03150345

# ------------------------------------------------------------------------------
# 12 month schedule
# ------------------------------------------------------------------------------

params <- simulate_parameters(dose_schedule = "12mo", 
                              incidence_shigella_0_6 = 0.03731, 
                              incidence_shigella_6_12 = 0.05796, 
                              incidence_severe_shigella_0_6 = 0.02025,
                              incidence_severe_shigella_6_12 = 0.02548)
data <- simulate_data(parameters = params, n = 1e7, type = "counterfactual")

long_term <- truth_nat_inf_12mo(data = data)
pop <- truth_pop_12mo(data = data)
short_term <- truth_short_term(data = data)

# > long_term
# [1] 0.0359065
# > pop
# [1] 0.003459764
# > short_term
# [1] 0.02489571

short_term <- truth_short_term(data = data, 
                               V_u_months = c(3,6,9,12),
                               V_u_week_interval = c(0, 4.3, 17.3, 30.4, 43.5))

# > short_term
# [1] 0.02346045

short_term <- truth_short_term(data = data, 
                               V_u_months = c(6,12),
                               V_u_week_interval = c(0, 17.3, 52))

# > short_term
# [1] 0.03121037

short_term <- truth_short_term(data = data, 
                               V_u_months = c(6,12),
                               V_u_week_interval = c(0, 26.5, 52))

# > short_term
# [1] 0.02378575

short_term <- truth_short_term(data = data, 
                               V_u_months = c(6,12),
                               V_u_week_interval = c(0, 17.3, 43.5))

# > short_term
# [1] 0.02977512

short_term <- truth_short_term(data = data, 
                               V_u_months = c(12),
                               V_u_week_interval = c(0, 52))
# short_term
# [1] 0.0359065 == long_term yayyyy
