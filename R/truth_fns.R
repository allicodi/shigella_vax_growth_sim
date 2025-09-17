# ------------------------------------------------------------------------------
# Function to get truth long term growth effects
# ------------------------------------------------------------------------------

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
                        E_hat = rep(NA, umax),
                        E_hat_Z0 = rep(NA, umax),
                        E_hat_Z1 = rep(NA, umax))
  
  for(u in 1:umax){
    E_wt_df$wt[u] <- mean(as.numeric(data$S_inf_time_Z0[data$S_inf_Z0 == 1] == u))
    
    idx <- which(data$S_inf_Z0 == 1 & data$S_inf_time_Z0 == u)
    
    Y_V_u_name_Z0 <- paste0("Y_", V_u_months[findInterval(u, V_u_week_interval, rightmost.closed = FALSE, left.open = TRUE)], "_Z0")
    Y_V_u_name_Z1 <- paste0("Y_", V_u_months[findInterval(u, V_u_week_interval, rightmost.closed = FALSE, left.open = TRUE)], "_Z1")
    
    E_wt_df$E_hat_Z1[u] <- mean(data[[Y_V_u_name_Z1]][idx])
    E_wt_df$E_hat_Z0[u] <- mean(data[[Y_V_u_name_Z0]][idx])
    E_wt_df$E_hat[u] <- mean(data[[Y_V_u_name_Z1]][idx] - data[[Y_V_u_name_Z0]][idx]) 
    
  }
  
  E_wt_df$wt_x_E_hat <- E_wt_df$wt * E_wt_df$E_hat
  E_wt_df$wt_x_E_hat_Z0 <- E_wt_df$wt * E_wt_df$E_hat_Z0
  E_wt_df$wt_x_E_hat_Z1 <- E_wt_df$wt * E_wt_df$E_hat_Z1
  
  return(list(short_term = sum(E_wt_df$wt_x_E_hat),
                    short_term_Z0 = sum(E_wt_df$wt_x_E_hat_Z0),
                    short_term_Z1 = sum(E_wt_df$wt_x_E_hat_Z1)))
  
}

