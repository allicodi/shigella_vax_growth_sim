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

# ------------------------------------------------------------------------------

### Plot of infection times 
hist(data$S_inf_time_Z0[data$S_inf_time_Z0 != 52], breaks = 52, 
     main = "Histogram of Shigella infection time in weeks", 
     xlab = "Shigella infection time (weeks)")

hist(data$S_inf_time_Z1[data$S_inf_time_Z1 != 52], breaks = 52, col = 'red', add = TRUE)

data$S_inf_time_Z0_month <- data$S_inf_time_Z0 / (52 / 12)
data$S_inf_time_Z1_month <- data$S_inf_time_Z1 / (52 / 12)

hist(data$S_inf_time_Z0_month[data$S_inf_time_Z0 != 52], breaks = 12, 
     main = "Histogram of Shigella infection time in months", 
     xlab = "Shigella infection time (months)")

hist(data$S_inf_time_Z1_month[data$S_inf_time_Z1 != 52], breaks = 12, col = 'red', add = TRUE)
# ^^why does it look like more on the 3 month marks??

### Month of min growth measurement

# dataframe of infected
matrix_range <- 1:12

data_S_inf_Z0 <- data[data$S_inf_Z0 == 1,]
data_S_inf_Z1 <- data[data$S_inf_Z1 == 1,]

# Get months past infection for spline model based on infection time variable
Y_t_df_Z0 <- matrix(0, nrow = nrow(data_S_inf_Z0), ncol = 12)
Y_t_df_Z1 <- matrix(0, nrow = nrow(data_S_inf_Z1), ncol = 12)

# skip for 0 and 12 (0 = no infection, 12 = no follow up month)
start_months_Z0 <- data_S_inf_Z0$S_inf_time_Z0_month
start_months_Z1 <- data_S_inf_Z1$S_inf_time_Z1_month

# fill in 0s up until/including infection time
# fill rest with sequence starting at 1 go until total n 0s + other = 12
# ex. if infection time = 4.3, ceiling(start) = 5, fill in 0 0 0 0 0.7 1.7 2.7 3.7 4.7 5.7 6.7 7.7
filled_rows_Z0 <- lapply(start_months_Z0, function(start) {
  n_mnth_before_inf <- ceiling(start) - 1
  n_mnth_after_inf <- length(matrix_range) - ceiling(start) + 1
  starting_inf_adj <- start - floor(start)
  c(rep(0, n_mnth_before_inf), seq(starting_inf_adj, starting_inf_adj + n_mnth_after_inf - 1, by = 1))
})

filled_rows_Z1 <- lapply(start_months_Z1, function(start) {
  n_mnth_before_inf <- ceiling(start) - 1
  n_mnth_after_inf <- length(matrix_range) - ceiling(start) + 1
  starting_inf_adj <- start - floor(start)
  c(rep(0, n_mnth_before_inf), seq(starting_inf_adj, starting_inf_adj + n_mnth_after_inf - 1, by = 1))
})

# Assign each row into Y_t_df
Y_t_df_Z0 <- do.call(rbind, filled_rows_Z0)
Y_t_df_Z0 <- as.data.frame(Y_t_df_Z0)
colnames(Y_t_df_Z0) <- paste0("Y_", matrix_range, "_Z0")

Y_t_df_Z1 <- do.call(rbind, filled_rows_Z1)
Y_t_df_Z1 <- as.data.frame(Y_t_df_Z1)
colnames(Y_t_df_Z1) <- paste0("Y_", matrix_range, "_Z1")

# Get column with minimum value
cols_Z0 <- paste0("Y_", 1:12, "_Z0")  # column names from Y_1_Z0 to Y_12_Z0
min_col_Z0 <- apply(data_S_inf_Z0[cols_Z0], 1, which.min)
Y_t_df_Z0$min_col_Z0 <- min_col_Z0

Y_t_df_Z0$min_growth_time <- Y_t_df_Z0[cbind(1:nrow(Y_t_df_Z0), Y_t_df_Z0$min_col_Z0)]
Y_t_df_Z0$I_sev <- data_S_inf_Z0$S_sev_Z0

cols_Z1 <- paste0("Y_", 1:12, "_Z1")  # column names from Y_1_Z0 to Y_12_Z0
min_col_Z1 <- apply(data_S_inf_Z1[cols_Z1], 1, which.min)
Y_t_df_Z1$min_col_Z1 <- min_col_Z1

Y_t_df_Z1$min_growth_time <- Y_t_df_Z1[cbind(1:nrow(Y_t_df_Z1), Y_t_df_Z1$min_col_Z1)]
Y_t_df_Z1$I_sev <- data_S_inf_Z1$S_sev_Z1

# min growth time often before infection?? 

# do again only considering post infection

# Get column with minimum value

# Replace 0 with NA for the which.min search
cols_Z0 <- paste0("Y_", 1:12, "_Z0")
cols_Z1 <- paste0("Y_", 1:12, "_Z1")

# Step 1: logical mask where TRUE = post infection
post_inf_mask_Z0 <- Y_t_df_Z0[cols_Z0] != 0
post_inf_mask_Z1 <- Y_t_df_Z1[cols_Z1] != 0

# Step 2: copy data and set pre-infection values to NA
data_S_inf_Z0_mask <- data_S_inf_Z0[cols_Z0]; data_S_inf_Z0_mask[!post_inf_mask_Z0] <- NA
data_S_inf_Z1_mask <- data_S_inf_Z1[cols_Z1]; data_S_inf_Z1_mask[!post_inf_mask_Z1] <- NA

min_col_Z0_post_inf <- apply(data_S_inf_Z0_mask[cols_Z0], 1, which.min)
Y_t_df_Z0$min_col_Z0_post_inf <- as.numeric(min_col_Z0_post_inf)

Y_t_df_Z0$min_growth_time_post_inf <- Y_t_df_Z0[cbind(1:nrow(Y_t_df_Z0), Y_t_df_Z0$min_col_Z0_post_inf)]

summary(Y_t_df_Z0$min_growth_time_post_inf[Y_t_df_Z0$I_sev == 0])
summary(Y_t_df_Z0$min_growth_time_post_inf[Y_t_df_Z0$I_sev == 1])

# ^ but also this is not quite what we want either? bc this isn't effect this is just outcome? 