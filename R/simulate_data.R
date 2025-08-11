# ------------------------------------------------------------------------------
# Function to simulate data based on geometric distribution, other misc changes
# ------------------------------------------------------------------------------

here::i_am("R/simulate_data.R")

source(here::here("R/simulate_parameters.R"))

# Data structure: 

# X = covariate (baseline HAZ)
# Z = treatment (vaccine)
# S_inf = infection (Shigella)
# S_sev = severe infection (moderate to severe Shigella)
# Y_t = outcome (lfazd90)

#' Function to simulate data for short-term growth simulation
#' 
#' @param parameters parameters object from simulate_parameters
#' @param n number of observations to simulate
#' @param VE_mild vaccine efficacy for less-severe disease
#' @param VE_severe vaccine efficacy for GEMS defined MSD
#' @param seed seed to set for replicability
#' @param type type of data to return; observed, counterfactual, or both
#' 
#' @returns dataframe of simulation data with baseline HAZ (X), vaccine (Z), Shigella infection (S_inf), 
#' severe Shigella infection (S_sev), Shigella infection time (S_inf_time), monthly growth measurements Y_t
simulate_data <- function(parameters,
                          n = 1e6,
                          VE_mild = 0.40,
                          VE_severe = 0.60,
                          seed = 12345,
                          type = "counterfactual"){
  set.seed(seed)
  
  # ---------------------------------------------------------------------------
  # Z: Vaccination (RCT 1:1) --------------------------------------------------
  # ---------------------------------------------------------------------------
  Z <- rbinom(n, size=1, prob=0.5)
  
  # ---------------------------------------------------------------------------
  # X: Baseline HAZ -----------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  # Simulate baseline growth 
  X <- rnorm(n, mean = parameters$mean_X, sd = parameters$sd_X)
  
  # ---------------------------------------------------------------------------
  # S: Shigella infection (S_inf / S_sev for regular and severe; S_inf_time) --
  # ---------------------------------------------------------------------------
  
  # 4) Use the conditional hazard with intercept + coef found above to model infection / severe infection time in absence of vaccine
  
  # Shigella infection 
  
  # Shigella months 0-6 post-baseline
  hazard_inf_0_6 <- hazard(intercept = parameters$hazard_S__X_int_0_6, 
                           haz_coef = parameters$hazard_S__X_coef_0_6,
                           haz = X)
  
  # SHOULD THIS BE +1 ?? bc could be 0? then those are getting lumped in with the uninfecteds
  S_inf_time_Z0 <- rgeom(n, prob = hazard_inf_0_6) + 1
  S_inf_time_Z0 <- ifelse(S_inf_time_Z0 <= 26, S_inf_time_Z0, 0) # if infected after week 26, simulate second half infection (fill in 0 temp, add 26 later)
  
  # Shigella months 6-12 post-baseline
  hazard_inf_6_12 <- hazard(intercept = parameters$hazard_S__X_int_6_12,
                            haz_coef = parameters$hazard_S__X_coef_6_12,
                            haz = X)
  
  S_inf_time_Z0[S_inf_time_Z0 == 0] <- 26 + rgeom(length(S_inf_time_Z0[S_inf_time_Z0 == 0]), hazard_inf_6_12[S_inf_time_Z0 == 0]) + 1
  S_inf_time_Z0 <- ifelse(S_inf_time_Z0 <= 52, S_inf_time_Z0, 0) # if infected after week 52, censor at 52 (0 temp, fill in 52 later)
  
  S_inf_Z0 <- ifelse(S_inf_time_Z0 != 0, 1, 0)
  
  # Severe Shigella infection
  
  # Severe Shigella months 0-6 post-baseline
  hazard_sev_0_6 <- rep(0,n)
  hazard_sev_0_6[S_inf_Z0 == 1 & S_inf_time_Z0 <= 26] <- hazard(intercept = parameters$hazard_S_sev__X_int_0_6,
                                                                haz_coef = parameters$hazard_S_sev__X_coef_0_6,
                                                                haz = X[S_inf_Z0 == 1 & S_inf_time_Z0 <= 26])
  
  S_sev_Z0_0_6 <- rbinom(n, 1, prob = hazard_sev_0_6) 
  
  # Severe Shigella months 6-12 post-baseline
  hazard_sev_6_12 <- rep(0, n)
  hazard_sev_6_12[S_inf_Z0 == 1 & S_inf_time_Z0 > 26] <- hazard(intercept = parameters$hazard_S_sev__X_int_6_12,
                                                                 haz_coef = parameters$hazard_S_sev__X_coef_6_12,
                                                                 haz = X[S_inf_Z0 == 1 & S_inf_time_Z0 > 26])
  # save separately then combine all the 1s
  S_sev_Z0_6_12 <- rbinom(n, 1, prob = hazard_sev_6_12)
  S_sev_Z0 <- S_sev_Z0_0_6 + S_sev_Z0_6_12
  #S_sev_Z0 <- ifelse(S_sev_Z0 == 0, S_sev_Z0_6_12, S_sev_Z0)
  
  # 5) Add in vaccine
  prob_S_inf_Z1 <- rep(0,n)
  prob_S_inf_Z1[S_inf_Z0 == 1] <- 1 - VE_mild
  S_inf_Z1 <- rbinom(n, 1, prob_S_inf_Z1)
  
  S_inf_time_Z1 <- ifelse(S_inf_Z0 == 1 & S_inf_Z1 == 0, 0, S_inf_time_Z0)
  
  prob_S_sev_Z1 <- rep(0, n)
  prob_S_sev_Z1[S_sev_Z0 == 1 & S_inf_Z1 == 1] <- (1 - VE_severe) / (1 - VE_mild)
  S_sev_Z1 <- rbinom(n, 1, prob_S_sev_Z1)
  
  # ---------------------------------------------------------------------------
  # Y_t1 to Y_t12: Monthly HAZ outcome 
  # ---------------------------------------------------------------------------
  
  # Changed back to 1:12 for each setting for ease of estimation later
  matrix_range <- 1:12
  
  # if(parameters$dose_schedule == "6mo"){
  #   matrix_range <- 7:18
  # } else{
  #   matrix_range <- 13:24
  # }
  
  # noise first so same for each Z0 and Z1
  noise_matrix <- sapply(matrix_range, function(i) rnorm(n, mean = 0, sd = parameters$monthly_growth_model$sd[i]))
  
  ### Simulate growth in absence of infection ----------------------------------
  Y_vec_no_inf <- vector(mode = "list", length = 12)
  Y_vec_no_inf[[1]] <- parameters$monthly_growth_model$beta_0[matrix_range[1]] +
    parameters$monthly_growth_model$beta_1[matrix_range[1]] * X + # baseline growth
    noise_matrix[,1]
  
  for(i in 2:length(Y_vec_no_inf)){
    Y_vec_no_inf[[i]] <- parameters$monthly_growth_model$beta_0[matrix_range[i]] +
      parameters$monthly_growth_model$beta_1[matrix_range[i]] * Y_vec_no_inf[[i-1]] +  # previous month's growth
      noise_matrix[,i]
  }
  
  Y_df_no_inf <- data.frame(Y_vec_no_inf)
  colnames(Y_df_no_inf) <- paste0("Y_", matrix_range)
  
  ### Simulate growth with infection under Z0 and under Z1 ---------------------
  
  #' Helper function for simulation of growth with infection
  #' 
  #' @param S_inf_time time of infection in months
  #' @param S_inf indicator of infection
  #' @param S_sev indicator of severe infection
  #' @param n number of observations
  #' @param matrix_range depending on dose schedule, 7:12 or 13:24
  #' @param Y_df_no_inf monthly growth estimates in absence of infection
  #' @param parameters parameters object from simulate_parameters
  simulate_growth_with_infection <- function(S_inf_time,
                                             S_inf,
                                             S_sev,
                                             n,
                                             matrix_range,
                                             Y_df_no_inf,
                                             parameters){
    # Convert infection time to months
    # maybe no round
    S_inf_time_month <- S_inf_time / (52 / 12)
    
    # Get months past infection for spline model based on infection time variable
    Y_t_df <- matrix(0, nrow = n, ncol = length(matrix_range))
    
    # skip for 0 and 12 (0 = no infection, 12 = no follow up month)
    valid_idx <- which(S_inf_time_month > 0 & S_inf_time_month < 12)
    start_months <- S_inf_time_month[valid_idx]
    
    # fill in 0s up until/including infection time
    # fill rest with sequence starting at 1 go until total n 0s + other = 12
    # ex. if infection time = 4.3, ceiling(start) = 5, fill in 0 0 0 0 0.7 1.7 2.7 3.7 4.7 5.7 6.7 7.7
    filled_rows <- lapply(start_months, function(start) {
      n_mnth_before_inf <- ceiling(start) - 1
      n_mnth_after_inf <- length(matrix_range) - ceiling(start) + 1
      starting_inf_adj <- start - floor(start)
      c(rep(0, n_mnth_before_inf), seq(starting_inf_adj, starting_inf_adj + n_mnth_after_inf - 1, by = 1))
    })
    
    # Assign each row into Y_t_df
    Y_t_df[valid_idx, ] <- do.call(rbind, filled_rows)
    Y_t_df <- as.data.frame(Y_t_df)
    colnames(Y_t_df) <- paste0("Y_", matrix_range)
    
    # Infection type: 0 = none, 1 = LSD, 2 = MSD
    S_inf_type <- ifelse(S_sev == 1, 2, ifelse(S_inf == 1, 1, 0))
    
    # adjustment dataframes
    # Easiest to get preds columnwise in both cases then rebuild??
    adj_df_lsd <- matrix(nrow = nrow(Y_t_df), ncol = ncol(Y_t_df)) # adjustment given LSD at time t
    adj_df_msd <- matrix(nrow = nrow(Y_t_df), ncol = ncol(Y_t_df)) # adjustment given MSD at time t
    # make 0 if 0
    for (col in 1:ncol(Y_t_df)) {
      adj_df_lsd[, col] <- predict(parameters$effect_shigella_growth_fits$lsd, 
                                   newdata = data.frame(x = Y_t_df[, col]))
      adj_df_msd[, col] <- predict(parameters$effect_shigella_growth_fits$msd, 
                                   newdata = data.frame(x = Y_t_df[, col]))
    }
    
    # Return row of 0s or from lsd or msd matrices depending on infection type
    adj_df_S <- do.call(rbind, lapply(1:nrow(Y_t_df), function(i){
      if (S_inf_type[i] == 0) {
        rep(0, length(matrix_range))
      } else if (S_inf_type[i] == 1) {
        adj_df_lsd[i, ]
      } else {
        adj_df_msd[i, ]
      }
    }))
    
    adj_df_S <- data.frame(adj_df_S)
    colnames(adj_df_S) <- paste0("adj_", matrix_range)
    
    # Add to no-infection growth model
    Y_df <- Y_df_no_inf + adj_df_S
    
    return(Y_df)
  }
  
  # Simulate growth under Z0
  Y_df_Z0 <- simulate_growth_with_infection(
    S_inf_time = S_inf_time_Z0,
    S_inf = S_inf_Z0,
    S_sev = S_sev_Z0,
    n = n,
    matrix_range = matrix_range,
    Y_df_no_inf = Y_df_no_inf,
    parameters = parameters
  )
  
  # Simulate growth under Z1
  Y_df_Z1 <- simulate_growth_with_infection(
    S_inf_time = S_inf_time_Z1,
    S_inf = S_inf_Z1,
    S_sev = S_sev_Z1,
    n = n,
    matrix_range = matrix_range,
    Y_df_no_inf = Y_df_no_inf,
    parameters = parameters
  )
  
  # Finally convert inf_time_0 to censor at 52 (vs 0 as is up to this point)
  S_inf_time_Z0 <- ifelse(S_inf_time_Z0 == 0, 52, S_inf_time_Z0)
  S_inf_time_Z1 <- ifelse(S_inf_time_Z1 == 0, 52, S_inf_time_Z1)
  
  # ---------------------------------------------------------------------------
  # Construct final dataset, paramaters, effects 
  # ---------------------------------------------------------------------------
  
  # Counterfactual
  if(type == "counterfactual"){
    data <- data.frame(id = 1:n,
                       X = X,
                       Z = Z,
                       S_inf_Z0 = S_inf_Z0,
                       S_sev_Z0 = S_sev_Z0,
                       S_inf_time_Z0 = S_inf_time_Z0,
                       S_inf_Z1 = S_inf_Z1,
                       S_sev_Z1 = S_sev_Z1,
                       S_inf_time_Z1 = S_inf_time_Z1)
    
    colnames(Y_df_Z0) <- paste0(colnames(Y_df_Z0), "_Z0")
    colnames(Y_df_Z1) <- paste0(colnames(Y_df_Z1), "_Z1")
    
    data <- cbind(data, Y_df_Z0, Y_df_Z1)
  } else if(type == "observed"){
    # Observed only
    data <- data.frame(id = 1:n,
                       X = X,
                       Z = Z)
    
    data$S_inf <- ifelse(Z == 0, S_inf_Z0, S_inf_Z1)
    data$S_sev <- ifelse(Z == 0, S_sev_Z0, S_sev_Z1)
    data$S_inf_time <- ifelse(Z == 0, S_inf_time_Z0, S_inf_time_Z1)
    
    for (i in matrix_range) {
      colname <- paste0("Y_", i)
      data[[colname]] <- ifelse(Z == 0, 
                                Y_df_Z0[[colname]], 
                                Y_df_Z1[[colname]])
    }
  } else{
    # Both
    data <- data.frame(id = 1:n,
                       X = X,
                       Z = Z,
                       S_inf_Z0 = S_inf_Z0,
                       S_sev_Z0 = S_sev_Z0,
                       S_inf_time_Z0 = S_inf_time_Z0,
                       S_inf_Z1 = S_inf_Z1,
                       S_sev_Z1 = S_sev_Z1,
                       S_inf_time_Z1 = S_inf_time_Z1)
    
    colnames(Y_df_Z0) <- paste0(colnames(Y_df_Z0), "_Z0")
    colnames(Y_df_Z1) <- paste0(colnames(Y_df_Z1), "_Z1")
    
    data <- cbind(data, Y_df_Z0, Y_df_Z1)
    
    data$S_inf <- ifelse(Z == 0, S_inf_Z0, S_inf_Z1)
    data$S_sev <- ifelse(Z == 0, S_sev_Z0, S_sev_Z1)
    data$S_inf_time <- ifelse(Z == 0, S_inf_time_Z0, S_inf_time_Z1)
    
    for (i in matrix_range) {
      colname <- paste0("Y_", i)
      data[[colname]] <- ifelse(Z == 0, 
                                Y_df_Z0[[paste0("Y_", i,"_Z0")]], 
                                Y_df_Z1[[paste0("Y_", i,"_Z1")]])
    }
    
  }
  
  return(data)
  
}
