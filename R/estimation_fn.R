#' Function for estimation of short term growth effect
#' 
#' @param data dataset of observed data from simulate_data_short_term 
#' @param parameters parameters from simulate_parameters
#' @param V_u_months months to measure growth 
#' @param V_u_week_interval intervals in weeks corresponding to the infection times that fall into each V_u_month measurement. Length should be length(V_u_months) + 1
#' @param pool boolean indicating pool over visits or not, default FALSE
#' @param I_S__X_Z0_0_6_formula formula for shigella infection in unvaccinated first six months of trial
#' @param I_S__X_Z0_6_12_formula formula for shigella infection in unvaccinated last six months of trial
#' @param Y_out__Z1_Suminus10_X_formula formula for growth outcome in vaccinated
#' @param Y_out__Z0_Suminus11_X_T_formula formula for growth outcome in vaccinated, not pooled 
#' @param Y_out__Z0_Suminus11_X_T_V_formula formula for growth outcome in vaccinated, pooled 
estimate_short_term <- function(data,
                                parameters,
                                V_u_months = c(3, 6, 9, 12),
                                V_u_week_interval = c(0, 4.3, 17.3, 30.4, 52),
                                pool = FALSE,
                                I_S__X_Z0_0_6_formula = "I_shig_all ~ X",
                                I_S__X_Z0_6_12_formula = "I_shig_all ~ X",
                                Y_out__Z1_Suminus10_X_formula = "Y_out ~ X",
                                Y_out__Z0_Suminus11_X_T_formula = "Y_out ~ X + Ti",
                                Y_out__Z0_Suminus11_X_T_V_formula = "Y_out ~ X + Ti + V"){
  
  # umax = maximum week that can be included (last obs in V_u_week_max)
  umax <- V_u_week_interval[length(V_u_week_interval)]

  # -------------------------------------------------------
  # Dataset formatting
  # -------------------------------------------------------
  
  # Add column shig where = 2 severe, = 1 mild, = 0 no
  data <- data %>%
    mutate(shig_all = if_else(S_sev == 1, 2, 
                          if_else(S_inf == 1, 1, 0)))
  # Find which month each S_inf_time falls into
  month_match <- rep(NA, nrow(data))
  
  # If not infected, final_month
  month_match[which(data$S_inf == 0)] <- V_u_months[length(V_u_months)]
  
  # If infected, find which interval fall into- if does not fall into interval, remains NA
  month_match[which(data$S_inf == 1)] <- V_u_months[
    findInterval(data$S_inf_time[which(data$S_inf == 1)], V_u_week_interval, rightmost.closed = FALSE, left.open = TRUE) 
  ]
  
  # Make column names
  Y_col_name <- paste0("Y_", month_match)
  
  # Pull results into Y_out column
  data$Y_out <- data[cbind(seq_len(nrow(data)), match(Y_col_name, names(data)))]
  
  # make long dataset
  long_data <- expand.grid(id = 1:nrow(data), week = 1:52) %>%
    arrange(id) %>%
    left_join(data, by = "id") %>%
    filter(!(shig_all != 0 & week > S_inf_time)) %>%
    mutate(shig_all = if_else(S_inf_time != 0 & week != S_inf_time, 0, shig_all)) %>%
    mutate(I_shig_all = if_else(shig_all != 0, 1, 0)) %>%
    select(-S_inf_time)
  
  # two separate models for first and second half of trial
  long_data_0_6 <- long_data[long_data$week <= 26,]
  long_data_6_12 <- long_data[long_data$week >= 27,] # remove people who were infected in the first half (they won't have any measurements 27 or later)
  
  # -------------------------------------------------------
  # Get densities
  # -------------------------------------------------------
  
  fit_0_6_Z0 <- glm(as.formula(I_S__X_Z0_0_6_formula), 
                    data = long_data_0_6[long_data_0_6$Z == 0, ], 
                    family = binomial())
  
  fit_6_12_Z0 <- glm(as.formula(I_S__X_Z0_6_12_formula), 
                     data = long_data_6_12[long_data_6_12$Z == 0, ], 
                     family = binomial())
  
  haz_Z0_0_6 <- predict(fit_0_6_Z0, newdata = data, type = 'response')
  density_Z0_0_6 <- do.call(rbind, sapply(haz_Z0_0_6, function(haz){
    haz * (1 - haz)^((1:26) - 1)
  }, simplify = FALSE))
  
  haz_Z0_6_12 <- predict(fit_6_12_Z0, newdata = data, type = 'response')
  cuminc_Z0_0_6 <- rowSums(density_Z0_0_6)
  
  density_Z0_6_12 <- do.call(rbind, 
                             mapply(
                               haz_6_12 = haz_Z0_6_12, cum_inc_0_6 = cuminc_Z0_0_6, 
                               FUN = function(haz_6_12, cum_inc_0_6){
                                 haz_6_12 * (1 - haz_6_12)^((1:26) - 1) * (1 - cum_inc_0_6)
                               }, SIMPLIFY = FALSE
                             )
  )
  
  density_Z0_0_12 <- cbind(density_Z0_0_6, density_Z0_6_12)
  
  # ----------------------------------------------------------------------------
  # Growth models
  # ----------------------------------------------------------------------------
  
  # E[Y_Vu | Z = 1, S_u-1 = 0, X] ----------------------------------------------
  
  # 1. subset to vaccinated
  data_Z1 <- data[data$Z == 1,]
  
  pred_Y_out__Z1_Suminus1_0_X <- data.frame(setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = umax)),
                                                     paste0("pred_", 1:umax)))
  
  for(u in 1:floor(umax)){
    
    # 2. determine which growth measurement we are considering for outcome of the regression 
    #.   based on the value of u. ex. if u = 1 we select the three month growth outcome
    Y_V_u_name <- paste0("Y_", V_u_months[findInterval(u, V_u_week_interval, rightmost.closed = FALSE, left.open = TRUE)] )
    
    # 3. Further subset the data to individuals who have not yet become ill by u (i.e. 
    #.   remain illness free through u - 1). This subset of data will potentially include individuals
    #.   who become infected either at u or between date and V(u)th visit. That's ok
    data_Z1_u <- data_Z1[which(data_Z1$S_inf_time >= u), ]
    data_Z1_u$Y_out <- data_Z1_u[[Y_V_u_name]]
    
    # 4. Regress selected growth outcome on X in this subset of data. This gives your estimate of E[Y_Vu | Z = 1, S_u-1 = 0, X]
    # TODO could be smart and look to see if ppl infected in between
    fit_Y_out__Z1_Suminus1_0_X <- glm(as.formula(Y_out__Z1_Suminus10_X_formula), 
                                      data = data_Z1_u,
                                      family = gaussian())
    
    # predict on full data
    pred_Y_out__Z1_Suminus1_0_X[,u] <- predict(fit_Y_out__Z1_Suminus1_0_X, newdata = data)
  }
  
  # E[Y_Vu | Z = 0, S_u-1 = 1, X] ----------------------------------------------
  
  ### If not dense in u (version 1): ###
  
  if(!pool){
    # 1. subset to unvaccinated
    data_Z0 <- data[data$Z == 0,]
    
    pred_Y_out__Z0_Suminus1_1_X <- data.frame(setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = floor(umax))),
                                                       paste0("pred_", 1:floor(umax))))
    
    # # 2. further subset to all individuals who fall ill on any u in U-1(v)
    for(i in 1:length(V_u_months)){
      
      # this should be min, max ranges for each visit
      v_umin <- V_u_week_interval[i]
      v_umax <- V_u_week_interval[i+1]
      
      # original but that doesn't require ill? (if at censoring time will be included still so wrong i think?)
      # data_Z0_Uminus1_v <- data_Z0[which(data_Z0$S_inf_time > v_umin & data_Z0$S_inf_time <= v_umax),] # bigger than min and less than max
      
      data_Z0_Uminus1_v <- data_Z0[which(data_Z0$S_inf == 1 & data_Z0$S_inf_time > v_umin & data_Z0$S_inf_time <= v_umax),] # bigger than min and less than max
      
      # 3. create a variable for each individual, say Ti, indicating the time period at which each individual
      #    falls ill, ie set Ti equal to the u such that dSu,i = 1
      
      # same as S_inf_time, rename for clarity notation
      data_Z0_Uminus1_v$Ti <- data_Z0_Uminus1_v$S_inf_time
      
      # 4. regress the selected growth outcome on X and T in this subset of data
      fit_Y_out__Z0_Suminus1_1_X_T <- glm(as.formula(Y_out__Z0_Suminus11_X_T_formula),
                                          data = data_Z0_Uminus1_v,
                                          family = gaussian())
      
      # 5. predicting from this model setting T = u for all individuals provides an estimate of E[Y_Vu | Z = 0, dSu = 1, X]
      if(v_umin == 0) v_umin <- 1
      for(u in ceiling(v_umin):floor(v_umax)){
        data_u <- data; data_u$Ti <- u
        pred_Y_out__Z0_Suminus1_1_X[,u] <- predict(fit_Y_out__Z0_Suminus1_1_X_T,
                                                   newdata = data_u,
                                                   type = 'response')
      }
    }
  }
  
  ### If not dense in u (version 2- pool across visits): ###
  
  if(pool){
    # 1. subset to unvaccinated
    data_Z0 <- data[data$Z == 0,]
    
    pred_Y_out__Z0_Suminus1_1_X <- data.frame(setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = floor(umax))),
                                                       paste0("pred_", 1:floor(umax))))
    
    pred_data <- data.frame()
    
    # # 2. further subset to all individuals who fall ill on any u in U-1(v)
    for(i in 1:length(V_u_months)){
      
      # this should be min, max ranges for each visit
      v_umin <- V_u_week_interval[i]
      v_umax <- V_u_week_interval[i+1]
      
      data_Z0_Uminus1_v <- data_Z0[which(data_Z0$S_inf_time > v_umin & data_Z0$S_inf_time <= v_umax),] # bigger than min and less than max
      
      # 3. create a variable for each individual, say Ti, indicating the time period at which each individual
      #    falls ill, ie set Ti equal to the u such that dSu,i = 1
      
      # same as S_inf_time
      data_Z0_Uminus1_v$Ti <- data_Z0_Uminus1_v$S_inf_time
      
      # QUESTION also assume need to get Y_out same way as the other one?
      
      data_Z0_Uminus1_v$Y_out <- data_Z0_Uminus1_v[[paste0("Y_", V_u_months[i])]]
      data_Z0_Uminus1_v$V <- rep(paste0("v_", V_u_months[i]), nrow(data_Z0_Uminus1_v))
      
      pred_data <- rbind(pred_data, data_Z0_Uminus1_v)
    }
    
    pred_data$V <- factor(pred_data$V, levels = paste0("v_", V_u_months))
    
    # 6. Regress the pooled growth outcome createde on X, V, and T
    fit_Y_out__Z0_Suminus1_1_X_T_V <- glm(as.formula(Y_out__Z0_Suminus11_X_T_V_formula),
                                          data = pred_data,
                                          family = gaussian())
    
    for(u in 1:umax){
      
      temp_data <- data
      temp_data$Ti <- u
      
      v_col_val <- paste0("v_", V_u_months[findInterval(u, V_u_week_interval, rightmost.closed = FALSE, left.open = TRUE)])
      temp_data$V <- v_col_val
    
      pred_Y_out__Z0_Suminus1_1_X[,u] <- predict(fit_Y_out__Z0_Suminus1_1_X_T_V,
                                                   newdata = temp_data,
                                                   type = 'response')
      
    }
  }
  
  # Put the final estimand together ----------------------------------------------
  
  # P(S_umax = 1 | Z = 0) == incidence of infection in the unvaccinated? pre week 43?
  
  #make just preds not 4:end
  
  # cuminc_0_12 <- rowSums(density_Z0_0_12)
  # mean(cuminc_0_12)
  
  P_Sumax_1__Z0_X <- rowSums(density_Z0_0_12)
  #colnames(P_Sumax_1__Z0_X) <- paste0("lambda_",1:52)
  
  P_Sumax_1__Z0 <- mean(P_Sumax_1__Z0_X )
  
  # P_S_umax_1__Z_0 <- mean(data$S_inf[data$Z == 0 & data$S_inf_time<=umax])
  
  estimate <- vector("numeric", length = umax)
  estimate_Z0 <- vector("numeric", length = umax)
  estimate_Z1 <- vector("numeric", length = umax)
  
  for(u in 1:umax){
    P_dSu_1__Z_0_X <- density_Z0_0_12[,u]
    E_Y_Vu__Z1_Suminus1_0_X <- pred_Y_out__Z1_Suminus1_0_X[,paste0("pred_",u)]
    E_Y_Vu__Z0_Suminus1_1_X <- pred_Y_out__Z0_Suminus1_1_X[,paste0("pred_",u)]
    
    estimate[u] <- mean((P_dSu_1__Z_0_X / P_Sumax_1__Z0) * (E_Y_Vu__Z1_Suminus1_0_X - E_Y_Vu__Z0_Suminus1_1_X))
    
    estimate_Z1[u] <- mean((P_dSu_1__Z_0_X / P_Sumax_1__Z0) * (E_Y_Vu__Z1_Suminus1_0_X))
    estimate_Z0[u] <- mean((P_dSu_1__Z_0_X / P_Sumax_1__Z0) * (E_Y_Vu__Z0_Suminus1_1_X))
    
  }
  
  final_growth_effect <- sum(estimate)
  final_estimate_Z0 <- sum(estimate_Z0)
  final_estimate_Z1 <- sum(estimate_Z1)
  
  return(list(growth_effect = final_growth_effect, 
              estimate_Z0 = final_estimate_Z0,
              estimate_Z1 = final_estimate_Z1))
}

#' Function for nat inf estimation - call gcomp or aipw with appropriate args
#' 
#' @returns additive effect estimate
est_nat_inf <- function(data,
                        estimator,
                        pkg_models,
                        Y_name,
                        exclusion_restriction,
                        cross_world,
                        two_part_model, 
                        Z_name = "Z", 
                        X_name = "X",
                        S_name = "S_inf"){
  
  if(estimator == "gcomp"){
    all_est <- vegrowth::do_gcomp_nat_inf(data = data, 
                                      models = pkg_models,
                                      Z_name = Z_name,
                                      X_name = X_name, 
                                      exclusion_restriction = exclusion_restriction,
                                      cross_world = cross_world,
                                      two_part_model = two_part_model)
    
    est <- list(additive_effect = all_est['additive_effect'],
                additive_se = NA,
                if_matrix = NA) 
    
  } else if(estimator == "aipw"){
    all_est <- vegrowth::do_aipw_nat_inf(data = data, 
                                     models = pkg_models, 
                                     Y_name = Y_name,
                                     exclusion_restriction = exclusion_restriction,
                                     cross_world = cross_world,
                                     Z_name = Z_name, 
                                     X_name = X_name, 
                                     S_name = S_name,
                                     return_se = TRUE,
                                     two_part_model = two_part_model)
    
    est <- list(additive_effect = all_est['additive_effect'],
                additive_se = all_est['additive_se'],
                if_matrix = attr(all_est, "if_matrix")) # workaround for package typing issues

  } else if (estimator == "unadj"){
    all_est <- vegrowth::do_unadj_nat_inf(data = data,
                                      Z_name = "Z",
                                      Y_name = Y_name,
                                      S_name = "S_inf")
    
    est <- list(additive_effect = all_est['additive_effect'],
                additive_se = NA,
                if_matrix = NA) 
  } else{
    stop("Unknown estimator: ", estimator)
  }
  
  return(est)
  
}


#' Function for population effect estimation - call gcomp or aipw with appropriate args
#' 
#' @returns additive effect estimate
est_pop <- function(data,
                    estimator,
                    pkg_models,
                    Y_name,
                    two_part_model, 
                    Z_name = "Z", 
                    X_name = "X",
                    S_name = "S_inf"){
  
  if(estimator == "gcomp"){
    all_est <- vegrowth::do_gcomp_pop(data = data, 
                                  models = pkg_models,
                                  Z_name = Z_name,
                                  X_name = X_name,
                                  two_part_model = two_part_model)
    
    est <- list(additive_effect = all_est['additive_effect'],
                additive_se = NA,
                if_matrix = NA) 
    
  } else if(estimator == "aipw"){
    all_est <- vegrowth::do_aipw_pop(data = data, 
                                 models = pkg_models, 
                                 Y_name = Y_name,
                                 Z_name = Z_name, 
                                 X_name = X_name, 
                                 return_se = TRUE,
                                 two_part_model = two_part_model) 
    
    est <- list(additive_effect = all_est['additive_effect'],
                additive_se = all_est['additive_se'],
                if_matrix = attr(all_est, "if_matrix")) # workaround for package typing issues
    
  } else{
    stop("Unknown estimator: ", estimator)
  }
  
  return(est)
  
}

