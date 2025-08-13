
one_boot <- function(data, config, parameters){
  
  # Create bootstrap data
  n <- dim(data)[1]
  boot_row_idx <- sample(1:n, replace=TRUE)
  boot_data <- data[boot_row_idx,]
  
  # Long term & population effect estimation -----------------------------
  
  # 1. Fit Models
  if(config$long_term | config$population){
    estimand <- c()
    if(config$long_term){
      estimand <- c(estimand, "nat_inf")
    } 
    
    if(config$population){
      estimand <- c(estimand, "pop")
    }
    
    pkg_models <- vegrowth::fit_models(data = boot_data,
                                       Y_name = "Y_12",
                                       Z_name = "Z", 
                                       X_name = "X", 
                                       S_name = "S_inf", 
                                       estimand = estimand, 
                                       method = "gcomp",
                                       family = "gaussian")
  } else{
    pkg_models <- NULL
  }
  
  # 2. Call estimation functions from package
  if(config$long_term){
    est_long_term <- vegrowth::do_gcomp_nat_inf(data = boot_data, models = pkg_models)['additive_effect']
  } else{
    est_long_term <- NULL
  }
  
  if(config$population){
    est_pop <- vegrowth::do_gcomp_pop(data = boot_data, models = pkg_models, Z_name = "Z", X_name = "X")['additive_effect']
  } else{
    est_pop <- NULL
  }
  
  # Short-term effect estimation -----------------------------
  
  # 1. Call effect estimation function (models fit within)
  if(config$short_term){
    est_short_term <- estimate_short_term(data = boot_data, 
                                          parameters = parameters, 
                                          V_u_months = as.numeric(config$V_u_months), 
                                          V_u_week_interval = as.numeric(config$V_u_week_interval), 
                                          pool = config$pool, 
                                          I_S__X_Z0_0_6_formula = config$I_S__X_Z0_0_6_formula,
                                          I_S__X_Z0_6_12_formula = config$I_S__X_Z0_6_12_formula, 
                                          Y_out__Z1_Suminus10_X_formula = config$Y_out__Z1_Suminus10_X_formula,
                                          Y_out__Z0_Suminus11_X_T_formula = config$Y_out__Z0_Suminus11_X_T_formula, 
                                          Y_out__Z0_Suminus11_X_T_V_formula = config$Y_out__Z0_Suminus11_X_T_V_formula)
  } else{
    est_short_term <- NULL
  }
  
  return(c(est_short_term = as.numeric(est_short_term),
                    est_pop = as.numeric(est_pop),
                    est_long_term = as.numeric(est_long_term)))
  
}