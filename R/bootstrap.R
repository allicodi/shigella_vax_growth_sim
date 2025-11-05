
one_boot <- function(data, config, parameters){
  
  # Create bootstrap data
  n <- dim(data)[1]
  boot_row_idx <- sample(1:n, replace=TRUE)
  boot_data <- data[boot_row_idx,]
  #boot_data$og_id <- boot_data$id
  #boot_data$id <- 1:nrow(boot_data)
  
  # 1. Fit Models
  estimand <- c()
  
  if(any(c(config$nat_inf_ER, config$nat_inf_no_ER) == TRUE)){
    estimand <- c(estimand, "nat_inf")
  } 
  
  if(config$population){
    estimand <- c(estimand, "pop")
  }
  
  # Get unique timepoints needed in config$intervals
  Y_out <- unique(do.call(c, config$intervals))
  
  results <- data.frame(Y_out = paste0("Y_", Y_out),
                        nat_inf_ER = NA,
                        nat_inf_no_ER = NA,
                        nat_inf_unadj = NA,
                        pop = NA)
  
  # Get effect for all individual Y_outs
  for(i in 1:length(Y_out)){
    
    # Name of outcome variable
    Y_name <- paste0("Y_", Y_out[i])
    
    # Fit models 
    if(any(c(config$nat_inf_ER, config$nat_inf_no_ER, config$population) == TRUE)){
      pkg_models <- vegrowth::fit_models(data = boot_data,
                                         Y_name = Y_name, 
                                         Z_name = "Z", 
                                         X_name = "X", 
                                         S_name = "S_inf", 
                                         estimand = estimand, 
                                         method = "gcomp", 
                                         exclusion_restriction = TRUE, 
                                         family = "gaussian")
    }
    
    # Call vegrowth functions for given outcome, nat inf and pop estimators
    if(config$nat_inf_ER){
      results$nat_inf_ER[i] <-  vegrowth::do_gcomp_nat_inf(data = boot_data, 
                                                           models = pkg_models,
                                                           Z_name = "Z",
                                                           X_name = "X", 
                                                           exclusion_restriction = TRUE)['additive_effect']
    }
    
    if(config$nat_inf_no_ER){
      results$nat_inf_no_ER[i] <-  vegrowth::do_gcomp_nat_inf(data = boot_data, 
                                                              models = pkg_models,
                                                              Z_name = "Z",
                                                              X_name = "X", 
                                                              exclusion_restriction = FALSE)['additive_effect']
    }
    
    if(config$nat_inf_unadj){
      results$nat_inf_unadj[i] <- vegrowth::do_unadj_nat_inf(data = boot_data,
                                                          Z_name = "Z",
                                                          Y_name = Y_name,
                                                          S_name = "S_inf")['additive_effect']
    }
    
    if(config$population){
      results$pop[i] <- vegrowth::do_gcomp_pop(data = boot_data, 
                                               models = pkg_models,
                                               Z_name = "Z",
                                               X_name = "X")['additive_effect']
    } 
    
    
  }
  
  # Get effect for averaged Y_outs
  which_intervals <- do.call(c, lapply(config$intervals, function(x) length(x) > 1))
  avg_intervals <- config$intervals[which_intervals]
  
  for(i in avg_intervals){
    # Create a name for the averaged interval (e.g., "Y_6_9")
    avg_name <- paste0("Y_", paste(i, collapse = "_"))
    
    # Get the corresponding Y_ variable names
    Y_names <- paste0("Y_", i)
    
    # Get the subset of results corresponding to those outcomes
    sub_df <- results[results$Y_out %in% Y_names, ]
    
    # Create a new row with the averaged estimates
    new_row <- data.frame(
      Y_out = avg_name,
      nat_inf_ER = if ("nat_inf_ER" %in% names(sub_df)) mean(sub_df$nat_inf_ER, na.rm = TRUE) else NA,
      nat_inf_no_ER = if ("nat_inf_no_ER" %in% names(sub_df)) mean(sub_df$nat_inf_no_ER, na.rm = TRUE) else NA,
      nat_inf_unadj = if ("nat_inf_unadj" %in% names(sub_df)) mean(sub_df$nat_inf_unadj, na.rm = TRUE) else NA,
      pop = if ("pop" %in% names(sub_df)) mean(sub_df$pop, na.rm = TRUE) else NA
    )
    
    # Append to results
    results <- rbind(results, new_row)
  }
  
  return(results)
  
}
