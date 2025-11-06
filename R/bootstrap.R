
one_boot <- function(data, config, parameters){
  
  # Create bootstrap data
  n <- dim(data)[1]
  boot_row_idx <- sample(1:n, replace=TRUE)
  boot_data <- data[boot_row_idx,]

  # 1. Fit Models
  estimand <- c()
  
  if(any(c(config$nat_inf_ER_1, config$nat_inf_ER_2, config$nat_inf_no_ER) == TRUE)){
    estimand <- c(estimand, "nat_inf")
  } 
  
  if(any(config$population_1, config$population_2) == TRUE){
    estimand <- c(estimand, "pop")
  }
  
  # Get unique timepoints needed in config$intervals
  Y_out <- unique(do.call(c, config$intervals))
  
  results <- expand.grid(Y_out = paste0("Y_", Y_out),
                         estimator = config$estimators,
                         nat_inf_ER_1 = NA,
                         nat_inf_ER_2 = NA,
                         nat_inf_no_ER = NA,
                         nat_inf_unadj = NA,
                         pop_1 = NA,
                         pop_2 = NA)
  
  # Get effect for all individual Y_outs
  for(i in 1:length(Y_out)){
    
    # Name of outcome variable
    Y_name <- paste0("Y_", Y_out[i])
    
    # Fit models
    if(any(c(config$nat_inf_ER_1, config$nat_inf_ER_2, config$nat_inf_no_ER, config$population_1, config$population_2) == TRUE)){
      pkg_models <- vegrowth::fit_models(data = boot_data,
                                         Y_name = Y_name, 
                                         Z_name = "Z", 
                                         X_name = "X", 
                                         S_name = "S_inf", 
                                         estimand = estimand, 
                                         method = config$estimators, 
                                         exclusion_restriction = TRUE, 
                                         family = "gaussian")
    }
    
    for(j in 1:length(config$estimators)){
      
      estimator <- config$estimators[j]
      
      if(config$nat_inf_ER_1){
        results$nat_inf_ER_1[results$Y_out == Y_name & results$estimator == estimator] <- est_nat_inf(data = boot_data, 
                                                                                                      pkg_models = pkg_models,
                                                                                                      Y_name = Y_name, 
                                                                                                      exclusion_restriction = TRUE, 
                                                                                                      two_part_model = FALSE, 
                                                                                                      estimator = estimator)
      }
      
      if(config$nat_inf_ER_2){
        results$nat_inf_ER_2[results$Y_out == Y_name & results$estimator == estimator] <- est_nat_inf(data = boot_data, 
                                                                                                      pkg_models = pkg_models,
                                                                                                      Y_name = Y_name, 
                                                                                                      exclusion_restriction = TRUE, 
                                                                                                      two_part_model = TRUE, 
                                                                                                      estimator = estimator)
      }
      
      if(config$nat_inf_no_ER){
        results$nat_inf_no_ER[results$Y_out == Y_name & results$estimator == estimator] <- est_nat_inf(data = boot_data, 
                                                                                                       pkg_models = pkg_models,
                                                                                                       Y_name = Y_name, 
                                                                                                       exclusion_restriction = FALSE, 
                                                                                                       two_part_model = FALSE, 
                                                                                                       estimator = estimator)
      }
      
      if(config$nat_inf_unadj){
        # same regardless of estimator; just do for j == 1
        if(j == 1){
          results$nat_inf_unadj[results$Y_out == Y_name & results$estimator == estimator] <-  vegrowth::do_unadj_nat_inf(data = boot_data,
                                                                                                                         Z_name = "Z",
                                                                                                                         Y_name = Y_name,
                                                                                                                         S_name = "S_inf")['additive_effect']
        } 
        
      }
      
      if(config$population_1){
        results$pop_1[results$Y_out == Y_name & results$estimator == estimator] <- est_pop(data = boot_data,
                                                                                           pkg_models = pkg_models,
                                                                                           Y_name = Y_name, 
                                                                                           estimator = estimator, 
                                                                                           two_part_model = FALSE)
      }
      
      if(config$population_2){
        results$pop_2[results$Y_out == Y_name & results$estimator == estimator] <- est_pop(data = boot_data,
                                                                                           pkg_models = pkg_models,
                                                                                           Y_name = Y_name, 
                                                                                           estimator = estimator, 
                                                                                           two_part_model = TRUE)
      }
      
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
    
    for(estimator in config$estimators){
      sub_df <- results[results$Y_out %in% Y_names &
                          results$estimator == estimator, ]
      
      # Create a new row with the averaged estimates
      new_row <- data.frame(
        Y_out = avg_name,
        estimator = estimator, 
        nat_inf_ER_1 =  mean(sub_df$nat_inf_ER_1),
        nat_inf_ER_2 =  mean(sub_df$nat_inf_ER_2),
        nat_inf_no_ER = mean(sub_df$nat_inf_no_ER),
        nat_inf_unadj = mean(sub_df$nat_inf_unadj),
        pop_1 = mean(sub_df$pop_1),
        pop_2 = mean(sub_df$pop_2)
      )
      
      # Append to results
      results <- rbind(results, new_row)
    }
    
  }
  
  return(results)
  
}
