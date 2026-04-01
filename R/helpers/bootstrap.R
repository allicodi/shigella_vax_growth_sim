
one_boot <- function(data, config, setting_grid, parameters){
  
  # Create bootstrap data
  n <- dim(data)[1]
  boot_row_idx <- sample(1:n, replace=TRUE)
  boot_data <- data[boot_row_idx,]

  # Long term & population effect estimation -----------------------------
  
  res_df <- data.frame()
  
  # 1. Fit Models
  
  # Get unique timepoints needed in config$intervals
  Y_out <- unique(do.call(c, config$intervals))
  
  # Get effect for all individual Y_outs
  for(i in 1:length(Y_out)){
    
    # Name of outcome variable
    Y_name <- paste0("Y_", Y_out[i])
    
    # Only fit models for gcomp (aipw uses closed form se, unadj doesn't need models)
    if(nrow(setting_grid > 0) & c("gcomp" %in% setting_grid$estimator)){
      pkg_models <- vaxstrat::fit_models(data = boot_data,
                                         Y_name = Y_name, 
                                         Z_name = "Z", 
                                         X_name = "X", 
                                         S_name = "S_inf", 
                                         estimand = config$estimand, 
                                         method = config$estimators, 
                                         exclusion_restriction = TRUE, 
                                         family = "gaussian")
    } else{
      pkg_models <- NULL
    }
    
    for(j in 1:nrow(setting_grid)){
      setting <- setting_grid[j,]
      
      if(setting$estimand == "nat_inf"){
        res <- est_nat_inf(data = boot_data,
                           estimator = setting$estimator,
                           pkg_models = pkg_models,
                           Y_name = Y_name, 
                           exclusion_restriction = setting$er,
                           cross_world = setting$cw,
                           two_part_model = setting$two_stage)
        
        
      } else{
        res <- est_pop(data = boot_data,
                       estimator = setting$estimator,
                       pkg_models = pkg_models,
                       Y_name = Y_name, 
                       two_part_model = setting$two_stage)
      }
          
      row <- data.frame(Y_out = Y_name, estimate = as.numeric(res['additive_effect']), se = NA, setting)
      res_df <- rbind(res_df, row)
      
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
    for(j in 1:nrow(setting_grid)){
      setting <- setting_grid[j,]
      
      if(setting$estimator != "unadj"){
        sub_df <- res_df[res_df$Y_out %in% Y_names &
                           res_df$estimator == setting$estimator &
                           res_df$estimand == setting$estimand &
                           res_df$er == setting$er &
                           res_df$cw == setting$cw &
                           res_df$two_stage == setting$two_stage,]
      } else{
        sub_df <- res_df[res_df$Y_out %in% Y_names &
                           res_df$estimator == setting$estimator &
                           res_df$estimand == setting$estimand,]
      }
      
      row <- data.frame(Y_out = avg_name, 
                        estimate = mean(sub_df$estimate),
                        se = NA,
                        setting)
      
      res_df <- rbind(res_df, row)
    }
    
  }
  
  return(res_df)
  
}
