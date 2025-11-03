# ----------------------------------------------------------------------------
# Functions for evaluating performance for growth effect estimation
# ----------------------------------------------------------------------------

get_bias <- function(results, truth, config, n){
  
  # get subset of results with matching sample size
  results_n <- results[results$n == n,]
  
  out <- list(n = n)
  
  for(i in 1:length(config$intervals)){
    int <- config$intervals[[i]]
    Y_out <- paste0("Y_", paste0(config$intervals[[i]], collapse = "_"))
    
    if(config$nat_inf){
      out[[Y_out]]$nat_inf <- mean(as.numeric(results_n$nat_inf[results_n$Y_out == Y_out]) - 
                                     truth[[paste0("nat_inf_", Y_out)]])
    }
    
    if(config$population){
      out[[Y_out]]$pop <- mean(as.numeric(results_n$pop[results_n$Y_out == Y_out]) - 
                                     truth[[paste0("pop_", Y_out)]])
    }
    
  }
  
  return(out)
  
}

get_coverage <- function(results, truth, config, n){
  
  # get subset of results with matching sample size
  results_n <- results[results$n == n,]
  
  out <- list(n = n)
  
  for(i in 1:length(config$intervals)){
    int <- config$intervals[[i]]
    Y_out <- paste0("Y_", paste0(int, collapse = "_"))
    
    if(config$nat_inf){
      coverage_vec_nat_inf <- ifelse(as.numeric(results_n$nat_inf_lower[results_n$Y_out == Y_out]) < truth[[paste0("nat_inf_", Y_out)]] &
                                       as.numeric(results_n$nat_inf_upper[results_n$Y_out == Y_out]) > truth[[paste0("nat_inf_", Y_out)]], 1, 0)
      out[[Y_out]]$nat_inf <- mean(coverage_vec_nat_inf)
    }
    
    if(config$pop){
      coverage_vec_pop <- ifelse(as.numeric(results_n$pop_lower[results_n$Y_out == Y_out]) < truth[[paste0("pop_", Y_out)]] &
                                       as.numeric(results_n$pop_upper[results_n$Y_out == Y_out]) > truth[[paste0("pop_", Y_out)]], 1, 0)
      out[[Y_out]]$pop <- mean(coverage_vec_pop)
    }
    
  }
  
  return(out)
  
}


get_power <- function(results, truth, config, n){
  
  # get subset of results with matching sample size
  results_n <- results[results$n == n,]
  
  out <- list(n = n)
  
  for(i in 1:length(config$intervals)){
    int <- config$intervals[[i]]
    Y_out <- paste0("Y_", paste0(int, collapse = "_"))
    
    if(config$nat_inf){
      power_vec_nat_inf <- ifelse(results_n$nat_inf_reject[results_n$Y_out == Y_out] == TRUE, 1, 0)
      out[[Y_out]]$nat_inf <- mean(power_vec_nat_inf)
    }
    
    if(config$nat_inf){
      power_vec_pop <- ifelse(results_n$pop_reject[results_n$Y_out == Y_out] == TRUE, 1, 0)
      out[[Y_out]]$pop <- mean(power_vec_pop)
    }
    
  }
  
  return(out)
  
}

get_neg_pt_est <- function(results, truth, config, n){
  
    # get subset of results with matching sample size
  results_n <- results[results$n == n,]
  
  out <- list(n = n)
  
  for(i in 1:length(config$intervals)){
    int <- config$intervals[[i]]
    Y_out <- paste0("Y_", paste0(config$intervals[[i]], collapse = "_"))
    
    if(config$nat_inf){
      out[[Y_out]]$nat_inf <- mean(ifelse(results_n$nat_inf[results_n$Y_out == Y_out] < 0, 1, 0))
    }
    
    if(config$population){
      out[[Y_out]]$pop <- mean(ifelse(results_n$pop[results_n$Y_out == Y_out] < 0, 1, 0))
    }
    
  }

  return(out)
  
}



