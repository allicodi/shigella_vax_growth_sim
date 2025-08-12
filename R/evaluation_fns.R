# ----------------------------------------------------------------------------
# Functions for evaluating performance for growth effect estimation
# ----------------------------------------------------------------------------

get_bias <- function(results, truth, n){
  
  # get subset of results with matching sample size
  results_n <- results[results$n == n,]
  
  out <- list(n = n)
  
  if(!is.null(results_n$est_long_term)){
    out$bias_long_term <- mean(results_n$est_long_term - truth$long_term)
  }
  
  if(!is.null(results_n$est_short_term)){
    out$bias_short_term <- mean(results_n$est_short_term - truth$short_term)
  }
  
  if(!is.null(results_n$est_pop)){
    out$bias_pop <- mean(results_n$est_pop - truth$population)
  }
  
  return(out)
  
}

get_coverage <- function(results, truth, n){
  
  # get subset of results with matching sample size
  results_n <- results[results$n == n,]
  
  out <- list(n = n)
  
  if(!is.null(results_n$est_long_term)){
    coverage_vec_long_term <- ifelse(results_n$long_term_lower_ci < truth$long_term &
                                       results_n$long_term_upper_ci > truth$long_term, 1, 0)
    out$coverage_long_term <- mean(coverage_vec_long_term)
  }
  
  if(!is.null(results_n$est_short_term)){
    coverage_vec_short_term <- ifelse(results_n$short_term_lower_ci < truth$short_term &
                                       results_n$short_term_upper_ci > truth$short_term, 1, 0)
    out$coverage_short_term <- mean(coverage_vec_short_term)
  }
  
  if(!is.null(results_n$est_pop)){
    coverage_vec_pop <- ifelse(results_n$pop_lower_ci < truth$population &
                                       results_n$pop_upper_ci > truth$population, 1, 0)
    out$coverage_pop <- mean(coverage_vec_pop)
  }
  
  return(out)
  
}


get_power <- function(results, truth, n){
  
  # get subset of results with matching sample size
  results_n <- results[results$n == n,]
  
  out <- list(n = n)
  
  if(!is.null(results_n$est_long_term)){
    power_vec_long_term <- ifelse(results_n$long_term_reject == TRUE, 1, 0 )
    out$power_long_term <- mean(power_vec_long_term)
  }
  
  if(!is.null(results_n$est_short_term)){
    power_vec_short_term <- ifelse(results_n$short_term_reject == TRUE, 1, 0 )
    out$power_short_term <- mean(power_vec_short_term)
  }
  
  if(!is.null(results_n$est_pop)){
    power_vec_pop <- ifelse(results_n$pop_reject == TRUE, 1, 0 )
    out$power_pop <- mean(power_vec_pop)
  }
  
  return(out)
  
}

get_neg_pt_est <- function(results, truth, n){
  
  # get subset of results with matching sample size
  results_n <- results[results$n == n,]
  
  out <- list(n = n)
  
  if(!is.null(results_n$est_long_term)){
    out$prop_neg_long_term <- mean(ifelse(results_n$est_long_term < 0, 1, 0 ))
  }
  
  if(!is.null(results_n$est_short_term)){
    out$prop_neg_short_term <- mean(ifelse(results_n$est_short_term < 0, 1, 0 ))
  }
  
  if(!is.null(results_n$est_pop)){
    out$prop_neg_pop <- mean(ifelse(results_n$est_pop < 0, 1, 0 ))
  }
  
  return(out)
  
}



