# ----------------------------------------------------------------------------
# Functions for evaluating performance for growth effect estimation
# ----------------------------------------------------------------------------

get_bias <- function(results, truth, config, n) {
  results_n <- results[results$n == n, ]
  out <- data.frame()
  
  for (i in seq_along(config$intervals)) {
    int <- config$intervals[[i]]
    Y_out <- paste0("Y_", paste0(int, collapse = "_"))
    
    # nat_inf
    if (config$nat_inf) {
      bias_nat_inf <- mean(
        as.numeric(results_n$nat_inf[results_n$Y_out == Y_out]) -
          truth[[paste0("nat_inf_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf",
        bias = bias_nat_inf
      ))
    }
    
    # population
    if (config$population) {
      bias_pop <- mean(
        as.numeric(results_n$pop[results_n$Y_out == Y_out]) -
          truth[[paste0("pop_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "pop",
        bias = bias_pop
      ))
    }
  }
  
  return(out)
}


get_coverage <- function(results, truth, config, n) {
  results_n <- results[results$n == n, ]
  out <- data.frame()
  
  for (i in seq_along(config$intervals)) {
    int <- config$intervals[[i]]
    Y_out <- paste0("Y_", paste0(int, collapse = "_"))
    
    if (config$nat_inf) {
      cov_nat_inf <- mean(
        as.numeric(results_n$nat_inf_lower[results_n$Y_out == Y_out]) < truth[[paste0("nat_inf_", Y_out)]] &
          as.numeric(results_n$nat_inf_upper[results_n$Y_out == Y_out]) > truth[[paste0("nat_inf_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf",
        coverage = cov_nat_inf
      ))
    }
    
    if (config$population) {
      cov_pop <- mean(
        as.numeric(results_n$pop_lower[results_n$Y_out == Y_out]) < truth[[paste0("pop_", Y_out)]] &
          as.numeric(results_n$pop_upper[results_n$Y_out == Y_out]) > truth[[paste0("pop_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "pop",
        coverage = cov_pop
      ))
    }
  }
  
  return(out)
}


get_power <- function(results, truth, config, n) {
  results_n <- results[results$n == n, ]
  out <- data.frame()
  
  for (i in seq_along(config$intervals)) {
    int <- config$intervals[[i]]
    Y_out <- paste0("Y_", paste0(int, collapse = "_"))
    
    if (config$nat_inf) {
      power_nat_inf <- mean(results_n$nat_inf_reject[results_n$Y_out == Y_out])
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf",
        power = power_nat_inf
      ))
    }
    
    if (config$population) {
      power_pop <- mean(results_n$pop_reject[results_n$Y_out == Y_out])
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "pop",
        power = power_pop
      ))
    }
  }
  
  return(out)
}


get_neg_pt_est <- function(results, truth, config, n) {
  results_n <- results[results$n == n, ]
  out <- data.frame()
  
  for (i in seq_along(config$intervals)) {
    int <- config$intervals[[i]]
    Y_out <- paste0("Y_", paste0(int, collapse = "_"))
    
    if (config$nat_inf) {
      neg_nat_inf <- mean(results_n$nat_inf[results_n$Y_out == Y_out] < 0)
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf",
        prop_neg = neg_nat_inf
      ))
    }
    
    if (config$population) {
      neg_pop <- mean(results_n$pop[results_n$Y_out == Y_out] < 0)
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "pop",
        prop_neg = neg_pop
      ))
    }
  }
  
  return(out)
}
