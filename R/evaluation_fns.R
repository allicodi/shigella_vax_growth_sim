# ----------------------------------------------------------------------------
# Functions for evaluating performance for growth effect estimation
# ----------------------------------------------------------------------------

get_bias <- function(results, truth, config, n) {
  results_n <- results[results$n == n, ]
  out <- data.frame()
  
  for (i in seq_along(config$intervals)) {
    int <- config$intervals[[i]]
    Y_out <- paste0("Y_", paste0(int, collapse = "_"))
    
    # nat_inf_ER
    if (config$nat_inf_ER) {
      bias_nat_inf <- mean(
        as.numeric(results_n$nat_inf_ER[results_n$Y_out == Y_out]) -
          truth[[paste0("nat_inf_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf_ER",
        bias = bias_nat_inf
      ))
    }
    
    # nat_inf_no_ER
    if (config$nat_inf_no_ER) {
      bias_nat_inf <- mean(
        as.numeric(results_n$nat_inf[results_n$Y_out == Y_out]) -
          truth[[paste0("nat_inf_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf_no_ER",
        bias = bias_nat_inf
      ))
    }
    
    # nat_inf_unadj
    if (config$nat_inf_unadj) {
      bias_nat_inf <- mean(
        as.numeric(results_n$nat_inf_unadj[results_n$Y_out == Y_out]) -
          truth[[paste0("nat_inf_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf_unadj",
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
    
    # No ER
    if (config$nat_inf_no_ER) {
      cov_nat_inf <- mean(
        as.numeric(results_n$nat_inf_no_ER_lower[results_n$Y_out == Y_out]) < truth[[paste0("nat_inf_", Y_out)]] &
          as.numeric(results_n$nat_inf_no_ER_upper[results_n$Y_out == Y_out]) > truth[[paste0("nat_inf_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf_no_ER",
        coverage = cov_nat_inf
      ))
    }
    
    # ER
    if (config$nat_inf_ER) {
      cov_nat_inf <- mean(
        as.numeric(results_n$nat_inf_ER_lower[results_n$Y_out == Y_out]) < truth[[paste0("nat_inf_", Y_out)]] &
          as.numeric(results_n$nat_inf_ER_upper[results_n$Y_out == Y_out]) > truth[[paste0("nat_inf_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf_ER",
        coverage = cov_nat_inf
      ))
    }
    
    # unadj
    if (config$nat_inf_unadj) {
      cov_nat_inf <- mean(
        as.numeric(results_n$nat_inf_unadj_lower[results_n$Y_out == Y_out]) < truth[[paste0("nat_inf_", Y_out)]] &
          as.numeric(results_n$nat_inf_unadj_upper[results_n$Y_out == Y_out]) > truth[[paste0("nat_inf_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf_unadj",
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
    
    # No ER
    if (config$nat_inf_no_ER) {
      power_nat_inf <- mean(results_n$nat_inf_no_ER_reject[results_n$Y_out == Y_out])
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf_no_ER",
        power = power_nat_inf
      ))
    }
    
    # ER
    if (config$nat_inf_ER) {
      power_nat_inf <- mean(results_n$nat_inf_ER_reject[results_n$Y_out == Y_out])
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf_ER",
        power = power_nat_inf
      ))
    }
    
    # unadj
    if (config$nat_inf_unadj) {
      power_nat_inf <- mean(results_n$nat_inf_unadj_reject[results_n$Y_out == Y_out])
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf_unadj",
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
    
    if (config$nat_inf_ER) {
      neg_nat_inf <- mean(results_n$nat_inf_ER[results_n$Y_out == Y_out] < 0)
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf_ER",
        prop_neg = neg_nat_inf
      ))
    }
    
    if (config$nat_inf_no_ER) {
      neg_nat_inf <- mean(results_n$nat_inf_no_ER[results_n$Y_out == Y_out] < 0)
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf_no_ER",
        prop_neg = neg_nat_inf
      ))
    }
    
    if (config$nat_inf_unadj) {
      neg_nat_inf <- mean(results_n$nat_inf_unadj[results_n$Y_out == Y_out] < 0)
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimand = "nat_inf_unadj",
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
