# ----------------------------------------------------------------------------
# Functions for evaluating performance for growth effect estimation
# ----------------------------------------------------------------------------

get_bias <- function(results, truth, config, n, estimator) {
  results_n <- results[results$n == n &
                         results$estimator == estimator, ]
  out <- data.frame()
  
  for (i in seq_along(config$intervals)) {
    int <- config$intervals[[i]]
    Y_out <- paste0("Y_", paste0(int, collapse = "_"))
    
    # nat_inf_ER_1
    if (config$nat_inf_ER_1) {
      bias_nat_inf <- mean(
        as.numeric(results_n$nat_inf_ER_1[results_n$Y_out == Y_out]) -
          truth[[paste0("nat_inf_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "nat_inf_ER_1",
        bias = bias_nat_inf
      ))
    }
    
    # nat_inf_ER_2
    if (config$nat_inf_ER_2) {
      bias_nat_inf <- mean(
        as.numeric(results_n$nat_inf_ER_1[results_n$Y_out == Y_out]) -
          truth[[paste0("nat_inf_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "nat_inf_ER_2",
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
        estimator = estimator,
        method = "nat_inf_no_ER",
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
        estimator = estimator,
        method = "nat_inf_unadj",
        bias = bias_nat_inf
      ))
    }
    
    # population_1
    if (config$population_1) {
      bias_pop <- mean(
        as.numeric(results_n$pop_1[results_n$Y_out == Y_out]) -
          truth[[paste0("pop_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        method = "pop_1",
        bias = bias_pop
      ))
    }
    
    # population_1
    if (config$population_2) {
      bias_pop <- mean(
        as.numeric(results_n$pop_2[results_n$Y_out == Y_out]) -
          truth[[paste0("pop_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "pop_2",
        bias = bias_pop
      ))
    }
  }
  
  return(out)
}


get_coverage <- function(results, truth, config, n, estimator) {
  results_n <- results[results$n == n &
                         results$estimator == estimator, ]
  out <- data.frame()
  
  for (i in seq_along(config$intervals)) {
    int <- config$intervals[[i]]
    Y_out <- paste0("Y_", paste0(int, collapse = "_"))
    
    # No ER
    if (config$nat_inf_no_ER_1) {
      cov_nat_inf <- mean(
        as.numeric(results_n$nat_inf_no_ER_lower[results_n$Y_out == Y_out]) < truth[[paste0("nat_inf_", Y_out)]] &
          as.numeric(results_n$nat_inf_no_ER_upper[results_n$Y_out == Y_out]) > truth[[paste0("nat_inf_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "nat_inf_no_ER",
        coverage = cov_nat_inf
      ))
    }
    
    # ER
    if (config$nat_inf_ER_1) {
      cov_nat_inf <- mean(
        as.numeric(results_n$nat_inf_ER_1_lower[results_n$Y_out == Y_out]) < truth[[paste0("nat_inf_", Y_out)]] &
          as.numeric(results_n$nat_inf_ER_1_upper[results_n$Y_out == Y_out]) > truth[[paste0("nat_inf_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "nat_inf_ER_1",
        coverage = cov_nat_inf
      ))
    }
    
    if (config$nat_inf_ER_2) {
      cov_nat_inf <- mean(
        as.numeric(results_n$nat_inf_ER_2_lower[results_n$Y_out == Y_out]) < truth[[paste0("nat_inf_", Y_out)]] &
          as.numeric(results_n$nat_inf_ER_2_upper[results_n$Y_out == Y_out]) > truth[[paste0("nat_inf_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "nat_inf_ER_2",
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
        estimator = estimator,
        method = "nat_inf_unadj",
        coverage = cov_nat_inf
      ))
    }
    
    if (config$population_1) {
      cov_pop <- mean(
        as.numeric(results_n$pop_1_lower[results_n$Y_out == Y_out]) < truth[[paste0("pop_", Y_out)]] &
          as.numeric(results_n$pop_1_upper[results_n$Y_out == Y_out]) > truth[[paste0("pop_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "pop_1",
        coverage = cov_pop
      ))
    }
    
    if (config$population_2) {
      cov_pop <- mean(
        as.numeric(results_n$pop_2_lower[results_n$Y_out == Y_out]) < truth[[paste0("pop_", Y_out)]] &
          as.numeric(results_n$pop_2_upper[results_n$Y_out == Y_out]) > truth[[paste0("pop_", Y_out)]]
      )
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "pop_2",
        coverage = cov_pop
      ))
    }
  }
  
  return(out)
}


get_power <- function(results, truth, config, n, estimator) {
  results_n <- results[results$n == n &
                         results$estimator == estimator, ]
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
        estimator = estimator,
        method = "nat_inf_no_ER",
        power = power_nat_inf
      ))
    }
    
    # ER
    if (config$nat_inf_ER_1) {
      power_nat_inf <- mean(results_n$nat_inf_ER_1_reject[results_n$Y_out == Y_out])
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "nat_inf_ER_1",
        power = power_nat_inf
      ))
    }
    
    # ER
    if (config$nat_inf_ER_2) {
      power_nat_inf <- mean(results_n$nat_inf_ER_2_reject[results_n$Y_out == Y_out])
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        method = "nat_inf_ER_2",
        power = power_nat_inf
      ))
    }
    
    # unadj
    if (config$nat_inf_unadj) {
      power_nat_inf <- mean(results_n$nat_inf_unadj_reject[results_n$Y_out == Y_out])
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "nat_inf_unadj",
        power = power_nat_inf
      ))
    }
    
    if (config$population_1) {
      power_pop <- mean(results_n$pop_1_reject[results_n$Y_out == Y_out])
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "pop_1",
        power = power_pop
      ))
    }
    
    if (config$population_2) {
      power_pop <- mean(results_n$pop_2_reject[results_n$Y_out == Y_out])
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "pop_2",
        power = power_pop
      ))
    }
  }
  
  return(out)
}


get_neg_pt_est <- function(results, truth, config, n, estimator) {
  results_n <- results[results$n == n &
                         results$estimator == estimator, ]
  out <- data.frame()
  
  for (i in seq_along(config$intervals)) {
    int <- config$intervals[[i]]
    Y_out <- paste0("Y_", paste0(int, collapse = "_"))
    
    if (config$nat_inf_ER_1) {
      neg_nat_inf <- mean(results_n$nat_inf_ER_1[results_n$Y_out == Y_out] < 0)
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        method = "nat_inf_ER_1",
        prop_neg = neg_nat_inf
      ))
    }
    
    if (config$nat_inf_ER_2) {
      neg_nat_inf <- mean(results_n$nat_inf_ER_2[results_n$Y_out == Y_out] < 0)
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "nat_inf_ER_2",
        prop_neg = neg_nat_inf
      ))
    }
    
    if (config$nat_inf_no_ER) {
      neg_nat_inf <- mean(results_n$nat_inf_no_ER[results_n$Y_out == Y_out] < 0)
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "nat_inf_no_ER",
        prop_neg = neg_nat_inf
      ))
    }
    
    if (config$nat_inf_unadj) {
      neg_nat_inf <- mean(results_n$nat_inf_unadj[results_n$Y_out == Y_out] < 0)
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "nat_inf_unadj",
        prop_neg = neg_nat_inf
      ))
    }
    
    if (config$population_1) {
      neg_pop <- mean(results_n$pop_1[results_n$Y_out == Y_out] < 0)
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "pop_1",
        prop_neg = neg_pop
      ))
    }
    
    if (config$population_2) {
      neg_pop <- mean(results_n$pop_2[results_n$Y_out == Y_out] < 0)
      out <- rbind(out, data.frame(
        n = n,
        Y_out = Y_out,
        estimator = estimator,
        method = "pop_2",
        prop_neg = neg_pop
      ))
    }
  }
  
  return(out)
}
