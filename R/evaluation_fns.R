# ----------------------------------------------------------------------------
# Functions for evaluating performance for growth effect estimation
# ----------------------------------------------------------------------------

#' Function to compute bias for vaccine trial simulations
#' 
#' @param results dataframe containing results from simulations
#' @param truth list containing truth for each estimand and endpoint combination
#' @param config config for given setting
get_bias <- function(results, truth, config) {
  # compute true value for each row based on estimand and Y_out
  results <- results %>%
    dplyr::mutate(
      truth_val = purrr::map2_dbl(
        estimand, Y_out,
        ~ truth[[paste0(.x, "_", .y)]]
      ),
      bias = estimate - truth_val
    )
  
  # summarize bias by relevant grouping structure
  results_summary <- results %>%
    dplyr::group_by(Y_out, estimand, estimator, er, cw, two_stage, n) %>%
    dplyr::summarise(
      bias = mean(bias, na.rm = TRUE),
      sd_bias = sd(bias, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(results_summary)
}

#' Function to compute CI coverage
#' 
#' @param results dataframe containing results from simulations
#' @param truth list containing truth for each estimand and endpoint combination
#' @param config config for given setting
get_coverage <- function(results, truth, config) {
  results <- results %>%
    dplyr::mutate(
      truth_val = purrr::map2_dbl(
        estimand, Y_out,
        ~ truth[[paste0(.x, "_", .y)]]
      ),
      covered = (lower_ci < truth_val) & (upper_ci > truth_val)
    )
  
  results_summary <- results %>%
    dplyr::group_by(n, estimator, Y_out, estimand, er, cw, two_stage) %>%
    dplyr::summarise(
      coverage = mean(covered, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(results_summary)
}

#' Function to compute power for vaccine trial simulations
#' 
#' @param results dataframe containing results from simulations
#' @param truth list containing truth for each estimand and endpoint combination
#' @param config config for given setting
get_power <- function(results, truth, config) {
  results_summary <- results %>%
    dplyr::group_by(n, estimator, Y_out, estimand, er, cw, two_stage) %>%
    dplyr::summarise(
      power = mean(as.numeric(reject), na.rm = TRUE),
      .groups = "drop"
    )
  
  return(results_summary)
}

#' Function to compute proportion of negative point estimates for vaccine trial simulations
#' 
#' @param results dataframe containing results from simulations
#' @param truth list containing truth for each estimand and endpoint combination
#' @param config config for given setting
get_neg_pt_est <- function(results, truth, config) {
  results_summary <- results %>%
    dplyr::group_by(n, estimator, Y_out, estimand, er, cw, two_stage) %>%
    dplyr::summarise(
      prop_neg = mean(estimate < 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(results_summary)
}
