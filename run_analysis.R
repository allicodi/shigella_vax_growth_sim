# ------------------------------------------------------------------------------
# Script to run analysis for given configuration settings 
# ------------------------------------------------------------------------------

.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

options(echo = TRUE)

here::i_am("run_analysis.R")

source(here::here("R/parameter_generation_fns.R"))
source(here::here("R/simulate_data.R"))
source(here::here("R/estimation_fn.R"))
source(here::here("R/bootstrap.R"))

# get seed & config settings from bash script
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
setting <- Sys.getenv("SETTING")
parameters <- readRDS(Sys.getenv("PARAMETERS_FILE"))

# issues with weeks not being overwritten
#config <- config::get(file = "config.yml", config = setting)

cfg <- yaml::read_yaml("config.yml")
config <- cfg[[setting]]

results <- lapply(config$n_sample_size, function(n){
  
  # Simulate data
  data <- simulate_data(parameters = parameters, 
                        n = n, 
                        VE_mild = config$VE_mild, 
                        VE_severe = config$VE_severe, 
                        seed = seed, 
                        type = "observed")
  
  # Long term & population effect estimation -----------------------------
  
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
      pkg_models <- vegrowth::fit_models(data = data,
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
        results$nat_inf_ER_1[results$Y_out == Y_name & results$estimator == estimator] <- est_nat_inf(data = data, 
                                                                                                      pkg_models = pkg_models,
                                                                                                      Y_name = Y_name, 
                                                                                                      exclusion_restriction = TRUE, 
                                                                                                      two_part_model = FALSE, 
                                                                                                      estimator = estimator)
      }
      
      if(config$nat_inf_ER_2){
        results$nat_inf_ER_2[results$Y_out == Y_name & results$estimator == estimator] <- est_nat_inf(data = data, 
                                                                                                     pkg_models = pkg_models,
                                                                                                     Y_name = Y_name, 
                                                                                                     exclusion_restriction = TRUE, 
                                                                                                     two_part_model = TRUE, 
                                                                                                     estimator = estimator)
      }
      
      if(config$nat_inf_no_ER){
        results$nat_inf_no_ER[results$Y_out == Y_name & results$estimator == estimator] <- est_nat_inf(data = data, 
                                                                                                       pkg_models = pkg_models,
                                                                                                       Y_name = Y_name, 
                                                                                                       exclusion_restriction = FALSE, 
                                                                                                       two_part_model = FALSE, 
                                                                                                       estimator = estimator)
      }
      
      if(config$nat_inf_unadj){
        # same regardless of estimator; just do for j == 1
        if(j == 1){
          results$nat_inf_unadj[results$Y_out == Y_name & results$estimator == estimator] <-  vegrowth::do_unadj_nat_inf(data = data,
                                                                                                                        Z_name = "Z",
                                                                                                                        Y_name = Y_name,
                                                                                                                        S_name = "S_inf")['additive_effect']
        } 
        
      }
      
      if(config$population_1){
        results$pop_1[results$Y_out == Y_name & results$estimator == estimator] <- est_pop(data = data,
                                                                                            pkg_models = pkg_models,
                                                                                            Y_name = Y_name, 
                                                                                            estimator = estimator, 
                                                                                            two_part_model = FALSE)
      }
      
      if(config$population_2){
        results$pop_2[results$Y_out == Y_name & results$estimator == estimator] <- est_pop(data = data,
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
  
  # Bootstrap Estimates ----------------------------------------
  
  # 1. Do n_boot bootstrap replicates
  boot_res_list <- replicate(config$n_boot, one_boot(data, config, parameters), simplify = FALSE)
  boot_res_df <- dplyr::bind_rows(boot_res_list, .id = "boot_id")
  
  # 2. Compute SEs, CIs, and rejection indicators per Y_out and estimand
  boot_summary <- boot_res_df %>%
    dplyr::group_by(Y_out, estimator) %>%
    dplyr::summarise(
      nat_inf_ER_1_se = sd(nat_inf_ER_1, na.rm = TRUE),
      nat_inf_ER_1_lower = quantile(nat_inf_ER_1, 0.025, na.rm = TRUE),
      nat_inf_ER_1_upper = quantile(nat_inf_ER_1, 0.975, na.rm = TRUE),
      
      nat_inf_ER_2_se = sd(nat_inf_ER_2, na.rm = TRUE),
      nat_inf_ER_2_lower = quantile(nat_inf_ER_2, 0.025, na.rm = TRUE),
      nat_inf_ER_2_upper = quantile(nat_inf_ER_2, 0.975, na.rm = TRUE),
      
      nat_inf_no_ER_se = sd(nat_inf_no_ER, na.rm = TRUE),
      nat_inf_no_ER_lower = quantile(nat_inf_no_ER, 0.025, na.rm = TRUE),
      nat_inf_no_ER_upper = quantile(nat_inf_no_ER, 0.975, na.rm = TRUE),
      
      nat_inf_unadj_se = sd(nat_inf_unadj, na.rm = TRUE),
      nat_inf_unadj_lower = quantile(nat_inf_unadj, 0.025, na.rm = TRUE),
      nat_inf_unadj_upper = quantile(nat_inf_unadj, 0.975, na.rm = TRUE),
      
      pop_1_se = sd(pop_1, na.rm = TRUE),
      pop_1_lower = quantile(pop_1, 0.025, na.rm = TRUE) ,
      pop_1_upper = quantile(pop_1, 0.975, na.rm = TRUE) ,
      
      pop_2_se = sd(pop_2, na.rm = TRUE),
      pop_2_lower = quantile(pop_2, 0.025, na.rm = TRUE) ,
      pop_2_upper = quantile(pop_2, 0.975, na.rm = TRUE) ,
      
      .groups = "drop"
    )
  
  # 3. Merge bootstrap summaries with point estimates
  results_full <- dplyr::left_join(results, boot_summary, by = c("Y_out", "estimator"))
  
  # 4. Add reject indicator columns
  if (config$nat_inf_ER_1) {
    results_full <- results_full %>%
      dplyr::mutate(
        nat_inf_ER_1_reject = (abs(nat_inf_ER_1 - config$null_hypothesis_value) / nat_inf_ER_1_se) > qnorm(1 - config$alpha_level / 2)
      )
  }
  
  if (config$nat_inf_ER_2) {
    results_full <- results_full %>%
      dplyr::mutate(
        nat_inf_ER_2_reject = (abs(nat_inf_ER_2 - config$null_hypothesis_value) / nat_inf_ER_2_se) > qnorm(1 - config$alpha_level / 2)
      )
  }
  
  if (config$nat_inf_no_ER) {
    results_full <- results_full %>%
      dplyr::mutate(
        nat_inf_no_ER_reject = (abs(nat_inf_no_ER - config$null_hypothesis_value) / nat_inf_no_ER_se) > qnorm(1 - config$alpha_level / 2)
      )
  }
  
  if (config$nat_inf_unadj) {
    results_full <- results_full %>%
      dplyr::mutate(
        nat_inf_unadj_reject = (abs(nat_inf_unadj - config$null_hypothesis_value) / nat_inf_unadj_se) > qnorm(1 - config$alpha_level / 2)
      )
  }
  
  if (config$population_1) {
    results_full <- results_full %>%
      dplyr::mutate(
        pop_1_reject = (abs(pop_1 - config$null_hypothesis_value) / pop_1_se) > qnorm(1 - config$alpha_level / 2)
      )
  }
  
  if (config$population_2) {
    results_full <- results_full %>%
      dplyr::mutate(
        pop_2_reject = (abs(pop_2 - config$null_hypothesis_value) / pop_2_se) > qnorm(1 - config$alpha_level / 2)
      )
  }
  
  # 5. Final combined result
  result <- data.frame(
    seed = seed,
    n = n,
    results_full
  )
  
  return(result)
  
})

results <- as.data.frame(do.call(rbind, results))

saveRDS(results, paste0("/projects/dbenkes/allison/shigella_vaccine_trial/results/", setting, "_seed_", seed, ".Rds"))
