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
  
  if(config$nat_inf){
    estimand <- c(estimand, "nat_inf")
  } 
  
  if(config$population){
    estimand <- c(estimand, "pop")
  }
  
  # Get unique timepoints needed in config$intervals
  Y_out <- unique(do.call(c, config$intervals))
  
  results <- data.frame(Y_out = paste0("Y_", Y_out),
                        nat_inf = NA,
                        pop = NA)
  
  # Get effect for all individual Y_outs
  for(i in 1:length(Y_out)){
    
    # Name of outcome variable
    Y_name <- paste0("Y_", Y_out[i])
    
    # Fit models 
    pkg_models <- vegrowth::fit_models(data = data,
                                       Y_name = Y_name, 
                                       Z_name = "Z", 
                                       X_name = "X", 
                                       S_name = "S_inf", 
                                       estimand = estimand, 
                                       method = "gcomp", 
                                       exclusion_restriction = TRUE, 
                                       family = "gaussian")
    
    # Call vegrowth functions for given outcome, nat inf and pop estimators
    if(config$nat_inf){
      results$nat_inf[i] <-  vegrowth::do_gcomp_nat_inf(data = data, 
                                                        models = pkg_models,
                                                        Z_name = "Z",
                                                        X_name = "X", 
                                                        exclusion_restriction = TRUE)['additive_effect']
    }
    
    if(config$population){
      results$pop[i] <- vegrowth::do_gcomp_pop(data = data, 
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
      nat_inf = if ("nat_inf" %in% names(sub_df)) mean(sub_df$nat_inf, na.rm = TRUE) else NA,
      pop = if ("pop" %in% names(sub_df)) mean(sub_df$pop, na.rm = TRUE) else NA
    )
    
    # Append to results
    results <- rbind(results, new_row)
  }
  
  # Bootstrap Estimates ----------------------------------------
  
  # 1. Do n_boot bootstrap replicates
  boot_res_list <- replicate(config$n_boot, one_boot(data, config, parameters), simplify = FALSE)
  boot_res_df <- dplyr::bind_rows(boot_res_list, .id = "boot_id")
  
  # 2. Compute SEs, CIs, and rejection indicators per Y_out and estimand
  boot_summary <- boot_res_df %>%
    dplyr::group_by(Y_out) %>%
    dplyr::summarise(
      nat_inf_se = if ("nat_inf" %in% names(.)) sd(nat_inf, na.rm = TRUE) else NA_real_,
      nat_inf_lower = if ("nat_inf" %in% names(.)) quantile(nat_inf, 0.025, na.rm = TRUE) else NA_real_,
      nat_inf_upper = if ("nat_inf" %in% names(.)) quantile(nat_inf, 0.975, na.rm = TRUE) else NA_real_,
      pop_se = if ("pop" %in% names(.)) sd(pop, na.rm = TRUE) else NA_real_,
      pop_lower = if ("pop" %in% names(.)) quantile(pop, 0.025, na.rm = TRUE) else NA_real_,
      pop_upper = if ("pop" %in% names(.)) quantile(pop, 0.975, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  
  # 3. Merge bootstrap summaries with point estimates
  results_full <- dplyr::left_join(results, boot_summary, by = "Y_out")
  
  # 4. Add reject indicator columns
  if (config$nat_inf) {
    results_full <- results_full %>%
      dplyr::mutate(
        nat_inf_reject = (abs(nat_inf - config$null_hypothesis_value) / nat_inf_se) > qnorm(1 - config$alpha_level / 2)
      )
  }
  
  if (config$population) {
    results_full <- results_full %>%
      dplyr::mutate(
        pop_reject = (abs(pop - config$null_hypothesis_value) / pop_se) > qnorm(1 - config$alpha_level / 2)
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
