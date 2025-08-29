# ------------------------------------------------------------------------------
# Script to run analysis for given configuration settings 
# ------------------------------------------------------------------------------

.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))
#.libPaths("/apps/R/4.4.0/lib64/R/site/library")
#.libPaths(c("/home/acodi/Rlibs", "/apps/R/4.4.0/lib64/R/site/library", .libPaths()))
#.libPaths(c("~/Rlibs", .libPaths()))

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

config <- config::get(file = "config.yml", config = setting)

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
  if(config$long_term | config$population){
    estimand <- c()
    if(config$long_term){
      estimand <- c(estimand, "nat_inf")
    } 
    
    if(config$population){
      estimand <- c(estimand, "pop")
    }
    
    pkg_models <- vegrowth::fit_models(data = data,
                                       Y_name = "Y_12",
                                       Z_name = "Z", 
                                       X_name = "X", 
                                       S_name = "S_inf", 
                                       estimand = estimand, 
                                       method = "gcomp",
                                       family = "gaussian")
  } else{
    pkg_models <- NULL
  }
  
  # 2. Call estimation functions from package
  if(config$long_term){
    est_long_term <- vegrowth::do_gcomp_nat_inf(data = data, models = pkg_models)['additive_effect']
  } else{
    est_long_term <- NULL
  }
  
  if(config$population){
    est_pop <- vegrowth::do_gcomp_pop(data = data, models = pkg_models, Z_name = "Z", X_name = "X")['additive_effect']
  } else{
    est_pop <- NULL
  }
  
  # Short-term effect estimation -----------------------------
  
  # 1. Call effect estimation function (models fit within)
  if(config$short_term){
    est_short_term <- estimate_short_term(data = data, 
                                          parameters = paramaters, 
                                          V_u_months = as.numeric(config$V_u_months), 
                                          V_u_week_interval = as.numeric(config$V_u_week_interval), 
                                          pool = config$pool, 
                                          I_S__X_Z0_0_6_formula = config$I_S__X_Z0_0_6_formula,
                                          I_S__X_Z0_6_12_formula = config$I_S__X_Z0_6_12_formula, 
                                          Y_out__Z1_Suminus10_X_formula = config$Y_out__Z1_Suminus10_X_formula,
                                          Y_out__Z0_Suminus11_X_T_formula = config$Y_out__Z0_Suminus11_X_T_formula, 
                                          Y_out__Z0_Suminus11_X_T_V_formula = config$Y_out__Z0_Suminus11_X_T_V_formula)
  } else{
    est_short_term <- NULL
  }
  
  # Bootstrap Estimates ----------------------------------------
  
  # 1. Do n_boot bootstrap replicates
  
  # boot_res <- vector("list", config$n_boot)
  # 
  # for (i in seq_len(config$n_boot)) {
  #   message("Running bootstrap replicate ", i, " of ", config$n_boot)
  #   boot_res[[i]] <- tryCatch(
  #     one_boot(data, config, parameters),
  #     error = function(e) {
  #       message("Error in bootstrap replicate ", i, ": ", e$message)
  #       return(NULL)  # Or NA
  #     }
  #   )
  # }
  
  boot_res <- replicate(config$n_boot, one_boot(data, config, parameters))
  boot_res_df <- as.data.frame(t(boot_res))
  
  # ^^ if any are NA, should I repeat?
  any_NA <- any(is.na(boot_res_df$est_short_term))
  attempt <- 1
  
  # return sum(any_NA)
  while(any_NA & attempt <= 5){
    idx <- which(is.na(boot_res_df$est_short_term))
    
    # assuming won't happen again
    boot_res_2 <- replicate(length(idx), one_boot(data, config, parameters))
    boot_res_df_2 <- as.data.frame(t(boot_res_2))
    
    boot_res_df[idx,] <- boot_res_df_2
    
    any_NA <- any(is.na(boot_res_df$est_short_term))
    attempt <- attempt + 1
    
  }
  
  # 2. Get bootstrap CIs & hypothesis test
  if(config$long_term){
    long_term_se <- sd(boot_res_df$est_long_term)
    long_term_lower_ci <- quantile(boot_res_df$est_long_term, p = 0.025)
    long_term_upper_ci <- quantile(boot_res_df$est_long_term, p = 0.975)
    
    long_term_reject <- (abs(est_long_term - config$null_hypothesis_value) / long_term_se) > qnorm(1 - config$alpha_level/2)
    
  } else{
    long_term_se <- NULL
    long_term_lower_ci <- NULL
    long_term_upper_ci <- NULL
    
    long_term_reject <- NULL
  }
  
  if(config$pop){
    pop_se <- sd(boot_res_df$est_pop)
    pop_lower_ci <- quantile(boot_res_df$est_pop, p = 0.025)
    pop_upper_ci <- quantile(boot_res_df$est_pop, p = 0.975)
    
    pop_reject <- (abs(est_pop - config$null_hypothesis_value) / pop_se) > qnorm(1 - config$alpha_level/2)
  } else{
    pop_se <- NULL
    pop_lower_ci <- NULL
    pop_upper_ci <- NULL
  }
  
  if(config$short_term){
    short_term_se <- sd(boot_res_df$est_short_term)
    short_term_lower_ci <- quantile(boot_res_df$est_short_term, p = 0.025)
    short_term_upper_ci <- quantile(boot_res_df$est_short_term, p = 0.975)
    
    short_term_reject <- (abs(est_short_term - config$null_hypothesis_value) / short_term_se) > qnorm(1 - config$alpha_level/2)
  } else{
    short_term_se <- NULL
    short_term_lower_ci <- NULL
    short_term_upper_ci <- NULL
  }
  
  result <- list(seed = seed, 
                 n = n, 
                 est_long_term = est_long_term, 
                 long_term_se = long_term_se,
                 long_term_lower_ci = long_term_lower_ci, 
                 long_term_upper_ci = long_term_upper_ci, 
                 long_term_reject = long_term_reject, 
                 est_pop = est_pop,
                 pop_se = pop_se,
                 pop_lower_ci = pop_lower_ci, 
                 pop_upper_ci = pop_upper_ci, 
                 pop_reject = pop_reject,
                 est_short_term = est_short_term,
                 short_term_se = short_term_se,
                 short_term_lower_ci = short_term_lower_ci, 
                 short_term_upper_ci = short_term_upper_ci,
                 short_term_reject = short_term_reject)
  
  # save list incrementally too so if gets killed don't need to start from scratch
  saveRDS(results, paste0("/projects/dbenkes/allison/shigella_vaccine_trial/results/", setting, "_n_", n, "_seed_", seed, ".Rds"))

  return(result)
  
})

results <- as.data.frame(do.call(rbind, results))

saveRDS(results, paste0("/projects/dbenkes/allison/shigella_vaccine_trial/results/", setting, "_seed_", seed, ".Rds"))
