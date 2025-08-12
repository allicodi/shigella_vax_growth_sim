# ------------------------------------------------------------------------------
# Script to run analysis for given configuration settings 
# ------------------------------------------------------------------------------

here::i_am("run_analysis.R")

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
                                          V_u_months = config$V_u_months, 
                                          V_u_week_interval = config$V_u_week_interval, 
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
  boot_res <- replicate(config$n_boot, one_boot(data, config))
  boot_res <- data.frame(do.call(rbind, boot_res))
  
  # 2. Get bootstrap CIs & hypothesis test
  if(config$long_term){
    long_term_se <- sd(boot_res$est_long_term)
    long_term_lower_ci <- quantile(boot_res$est_long_term, p = 0.025)
    long_term_upper_ci <- quantile(boot_res$est_long_term, p = 0.975)
    
    long_term_reject <- (abs(est_long_term - null_hypothesis_value) / long_term_se) > qnorm(1 - alpha_level/2)
    
  } else{
    long_term_se <- NULL
    long_term_lower_ci <- NULL
    long_term_upper_ci <- NULL
    
    long_term_reject <- NULL
  }
  
  if(config$pop){
    pop_se <- sd(boot_res$est_pop)
    pop_lower_ci <- quantile(boot_res$est_pop, p = 0.025)
    pop_upper_ci <- quantile(boot_res$est_pop, p = 0.975)
    
    pop_reject <- (abs(est_pop - null_hypothesis_value) / pop_se) > qnorm(1 - alpha_level/2)
  } else{
    pop_se <- NULL
    pop_lower_ci <- NULL
    pop_upper_ci <- NULL
  }
  
  if(config$short_term){
    short_term_se <- sd(boot_res$est_short_term)
    short_term_lower_ci <- quantile(boot_res$est_short_term, p = 0.025)
    short_term_upper_ci <- quantile(boot_res$est_short_term, p = 0.975)
    
    short_term_reject <- (abs(est_short_term - null_hypothesis_value) / short_term_se) > qnorm(1 - alpha_level/2)
  } else{
    short_term_se <- NULL
    short_term_lower_ci <- NULL
    short_term_upper_ci <- NULL
  }
  
  return(list(n = n, 
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
              short_term_reject = short_term_reject))
  
})

results <- as.data.frame(do.call(rbind, results))
results$seed <- seed

saveRDS(results, "/projects/dbenkes/allison/shigella_vaccine_trial/results/", setting, "_seed_", seed, ".Rds")
