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

setting_grid <- expand.grid(estimand = config$estimand, 
                            estimator = config$estimators, 
                            er = config$exclusion_restriction, 
                            cw = config$cross_world,
                            two_stage = config$two_stage)

# elim any settings that do not exist (ex. where ER & CW both == FALSE, CW + 2 part, ER + CW + 2part)
elim <- which(((setting_grid$estimand == "nat_inf" & setting_grid$cw == FALSE & setting_grid$er == FALSE)) | # nat inf both false
                (setting_grid$estimand == "nat_inf" & setting_grid$cw == TRUE & setting_grid$two_stage == TRUE) | # only 1part for cross world
                (setting_grid$estimand == "pop" & (setting_grid$cw == TRUE | setting_grid$er == TRUE))) # only need to run pop once (ER + CW do not apply; no assumptions)

setting_grid <- setting_grid[-elim,]

# add unadj if applicable
if(config$nat_inf_unadj){
  setting_grid <- rbind(setting_grid, data.frame(estimand = "nat_inf",
                                                 estimator = "unadj",
                                                 er = NA, 
                                                 cw = NA,
                                                 two_stage = NA))
}

results <- lapply(config$n_sample_size, function(n){
  
  # Simulate data
  data <- simulate_data(parameters = parameters, 
                        n = n, 
                        VE_mild = config$VE_mild, 
                        VE_severe = config$VE_severe, 
                        seed = seed, 
                        type = "observed")
  
  # Long term & population effect estimation -----------------------------
  
  res_df <- data.frame()
  
  # 1. Fit Models
  
  # Get unique timepoints needed in config$intervals
  Y_out <- unique(do.call(c, config$intervals))

  # Get effect for all individual Y_outs
  for(i in 1:length(Y_out)){
    
    # Name of outcome variable
    Y_name <- paste0("Y_", Y_out[i])
    
    # Fit models
    if(nrow(setting_grid > 0)){
      pkg_models <- vegrowth::fit_models(data = data,
                                         Y_name = Y_name, 
                                         Z_name = "Z", 
                                         X_name = "X", 
                                         S_name = "S_inf", 
                                         estimand = config$estimand, 
                                         method = config$estimators, 
                                         exclusion_restriction = TRUE, 
                                         family = "gaussian")
    }
    
    for(j in 1:nrow(setting_grid)){
      setting <- setting_grid[j,]
      
      if(setting$estimand == "nat_inf"){
        res <- est_nat_inf(data = data,
                           estimator = setting$estimator,
                           pkg_models = pkg_models,
                           Y_name = Y_name, 
                           exclusion_restriction = setting$er,
                           cross_world = setting$cw,
                           two_part_model = setting$two_stage)
        
        
      } else{
        res <- est_pop(data = data,
                       estimator = setting$estimator,
                       pkg_models = pkg_models,
                       Y_name = Y_name, 
                       two_part_model = setting$two_stage)
        
      }
      
      row <- data.frame(Y_out = Y_name, estimate = res, setting)
      res_df <- rbind(res_df, row)
      
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
    for(j in 1:nrow(setting_grid)){
      setting <- setting_grid[j,]
      
      if(setting$estimator != "unadj"){
        sub_df <- res_df[res_df$Y_out %in% Y_names &
                           res_df$estimator == setting$estimator &
                           res_df$estimand == setting$estimand &
                           res_df$er == setting$er &
                           res_df$cw == setting$cw &
                           res_df$two_stage == setting$two_stage,]
      } else{
        sub_df <- res_df[res_df$Y_out %in% Y_names &
                           res_df$estimator == setting$estimator &
                           res_df$estimand == setting$estimand,]
      }
      
      row <- data.frame(Y_out = avg_name, 
                        estimate = mean(sub_df$estimate),
                        setting)
      
      res_df <- rbind(res_df, row)
    }
   
  }
  
  # Bootstrap Estimates ----------------------------------------
  
  # 1. Do n_boot bootstrap replicates
  boot_res_list <- replicate(config$n_boot, one_boot(data, config, setting_grid, parameters), simplify = FALSE)
  boot_res_df <- dplyr::bind_rows(boot_res_list, .id = "boot_id")
  
  # 2. Compute SEs, CIs, and rejection indicators per Y_out and estimand
  boot_summary <- boot_res_df %>%
    dplyr::group_by(Y_out, estimand, estimator, er, cw, two_stage) %>%
    dplyr::summarise(
      se = sd(estimate, na.rm = TRUE),
      lower_ci = quantile(estimate, 0.025, na.rm = TRUE),
      upper_ci = quantile(estimate, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 3. Merge bootstrap summaries with point estimates
  results_full <- dplyr::left_join(res_df, boot_summary, by = c("Y_out", "estimand", "estimator", "er", "cw", "two_stage"))
  
  # 4. Add reject indicator columns
  results_full <- results_full %>%
    dplyr::mutate(
      reject = (abs(estimate - config$null_hypothesis_value) / se) > qnorm(1 - config$alpha_level / 2)
    )
  
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
