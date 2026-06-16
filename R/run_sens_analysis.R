# ------------------------------------------------------------------------------
# Script to run analysis for given configuration settings 
# ------------------------------------------------------------------------------

.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

options(echo = TRUE)

here::i_am("R/run_sens_analysis.R")

library(dplyr)

source(here::here("R/helpers/parameter_generation_fns.R"))
source(here::here("R/helpers/simulate_data.R"))
source(here::here("R/helpers/estimation_fn.R"))
source(here::here("R/helpers/bootstrap.R"))

# get seed & config settings from bash script
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
setting <- Sys.getenv("SETTING")
parameters <- readRDS(Sys.getenv("PARAMETERS_FILE"))

# set seed
set.seed(seed)

# ALTERNATIVE
# get setting from bash script
# cargs <- commandArgs(TRUE)
# setting <- cargs[1]4
# parameter_dir <- cargs[2]

cfg <- yaml::read_yaml("config_unmeas_conf.yml")
config <- cfg[[setting]]

setting_grid <- expand.grid(estimand = config$estimand, 
                            estimator = config$estimators, 
                            er = config$exclusion_restriction, 
                            cw = config$cross_world,
                            two_stage = config$two_stage,
                            beta_U_S = config$beta_U_S,
                            beta_Y_Y = config$beta_U_Y,
                            p_U_is_1 = config$p_U_is_1)

# elim any settings that do not exist (ex. where ER & CW both == FALSE, CW + 2 part, ER + CW + 2part)
elim <- which(((setting_grid$estimand == "nat_inf" & setting_grid$cw == FALSE & setting_grid$er == FALSE)) | # nat inf both false
                (setting_grid$estimand == "nat_inf" & setting_grid$cw == TRUE & setting_grid$two_stage == TRUE)  ) #| # only 1part for cross world
#(setting_grid$estimand == "pop" & (setting_grid$cw == TRUE | setting_grid$er == TRUE))) # only need to run pop once (ER + CW do not apply; no assumptions)

if(length(elim) > 0){
  setting_grid <- setting_grid[-elim,]
}

# if pop in there multiple times for given estimator, eliminate (ER + CW do not apply; no assumptions)
setting_grid <- setting_grid %>%
  group_by(estimand, estimator) %>%
  # keep all non-pop rows; for pop keep only one per estimator
  filter(
    estimand != "pop" |
      row_number() == 1
  ) %>%
  ungroup() %>%
  # for the remaining pop rows, set er and cw to FALSE (no assumptions)
  mutate(
    er = if_else(estimand == "pop", FALSE, er),
    cw = if_else(estimand == "pop", FALSE, cw)
  )

# add unadj if applicable
if(config$nat_inf_unadj){
  setting_grid <- rbind(setting_grid, data.frame(estimand = "nat_inf",
                                                 estimator = "unadj",
                                                 er = NA, 
                                                 cw = NA,
                                                 two_stage = NA,
                                                 beta_U_S = config$beta_U_S,
                                                 beta_Y_Y = config$beta_U_Y,
                                                 p_U_is_1 = config$p_U_is_1))
}

results <- lapply(config$n_sample_size, function(n){
  
  # Simulate data
  data <- simulate_data(parameters = parameters, 
                        n = n, 
                        VE_mild = config$VE_mild, 
                        VE_severe = config$VE_severe, 
                        seed = seed, 
                        type = "observed",
                        beta_U_S = config$beta_U_S,
                        beta_Y_Y = config$beta_U_Y,
                        p_U_is_1 = config$p_U_is_1)
  
  # if dropout, remove dropout% of observations 
  if(config$dropout > 0){
    n_drop <- ceiling(nrow(data) * config$dropout)
    which_drop <- sample(1:nrow(data), n_drop, replace = FALSE)
    data <- data[-which_drop,]
  }
  
  # Long term & population effect estimation -----------------------------
  
  #res_df <- data.frame()
  
  # 1. Fit Models
  
  # Get unique timepoints needed in config$intervals
  Y_out <- unique(do.call(c, config$intervals))
  
  # list for if matrix, df for ease later on
  res_list <- vector("list", length = length(Y_out))
  res_df <- data.frame()
  
  # Get effect for all individual Y_outs
  for(i in 1:length(Y_out)){
    
    # Name of outcome variable
    Y_name <- paste0("Y_", Y_out[i])
    
    # Element in list for each outcome variable, those elements are lists with one item for each setting
    names(res_list)[i] <- Y_name
    res_list[[Y_name]] <- vector("list", length = nrow(setting_grid))
    
    # Fit models
    if(nrow(setting_grid) > 0){
      pkg_models <- vaxstrat::fit_models(data = data,
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
      
      # add settings to element in list
      res_list[[Y_name]][[j]]$setting <- setting
      
      if(setting$estimand == "nat_inf"){
        
        res <- est_nat_inf(data = data,
                           estimator = setting$estimator,
                           pkg_models = pkg_models,
                           Y_name = Y_name, 
                           exclusion_restriction = setting$er,
                           cross_world = setting$cw,
                           two_part_model = setting$two_stage)
        
        res_list[[Y_name]][[j]]$res <- res
        
      } else{
        res <- est_pop(data = data,
                       estimator = setting$estimator,
                       pkg_models = pkg_models,
                       Y_name = Y_name, 
                       two_part_model = setting$two_stage)
        
        res_list[[Y_name]][[j]]$res <- res
        
      }
      
      row <- data.frame(Y_out = Y_name, estimate = res$additive_effect, se = res$additive_se, setting)
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
      
      # use influence functions to get standard error
      if(setting$estimator == "aipw"){
        
        if_matrix <- matrix(data = NA, ncol = 0, nrow = nrow(data))
        
        for(Y_name in Y_names){
          Y_name_if_matrix <- res_list[[Y_name]][[j]]$res$if_matrix
          colnames(Y_name_if_matrix) <- paste0(Y_name, "_", colnames(Y_name_if_matrix))
          if_matrix <- cbind(if_matrix, Y_name_if_matrix)
        }
        
        # denominator == number of timpoints being averaged. 
        # ex. Y_6_12 denom = 2, so gradient is 1/2, -1/2, 1/2, -1/2 to average effects from Y_6 and Y_12
        denom <- length(Y_names)
        gradient <- matrix(rep(c(1/denom, -1/denom), length(Y_names)), ncol = 1)
        
        cov_matrix <- cov(if_matrix) / nrow(data)
        se_interval <- sqrt(t(gradient) %*% cov_matrix %*% gradient)
        
      } else{
        # otherwise get with bootstrap later
        se_interval <- NA
      }
      
      row <- data.frame(Y_out = avg_name, 
                        estimate = mean(sub_df$estimate),
                        se = se_interval,
                        setting)
      
      res_df <- rbind(res_df, row)
    }
    
  }
  
  # 1. No bootstrap, just fill in CIs
  results_full <- res_df %>%
    mutate(
      # fill in lower_ci and upper_ci
      lower_ci = estimate - 1.96 * se,
      upper_ci = estimate + 1.96 * se
    )
  
  # 2. Add reject indicator columns
  results_full <- results_full %>%
    dplyr::mutate(
      reject = (abs(estimate - config$null_hypothesis_value) / se) > qnorm(1 - config$alpha_level / 2)
    )
  
  # 3. Final combined result
  result <- data.frame(
    seed = seed,
    n = n,
    results_full
  )
  
  return(result)
  
})

results <- as.data.frame(do.call(rbind, results))

saveRDS(results, paste0("/projects/dbenkes/allison/shigella_vaccine_trial/sens_results/", setting, "_seed_", seed, ".Rds"))
