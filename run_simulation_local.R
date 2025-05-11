
here::i_am("run_simulation_local.R")

devtools::load_all("../packages/vegrowth/")

source(here::here("R/00_simulate_data.R"))
source(here::here("R/01_get_truth.R"))
source(here::here("R/02_do_one_ci_sim.R"))
source(here::here("R/03_evaluate_performance.R"))

library(parallel)
library(ggplot2)
library(dplyr)

run_simulation <- function(n_seeds = 1000,
                           n_boot = 1000,
                           n_sample_size = 4000,
                           incidence_shigella = 0.0618,
                           incidence_severe_shigella = 0.0240,
                           VE_mild = 0.4,
                           VE_severe = 0.6,
                           lfaz_bl_effect_shig = -0.3147,
                           lfaz_bl_effect_severe_shig = -0.1990,
                           shig_coef_mild = -0.052,
                           shig_coef_severe = -0.072,
                           sd_growth = "prev-1",
                           est = c("gcomp_pop_estimand", "gcomp",
                                   "efficient_gcomp", "efficient_aipw", "efficient_tmle",
                                   "choplump", "hudgens_lower", "hudgens_upper",
                                   "hudgens_lower_doomed", "hudgens_upper_doomed",
                                   "hudgens_adj_lower", "hudgens_adj_upper"),
                           zhifei = FALSE,
                           catch_up = TRUE,
                           G_X_model = "G ~ X",
                           Yinf_X_model = "Y_inf ~ X"){
 
  # Get combinations of parameters to iterate over
  param_combo_df <- expand.grid(incidence_shigella = incidence_shigella,
                                incidence_severe_shigella = incidence_severe_shigella,
                                VE_mild = VE_mild,
                                VE_severe = VE_severe,
                                lfaz_bl_effect_shig = lfaz_bl_effect_shig,
                                lfaz_bl_effect_severe_shig = lfaz_bl_effect_severe_shig,
                                shig_coef_mild = shig_coef_mild,
                                shig_coef_severe = shig_coef_severe,
                                sd = sd_growth,
                                catch_up = catch_up,
                                zhifei = zhifei, 
                                G_X_model = G_X_model,
                                Yinf_X_model = Yinf_X_model)
  
  param_combo_df$G_X_model <- as.character(param_combo_df$G_X_model)
  param_combo_df$Yinf_X_model <- as.character(param_combo_df$Yinf_X_model)
  
  # Vectors to hold truth
  true_growth_effect <- vector(mode = "numeric", length = nrow(param_combo_df))
  true_pop_growth_effect <- vector(mode = "numeric", length = nrow(param_combo_df))
  
  # For each combination, find coefficient needed to achieve desired vaccine efficacy
  # Also find true value for each combo
  for(i in 1:nrow(param_combo_df)){
    row <- param_combo_df[i,]
    
    # Get truth
    truths <- get_truth(incidence_shigella = row$incidence_shigella,
                        incidence_severe_shigella = row$incidence_severe_shigella,
                        VE_mild = row$VE_mild,
                        VE_severe = row$VE_severe,
                        lfaz_bl_effect_shig = row$lfaz_bl_effect_shig,
                        lfaz_bl_effect_severe_shig = row$lfaz_bl_effect_severe_shig,
                        shig_coef_mild = row$shig_coef_mild,
                        shig_coef_severe = row$shig_coef_severe,
                        sd_growth = row$sd,
                        zhifei = row$zhifei,
                        catch_up = row$catch_up)
    
    true_growth_effect[i] <- truths$true_ge
    true_pop_growth_effect[i] <- truths$true_pop_ge
    
  }
  
  # add truths to parameter grid
  param_combo_df$true_growth_effect <- true_growth_effect
  param_combo_df$true_pop_growth_effect <- true_pop_growth_effect
  
  # add seed and sample size to grid 
  param_grid <- expand.grid(
    seed = 1:n_seeds,
    n = n_sample_size
  )
  
  param_grid <- merge(param_grid, param_combo_df)
  
  # if/when this goes to cluster, change so these are all separate jobs?
  results <- mclapply(1:nrow(param_grid), function(x, n_boot, est){
    params <- param_grid[x,]
    do_one_ci_sim(params, n_boot, est)
  }, n_boot = n_boot, est = est, mc.cores = 5)
  
  #results <- do.call(rbind, results)
  #rownames(results) <- 1:nrow(results)
  results <- do.call(rbind, results)
  
  # Drop seed to get datafame of unique n / parameter combinations
  evaluation_df <- unique(param_grid[,!(names(param_grid) %in% c("seed"))])
  
  bias_df <- data.frame()
  coverage_df <- data.frame()
  power_df <- data.frame()
  compare_df <- data.frame()
  prop_neg_df <- data.frame()
  ci_width_df <- data.frame()
  
  for(i in 1:nrow(evaluation_df)){
    # bias_df <- rbind(bias_df, unlist(get_bias(results, evaluation_df[i,], est)))
    # coverage_df <- rbind(coverage_df, unlist(get_coverage(results, evaluation_df[i,], est)))
    # power_df <- rbind(power_df, unlist(get_power(results, evaluation_df[i,], est)))
    # compare_df <- rbind(compare_df, unlist(get_se_compare(results, evaluation_df[i,], est)))
    # prop_neg_df <- rbind(prop_neg_df, unlist(get_neg_pt_est(results, evaluation_df[i,], est)))
    # ci_width_df <- rbind(ci_width_df, unlist(get_ci_width(results, evaluation_df[i,], est)))
    
    bias_df <- rbind(bias_df, do.call(cbind, get_bias(results, evaluation_df[i,], est)))
    coverage_df <- rbind(coverage_df, do.call(cbind, get_coverage(results, evaluation_df[i,], est)))
    power_df <- rbind(power_df, do.call(cbind, get_power(results, evaluation_df[i,], est)))
    compare_df <- rbind(compare_df, do.call(cbind, get_se_compare(results, evaluation_df[i,], est)))
    prop_neg_df <- rbind(prop_neg_df, do.call(cbind, get_neg_pt_est(results, evaluation_df[i,], est)))
    ci_width_df <- rbind(ci_width_df, do.call(cbind, get_ci_width(results, evaluation_df[i,], est)))
  }
  
  # colnames(bias_df) <- paste0("bias_", est)
  # colnames(coverage_df) <- paste0("coverage_", est)
  # colnames(power_df) <- paste0("power_", est)
  # colnames(compare_df) <- paste0("se_ratio_", est)
  # colnames(prop_neg_df) <- paste0("prop_neg_", est)
  # colnames(ci_width_df) <- paste0("ci_width_", est)
  
  bias_df <- cbind(evaluation_df, bias_df)
  coverage_df <- cbind(evaluation_df, coverage_df)
  power_df <- cbind(evaluation_df, power_df)
  compare_df <- cbind(evaluation_df, compare_df)
  prop_neg_df <- cbind(evaluation_df, prop_neg_df)
  ci_width_df <- cbind(evaluation_df, ci_width_df)
  
  return(list(results = results,
              evaluation_df = evaluation_df,
              bias_df = bias_df,
              coverage_df = coverage_df,
              power_df = power_df,
              compare_df = compare_df,
              ci_width_df = ci_width_df)) 
  
}

# ------------------------------------------------------------------------
# debug(run_simulation)
test <- run_simulation(n_seeds = 1000,
                       n_boot = 1000,
                       n_sample_size = c(4000, 8000, 10000, 25000, 50000),
                       incidence_shigella = 0.0618,
                       incidence_severe_shigella = 0.0246,
                       VE_mild = 0.4,
                       VE_severe = 0.6,
                       lfaz_bl_effect_shig= -0.3147,
                       lfaz_bl_effect_severe_shig= -0.1990,
                       shig_coef_mild= -0.052,
                       shig_coef_severe= -0.072,
                       sd = "baseline",
                       G_X_model = "G ~ X",
                       Yinf_X_model = c("Y_inf ~ X"),
                       est = c("gcomp", 
                               "efficient_aipw", "efficient_tmle",
                               "hudgens_lower", "hudgens_upper", 
                               "hudgens_lower_doomed", "hudgens_upper_doomed"))

# TODO
# Check plotting functions
# Make cluster version + run
