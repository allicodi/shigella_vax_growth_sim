
here::i_am("get_true_values.R")

source(here::here("R/01_get_truth.R"))

# get setting from bash script
cargs <- commandArgs(TRUE)
setting <- cargs[1]

config <- config::get(file = "config.yml", config = setting)
if(is.null(config$zhifei)) config$zhifei <- FALSE
if(is.null(config$catch_up)) config$catch_up <- FALSE

# Get combinations of parameters 
param_combo_df <- expand.grid(incidence_shigella = config$incidence_shigella,
                              incidence_severe_shigella = config$incidence_severe_shigella,
                              VE_mild = config$VE_mild,
                              VE_severe = config$VE_severe,
                              lfaz_bl_effect_shig = config$lfaz_bl_effect_shig,
                              lfaz_bl_effect_severe_shig = config$lfaz_bl_effect_severe_shig,
                              shig_coef_mild = config$shig_coef_mild,
                              shig_coef_severe = config$shig_coef_severe,
                              sd = config$sd,
                              catch_up = config$catch_up)

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
                      zhifei = config$zhifei,
                      catch_up = config$catch_up)
  
  # truths <- get_truth_quad(incidence_shigella = row$incidence_shigella,
  #                          incidence_severe_shigella = row$incidence_severe_shigella,
  #                          VE_mild = row$VE_mild,
  #                          VE_severe = row$VE_severe,
  #                          lfaz_bl_effect_shig = row$lfaz_bl_effect_shig,
  #                          lfaz_bl_effect_severe_shig = row$lfaz_bl_effect_severe_shig,
  #                          shig_coef_mild = row$shig_coef_mild,
  #                          shig_coef_severe = row$shig_coef_severe,
  #                          sd_growth = row$sd,
  #                          zhifei = config$zhifei,
  #                          catch_up = config$catch_up)
  
  true_growth_effect[i] <- truths$true_ge
  true_pop_growth_effect[i] <- truths$true_pop_ge
  
}

# add coefficients corresponding to desired efficacy to parameter grid
param_combo_df$true_growth_effect <- true_growth_effect
param_combo_df$true_pop_growth_effect <- true_pop_growth_effect

# Save param_combo_df to file to be loaded by run_simulation_cluster.R
#saveRDS(param_combo_df, file = paste0("/projects/dbenkes/allison/shigella/.temp_results/quad_precomputed_values_add_growth_",setting,".Rds"))
saveRDS(param_combo_df, file = paste0("/projects/dbenkes/allison/shigella/.temp_results/precomputed_values_add_growth_",setting,".Rds"))
