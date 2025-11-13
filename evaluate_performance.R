
options(echo = TRUE)

.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

library(dplyr)

here::i_am("evaluate_performance.R")

source(here::here("R/evaluation_fns.R"))

# get seed & config settings from bash script
cargs <- commandArgs(TRUE)
setting <- cargs[1]
truth_file <- cargs[2]

# read in config file
cfg <- yaml::read_yaml("config.yml")
config <- cfg[[setting]]

# read in truth file 
truth_df <- readRDS(truth_file)

# get list of results files matching pattern
dir <- "/projects/dbenkes/allison/shigella_vaccine_trial/results"
pattern <- paste0(setting, "_seed_.*\\.Rds$")
all_files <- list.files(dir, pattern = pattern, full.names = TRUE)

results_list <- lapply(all_files, readRDS)
results <- bind_rows(results_list)

bias_df     <- get_bias(results, truth_df, config)
coverage_df <- get_coverage(results, truth_df, config)
power_df    <- get_power(results, truth_df, config)
prop_neg_df <-  get_neg_pt_est(results, truth_df, config)

final_results <- list(
  results = results,
  truth_df = truth_df,
  bias_df = bias_df,
  coverage_df = coverage_df,
  power_df = power_df,
  prop_neg_df = prop_neg_df
)

saveRDS(final_results, file = paste0(dir, "/", setting, "_evaluation_results.Rds"))
