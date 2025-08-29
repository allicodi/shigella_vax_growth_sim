
options(echo = TRUE)

.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

here::i_am("evaluate_performance.R")

source(here::here("R/evaluation_fns.R"))

# get seed & config settings from bash script
cargs <- commandArgs(TRUE)
setting <- cargs[1]
truth_file <- cargs[2]

# read in config file
config <- config::get(file = "config.yml", config = setting)

# read in truth file 
truth_df <- readRDS(truth_file)

# get list of results files matching pattern
dir <- "/projects/dbenkes/allison/shigella_vaccine_trial/"
pattern <- paste0("results/", setting, "_seed_.*\\.Rds$")
all_files <- list.files(dir, pattern = pattern, full.names = TRUE)

# Load and combine results
results_list <- lapply(valid_files, readRDS)
results <- do.call(rbind, results_list)

bias_df <- data.frame()
coverage_df <- data.frame()
power_df <- data.frame()
prop_neg_df <- data.frame()

for(n in unique(results$n)){
  bias_df <- rbind(bias_df, do.call(cbind, get_bias(results, truth_df, n)))
  coverage_df <- rbind(coverage_df, do.call(cbind, get_coverage(results, truth_df, n)))
  power_df <- rbind(power_df, do.call(cbind, get_power(results, truth_df, n)))
  prop_neg_df <- rbind(prop_neg_df, do.call(cbind, get_neg_pt_est(results, truth_df, n)))
}

final_results <- list(
  results = results,
  truth_df = truth_df,
  bias_df = bias_df,
  coverage_df = coverage_df,
  power_df = power_df,
  prop_neg_df = prop_neg_df
)

saveRDS(final_results, file = paste0(dir, "/results/", setting, "_evaluation_results.Rds"))
