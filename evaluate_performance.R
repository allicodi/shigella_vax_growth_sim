here::i_am("evaluate_performance.R")

source(here::here("R/03_evaluate_performance.R"))

# get seed & config settings from bash script
cargs <- commandArgs(TRUE)
setting <- cargs[1]

# read in config file
config <- config::get(file = "config.yml", config = setting)
est <- config$est

# Read in results by seed and aggregate into one dataframe
dir <- "/projects/dbenkes/allison/shigella/.temp_results"
pattern <- paste0("test_dr_result_seed_.*_setting_", setting, "\\.rds$")
all_files <- list.files(dir, pattern = pattern, full.names = TRUE)

# Filter files by date (12/11/2024 or later)
date_threshold <- as.POSIXct("2024-12-11")
file_info <- file.info(all_files)
valid_files <- rownames(file_info[file_info$mtime >= date_threshold, ])

# Load and combine results
results_list <- lapply(valid_files, readRDS)
results <- do.call(rbind, results_list)

# read in truth + param combos
truth_df <- readRDS(paste0(dir, "/quad_precomputed_values_add_growth_",setting,".Rds"))
truth_df <- merge(truth_df, expand.grid(n = config$n_sample_size))

bias_df <- data.frame()
coverage_df <- data.frame()
power_df <- data.frame()
compare_df <- data.frame()
prop_neg_df <- data.frame()
ci_width_df <- data.frame()

for(i in 1:nrow(truth_df)){
  bias_df <- rbind(bias_df, do.call(cbind, get_bias(results, truth_df[i,], est)))
  coverage_df <- rbind(coverage_df, do.call(cbind, get_coverage(results, truth_df[i,], est)))
  power_df <- rbind(power_df, do.call(cbind, get_power(results, truth_df[i,], est)))
  compare_df <- rbind(compare_df, do.call(cbind, get_se_compare(results, truth_df[i,], est)))
  prop_neg_df <- rbind(prop_neg_df, do.call(cbind, get_neg_pt_est(results, truth_df[i,], est)))
  ci_width_df <- rbind(ci_width_df, do.call(cbind, get_ci_width(results, truth_df[i,], est)))
}

bias_df <- cbind(truth_df, bias_df)
coverage_df <- cbind(truth_df, coverage_df)
power_df <- cbind(truth_df, power_df)
compare_df <- cbind(truth_df, compare_df)
prop_neg_df <- cbind(truth_df, prop_neg_df)
ci_width_df <- cbind(truth_df, ci_width_df)

# Save the final results
final_results <- list(
  results = results,
  truth_df = truth_df,
  bias_df = bias_df,
  coverage_df = coverage_df,
  power_df = power_df,
  compare_df = compare_df,
  prop_neg_df = prop_neg_df,
  ci_width_df = ci_width_df
)

saveRDS(final_results, file = paste0("/projects/dbenkes/allison/shigella/results/test_dr_", setting, "_results.Rds"))
saveRDS(final_results, file = paste0("results/test_dr_", setting, "_results.Rds"))

