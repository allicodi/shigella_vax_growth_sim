# ---------------------------------------------------------------------------------------------
# Function to get truth for given config settings
# ---------------------------------------------------------------------------------------------

.libPaths(c("/apps/R/4.4.0/lib64/R/site/library","/apps/R/4.4.0/lib64/R/library", "~/Rlibs_ve_trial"))

here::i_am("get_truth.R")

source(here::here("R/parameter_generation_fns.R"))
source(here::here("R/simulate_data.R"))
source(here::here("R/truth_fns.R"))

# get setting from bash script
cargs <- commandArgs(TRUE)
setting <- cargs[1]
parameters_file <- cargs[2]
truth_dir <- cargs[3]

cfg <- yaml::read_yaml("config.yml")
config <- cfg[[setting]]

parameters <- readRDS(parameters_file)
data <- simulate_data(parameters = parameters,
                      n = 1e7, 
                      VE_mild = config$VE_mild, 
                      VE_severe = config$VE_severe, 
                      seed = 12345,
                      type = "counterfactual")

truth <- list()

# Get unique timepoints needed in config$intervals
Y_out <- unique(do.call(c, config$intervals))

# Get truth for all individual times
for(Y in Y_out){
  if(config$nat_inf){
    
    if(Y == 3){
      truth$nat_inf_Y_3 <- truth_nat_inf_3mo(data = data)
    } else if (Y == 6){
      truth$nat_inf_Y_6 <- truth_nat_inf_6mo(data = data)
    } else if (Y == 9){
      truth$nat_inf_Y_9 <- truth_nat_inf_9mo(data = data)
    } else{
      truth$nat_inf_Y_12 <- truth_nat_inf_12mo(data = data)
    }
    
  } 
  
  if(config$pop){
    
    if(Y == 3){
      truth$pop_Y_3 <- truth_pop_3mo(data = data)
    } else if (Y == 6){
      truth$pop_Y_6 <- truth_pop_6mo(data = data)
    } else if (Y == 9){
      truth$pop_Y_9 <- truth_pop_9mo(data = data)
    } else{
      truth$pop_Y_12 <- truth_pop_12mo(data = data)
    }
    
  }
}

# Get effect for averaged Y_outs
which_intervals <- do.call(c, lapply(config$intervals, function(x) length(x) > 1))
avg_intervals <- config$intervals[which_intervals]

for (i in avg_intervals) {
  # Create a name for the averaged interval (e.g., "Y_6_9")
  avg_name <- paste0("Y_", paste(i, collapse = "_"))
  
  # Construct the names of the individual timepoint entries
  Y_names <- paste0("Y_", i)
  
  # For each type of truth that applies
  if (config$nat_inf) {
    # Collect available nat_inf_Y_* entries for this interval
    nat_vals <- unlist(truth[paste0("nat_inf_", Y_names)])
    if (length(nat_vals) > 0) {
      truth[[paste0("nat_inf_", avg_name)]] <- mean(nat_vals)
    }
  }
  
  if (config$pop) {
    # Collect available pop_Y_* entries for this interval
    pop_vals <- unlist(truth[paste0("pop_", Y_names)])
    if (length(pop_vals) > 0) {
      truth[[paste0("pop_", avg_name)]] <- mean(pop_vals)
    }
  }
}

saveRDS(truth, paste0(truth_dir, "/truth_", setting, ".Rds"))
