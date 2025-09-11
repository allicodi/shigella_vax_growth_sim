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

cfg <- yaml::read_yaml("config_contour.yml")
config <- cfg[[setting]]

parameters <- readRDS(parameters_file)
data <- simulate_data(parameters = parameters,
                      n = 1e7, 
                      VE_mild = config$VE_mild, 
                      VE_severe = config$VE_severe, 
                      seed = 12345,
                      type = "counterfactual")

truth <- list()

if(config$long_term){
  truth$long_term <- truth_nat_inf_12mo(data = data)
} 

if(config$population){
  truth$population <- truth_pop_12mo(data = data)
}

if(config$short_term){
  truth$short_term <- truth_short_term(data = data, 
                                       V_u_months = config$V_u_months, 
                                       V_u_week_interval = as.numeric(config$V_u_week_interval))
}

saveRDS(truth, paste0(truth_dir, "/truth_", setting, ".Rds"))
