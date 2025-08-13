# ---------------------------------------------------------------------------------------------
# Function to get parameters for given config setting
# ---------------------------------------------------------------------------------------------

.libPaths("~/Rlibs")

options(echo = TRUE)

sessionInfo()
.libPaths()

here::i_am("get_parameters.R")

source(here::here("R/parameter_generation_fns.R"))
source(here::here("R/simulate_parameters.R"))

# get setting from bash script
cargs <- commandArgs(TRUE)
setting <- cargs[1]
parameter_dir <- cargs[2]

config <- config::get(file = "config.yml", config = setting)

parameters <- simulate_parameters(dose_schedule = config$dose_schedule,
                                  site = config$site,
                                  incidence_shigella_0_6 = config$incidence_shigella_0_6, 
                                  incidence_severe_shigella_0_6 = config$incidence_severe_shigella_0_6, 
                                  incidence_shigella_6_12 = config$incidence_shigella_6_12, 
                                  incidence_severe_shigella_6_12 = config$incidence_severe_shigella_0_6, 
                                  effect_shigella_growth_formula = config$effect_shigella_growth_formula)

saveRDS(parameters, paste0(parameter_dir, "/parameters_", setting, ".Rds"))
