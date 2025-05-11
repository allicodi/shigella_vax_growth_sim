here::i_am("run_simulation_cluster.R")

devtools::load_all(here::here("../vegrowth/"))

source(here::here("R/00_simulate_data.R"))
source(here::here("R/01_get_truth.R"))
source(here::here("R/02_do_one_ci_sim.R"))
source(here::here("R/03_evaluate_performance.R"))

# get seed & config settings from bash script
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
setting <- Sys.getenv("SETTING")

# read in config file
config <- config::get(file = "config.yml", config = setting)
if(is.null(config$zhifei)) config$zhifei <- FALSE
if(is.null(config$catch_up)) config$catch_up <- FALSE

if(is.null(config$G_X_model)) config$G_X_model <- FALSE
if(is.null(config$Yinf_X_model)) config$Yinf_X_model <- FALSE

n_boot <- config$n_boot
est <- config$est

# Load precomputed truths and vaccine efficacy from get_true_values.R
param_combo_df <- readRDS(paste0("/projects/dbenkes/allison/shigella/.temp_results/precomputed_values_add_growth_",setting,".Rds"))
#param_combo_df <- readRDS(paste0("/projects/dbenkes/allison/shigella/.temp_results/quad_precomputed_values_add_growth_",setting,".Rds"))

# generate full parameter grid
param_grid <- expand.grid(
  seed = seed,
  n = config$n_sample_size,
  zhifei = config$zhifei,
  catch_up = config$catch_up,
  G_X_model = config$G_X_model,
  Yinf_X_model = config$Yinf_X_model
)

param_grid <- merge(param_grid, param_combo_df)

results <- lapply(1:nrow(param_grid), function(x, n_boot, est){
  params <- param_grid[x,]
  do_one_ci_sim(params = params, n_boot = n_boot, est = est)
}, n_boot = n_boot, est = est)

results <- do.call(rbind, results)
rownames(results) <- 1:nrow(results)

saveRDS(results, file = paste0("/projects/dbenkes/allison/shigella/.temp_results/optimistic_result_seed_", seed, "_setting_", setting, ".rds" ))
