CONFIG_FILE := config.yml
SETTING := default
PRECOMPUTED_DIR := /projects/dbenkes/allison/shigella/.temp_results
TRUTH_FILE := $(PRECOMPUTED_DIR)/precomputed_values_$(SETTING)
PARTITION := empire

full_analysis: run_simulation 

$(TRUTH_FILE).Rds: $(CONFIG_FILE) get_true_values.R
	Rscript get_true_values.R $(SETTING)

run_simulation: $(TRUTH_FILE).Rds run_simulation.sh run_simulation_cluster.R
	./run_simulation.sh $(PARTITION) $(SETTING) 1000

evaluate_performance: evaluate_performance.R 
	Rscript evaluate_performance.R $(SETTING)
