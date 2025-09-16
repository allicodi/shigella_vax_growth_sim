CONFIG_FILE := config.yml
SETTING := debug_default
PARAMETERS_DIR := /projects/dbenkes/allison/shigella_vaccine_trial/parameters
PARAMETERS_FILE := $(PARAMETERS_DIR)/parameters_$(SETTING).Rds
TRUTH_DIR := /projects/dbenkes/allison/shigella_vaccine_trial/truth
TRUTH_FILE := $(TRUTH_DIR)/truth_$(SETTING).Rds
PARTITION := empire
NSEEDS := 100

full_analysis: run_analysis evaluate_performance

$(PARAMETERS_FILE): $(CONFIG_FILE) get_parameters.R
	/apps/R/4.4.0/bin/Rscript get_parameters.R $(SETTING) $(PARAMETERS_DIR)

$(TRUTH_FILE): $(CONFIG_FILE) $(PARAMETERS_FILE) get_truth.R
	/apps/R/4.4.0/bin/Rscript get_truth.R $(SETTING) $(PARAMETERS_FILE) $(TRUTH_DIR)

run_analysis: $(PARAMETERS_FILE) run_simulation.sh run_analysis.R
	./run_simulation.sh $(PARTITION) $(SETTING) $(PARAMETERS_FILE) $(NSEEDS)

evaluate_performance: $(TRUTH_FILE) evaluate_performance.R 
	/apps/R/4.4.0/bin/Rscript evaluate_performance.R $(SETTING) $(TRUTH_FILE)

