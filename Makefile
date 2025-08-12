CONFIG_FILE := config.yml
SETTING := default
PARAMETERS_DIR := /projects/dbenkes/allison/shigella_vaccine_trial/parameters
PARAMETERS_FILE := $(PARAMETER_DIR)/parameters_$(SETTING).Rds
TRUTH_DIR := /projects/dbenkes/allison/shigella_vaccine_trial/truth
TRUTH_FILE := $(TRUTH_DIR)/truth_$(SETTING).Rds
PARTITION := empire

#full_analysis: run_analysis evaluate_performance
full_analysis: run_analysis 

$(PARAMETERS_FILE).Rds: $(CONFIG_FILE) get_parameters.R
	Rscript get_parameters.R $(SETTING) $(PARAMETERS_DIR)
	
$(TRUTH_FILE).Rds: $(CONFIG_FILE) $(PARAMETERS_FILE).Rds get_truth.R
	Rscript get_truth.R $(SETTING) $(PARAMETERS_FILE) $(TRUTH_DIR)

run_analysis: $(PARAMETERS_FILE) run_simulation.sh run_analysis.R
	./run_simulation.sh $(PARTITION) $(SETTING) $(PARAMETERS_FILE) 1000

evaluate_performance: evaluate_performance.R 
	Rscript evaluate_performance.R $(SETTING) $(TRUTH_FILE)