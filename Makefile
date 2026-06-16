CONFIG_FILE := config.yml
SETTING := default
PARAMETERS_DIR := /projects/dbenkes/allison/shigella_vaccine_trial/parameters
PARAMETERS_FILE := $(PARAMETERS_DIR)/parameters_$(SETTING).Rds
TRUTH_DIR := /projects/dbenkes/allison/shigella_vaccine_trial/truth
TRUTH_FILE := $(TRUTH_DIR)/truth_$(SETTING).Rds
PARTITION := empire
NSEEDS := 1000

full_analysis: run_analysis evaluate_performance

full_sensitivity_analysis: run_unmeas_conf_truth run_unmeas_conf_analysis

$(PARAMETERS_FILE): $(CONFIG_FILE) get_parameters.R
	/apps/R/4.4.0/bin/Rscript get_parameters.R $(SETTING) $(PARAMETERS_DIR)

$(TRUTH_FILE): $(CONFIG_FILE) $(PARAMETERS_FILE) get_truth.R run_truth.sh
	./run_truth.sh $(PARTITION) $(SETTING) $(PARAMETERS_FILE) $(TRUTH_DIR)
	
run_analysis: $(PARAMETERS_FILE) run_simulation.sh run_analysis.R
	./run_simulation.sh $(PARTITION) $(SETTING) $(PARAMETERS_FILE) $(NSEEDS)

evaluate_performance: $(TRUTH_FILE) evaluate_performance.R 
	/apps/R/4.4.0/bin/Rscript evaluate_performance.R $(SETTING) $(TRUTH_FILE)

.PHONY: truth
truth: $(TRUTH_FILE)

# ------------------------------------------------------------------ 
# Unmeasured-confounding sensitivity analyses 
# ------------------------------------------------------------------ 
run_unmeas_conf_truth: run_unmeas_conf_truth.sh R/get_sens_truth.R 
	./run_unmeas_conf_truth.sh $(PARTITION) $(SETTING) $(TRUTH_DIR) 
	
run_unmeas_conf_analysis: run_unmeas_conf_simulation.sh R/run_sens_analysis.R 
	./run_unmeas_conf_simulation.sh $(PARTITION) $(SETTING) $(NSEEDS)
