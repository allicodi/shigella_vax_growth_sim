# Code for "Estimating the impact of Shigella vaccines on growth outcomes and implications for clinical trial design"

This repository contains code used to run realistic simulated trials of *Shigella* vaccines with a focus on estimating the effect of the vaccine on growth. The code is configured to allow easy changes to the trial design and epidemiologic setting.

Included are all pieces of code needed to reproduce results of the manuscript "Estimating the impact of Shigella vaccines on growth outcomes and implications for clinical trial design," though the framework can be used to extend beyond the settings considered therein.

Below we include a detailed description of the general workflow. The workflow is intended to be executed on an HPC system using a `slurm` scheduler. Note: As configurations of HPC computing environments differ across clusters, it is not expected that this code will run without modification on any HPC system. 

**Note:** All analyses use the **[vaxstrat](https://github.com/allicodi/vaxstrat)** R package, which implements principal stratification methods for estimating causal effects of vaccination on post-infection outcomes. The package supports multiple estimands (Naturally Infected, Doomed, and population-level effects) and estimators (G-computation, AIPW, TMLE, and nonparametric bounds).


---

## Workflow Overview

The simulation pipeline consists of four main steps, managed through the `Makefile`:

1. **Parameter Generation** (`get_parameters.R`): Creates the data-generating process parameters from settings in `config.yml`
2. **Simulation Execution** (`run_analysis.R`): Runs the main simulation across multiple seeds and sample sizes
3. **Truth Calculation** (`get_truth.R`): Computes true causal effects via large Monte Carlo simulation
4. **Performance Evaluation** (`evaluate_performance.R`): Aggregates results and computes operating characteristics

Steps 1 and 2 must be run in order. Step 3 (truth calculation) can be run anytime before step 4, including after the main simulation completes.

**To run the complete workflow:**

```bash
make full_analysis SETTING=<scenario_name>
```


**To run individual steps:**
```bash
make $(PARAMETERS_FILE) SETTING=  # Generate parameters
make run_analysis SETTING=         # Run simulations  
make truth SETTING=                # Compute truth (can run anytime before evaluate)
make evaluate_performance SETTING= # Compute performance metrics
```

---

## `config.yml`

The `config.yml` file is used to define the data-generating process, the scientific assumptions under evaluation, and the analysis strategy for the simulation study. The parameters in `config.yml` can be grouped into several conceptual components that reflect the structure of the simulation study.

### 1. Trial Design and Setting

These parameters define the high-level structure of the simulated clinical trial.

- `n_sample_size` defines the set of trial sample sizes considered. Each value corresponds to a separate design point, allowing the study to characterize how performance scales with trial size.
- `dose_schedule` determines the vaccine schedule (e.g., series completed by 6-month vs 12-month).
- `site` identifies which EFGH site should be used to estimate the baseline growth distribution.
- `incidence_shigella_0_6` and `incidence_shigella_6_12` define incidence rates of Shigella infection during months 0 to 6 and months 6 to 12 of followup.
- `incidence_severe_shigella_0_6` and `incidence_severe_shigella_6_12` define the corresponding incidence of severe disease during the same trial periods.
- `effect_shigella_growth_formula` specifies the functional form relating time since Shigella infection to subsequent growth measurements. This model is fit to the estimated post-infection growth effects estimated using MAL-ED.
- `scale_growth_effect_0_6` and `scale_growth_effect_6_12` allows the user to modify the magnitude of the growth effect. These numbers should be between -1 and 1. Point estimates are scaled towards their 95\% confidence limits prior to fitting the `effect_shigella_growth_formula` above.
- `dropout` introduces loss to follow-up, affecting both power and bias.

### 2. Vaccine Efficacy Assumptions

These parameters define how vaccination modifies the infection process.

- `VE_mild` specifies vaccine efficacy against mild Shigella infection.
- `VE_severe` specifies vaccine efficacy against severe Shigella infection.

---

### 3. Estimands and Identification Assumptions for Estimation

These parameters determine which causal effects are targeted and what assumptions are imposed.

- `estimand` specifies the causal estimands of interest (e.g., effects in the naturally infected population vs the full population).
- `exclusion_restriction`, `cross_world`, and related flags control structural assumptions required for identification of naturally infected effects.
- `two_stage` governs how growth regression models are estimated for the exclusion restriction only Naturally Infected effect estimator (not ultimately used in the simulation).
- `nat_inf_unadj` includes an unadjusted estimate of the Naturally Infected effect based on only the exclusion restriction (not ultimately used in the simulation).
- `estimators` specifies which statistical methods are applied. AIPW only was used in the final simulation.
- `n_boot` if requested estimators utilize bootstrap, specifies the number of bootstrap replicates used (not ultimately used by AIPW in the simulation)

---

### 4. Follow-up Windows and Outcome Definition

- `intervals` defines the time points or grouped time points at which growth effects are evaluated. These intervals allow the analysis to assess effects at specific follow-up times (e.g., 3, 6, 12 months) and/or average effects over multiple time points. 

---

### 5. Inference and Testing Parameters

These parameters control how statistical inference is conducted.

- `null_hypothesis_value` defines the null effect used for hypothesis testing.
- `alpha_level` specifies the significance threshold.

---

## Parameter Generation: `get_parameters.R`

Before running simulations, this script generates and saves the data-generating process parameters for a given scenario. It reads configuration settings from `config.yml` and calls `simulate_parameters()` to:

- Estimate baseline growth distributions from the specified EFGH site
- Fit hazard models relating baseline HAZ to infection risk
- Calibrate model intercepts to match target incidence rates
- Generate monthly growth trajectory parameters
- Fit models for growth decrement after Shigella infection

The output is saved as `parameters_<SETTING>.Rds` and is required by both the main simulation and truth calculation scripts.

---

## Main Simulation: `run_analysis.R` and `run_simulation.sh`

The script `run_analysis.R` is the primary driver of the simulation pipeline. It is designed to be executed in parallel on a SLURM cluster using array jobs, where each task corresponds to a different random seed.

### Overview

For a given simulation scenario defined in `config.yml`, this script:

1. Loads configuration settings and pre-generated parameters  
2. Simulates trial data for a range of sample sizes  
3. Applies specified estimators to evaluate vaccine effects on growth  
4. Computes standard errors and confidence intervals  
5. Optionally performs bootstrap inference  
6. Saves results to disk for downstream aggregation and analysis  

Each SLURM task produces results for a single seed across all sample sizes and estimator configurations.

### Inputs

The script expects the following inputs via environment variables (typically set in a SLURM submission script):

- `SLURM_ARRAY_TASK_ID`: Used as the random seed for reproducibility  
- `SETTING`: Name of the scenario in `config.yml` to run  
- `PARAMETERS_FILE`: Path to an `.rds` file containing pre-generated parameters for the data-generating process  

### Submission via `run_simulation.sh`

The `run_simulation.sh` shell script submits the main simulation jobs to a SLURM cluster. It sets up the R environment, configures library paths, and launches an array job where each task corresponds to a different random seed.

**Usage:**
```bash
./run_simulation.sh    
```

Key inputs:
- `PARTITION`: SLURM partition to submit jobs to  
- `SETTING`: scenario name from `config.yml`  
- `PARAMETERS_FILE`: path to pre-generated simulation parameters  
- `NSEEDS`: number of parallel simulation replicates (array size)

---

## Truth Calculation: `get_truth.R` and `run_truth.sh`

To evaluate simulation performance (bias, coverage, power), true causal effects must be computed. The `get_truth.R` script simulates a very large trial (n = 10 million) and computes true effects for all estimands and timepoints specified in the configuration. Truth values can be calculated anytime after parameters are generated, though they are required before running `evaluate_performance.R`.

The script is submitted to SLURM via `run_truth.sh`:

**Usage:**
```bash
./run_truth.sh    
```

Key inputs:
- `PARTITION`: SLURM partition to submit jobs to  
- `SETTING`: scenario name from `config.yml`  
- `PARAMETERS_FILE`: path to pre-generated simulation parameters  
- `TRUTH_DIR`: directory where truth file will be saved

Output is saved as `truth_<SETTING>.Rds`.

---


## Post-processing and Results Summaries

After simulation runs are completed, the repository uses a small set of post-processing scripts to combine results across seeds, compute performance metrics, and generate  figures and tables.

### `evaluate_performance.R`

This script aggregates simulation results for a given scenario across all saved seed-specific output files, compares the estimated effects to the corresponding truth values, and computes key operating characteristics. These include bias, coverage, power, and the proportion of negative point estimates. The script saves these combined evaluation objects as a single `.Rds` file for downstream plotting and reporting.

### `make_results_figures.R`

This script reads the evaluated results objects and produces the main manuscript and supplement figures showing power and the proportion of negative estimated effects across sample sizes, trial settings, estimands, and outcome definitions. It also creates summary tables of true effect sizes across scenarios that is included in the supplement.

---

## Miscellaneous Scripts

### `make_simulation_design_figure.R`

This script generates a multi-panel figure summarizing key components of the simulation design for a given scenario. Panels include incidence rates by age and severity, the baseline HAZ distribution, the relationship between baseline growth and infection risk, the underlying growth trajectory, and the assumed effect of Shigella on growth.

---

### `make_simulation_design_tables_figs_supp.R`

This script produces supplementary tables and figures describing the simulation design across all scenarios. Outputs include tables of baseline HAZ distributions and incidence rates, as well as figures for growth trajectories and incidence rate ratios by age and recruitment strategy. Results are formatted for manuscript-ready inclusion and saved to `results/figures/`. 

## Sensitivity analysis

We also include a sensitivity analysis evaluating violations of the partial principal ignorability assumption (i.e., unmeasured confounding of infection status and post-infection growth outcomes). The analysis follows the same workflow as the primary simulation. The primary configuration file is replaced by `config_unmeas_conf.yml`, and primary bash scripts are replaced by `run_unmeasured_conf_simulation.sh` and `run_unmeas_conf_truth.sh`, with the remaining workflow unchanged.

## Questions and Citation

For questions about this code, please contact Allison Codi at [allison.codi@emory.edu].

If you use this code in your research, please cite the accompanying manuscript (under review, **[medRxiv](https://doi.org/10.64898/2026.04.03.26350105)**).
