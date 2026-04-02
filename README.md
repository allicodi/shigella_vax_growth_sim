# Code for "Estimating the impact of Shigella vaccines on growth outcomes and implications for clinical trial design"

This repository contains code used to run realistic simulated trials of *Shigella* vaccines with a focus on estimating the effect of the vaccine on growth. The code is configured to allow easy changes to the trial design and epidemiologic setting.

Included are all pieces of code needed to reproduce results of the manuscript "Estimating the impact of Shigella vaccines on growth outcomes and implications for clinical trial design," though the framework can be used to extend beyond the settings considered therein.

Below we include a detailed description of the general workflow. The workflow is intended to be executed on an HPC system using a `slurm` scheduler. Note: As configurations of HPC computing environments differ across clusters, it is not expected that this code will run without modification on any HPC system. 

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
- `n_boot` if requested estimators utilize bootstrap, specifies the number of bootstrap replicates used.

---

### 4. Follow-up Windows and Outcome Definition

- `intervals` defines the time points or grouped time points at which growth effects are evaluated. These intervals allow the analysis to assess effects at specific follow-up times (e.g., 3, 6, 12 months) and/or average effects over multiple time points. 

---

### 5. Inference and Testing Parameters

These parameters control how statistical inference is conducted.

- `null_hypothesis_value` defines the null effect used for hypothesis testing.
- `alpha_level` specifies the significance threshold.


## Main Simulation Script: `run_analysis.R`

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

---

### Inputs

The script expects the following inputs via environment variables (typically set in a SLURM submission script):

- `SLURM_ARRAY_TASK_ID`: Used as the random seed for reproducibility  
- `SETTING`: Name of the scenario in `config.yml` to run  
- `PARAMETERS_FILE`: Path to an `.rds` file containing pre-generated parameters for the data-generating process  


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

---

### `run_simulation.sh`

This shell script submits the main simulation jobs to a SLURM cluster. It sets up the R environment, configures library paths, and launches an array job where each task corresponds to a different random seed. The script passes the selected configuration (`SETTING`) and parameter file (`PARAMETERS_FILE`) to `run_analysis.R`, enabling parallel execution of simulation replicates across compute nodes.

Key inputs:
- `PARTITION`: SLURM partition to submit jobs to  
- `SETTING`: scenario name from `config.yml`  
- `PARAMETERS_FILE`: path to pre-generated simulation parameters  
- `NSEEDS`: number of parallel simulation replicates (array size)  

