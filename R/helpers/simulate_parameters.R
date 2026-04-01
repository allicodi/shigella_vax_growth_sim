
# here::i_am("R/simulate_parameters.R")

# source(here::here("R/parameter_generation_fns.R"))

#' Function to simulate parameters for Shigella VE growth simulation
#'
#' This function assembles and returns a list of parameters needed for simulating the impact of Shigella infection 
#' (and vaccine efficacy against it) on child growth. It incorporates baseline growth distributions, Shigella 
#' incidence rates, hazard models for Shigella and severe Shigella infection, and models for the growth effect of infection.
#' 
#' @param dose_schedule Character. Either `"6mo"` or `"12mo"`, indicating whether the vaccine schedule ends at 6 or 12 months.
#' @param site Character. EFGH site name to use when estimating the baseline growth distribution. Default is `"Overall"`.
#' @param incidence_shigella_0_6 Numeric. Cumulative 6-month incidence of Shigella infection during the first half of the trial (months 0–6). Typically from `get_incidence()`.
#' @param incidence_severe_shigella_0_6 Numeric. Cumulative 6-month incidence of *moderate-to-severe* Shigella infection during the first half of the trial (months 0–6). From `get_incidence()`.
#' @param incidence_shigella_6_12 Numeric. Cumulative 6-month incidence of Shigella infection during the second half of the trial (months 6–12). From `get_incidence()`.
#' @param incidence_severe_shigella_6_12 Numeric. Cumulative 6-month incidence of *moderate-to-severe* Shigella infection during the second half of the trial (months 6–12). From `get_incidence()`.
#' @param effect_shigella_growth_formula Formula. A spline-based formula (e.g., `y ~ -1 + x + I(pmax(0, x - 4)) + I(pmax(0, x - 8))`) used to model the effect of Shigella infection on monthly growth outcomes.
#' 
#' @returns list of parameters containing:
#' \describe{
#'   \item{dose_schedule}{Vaccine schedule used ("6mo" or "12mo").}
#'   \item{site}{Site name used to determine baseline growth parameters.}
#'   \item{mean_X}{Mean baseline LAZ at enrollment.}
#'   \item{sd_X}{Standard deviation of baseline LAZ at enrollment.}
#'   \item{hazard_S__X_int_0_6}{Intercept from hazard model for Shigella infection during months 0–6.}
#'   \item{hazard_S__X_coef_0_6}{Covariate coefficient for baseline LAZ in Shigella hazard model (months 0–6).}
#'   \item{incidence_S_mean_X_0_6}{Expected incidence of Shigella in months 0-6 at mean_X}
#'   \item{hazard_S__X_int_6_12}{Intercept from hazard model for Shigella infection during months 6–12.}
#'   \item{hazard_S__X_coef_6_12}{Covariate coefficient for baseline LAZ in Shigella hazard model (months 6–12).}
#'   \item{incidence_S_mean_X_6_12}{Expected incidence of Shigella in months 6–12 at mean_X}
#'   \item{hazard_S_sev__X_int_0_6}{Intercept from hazard model for severe Shigella infection during months 0–6.}
#'   \item{hazard_S_sev__X_coef_0_6}{Covariate coefficient in severe Shigella hazard model (months 0–6).}
#'   \item{incidence_S_sev_mean_X_0_6}{Expected incidence of severe Shigella in months 0-6 at mean_X.}
#'   \item{hazard_S_sev__X_int_6_12}{Intercept from hazard model for severe Shigella infection during months 6–12.}
#'   \item{hazard_S_sev__X_coef_6_12}{Covariate coefficient in severe Shigella hazard model (months 6–12).}
#'   \item{incidence_S_sev_mean_X_6_12}{Expected incidence of severe Shigella in months 6-12 at mean_X.}
#'   \item{effect_shigella_growth_fits}{List of fitted spline models estimating the growth impact of Shigella.}
#'   \item{monthly_growth_model}{Monthly LAZ trajectory model in the absence of Shigella infection.}
#' }
simulate_parameters <- function(dose_schedule = "6mo",                      # or 12mo
                                site = "Overall",
                                incidence_shigella_0_6 = 0.0253,            # from Maria; 6 month incidence
                                incidence_severe_shigella_0_6 = 0.0147,     # from Maria; 6 month incidence
                                incidence_shigella_6_12 = 0.0373,           # from Maria; 6 month incidence
                                incidence_severe_shigella_6_12 = 0.0203,    # from Maria; 6 month incidence
                                effect_shigella_growth_formula = y ~ -1 + x + I(pmax(0, x - 4)) + I(pmax(0, x - 8)),
                                scale_growth_effect_0_6 = 0, # factor to scale growth effects by between CIs in MALED model aka figure of the year
                                scale_growth_effect_6_12 = 0  # factor to scale growth effects by between CIs in MALED model aka figure of the year
                                ){
  
  # ---------------------------------------------------------------------------
  # Baseline HAZ: Mean and standard deviation for given age & site combo ------
  # ---------------------------------------------------------------------------
  
  # this will be a 6mo or 12mo distribution
  # 6mo = age 6-8 (don't have efgh before 6), # 12mo = age 9-11
  bl_growth_param <- get_baseline_growth(dose_schedule = dose_schedule,
                                         enroll_site = site)
  
  # ---------------------------------------------------------------------------
  # Effect of baseline HAZ on Shigella infection ------------------------------
  # ---------------------------------------------------------------------------
  
  hr_bl_haz_shig_inf <- hr_maled(dose_schedule = dose_schedule)
  
  # Find intercept that matches cumulative incidence and hazard for first period
  bl_haz_shig_combos_0_6 <- expand.grid(intercept = seq(-10, -4, by = 0.01),
                                    haz_coef = hr_bl_haz_shig_inf$coef_0_6,
                                    mean_X = bl_growth_param$mean_enr_haz)
  
  bl_haz_shig_cum_inc_0_6 <- cbind(bl_haz_shig_combos_0_6, data.frame(cum_inc_haz_mean = apply(bl_haz_shig_combos_0_6, 1, cum_inc_by_row)))
  bl_haz_shig_parameters_0_6 <- bl_haz_shig_cum_inc_0_6 %>%
    mutate(dist = sqrt(abs(cum_inc_haz_mean - incidence_shigella_0_6)^2)) %>%
    filter(dist == min(dist))
  
  # Repeat for cumulative incidence in second period
  bl_haz_shig_combos_6_12 <- expand.grid(intercept = seq(-10, -4, by = 0.01),
                                        haz_coef = hr_bl_haz_shig_inf$coef_6_12,
                                        mean_X = bl_growth_param$mean_enr_haz)
  
  bl_haz_shig_cum_inc_6_12 <- cbind(bl_haz_shig_combos_6_12,  data.frame(cum_inc_haz_mean = apply(bl_haz_shig_combos_6_12, 1, cum_inc_by_row)))
  bl_haz_shig_parameters_6_12 <- bl_haz_shig_cum_inc_6_12 %>%
    mutate(dist = sqrt(abs(cum_inc_haz_mean - incidence_shigella_6_12)^2)) %>%
    filter(dist == min(dist))
  
  # ---------------------------------------------------------------------------
  # Effect of baseline HAZ on severe infection --------------------------------
  # ---------------------------------------------------------------------------
  
  # Get coefficient and intercept for severe infection from model of MSD ~ BL Growth | Infected (using baseline growth in age stratum)
  effect_bl_haz_sev_shig <- get_effect_baseline_growth_sev_inf(dose_schedule = dose_schedule,
                                                               enroll_site = site)
  
  bl_haz_sev_shig_combos_0_6 <- expand.grid(intercept = seq(effect_bl_haz_sev_shig$sev_int_0_6 - 3, 
                                                        effect_bl_haz_sev_shig$sev_int_0_6 + 3, 
                                                        by = 0.01),
                                            haz_coef = effect_bl_haz_sev_shig$sev_coef_0_6)
  
  # get cumulative incidence of severe shigella at the mean HAZ based on coefficient + intercept combos from model, and cumulative inicidence of regular shigella from previous step
  bl_haz_sev_shig_combos_0_6$cum_inc_sev_inf_haz_mean_0_6 <- bl_haz_shig_parameters_0_6$cum_inc_haz_mean * plogis(bl_haz_sev_shig_combos_0_6$intercept + bl_haz_sev_shig_combos_0_6$haz_coef * bl_growth_param$mean_enr_haz)
  
  # get intercept for first interval conversion probability
  bl_haz_sev_shig_parameters_0_6 <- bl_haz_sev_shig_combos_0_6 %>%
    mutate(dist = sqrt((cum_inc_sev_inf_haz_mean_0_6 - incidence_severe_shigella_0_6)^2)) %>%
    filter(dist == min(dist)) 
  
  # repeat for second model
  bl_haz_sev_shig_combos_6_12 <- expand.grid(intercept = seq(effect_bl_haz_sev_shig$sev_int_6_12 - 3, 
                                                            effect_bl_haz_sev_shig$sev_int_6_12 + 3, 
                                                            by = 0.01),
                                            haz_coef = effect_bl_haz_sev_shig$sev_coef_6_12)
  
  bl_haz_sev_shig_combos_6_12$cum_inc_sev_inf_haz_mean_6_12 <- bl_haz_shig_parameters_6_12$cum_inc_haz_mean * plogis(bl_haz_sev_shig_combos_6_12$intercept + bl_haz_sev_shig_combos_6_12$haz_coef * bl_growth_param$mean_enr_haz)
  
  # get intercept for first interval conversion probability
  bl_haz_sev_shig_parameters_6_12 <- bl_haz_sev_shig_combos_6_12 %>%
    mutate(dist = sqrt((cum_inc_sev_inf_haz_mean_6_12 - incidence_severe_shigella_6_12)^2)) %>%
    filter(dist == min(dist)) %>%
    select(intercept, haz_coef, cum_inc_sev_inf_haz_mean_6_12, dist)
  
  # ---------------------------------------------------------------------------
  # Effect of Shigella on growth (S_inf --> Y_m; S_sev --> Y_m)
  # ---------------------------------------------------------------------------

  # this is model for all-ages, could change to age-stratified models? or do something with MSM? hard coded for now
  shigella_growth_meta_analysis_results <- readRDS(here::here("misc/results/case_control/shigella_growth_effect_data.Rds"))
  
  shigella_growth_fits <- fit_effect_shigella_growth_models(plot_df = shigella_growth_meta_analysis_results,
                                                            spline_formula = effect_shigella_growth_formula,
                                                            scale_growth_effect_0_6 = scale_growth_effect_0_6, 
                                                            scale_growth_effect_6_12 = scale_growth_effect_6_12,
                                                            dose_schedule = config$dose_schedule)
  
  # ---------------------------------------------------------------------------
  # Monthly HAZ outcome (Y_m) in absence of Shigella 
  # ---------------------------------------------------------------------------
  
  monthly_growth <- get_monthly_growth(dose_schedule = dose_schedule)
  
  # Parameters
  params <- list(dose_schedule = dose_schedule,
                 site = site,
                 # Distribution of baseline growth (X)
                 mean_X = bl_growth_param$mean_enr_haz, 
                 sd_X = bl_growth_param$sd_enr_haz,
                 # Hazard for Shigella months 0-6 | X
                 hazard_S__X_int_0_6 = bl_haz_shig_parameters_0_6$intercept,
                 hazard_S__X_coef_0_6 = bl_haz_shig_parameters_0_6$haz_coef,
                 incidence_S_mean_X_0_6 = bl_haz_shig_parameters_0_6$cum_inc_haz_mean,
                 # Hazard for Shigella months 6-12 | X
                 hazard_S__X_int_6_12 = bl_haz_shig_parameters_6_12$intercept,
                 hazard_S__X_coef_6_12 = bl_haz_shig_parameters_6_12$haz_coef,
                 incidence_S_mean_X_6_12 = bl_haz_shig_parameters_6_12$cum_inc_haz_mean,
                 # Hazard for severe Shigella months 0-6 | X
                 hazard_S_sev__X_int_0_6 = bl_haz_sev_shig_parameters_0_6$intercept,
                 hazard_S_sev__X_coef_0_6 = bl_haz_sev_shig_parameters_0_6$haz_coef,
                 incidence_S_sev_mean_X_0_6 = bl_haz_sev_shig_parameters_0_6$cum_inc_sev_inf_haz_mean_0_6,
                 # Hazard for severe Shigella months 6-12 | X
                 hazard_S_sev__X_int_6_12 = bl_haz_sev_shig_parameters_6_12$intercept,
                 hazard_S_sev__X_coef_6_12 = bl_haz_sev_shig_parameters_6_12$haz_coef,
                 incidence_S_sev_mean_X_6_12 = bl_haz_sev_shig_parameters_6_12$cum_inc_sev_inf_haz_mean_6_12,
                 # Effect of Shigella on growth 
                 effect_shigella_growth_fits = shigella_growth_fits,
                 # Monthly growth model in absence of infection
                 monthly_growth_model = monthly_growth)
  
  class(params) <- "short_term_growth_parameters"
  
  return(params)

}
