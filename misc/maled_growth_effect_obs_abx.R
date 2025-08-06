# -----------------------------------------------------------------------------
# Script for Shigella effect on growth regardless of abx
# -----------------------------------------------------------------------------

here::i_am("misc/maled_growth_effect_obs_abx.R")

devtools::load_all("~/Documents/shigella_projects/packages/abxGrowth/")

library(SuperLearner)

source(here::here("misc/SL.wrappers.R"))

# change option for age stratified
age_stratified <- TRUE

age_6_9 <- TRUE
age_9_12 <- FALSE
age_12_15 <- FALSE

lsd <- FALSE
msd <- TRUE

# to save msm results
msm_age_df_all <- data.frame()

for(month in 1:12){
  
  data <- readRDS(here::here(paste0("misc/maled_data/maled_",month,"mo.Rds")))
  
  if(age_stratified){
    if(age_6_9){
      data <- data[data$agemonths >= 6 & data$agemonths < 9,]
    } else if(age_9_12){
      data <- data[data$agemonths >= 9 & data$agemonths < 12,]
    } else{
      data <- data[data$agemonths >= 12 & data$agemonths < 15,]
    }
  }
  
  if(msd){
    # Get IDs corresponding to MSD cases
    maled_data_msd_caseids <- unique(data$case_id[which(data$MSD == 1)])
    
    # Subset to MSD cases and their matched controls
    data <- data[which(data$case_id %in% maled_data_msd_caseids),]
  } 
  
  if(lsd){
    maled_data_lsd_caseids <- unique(data$case_id[which(data$case == 1 & data$MSD == 0)])
    
    # Subset to LSD cases and their matched controls
    data <- data[which(data$case_id %in% maled_data_lsd_caseids),]
  }
  
  # Outcome 1, Missingness
  sl.library.with.abx <- list(c("SL.glm", "SL.screen.abx.lt.min_prop"),
                              c("SL.glm", "SL.screen.abx.glmnet"),
                              c("SL.glm.spline.age.haz", "SL.screen.abx.lt.min_prop"),
                              "SL.step.forward.spline.age.haz.abx",
                              "SL.glmnet",
                              "SL.ranger",
                              "SL.earth",
                              "SL.xgboost")
  
  # Propensity, outcome 2
  sl.library.without.abx <- list(c("SL.glm", "SL.screen.abx.lt.min_prop"),
                                 c("SL.glm", "screen.glmnet"),
                                 c("SL.glm.spline.age.haz", "SL.screen.abx.lt.min_prop"),
                                 "SL.step.forward.spline.age.haz",
                                 "SL.glmnet",
                                 "SL.ranger",
                                 "SL.earth",
                                 "SL.xgboost")
  
  # Impute missing values, cases and controls separate
  maled_case_data <- data[data$case == 1,]
  maled_control_data <- data[data$case == 0,]
  
  maled_case_data_imp <- impute_covariates(data = maled_case_data,
                                           site_var_name = "site",
                                           imp_by_site = TRUE,
                                           imp_covariates = c("site",
                                                              "sex",
                                                              "agemonths",
                                                              "baseline_haz",
                                                              "mated_bin", 
                                                              "drinkimp",
                                                              "sanitimp",
                                                              "adenovirus_40_41_new",                  
                                                              "aeromonas_new",                         
                                                              "astrovirus_new",                        
                                                              "campylobacter_pan_new",                 
                                                              "cryptosporidium_new",                   
                                                              "cyclospora_new",                        
                                                              "e_histolytica_new",                     
                                                              "isospora_new",                          
                                                              "norovirus_new",                         
                                                              "rotavirus_new",                         
                                                              "salmonella_new",                        
                                                              "sapovirus_new",  
                                                              "shigella_new",
                                                              "st_etec_new",                           
                                                              "tEPEC_new",                             
                                                              "v_cholerae_new",                        
                                                              "ETEC_new",                              
                                                              "e_bieneusi_new",                        
                                                              "eaec_new",
                                                              "dysentery",
                                                              "fever",
                                                              "fever_days",
                                                              "dehyd",
                                                              "lsstools",
                                                              "daysvomit"))
  
  maled_control_data_imp <- impute_covariates(data = maled_control_data,
                                              site_var_name = "site",
                                              imp_by_site = TRUE,
                                              imp_covariates = c("site",
                                                                 "sex",
                                                                 "agemonths",
                                                                 "baseline_haz",
                                                                 "mated_bin", 
                                                                 "drinkimp",
                                                                 "sanitimp"))
  
  maled_case_data_imp <- maled_case_data_imp[,colnames(maled_case_data_imp) %in% colnames(data), drop = FALSE]
  maled_control_data_imp <- maled_control_data_imp[,colnames(maled_control_data_imp) %in% colnames(data), drop = FALSE]
  maled_data_imp <- rbind(maled_case_data_imp, maled_control_data_imp)
  
  # Add "age" and "enr_haz" to match wrappers (go back and rename in raw data later)
  maled_data_imp$age <- maled_data_imp$agemonths
  maled_data_imp$enr_haz <- maled_data_imp$baseline_haz
  
  one_hot_maled <- one_hot_encode(data = maled_data_imp,
                                  laz_var_name = "monthx_haz",
                                  covariate_list = c("site",
                                                     "sex",
                                                     "age",
                                                     "enr_haz",
                                                     "mated_bin", 
                                                     "drinkimp",
                                                     "sanitimp",
                                                     "I_followup_days",
                                                     "I_followup_days_x_followup_days"),
                                  abx_var_name = "all_abx",
                                  site_var_name = "site",
                                  site_interaction = TRUE,
                                  severity_list = c("dysentery",
                                                    "fever",
                                                    "fever_days",
                                                    "dehyd",
                                                    "lsstools",
                                                    "daysvomit",
                                                    "duration_pre_abx"),
                                  age_var_name = "age")
  
  # SET PARAMETERS
  data = one_hot_maled$data
  laz_var_name = "monthx_haz"
  abx_var_name = "all_abx"
  case_var_name = "case"
  followup_var_names = c("I_followup_days", "I_followup_days_x_followup_days")
  site_var_name = one_hot_maled$site_var_names
  covariate_list = one_hot_maled$covariate_list
  pathogen_quantity_list = c(  "adenovirus_40_41_new",                  
                               "aeromonas_new",                         
                               "astrovirus_new",                        
                               "campylobacter_pan_new",                 
                               "cryptosporidium_new",                   
                               "cyclospora_new",                        
                               "e_histolytica_new",                     
                               "isospora_new",                          
                               "norovirus_new",                         
                               "rotavirus_new",                         
                               "salmonella_new",                        
                               "sapovirus_new",  
                               "shigella_new",
                               "st_etec_new",                           
                               "tEPEC_new",                             
                               "v_cholerae_new",                        
                               "ETEC_new",                              
                               "e_bieneusi_new",                        
                               "eaec_new")
  severity_list = one_hot_maled$severity_list
  first_id_var_name = "first_id"
  outcome_type = "gaussian"
  sl.library.outcome.case = sl.library.with.abx
  sl.library.outcome.control = sl.library.without.abx
  sl.library.treatment = sl.library.without.abx
  sl.library.infection = sl.library.without.abx
  sl.library.missingness.case = sl.library.with.abx
  sl.library.missingness.control = sl.library.without.abx
  v_folds = 3
  return_models = TRUE
  msm = TRUE
  msm_var_name = "age"
  msm_formula = "age"
  case_control = TRUE
  seed = 12345
  
  msm_age_df <- data.frame(child_id = data$child_id,
                           case_id = data$case_id,
                           first_id = data$first_id,
                           case = data$case,
                           age = data$agemonths,
                           outcome_month = month)
  
  # ------------------------------------------------------------
  # STEP 0: Create subsets of data for model fitting
  # ------------------------------------------------------------
  
  set.seed(seed)
  
  # Rename abx levels if spaces in them
  abx_levels <- levels(factor(data[[abx_var_name]]))
  abx_levels_new <- ifelse(is.na(abx_levels), NA, gsub("[ /]", "_", abx_levels))
  data[[abx_var_name]] <- factor(data[[abx_var_name]], levels = abx_levels, labels = abx_levels_new)
  
  # get case vs control data
  case_data <- data[data[[case_var_name]] == 1,]
  control_data <- data[data[[case_var_name]] == 0,]
  
  case_data_idx <- which(data[[case_var_name]] == 1)
  control_data_idx <- which(data[[case_var_name]] == 0)
  
  # Case data prep
  I_Y_case <- ifelse(is.na(case_data[[laz_var_name]]), 1, 0)
  Y_case <- case_data[[laz_var_name]]
  covariates_case <- case_data[, covariate_list, drop = FALSE]
  severity_case <- case_data[, severity_list, drop = FALSE]
  pathogen_q_case <- case_data[, pathogen_quantity_list, drop = FALSE]
  abx_case <- case_data[, abx_var_name, drop = FALSE]
  
  # Complete case data (excluding missing outcome)
  case_data_complete <- case_data[!is.na(case_data[[laz_var_name]]), ]
  
  Y_case_complete <- case_data_complete[[laz_var_name]]
  covariates_case_complete <- case_data_complete[, covariate_list, drop = FALSE]
  severity_case_complete <- case_data_complete[, severity_list, drop = FALSE]
  pathogen_q_case_complete <- case_data_complete[, pathogen_quantity_list, drop = FALSE]
  abx_case_complete <- case_data_complete[, abx_var_name, drop = FALSE]
  
  # Control data prep
  I_Y_control <- ifelse(is.na(control_data[[laz_var_name]]), 1, 0)
  Y_control <- control_data[[laz_var_name]]
  if(is.null(covariate_list)){
    covariate_list <- covariate_list
  }
  covariates_control <- control_data[, covariate_list, drop = FALSE]
  
  # Complete control data (excluding missing outcome)
  control_data_complete <- control_data[!is.na(control_data[[laz_var_name]]), ]
  
  Y_control_complete <- control_data_complete[[laz_var_name]]
  covariates_control_complete <- control_data_complete[, covariate_list, drop = FALSE]
  
  # ------------------------------------------------------------
  # STEP 1: Fit & predict from outcome models
  # ------------------------------------------------------------
  
  ## Model 1a: Outcome model in cases
  outcome_model_1a <- SuperLearner::SuperLearner(Y = Y_case_complete, 
                                                 X = data.frame(abx_case_complete,
                                                                covariates_case_complete, 
                                                                severity_case_complete, 
                                                                pathogen_q_case_complete),
                                                 family = outcome_type, 
                                                 SL.library = sl.library.outcome.case,
                                                 cvControl = list(V = v_folds))
  
  ## Model 1b: Outcome model in controls
  outcome_model_1b <- SuperLearner::SuperLearner(Y = Y_control_complete, 
                                                 X = data.frame(covariates_control_complete),
                                                 family = outcome_type, 
                                                 SL.library = sl.library.outcome.control,
                                                 cvControl = list(V = v_folds))
  
  ## Predict outcomes
  
  # For each level of abx, predict setting abx = x, infection = 1 & abx = x, infection = 0
  abx_levels <- unique(case_data[[abx_var_name]])[!is.na(unique(case_data[[abx_var_name]]))]
  
  outcome_vectors_1a <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
  outcome_vectors_1b <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
  
  colnames(outcome_vectors_1a) <- paste0("abx_observed")
  colnames(outcome_vectors_1b) <- paste0("control")
  
  outcome_vectors_1a[case_data_idx, 1] <- stats::predict(outcome_model_1a, newdata = case_data[, c(abx_var_name,
                                                                                                   covariate_list,
                                                                                                   severity_list,
                                                                                                   pathogen_quantity_list)], type = "response")$pred
  
  # Replace NA controls with 0
  outcome_vectors_1a[is.na(outcome_vectors_1a)] <- 0
  
  outcome_vectors_1b[, 1] <- stats::predict(outcome_model_1b, newdata = data[,covariate_list], type = "response")$pred
  
  # ------------------------------------------------------------
  # STEP 2: Fit & predict from propensity models
  # ------------------------------------------------------------
  
  # If followup_days related variables, remove from missingness models
  if(!any(is.na(followup_var_names))){
    covariate_list_no_followup <- covariate_list[!(covariate_list %in% followup_var_names)]
  } else{
    covariate_list_no_followup <- covariate_list
  }
  
  prop_vectors_2a <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
  prop_vectors_3a <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
  prop_vectors_3b <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
  
  colnames(prop_vectors_2a) <- c("case")
  colnames(prop_vectors_3a) <- paste0("abx_observed")
  colnames(prop_vectors_3b) <- c("control")
  
  ## Part 2: Propensity models for shigella (or other infection) attribution
  
  # 2a_1 = Shigella Attributable ~ BL Cov
  prop_model_2a <- SuperLearner::SuperLearner(Y = data[[case_var_name]],
                                              X = data[, covariate_list, drop = FALSE],
                                              family = stats::binomial(),
                                              SL.library = sl.library.infection, 
                                              cvControl = list(V = v_folds))
  tmp_pred_2a <- prop_model_2a$SL.pred
  prop_vectors_2a[,1] <- tmp_pred_2a
  
  ## Part 3: Propensity models for missingness
  
  covariates_case_no_site <- case_data[, covariate_list_no_followup, drop = FALSE]
  covariates_control_no_site <- control_data[, covariate_list_no_followup, drop = FALSE]
  
  ## Missingness model in cases
  prop_model_3a <- SuperLearner::SuperLearner(Y = I_Y_case,
                                              X = data.frame(abx_case,
                                                             covariates_case_no_site,
                                                             severity_case,
                                                             pathogen_q_case),
                                              family = stats::binomial(),
                                              SL.library = sl.library.missingness.case,
                                              cvControl = list(V = v_folds))
  
  ## Missingness model in controls
  prop_model_3b <- SuperLearner::SuperLearner(Y = I_Y_control,
                                              X = data.frame(covariates_control_no_site),
                                              family = stats::binomial(),
                                              SL.library = sl.library.missingness.control,
                                              cvControl = list(V = v_folds))
  
  # Predict setting each abx level
  
  prop_vectors_3a[case_data_idx,1] <- prop_model_3a$SL.pred
  # Fill in NAs with 0
  prop_vectors_3a[is.na(prop_vectors_3a)] <- 0
  
  prop_vectors_3b[control_data_idx,1] <- prop_model_3b$SL.pred
  prop_vectors_3b[is.na(prop_vectors_3b)] <- 0
  
  # -------------------------------------------------
  # STEP 3: AIPW estimates and confidence intervals
  # -------------------------------------------------
  
  ## Plug-in estimates
  plug_ins_case <- mean(outcome_vectors_1a[case_data_idx,])
  plug_ins_control <- mean(outcome_vectors_1b[case_data_idx,])
  
  # 'individual level effect' save to do MSM stuff later
  outcome_vectors_difference <- outcome_vectors_1a - outcome_vectors_1b
  
  msm_age_df$outcome_difference <- outcome_vectors_difference$abx_observed
  msm_age_df_all <- rbind(msm_age_df_all, msm_age_df)
  
  ## EIF for bias corrections
  case_eifs <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
  control_eifs <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
  
  colnames(case_eifs) <- paste0("case_eif_observed")
  colnames(control_eifs) <- paste0("control_eif_observed")
  
  # Truncate any large propensity scores
  ps_trunc_level <- 0.01
  if(!is.na(ps_trunc_level)){
    #  any observation that is < ps_trunc_level or > 1-ps_trunc_level should be changed to ps_trunc_level or 1 - ps_trunc_level
    
    truncate_ps <- function(mat, ps_trunc_level){
      mat[mat < ps_trunc_level] <- ps_trunc_level
      mat[mat > 1- ps_trunc_level] <- 1- ps_trunc_level
      return(mat)
    }
    
    prop_vectors_2a <- truncate_ps(prop_vectors_2a, ps_trunc_level)
    prop_vectors_3a <- truncate_ps(prop_vectors_3a, ps_trunc_level)
    prop_vectors_3b <- truncate_ps(prop_vectors_3b, ps_trunc_level)
  }
  
  # 1 - Bias correction for case, abx level = a
  
  I_Case <- data[[case_var_name]]
  P_Case <- mean(prop_vectors_2a[,1])
  
  I_Delta_0 <- as.numeric(!is.na(data[[laz_var_name]])) # Indicator NOT missing
  #P_Delta_0__Case_all <- 1 - prop_vectors_3a[,i]
  P_Delta_0__Case_all <- 1 - prop_vectors_3a[,1]
  
  obs_outcome <- ifelse(is.na(data[[laz_var_name]]), 0, data[[laz_var_name]])  
  Qbar_Case_Abx_a_Covariates <- outcome_vectors_1a[,1]
  
  eif_case_vec <- (I_Case / P_Case) * (I_Delta_0 / P_Delta_0__Case_all) * (obs_outcome - Qbar_Case_Abx_a_Covariates) +
    (I_Case / P_Case) * (Qbar_Case_Abx_a_Covariates - plug_ins_case[1])
  
  # correct for / 0 with P_Abx_a__Case_Covariates, should be 0'd out
  eif_case_vec <- ifelse(is.nan(eif_case_vec), 0, eif_case_vec)
  
  case_eifs[,1] <- eif_case_vec
  
  
  # 2 - Bias correction for control
  I_control <- ifelse(data[[case_var_name]] == 1, 0, 1)
  P_control__Covariates <- 1 - prop_vectors_2a[,1]
  
  I_Delta_0 <- as.numeric(!is.na(data[[laz_var_name]])) # Indicator NOT missing
  I_Delta_0__Control_all <- 1 - prop_vectors_3b[,1]
  
  P_Case__all <- prop_vectors_2a[,1]
  
  Qbar_Control_Covariates <- outcome_vectors_1b[,1]
  
  eif_control_vec <- (I_control / P_control__Covariates) * (I_Delta_0 / I_Delta_0__Control_all) * (P_Case__all / P_Case) * (obs_outcome - Qbar_Control_Covariates) +
    (I_Case / P_Case) * (Qbar_Control_Covariates - plug_ins_control)
  
  control_eif <- eif_control_vec
  
  # Get AIPWs
  aipw_case <- plug_ins_case + colMeans(case_eifs)
  aipw_control <- plug_ins_control + colMeans(control_eif)  
  
  eif_matrix <- cbind(case_eifs, control_eif)
  
  # Get id for each participant and recreate EIFs based on this if present
  if(!is.null(first_id_var_name)){
    first_id_eif_matrix <- cbind(data.frame(first_id = data[[first_id_var_name]]), eif_matrix)
    first_id_eif_matrix <- aggregate(. ~ first_id, data = first_id_eif_matrix, FUN = sum)
    scaled_matrix <- first_id_eif_matrix[,-c(1)] * (nrow(first_id_eif_matrix) / nrow(eif_matrix))
    
  }else{
    scaled_matrix <- eif_matrix
  }
  
  aipws_effect <- vector("numeric", length = 1)
  eifs_effect <- data.frame(matrix(ncol = 1, nrow = nrow(scaled_matrix)))
  
  names(aipws_effect) <- paste0("effect_observed")
  colnames(eifs_effect) <- paste0("effect_observed")
  
  aipw_effect <- aipw_case[1] - aipw_control[1]
  aipws_effect[1] <- aipw_effect
  
  idx_1 <- 1
  idx_2 <- 1 + 1
  
  gradient <- rep(0, 2)
  
  gradient[idx_1] <- 1
  gradient[idx_2] <- -1
  
  gradient <- matrix(gradient, ncol = 1)
  eif_effect <- as.numeric(as.matrix(scaled_matrix) %*% gradient)
  eifs_effect[,1] <- eif_effect
  
  results_df <- data.frame(abx_levels = "observed",
                           abx_level_case = aipw_case,
                           abx_level_control = rep(aipw_control, 1),
                           effect_inf_abx_level = aipws_effect)
  
  eif_matrix_scaled <- cbind(scaled_matrix, eifs_effect)
  cov_matrix <- stats::cov(eif_matrix_scaled)
  eif_hat <- sqrt( diag(cov_matrix) / nrow(eif_matrix_scaled) )
  
  results_object <- list(results_df = results_df,
                         plug_ins_case = plug_ins_case,
                         plug_ins_control = plug_ins_control,
                         eif_matrix = eif_matrix_scaled,
                         se = eif_hat)
  
  class(results_object) <- "aipw_case_control_observed"
  
  
  aipw_models <- list(outcome_model_1a = outcome_model_1a,
                      outcome_model_1b = outcome_model_1b,
                      prop_model_2a = prop_model_2a,
                      prop_model_3a = prop_model_3a,
                      prop_model_3b = prop_model_3b)
  
  results <- list(results_object = results_object,
                  aipw_models = aipw_models)
  
  class(results) <- "agaipw_res"
  
  results$aipw_models <- NULL
  
  if(!age_stratified){
      
    if(msd){
      saveRDS(results, here::here(paste0("misc/results/case_control/ipd_case_control_no_abx_msd_", month, ".Rds")))
    } else if(lsd){
      saveRDS(results, here::here(paste0("misc/results/case_control/ipd_case_control_no_abx_lsd_", month, ".Rds")))
    } else{
      saveRDS(results, here::here(paste0("misc/results/case_control/ipd_case_control_no_abx_", month, ".Rds")))
    }
      

  } else{
    
    # ADD LSD HERE
    
    if(age_6_9){
      if(msd){
        saveRDS(results, here::here(paste0("misc/results/case_control/ipd_case_control_no_abx_6_9_msd_month_", month, ".Rds")))
      } else if(lsd){
        saveRDS(results, here::here(paste0("misc/results/case_control/ipd_case_control_no_abx_6_9_lsd_month_", month, ".Rds")))
      } else{
        saveRDS(results, here::here(paste0("misc/results/case_control/ipd_case_control_no_abx_6_9_month_", month, ".Rds")))
      }
      
    } else if (age_9_12){
      
      if(msd){
        saveRDS(results, here::here(paste0("misc/results/case_control/ipd_case_control_no_abx_9_12_msd_month_", month, ".Rds")))
      } else if(lsd){
        saveRDS(results, here::here(paste0("misc/results/case_control/ipd_case_control_no_abx_9_12_lsd_month_", month, ".Rds")))
      } else{
        saveRDS(results, here::here(paste0("misc/results/case_control/ipd_case_control_no_abx_9_12_month_", month, ".Rds")))
      }
      
    } else{
      if(msd){
        saveRDS(results, here::here(paste0("misc/results/case_control/ipd_case_control_no_abx_12_15_msd_month_", month, ".Rds")))
      } else if(lsd){
        saveRDS(results, here::here(paste0("misc/results/case_control/ipd_case_control_no_abx_12_15_lsd_month_", month, ".Rds")))
      } else{
        saveRDS(results, here::here(paste0("misc/results/case_control/ipd_case_control_no_abx_12_15_month_", month, ".Rds")))
      }
      
    }
  }
}

if(!age_stratified){
  if(msd){
    saveRDS(msm_age_df_all, here::here("misc/results/case_control/msm_age_msd_df.Rds"))
  } else if(lsd){
    saveRDS(msm_age_df_all, here::here("misc/results/case_control/msm_age_lsd_df.Rds"))
  } else{
    saveRDS(msm_age_df_all, here::here("misc/results/case_control/msm_age_df.Rds"))
  }
}


# --------------------------------------------------------------------------------------------------
  
  # Vector of months
  months <- 1:12
  
  # Load all the RDS files into a list
  mo_list <- lapply(months, function(mo) {
    readRDS(here::here(sprintf("misc/results/case_control/ipd_case_control_no_abx_%d.Rds", mo)))
  })
  
  mo_list_severe <- lapply(months, function(mo) {
    readRDS(here::here(sprintf("misc/results/case_control/ipd_case_control_no_abx_msd_%d.Rds", mo)))
  })
  
  mo_list_mild <- lapply(months, function(mo) {
    readRDS(here::here(sprintf("misc/results/case_control/ipd_case_control_no_abx_lsd_%d.Rds", mo)))
  })
  
  # Create a data frame with point estimates and confidence intervals
  plot_data <- data.frame(
    subgroup = paste0(months, " month"),
    pt_est = sapply(mo_list, function(mo) mo$results_object$results_df$effect_inf_abx_level),
    lower_ci = sapply(mo_list, function(mo) {
      est <- mo$results_object$results_df$effect_inf_abx_level
      se <- mo$results_object$se['effect_observed']
      est - 1.96 * se
    }),
    upper_ci = sapply(mo_list, function(mo) {
      est <- mo$results_object$results_df$effect_inf_abx_level
      se <- mo$results_object$se['effect_observed']
      est + 1.96 * se
    })
  )
  
  plot_data_severe <- data.frame(
    subgroup = paste0(months, " month"),
    pt_est = sapply(mo_list_severe, function(mo) mo$results_object$results_df$effect_inf_abx_level),
    lower_ci = sapply(mo_list_severe, function(mo) {
      est <- mo$results_object$results_df$effect_inf_abx_level
      se <- mo$results_object$se['effect_observed']
      est - 1.96 * se
    }),
    upper_ci = sapply(mo_list_severe, function(mo) {
      est <- mo$results_object$results_df$effect_inf_abx_level
      se <- mo$results_object$se['effect_observed']
      est + 1.96 * se
    })
  )
  
  plot_data_mild <- data.frame(
    subgroup = paste0(months, " month"),
    pt_est = sapply(mo_list_mild, function(mo) mo$results_object$results_df$effect_inf_abx_level),
    lower_ci = sapply(mo_list_mild, function(mo) {
      est <- mo$results_object$results_df$effect_inf_abx_level
      se <- mo$results_object$se['effect_observed']
      est - 1.96 * se
    }),
    upper_ci = sapply(mo_list_mild, function(mo) {
      est <- mo$results_object$results_df$effect_inf_abx_level
      se <- mo$results_object$se['effect_observed']
      est + 1.96 * se
    })
  )
  
  # Add numeric month column for plotting
  plot_data$month_num <- months
  plot_data_severe$month_num <- months
  plot_data_mild$month_num <- months
  
  # Plot
  # ggplot(plot_data, aes(x = month_num, y = pt_est)) +
  #   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "steelblue", alpha = 0.4) +
  #   geom_line(color = "steelblue", size = 1.2) +
  #   geom_point(color = "steelblue", size = 2.5) +
  #   scale_x_continuous(breaks = months) +
  #   labs(x = "Month", y = "Growth Effect",
  #        title = "Shigella growth effect in MAL-ED with observed abx use") + 
  #   theme_minimal(base_size = 14) 
  # 
  # ggplot(plot_data_severe, aes(x = month_num, y = pt_est)) +
  #   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "red3", alpha = 0.4) +
  #   geom_line(color = "red3", size = 1.2) +
  #   geom_point(color = "red3", size = 2.5) +
  #   scale_x_continuous(breaks = months) +
  #   labs(x = "Month", y = "Growth Effect",
  #        title = "Severe Shigella growth effect in MAL-ED with observed abx use") + 
  #   theme_minimal(base_size = 14) 
  # 
  # Add a column to identify the group
  plot_data$group <- "Any Shigella"
  plot_data_severe$group <- "Moderate-to-Severe Shigella"
  plot_data_mild$group <- "Less-severe Shigella"
  
  # Combine the datasets
  combined_plot_data <- rbind(plot_data_mild, plot_data, plot_data_severe)
  combined_plot_data$group <- factor(combined_plot_data$group, levels = c("Less-severe Shigella", "Any Shigella", "Moderate-to-Severe Shigella"))
  
  final_plot <- ggplot(combined_plot_data, aes(x = month_num, y = pt_est, color = group, fill = group)) +
    # CI ribbon for raw estimates
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3, color = NA) +
    # Raw line and points
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    # Smoothed piecewise linear fit
    stat_smooth(
      method = "glm",
      formula = y ~ -1 + x + I(pmax(0, x - 3)) + I(pmax(0, x - 9)),
      size = 1.2,
      linetype = "dotted",
      se = FALSE
    ) +
    # Axis and labels
    scale_x_continuous(breaks = months) +
    labs(
      x = "Month",
      y = "Growth Effect",
      title = "Longitudinal Shigella growth effects in MAL-ED with observed antibiotic use",
      color = "Shigella severity",
      fill = "Shigella severity"
    ) +
    scale_color_manual(values = c(
      "Any Shigella" = "#00468BFF",
      "Moderate-to-Severe Shigella" = "#ED0000FF",
      "Less-severe Shigella" = "#42B540FF"
    )) +
    scale_fill_manual(values = c(
      "Any Shigella" = "#00468BFF",
      "Moderate-to-Severe Shigella" = "#ED0000FF",
      "Less-severe Shigella" = "#42B540FF"
    )) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
  
  ggsave(here::here("misc/results/case_control/growth_effect_plot.jpeg"),final_plot,  dpi = 300, width = 12, height = 6)
  
  
  # Define the formula
  my_formula <- pt_est ~ -1 + x + I(pmax(0, x - 3)) + I(pmax(0, x - 9))
  
  # Fit the model for each group
  model_fits <- combined_plot_data %>%
    group_by(group) %>%
    group_map(~ glm(my_formula, data = .x), .keep = TRUE)
  
  # Optional: assign names to list of models
  names(model_fits) <- unique(combined_plot_data$group)
  
# --------------------------------------------------------------------------------------------------
# age stratified plots
  
  # Vector of months
  months <- 1:12
  
  # Load all the RDS files into a list
  mo_list_6_9 <- lapply(months, function(mo) {
    readRDS(here::here(sprintf("misc/results/case_control/ipd_case_control_no_abx_6_9_month_%d.Rds", mo)))
  })
  
  mo_list_severe_6_9 <- lapply(months, function(mo) {
    readRDS(here::here(sprintf("misc/results/case_control/ipd_case_control_no_abx_6_9_msd_month_%d.Rds", mo)))
  })
  
  mo_list_mild_6_9 <- lapply(months, function(mo) {
    readRDS(here::here(sprintf("misc/results/case_control/ipd_case_control_no_abx_6_9_lsd_month_%d.Rds", mo)))
  })
  
  # Load all the RDS files into a list
  mo_list_9_12 <- lapply(months, function(mo) {
    readRDS(here::here(sprintf("misc/results/case_control/ipd_case_control_no_abx_9_12_month_%d.Rds", mo)))
  })
  
  mo_list_severe_9_12 <- lapply(months, function(mo) {
    readRDS(here::here(sprintf("misc/results/case_control/ipd_case_control_no_abx_9_12_msd_month_%d.Rds", mo)))
  })
  
  mo_list_mild_9_12 <- lapply(months, function(mo) {
    readRDS(here::here(sprintf("misc/results/case_control/ipd_case_control_no_abx_9_12_lsd_month_%d.Rds", mo)))
  })
  
  # Load all the RDS files into a list
  mo_list_12_15 <- lapply(months, function(mo) {
    readRDS(here::here(sprintf("misc/results/case_control/ipd_case_control_no_abx_12_15_month_%d.Rds", mo)))
  })
  
  mo_list_severe_12_15 <- lapply(months, function(mo) {
    readRDS(here::here(sprintf("misc/results/case_control/ipd_case_control_no_abx_12_15_msd_month_%d.Rds", mo)))
  })
  
  mo_list_mild_12_15 <- lapply(months, function(mo) {
    readRDS(here::here(sprintf("misc/results/case_control/ipd_case_control_no_abx_12_15_lsd_month_%d.Rds", mo)))
  })
  
  # Helper to extract effect estimates and CIs
  extract_plot_data <- function(mo_list, label) {
    data.frame(
      pt_est = sapply(mo_list, function(mo) mo$results_object$results_df$effect_inf_abx_level),
      lower_ci = sapply(mo_list, function(mo) {
        est <- mo$results_object$results_df$effect_inf_abx_level
        se <- mo$results_object$se['effect_observed']
        est - 1.96 * se
      }),
      upper_ci = sapply(mo_list, function(mo) {
        est <- mo$results_object$results_df$effect_inf_abx_level
        se <- mo$results_object$se['effect_observed']
        est + 1.96 * se
      }),
      month_num = months,
      group = label
    )
  }
  
  # Process data for 6–9 months
  plot_6_9_any   <- extract_plot_data(mo_list_6_9, "Any Shigella")
  plot_6_9_sev   <- extract_plot_data(mo_list_severe_6_9, "Moderate-to-Severe Shigella")
  plot_6_9_mild  <- extract_plot_data(mo_list_mild_6_9, "Less-severe Shigella")
  combined_6_9   <- bind_rows(plot_6_9_mild, plot_6_9_any, plot_6_9_sev)
  combined_6_9$age_strata <- "6–9 months"
  
  # Process data for 9–12 months
  plot_9_12_any   <- extract_plot_data(mo_list_9_12, "Any Shigella")
  plot_9_12_sev   <- extract_plot_data(mo_list_severe_9_12, "Moderate-to-Severe Shigella")
  plot_9_12_mild  <- extract_plot_data(mo_list_mild_9_12, "Less-severe Shigella")
  combined_9_12   <- bind_rows(plot_9_12_mild, plot_9_12_any, plot_9_12_sev)
  combined_9_12$age_strata <- "9–12 months"
  
  # Process data for 12–15 months
  plot_12_15_any   <- extract_plot_data(mo_list_12_15, "Any Shigella")
  plot_12_15_sev   <- extract_plot_data(mo_list_severe_12_15, "Moderate-to-Severe Shigella")
  plot_12_15_mild  <- extract_plot_data(mo_list_mild_12_15, "Less-severe Shigella")
  combined_12_15   <- bind_rows(plot_12_15_mild, plot_12_15_any, plot_12_15_sev)
  combined_12_15$age_strata <- "12–15 months"
  
  # Combine all
  plot_df <- bind_rows(combined_6_9, combined_9_12, combined_12_15)
  plot_df$group <- factor(plot_df$group, levels = c("Less-severe Shigella", "Any Shigella", "Moderate-to-Severe Shigella"))
  plot_df$age_strata <- factor(plot_df$age_strata, levels = c("6–9 months", "9–12 months", "12–15 months"))
  
  # Create faceted plot
  facet_plot <- ggplot(plot_df, aes(x = month_num, y = pt_est, color = group, fill = group)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3, color = NA) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    stat_smooth(
      method = "glm",
      formula = y ~ -1 + x + I(pmax(0, x - 4)) + I(pmax(0, x - 8)),
      size = 1.2,
      linetype = "dotted",
      se = FALSE
    ) +
    scale_x_continuous(breaks = months) +
    facet_wrap(~ age_strata, nrow = 1) +
    labs(
      x = "Month",
      y = "Growth Effect",
      title = "Shigella growth effects in MAL-ED with observed antibiotic use by age",
      color = "Shigella severity",
      fill = "Shigella severity"
    ) +
    scale_color_manual(values = c(
      "Any Shigella" = "#00468BFF",
      "Moderate-to-Severe Shigella" = "#ED0000FF",
      "Less-severe Shigella" = "#42B540FF"
    )) +
    scale_fill_manual(values = c(
      "Any Shigella" = "#00468BFF",
      "Moderate-to-Severe Shigella" = "#ED0000FF",
      "Less-severe Shigella" = "#42B540FF"
    )) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
 
  