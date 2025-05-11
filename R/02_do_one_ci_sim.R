do_one_ci_sim <- function(
    params, n_boot = 1000, 
    est = c("gcomp_pop_estimand", "gcomp", "efficient_gcomp", 
            "efficient_aipw", "efficient_tmle", "choplump"),
    null_hypothesis_value = 0,
    alpha_level = 0.025
){
  set.seed(params$seed)
  
  data <- generate_EFGH_complex(n = params$n,
                                incidence_shigella = params$incidence_shigella,
                                incidence_severe_shigella = params$incidence_severe_shigella,
                                VE_mild = params$VE_mild,
                                VE_severe = params$VE_severe,
                                lfaz_bl_effect_shig = params$lfaz_bl_effect_shig,
                                lfaz_bl_effect_severe_shig = params$lfaz_bl_effect_severe_shig,
                                shig_coef_mild = params$shig_coef_mild,
                                shig_coef_severe = params$shig_coef_severe,
                                sd_growth = params$sd,
                                zhifei = params$zhifei,
                                catch_up = params$catch_up)
  
  # data <- generate_EFGH_complex_quad(n = params$n,
  #                                    incidence_shigella = params$incidence_shigella,
  #                                    incidence_severe_shigella = params$incidence_severe_shigella,
  #                                    VE_mild = params$VE_mild,
  #                                    VE_severe = params$VE_severe,
  #                                    lfaz_bl_effect_shig = params$lfaz_bl_effect_shig,
  #                                    lfaz_bl_effect_severe_shig = params$lfaz_bl_effect_severe_shig,
  #                                    shig_coef_mild = params$shig_coef_mild,
  #                                    shig_coef_severe = params$shig_coef_severe,
  #                                    sd_growth = params$sd,
  #                                    zhifei = params$zhifei,
  #                                    catch_up = params$catch_up)
  
  # results = object of class "vegrowth" = list with entry for each element of est
  results <- vegrowth(data =  data,
                      G_name = "G",
                      V_name = "V",
                      X_name = "X",
                      Y_name = "Y_inf",
                      est = est, 
                      n_boot = n_boot,
                      seed = params$seed, 
                      #G_X_model = params$G_X_model, 
                      #Y_X_model = params$Yinf_X_model, 
                      null_hypothesis_value = null_hypothesis_value, 
                      alpha_level = alpha_level, 
                      family = "gaussian", 
                      return_models = FALSE)
  
  # Unlist to get in old format?
  results_df <- lapply(results, function(x) {
    if(class(x) %in% c("choplump_res", "hudgens_lower_res", "hudgens_upper_res",
                       "hudgens_lower_res_doomed", "hudgens_upper_res_doomed")){
      df <- data.frame(
        pt_est = x$obs_diff,
        se = x$se,
        lower_ci = x$lower_ci,
        upper_ci = x$upper_ci,
        reject = x$reject,
        pval = x$pval,
        class = class(x),
        stringsAsFactors = FALSE
      )
    } else{
      df <- data.frame(
        pt_est = x$pt_est,
        se = x$se,
        lower_ci = unname(x$lower_ci),
        upper_ci = unname(x$upper_ci),
        reject = x$reject,
        pval = NA,
        class = class(x),
        stringsAsFactors = FALSE
      )
    }
    return(df)
  })
  
  results_df <- do.call(rbind, results_df)
  results_df <- cbind(results_df, params)
  
  return(results_df)
  
}