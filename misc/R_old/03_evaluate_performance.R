
# ---------------------------------------------------
# Function to get bias
# ---------------------------------------------------
get_bias <- function(results, row_eval, est){
  sub_res <- results[results$n == row_eval$n &
                       results$VE_mild == row_eval$VE_mild &
                       results$VE_severe == row_eval$VE_severe &
                       results$incidence_shigella == row_eval$incidence_shigella &
                       results$incidence_severe_shigella == row_eval$incidence_severe_shigella &
                       results$lfaz_bl_effect_shig == row_eval$lfaz_bl_effect_shig &
                       results$lfaz_bl_effect_severe_shig == row_eval$lfaz_bl_effect_severe_shig &
                       results$shig_coef_mild == row_eval$shig_coef_mild &
                       results$shig_coef_severe == row_eval$shig_coef_severe &
                       results$sd == row_eval$sd &
                       results$Yinf_X_model == row_eval$Yinf_X_model &
                       results$G_X_model == row_eval$G_X_model,]
  out <- list()
  
  for(method in est){
    if("gcomp" == method){
      out$bias_gcomp <- mean(sub_res$pt_est[sub_res$class == "gcomp_res"] - row_eval$true_growth_effect)
    }
    
    if("efficient_aipw" == method){
      out$bias_aipw <- mean(sub_res$pt_est[sub_res$class == "aipw_res"] - row_eval$true_growth_effect)
    }
    
    if("efficient_tmle" == method){
      out$bias_tmle <- mean(sub_res$pt_est[sub_res$class == "tmle_res"] - row_eval$true_growth_effect)
    }
    
    if("gcomp_pop_estimand" == method){
      out$bias_pop_comp <- mean(sub_res$pt_est[sub_res$class == "pop_gcomp_res"] - row_eval$true_pop_growth_effect)
    }
  }
  
  return(out)
}

# ---------------------------------------------------
# Function to get 95% CI coverage
# ---------------------------------------------------
get_coverage <- function(results, row_eval, est){
  sub_res <- results[results$n == row_eval$n &
                       results$VE_mild == row_eval$VE_mild &
                       results$VE_severe == row_eval$VE_severe &
                       results$incidence_shigella == row_eval$incidence_shigella &
                       results$incidence_severe_shigella == row_eval$incidence_severe_shigella &
                       results$lfaz_bl_effect_shig == row_eval$lfaz_bl_effect_shig &
                       results$lfaz_bl_effect_severe_shig == row_eval$lfaz_bl_effect_severe_shig &
                       results$shig_coef_mild == row_eval$shig_coef_mild &
                       results$shig_coef_severe == row_eval$shig_coef_severe &
                       results$sd == row_eval$sd &
                       results$Yinf_X_model == row_eval$Yinf_X_model &
                       results$G_X_model == row_eval$G_X_model,]
  
  out <- list()
  
  for(method in est){
    if("gcomp" == method){
      coverage_vec_ge <- ifelse(sub_res$lower_ci[sub_res$class == "gcomp_res"] < row_eval$true_growth_effect & sub_res$upper_ci[sub_res$class == "gcomp_res"] > row_eval$true_growth_effect, 1, 0)
      out$coverage_gcomp <- mean(coverage_vec_ge)
    }
    
    if("efficient_aipw" == method){
      coverage_vec_aipw <- ifelse(sub_res$lower_ci[sub_res$class == "aipw_res"] < row_eval$true_growth_effect & sub_res$upper_ci[sub_res$class == "aipw_res"] > row_eval$true_growth_effect, 1, 0)
      out$coverage_aipw <- mean(coverage_vec_aipw)
    }
    
    if("efficient_tmle" == method){
      coverage_vec_tmle <- ifelse(sub_res$lower_ci[sub_res$class == "tmle_res"]  < row_eval$true_growth_effect & sub_res$upper_ci[sub_res$class == "tmle_res"] > row_eval$true_growth_effect, 1, 0)
      out$coverage_tmle <- mean(coverage_vec_tmle)
    }
    
    if("gcomp_pop_estimand" == method){
      coverage_vec_pop <- ifelse(sub_res$lower_ci[sub_res$class == "pop_gcomp_res"] < row_eval$true_pop_growth_effect & sub_res$upper_ci[sub_res$class == "pop_gcomp_res"] > row_eval$true_pop_growth_effect, 1, 0)
      out$coverage_pop_comp <- mean(coverage_vec_pop)
    }
  }
  
  return(out)
}

# --------------------------------------------------
# Calculate power 
# --------------------------------------------------
get_power <- function(results, row_eval, est){
  sub_res <- results[results$n == row_eval$n &
                       results$VE_mild == row_eval$VE_mild &
                       results$VE_severe == row_eval$VE_severe &
                       results$incidence_shigella == row_eval$incidence_shigella &
                       results$incidence_severe_shigella == row_eval$incidence_severe_shigella &
                       results$lfaz_bl_effect_shig == row_eval$lfaz_bl_effect_shig &
                       results$lfaz_bl_effect_severe_shig == row_eval$lfaz_bl_effect_severe_shig &
                       results$shig_coef_mild == row_eval$shig_coef_mild &
                       results$shig_coef_severe == row_eval$shig_coef_severe &
                       results$sd == row_eval$sd &
                       results$Yinf_X_model == row_eval$Yinf_X_model &
                       results$G_X_model == row_eval$G_X_model,]
  
  out <- list()
  
  for(method in est){
    if("gcomp" == method){
      sub_res$reject_ge_bin <- ifelse(sub_res$reject[sub_res$class == "gcomp_res"] == TRUE, 1, 0)
      out$power_gcomp <- mean(sub_res$reject_ge_bin)
    }
    
    if("efficient_aipw" == method){
      sub_res$reject_ge_aipw_bin <- ifelse(sub_res$reject[sub_res$class == "aipw_res"] == TRUE, 1, 0)
      out$power_aipw <- mean(sub_res$reject_ge_aipw_bin)
    } 
    
    if("efficient_tmle" == method){
      sub_res$reject_ge_tmle_bin <- ifelse(sub_res$reject[sub_res$class == "tmle_res"] == TRUE, 1, 0)
      out$power_tmle <- mean(sub_res$reject_ge_tmle_bin)
    }
    
    if("gcomp_pop_estimand" == method){
      sub_res$reject_ge_pop_bin <- ifelse(sub_res$reject[sub_res$class == "pop_gcomp_res"] == TRUE, 1, 0)
      out$power_pop_gcomp <- mean(sub_res$reject_ge_pop_bin)
    }
    
    if("choplump" == method){
      sub_res$reject_ge_chop_lump_bin <- ifelse(sub_res$reject[sub_res$class == "choplump_res"] == TRUE, 1, 0)
      out$power_chop_lump <- mean(sub_res$reject_ge_chop_lump_bin)
    }
    
    if("hudgens_lower" == method){
      sub_res$reject_hudgens_lower_bin <- ifelse(sub_res$reject[sub_res$class == "hudgens_lower_res"] == TRUE, 1, 0)
      out$power_hudgens_lower <- mean(sub_res$reject_hudgens_lower_bin)
    }
    
    if("hudgens_upper" == method){
      sub_res$reject_hudgens_upper_bin <- ifelse(sub_res$reject[sub_res$class == "hudgens_upper_res"] == TRUE, 1, 0)
      out$power_hudgens_upper <- mean(sub_res$reject_hudgens_upper_bin)
    }
    
    if("hudgens_lower_doomed" == method){
      sub_res$reject_hudgens_lower_doomed_bin <- ifelse(sub_res$reject[sub_res$class == "hudgens_lower_res_doomed"] == TRUE, 1, 0)
      out$power_hudgens_lower_doomed <- mean(sub_res$reject_hudgens_lower_bin)
    }
    
    if("hudgens_upper_doomed" == method){
      sub_res$reject_hudgens_upper__doomedbin <- ifelse(sub_res$reject[sub_res$class == "hudgens_upper_res_doomed"] == TRUE, 1, 0)
      out$power_hudgens_upper_doomed <- mean(sub_res$reject_hudgens_upper_bin)
    }
    
    if("hudgens_adj_lower" == method){
      sub_res$reject_hudgens_adj_lower_bin <- ifelse(sub_res$reject[sub_res$class == "hudgens_adj_lower_res"] == TRUE, 1, 0)
      out$power_hudgens_lower <- mean(sub_res$reject_hudgens_adj_lower_bin)
    }
    
    if("hudgens_adj_upper" == method){
      sub_res$reject_hudgens_adj_upper_bin <- ifelse(sub_res$reject[sub_res$class == "hudgens_adj_upper_res"] == TRUE, 1, 0)
      out$power_hudgens_upper <- mean(sub_res$reject_hudgens_adj_upper_bin)
    }
    
  }
  
  return(out)
  
}

# --------------------------------------------------
# Function to compare mean of bootstrap se estimates to sd of pt_ests
# --------------------------------------------------
get_se_compare <- function(results, row_eval, est){
  sub_res <- results[results$n == row_eval$n &
                       results$VE_mild == row_eval$VE_mild &
                       results$VE_severe == row_eval$VE_severe &
                       results$incidence_shigella == row_eval$incidence_shigella &
                       results$incidence_severe_shigella == row_eval$incidence_severe_shigella &
                       results$lfaz_bl_effect_shig == row_eval$lfaz_bl_effect_shig &
                       results$lfaz_bl_effect_severe_shig == row_eval$lfaz_bl_effect_severe_shig &
                       results$shig_coef_mild == row_eval$shig_coef_mild &
                       results$shig_coef_severe == row_eval$shig_coef_severe &
                       results$sd == row_eval$sd &
                       results$Yinf_X_model == row_eval$Yinf_X_model &
                       results$G_X_model == row_eval$G_X_model,]
  
  out <- list()
  
  for(method in est){
    if("gcomp" == method){
      mean_bootstrap_se_ge <- mean(sub_res$se[sub_res$class == "gcomp_res"])
      sd_pt_est_ge <- sd(sub_res$pt_est[sub_res$class == "gcomp_res"])
      out$ratio_gcomp <- mean_bootstrap_se_ge / sd_pt_est_ge
    }
    
    if("efficient_aipw" == method){
      mean_bootstrap_se_aipw <- mean(sub_res$se[sub_res$class == "aipw_res"])
      sd_pt_est_ge_aipw <- sd(sub_res$pt_est[sub_res$class == "aipw_res"])
      out$ratio_aipw <- mean_bootstrap_se_aipw / sd_pt_est_ge_aipw
    }
    
    if("efficient_tmle" == method){
      mean_bootstrap_se_tmle <- mean(sub_res$se[sub_res$class == "tmle_res"])
      sd_pt_est_ge_tmle <- sd(sub_res$pt_est[sub_res$class == "tmle_res"])
      out$ratio_tmle <- mean_bootstrap_se_tmle / sd_pt_est_ge_tmle
    }
    
    if("gcomp_pop_estimand" == method){
      mean_bootstrap_se_pop <- mean(sub_res$se[sub_res$class == "pop_gcomp_res"])
      sd_pt_est_ge_pop <- sd(sub_res$pt_est[sub_res$class == "pop_gcomp_res"])
      out$ratio_pop_gcomp <- mean_bootstrap_se_pop / sd_pt_est_ge_pop
    }
    
    if("hudgens_adj_lower" == method){
      mean_bootstrap_se_hudgens_adj_lower <- mean(sub_res$se[sub_res$class == "hudgens_adj_lower_res"])
      sd_pt_est_hudgens_adj_lower <- sd(sub_res$pt_est[sub_res$class == "hudgens_adj_lower_res"])
      out$ratio_hudgens_adj_lower <- mean_bootstrap_se_hudgens_adj_lower / sd_pt_est_hudgens_adj_lower
    }
    
    if("hudgens_adj_upper" == method){
      mean_bootstrap_se_hudgens_adj_upper<- mean(sub_res$se[sub_res$class == "hudgens_adj_upper_res"])
      sd_pt_est_hudgens_adj_upper <- sd(sub_res$pt_est[sub_res$class == "hudgens_adj_upper_res"])
      out$ratio_hudgens_adj_upper <- mean_bootstrap_se_hudgens_adj_upper / sd_pt_est_hudgens_adj_upper
    }
  }
  
  return(out)
  
}

# -------------------------------------------------------
# Function to get proportion of negative point estimates
# -------------------------------------------------------
get_neg_pt_est <- function(results, row_eval, est){
  sub_res <- results[results$n == row_eval$n &
                       results$VE_mild == row_eval$VE_mild &
                       results$VE_severe == row_eval$VE_severe &
                       results$incidence_shigella == row_eval$incidence_shigella &
                       results$incidence_severe_shigella == row_eval$incidence_severe_shigella &
                       results$lfaz_bl_effect_shig == row_eval$lfaz_bl_effect_shig &
                       results$lfaz_bl_effect_severe_shig == row_eval$lfaz_bl_effect_severe_shig &
                       results$shig_coef_mild == row_eval$shig_coef_mild &
                       results$shig_coef_severe == row_eval$shig_coef_severe &
                       results$sd == row_eval$sd &
                       results$Yinf_X_model == row_eval$Yinf_X_model &
                       results$G_X_model == row_eval$G_X_model,]
  
  out <- list()
  
  for(method in est){
    if("gcomp" == method){
      out$prop_neg_gcomp <- mean(ifelse(sub_res$pt_est[sub_res$class == "gcomp_res"] < 0, 1, 0))
    }
    
    if("efficient_aipw" == method){
      out$prop_neg_aipw <- mean(ifelse(sub_res$pt_est[sub_res$class == "aipw_res"] < 0, 1, 0))
    }
    
    if("efficient_tmle" == method){
      out$prop_neg_tmle <- mean(ifelse(sub_res$pt_est[sub_res$class == "tmle_res"] < 0, 1, 0))
    }
    
    if("gcomp_pop_estimand" == method){
      out$prop_neg_pop_gcomp <- mean(ifelse(sub_res$pt_est[sub_res$class == "pop_gcomp_res"] < 0, 1, 0))
    }
  }
  
  return(out)
  
}

# -------------------------------------------------------
# Function to get average CI width
# -------------------------------------------------------
get_ci_width <- function(results, row_eval, est){
  sub_res <- results[results$n == row_eval$n &
                       results$VE_mild == row_eval$VE_mild &
                       results$VE_severe == row_eval$VE_severe &
                       results$incidence_shigella == row_eval$incidence_shigella &
                       results$incidence_severe_shigella == row_eval$incidence_severe_shigella &
                       results$lfaz_bl_effect_shig == row_eval$lfaz_bl_effect_shig &
                       results$lfaz_bl_effect_severe_shig == row_eval$lfaz_bl_effect_severe_shig &
                       results$shig_coef_mild == row_eval$shig_coef_mild &
                       results$shig_coef_severe == row_eval$shig_coef_severe &
                       results$sd == row_eval$sd &
                       results$Yinf_X_model == row_eval$Yinf_X_model &
                       results$G_X_model == row_eval$G_X_model,]
  
  out <- list()
  
  for(method in est){
    if("gcomp" == method){
      out$ci_width_gcomp <- mean(sub_res$upper_ci[sub_res$class == "gcomp_res"] - sub_res$lower_ci[sub_res$class == "gcomp_res"])
    }
    
    if("efficient_aipw" == method){
      out$ci_width_aipw <- mean(sub_res$upper_ci[sub_res$class == "aipw_res"] - sub_res$lower_ci[sub_res$class == "aipw_res"])
    }
    
    if("efficient_tmle" == method){
      out$ci_width_tmle <- mean(sub_res$upper_ci[sub_res$class == "tmle_res"] - sub_res$lower_ci[sub_res$class == "tmle_res"])
    }
    
    if("gcomp_pop_estimand" == method){
      out$ci_width_pop <- mean(sub_res$upper_ci[sub_res$class == "pop_gcomp_res"] - sub_res$lower_ci[sub_res$class == "pop_gcomp_res"])
    }
    
    if("hudgens_adj_lower" == method){
      out$ci_width_hudgens_adj_lower <- mean(sub_res$upper_ci[sub_res$class == "hudgens_adj_lower_res"] - sub_res$lower_ci[sub_res$class == "hudgens_adj_lower_res"])
    }
    
    if("hudgens_adj_upper" == method){
      out$ci_width_hudgens_adj_upper <- mean(sub_res$upper_ci[sub_res$class == "hudgens_adj_upper_res"] - sub_res$lower_ci[sub_res$class == "hudgens_adj_upper_res"])
    }
  }
  
  return(out)
  
}

# ---------------------------------------------------
# Function to plot power vs sample size curve
# ---------------------------------------------------
plot_power <- function(power_df){
  
  by_x <- ifelse(max(power_df$n) > 200000, 100000, 
                 ifelse(max(power_df$n) < 100000, 10000, 50000))
  
  # Get columns that start with `power_`
  power_cols <- colnames(power_df)[grepl("^power_", colnames(power_df))]
  
  for(i in 1:length(power_cols)){
    power_df[power_cols[i]] <- round(power_df[power_cols[i]], 5)
  }
  
  # Round truth values 
  true_cols <- colnames(power_df)[grepl("^true_", colnames(power_df))]
  for(i in 1:length(true_cols)){
    power_df[true_cols[i]] <- round(power_df[true_cols[i]], 5)
  }
  
  # Melt the dataframe to long format for ggplot
  power_long <- power_df %>%
    tidyr::pivot_longer(cols = all_of(power_cols), names_to = "method", values_to = "power")
  
  # Dynamically generate a color palette based on the number of methods
  color_palette <- scales::hue_pal()(length(power_cols))
  
  # Custom labels for methods in the legend
  method_labels <- c(
    "power_gcomp_pop_estimand" = "G-Comp: Population Estimand",
    "power_gcomp" = "G-Comp: VE Estimand",
    "power_choplump" = "Chop-Lump method",
    "power_hudgens_upper" = "Hudgens style method: upper bound",
    "power_hudgens_lower" = "Hudgens style method: lower bound"
  )
  
  # Custom labels for the facets
  facet_labels <- labeller(
    VE_mild = label_both,
    VE_severe = label_both,
    incidence_shigella = label_both,
    incidence_severe_shigella = label_both,
    lfaz_bl_effect_shig = label_both,
    lfaz_bl_effect_severe_shig = label_both,
    shig_coef_mild = label_both,
    shig_coef_severe = label_both,
    true_growth_effect = label_both,
    true_pop_growth_effect = label_both,
    sd = label_both
  )
  
  # Generate the plot with dynamic lines for each power_* column
  power_plot <- ggplot(power_long, aes(x = n, y = power, color = method)) +
    geom_line() +
    geom_point() +
    facet_wrap(vars(VE_mild,
                    VE_severe,
                    incidence_shigella,
                    incidence_severe_shigella,
                    lfaz_bl_effect_shig,
                    lfaz_bl_effect_severe_shig,
                    shig_coef_mild,
                    shig_coef_severe,
                    sd,
                    true_growth_effect,
                    true_pop_growth_effect),
               labeller = facet_labels) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
    scale_x_continuous(breaks = seq(0, max(power_df$n), by = by_x), 
                       labels = scales::comma) +  
    # Dynamically apply colors and custom labels
    scale_color_manual(values = setNames(color_palette, power_cols),
                       labels = method_labels) +
    labs(title = "Power vs sample size curve for Shigella growth effect",
         x = "Sample size",
         y = "Power",
         color = "Method")
  
  return(power_plot)
  
}

plot_width <- function(ci_width_df){
  
  by_x <- ifelse(max(ci_width_df$n) > 200000, 100000, 
                 ifelse(max(ci_width_df$n) < 100000, 10000, 50000))
  max_y <- ifelse(max(ci_width_df$ci_width_gcomp) > 0.6, 0.35, max(ci_width_df$ci_width_gcomp))
  
  ci_width_df$true_growth_effect <- round(ci_width_df$true_growth_effect, 5)
  ci_width_df$true_pop_growth_effect <- round(ci_width_df$true_pop_growth_effect, 5)
  
  # Custom labels for the facets
  facet_labels <- labeller(
    VE_mild = label_both,
    VE_severe = label_both,
    incidence_shigella = label_both,
    incidence_severe_shigella = label_both,
    lfaz_bl_effect_shig = label_both,
    lfaz_bl_effect_severe_shig = label_both,
    shig_coef_mild = label_both,
    shig_coef_severe = label_both,
    true_growth_effect = label_both,
    true_pop_growth_effect = label_both,
    sd = label_both
  )
  
  ci_plot <- ggplot(ci_width_df, aes(x = n)) +
    geom_line(aes(y = ci_width_gcomp_pop_estimand, color = "gcomp_pop_estimand")) +
    geom_line(aes(y = ci_width_gcomp, color = "gcomp")) +
    geom_point(aes(y = ci_width_gcomp_pop_estimand, color = "gcomp_pop_estimand")) +
    geom_point(aes(y = ci_width_gcomp, color = "gcomp")) +
    facet_wrap(vars(VE_mild,
                    VE_severe,
                    incidence_shigella,
                    incidence_severe_shigella,
                    lfaz_bl_effect_shig,
                    lfaz_bl_effect_severe_shig,
                    shig_coef_mild,
                    shig_coef_severe,
                    sd,
                    true_growth_effect,
                    true_pop_growth_effect),
               labeller = facet_labels) +
    scale_x_continuous(breaks = seq(0, max(ci_width_df$n), by = by_x), 
                       labels = scales::comma) +  
    coord_cartesian(ylim = c(0, max_y)) + 
    scale_color_manual(values = c("gcomp" = "steelblue", 
                                  "gcomp_pop_estimand" = "violetred1"),
                       labels = c("gcomp" = "G-Comp: VE Estimand", 
                                  "gcomp_pop_estimand" = "G-Comp: Population Estimand")) +
    labs(title = "Confidence interval width vs sample size curve for Shigella growth effect",
         x = "Sample size",
         y = "CI Width",
         color = "Method")
  
  return(ci_plot)
}

plot_prop_neg <- function(prop_neg_df){
  
  by_x <- ifelse(max(prop_neg_df$n) > 200000, 100000, 
                 ifelse(max(prop_neg_df$n) < 100000, 10000, 50000))
  
  prop_neg_df$true_growth_effect <- round(prop_neg_df$true_growth_effect, 5)
  prop_neg_df$true_pop_growth_effect <- round(prop_neg_df$true_pop_growth_effect, 5)
  
  # Custom labels for the facets
  facet_labels <- labeller(
    VE_mild = label_both,
    VE_severe = label_both,
    incidence_shigella = label_both,
    incidence_severe_shigella = label_both,
    lfaz_bl_effect_shig = label_both,
    lfaz_bl_effect_severe_shig = label_both,
    shig_coef_mild = label_both,
    shig_coef_severe = label_both,
    sd = label_both,
    true_growth_effect = label_both,
    true_pop_growth_effect = label_both
  )
  
  prop_neg_plot <- ggplot(prop_neg_df, aes(x = n)) +
    geom_line(aes(y = prop_neg_gcomp_pop_estimand, color = "gcomp_pop_estimand")) +
    geom_line(aes(y = prop_neg_gcomp, color = "gcomp")) +
    geom_point(aes(y = prop_neg_gcomp_pop_estimand, color = "gcomp_pop_estimand")) +
    geom_point(aes(y = prop_neg_gcomp, color = "gcomp")) +
    facet_wrap(vars(VE_mild,
                    VE_severe,
                    incidence_shigella,
                    incidence_severe_shigella,
                    lfaz_bl_effect_shig,
                    lfaz_bl_effect_severe_shig,
                    shig_coef_mild,
                    shig_coef_severe,
                    sd,
                    true_growth_effect,
                    true_pop_growth_effect),
               labeller = facet_labels) +
    scale_x_continuous(breaks = seq(0, max(prop_neg_df$n), by = by_x), 
                       labels = scales::comma) +  
    scale_color_manual(values = c("gcomp" = "steelblue", 
                                  "gcomp_pop_estimand" = "violetred1"),
                       labels = c("gcomp" = "G-Comp: VE Estimand", 
                                  "gcomp_pop_estimand" = "G-Comp: Population Estimand")) +
    labs(title = "Proportion of negative point estimates vs sample size curve for Shigella growth effect",
         x = "Sample size",
         y = "Proportion of negative point estimates",
         color = "Method")
  
  return(prop_neg_plot)
}



