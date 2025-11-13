# ------------------------------------------------------------------
# Shigella vaccine trial plots - power, proportion negative
# ------------------------------------------------------------------

library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(patchwork)

here::i_am("plot_results.R")

setting <- "default"

results <- readRDS(here::here(paste0("results/", setting, "_evaluation_results.Rds")))

power_df <- results$power_df 
prop_neg_df <- results$prop_neg_df 
bias_df <- results$bias_df 

# switch to make estimator = shape (dot for AIPW, triangle for G-Comp, square for unadjusted), and make the colors = to the different er / cw combinations
# -------------------------------------
# Helper to generate readable labels
# -------------------------------------
make_color_label <- function(estimator, estimand, er, cw, two_stage) {
  
  label <- ifelse(
    is.na(er) & is.na(cw) & is.na(two_stage), "(None)", # this is the unadjusted version
      ifelse(estimand == "pop", "(None)", 
        ifelse(er & cw & two_stage, "(ER, CW, two-part)",
               ifelse(er & cw & !two_stage, "(ER, CW)",
                      ifelse(er & !cw & two_stage, "(ER, two-part)",
                             ifelse(er & !cw & !two_stage, "(ER)",
                                    ifelse(!er & cw & two_stage, "(CW, two-part)",
                                           ifelse(!er & cw & !two_stage, "(CW)", ""))))))))
  
  # # Estimator names
  # estimator <- ifelse(estimator == "gcomp", "G-Comp",
  #                     ifelse(estimator == "aipw", "AIPW", "Unadjusted"))
  # 
  # # Suffix for adjustments
  # suffix <- ifelse(
  #   is.na(er) & is.na(cw) & is.na(two_stage), "",
  #   ifelse(er & cw & two_stage, "(ER, CW, two-part)",
  #          ifelse(er & cw & !two_stage, "(ER, CW)",
  #                 ifelse(er & !cw & two_stage, "(ER, two-part)",
  #                        ifelse(er & !cw & !two_stage, "(ER)",
  #                               ifelse(!er & cw & two_stage, "(CW, two-part)",
  #                                      ifelse(!er & cw & !two_stage, "(CW)", "")))))))
  
  return(label)
  #paste(estimator, suffix)
}

# -------------------------------------
# Add labels for plotting
# -------------------------------------
add_labels <- function(df) {
  df %>%
    mutate(
      color_label = make_color_label(estimator, estimand, er, cw, two_stage),
      line_type = ifelse(estimand == "nat_inf", "Naturally Infected", "Population"),
      shape_type = ifelse(estimator == "gcomp", "G-Comp", ifelse(estimator == "aipw", "AIPW", "Unadjusted")), 
      Y_out = recode(
        Y_out,
        Y_3 = "Y[3]",
        Y_6 = "Y[6]",
        Y_9 = "Y[9]",
        Y_12 = "Y[12]",
        Y_6_12 = "Y[6-12]",
        Y_3_6_9_12 = "Y[3-6-9-12]"
      ),
      Y_out = factor(
        Y_out,
        levels = c("Y[3]", "Y[6]", "Y[9]", "Y[12]", "Y[6-12]", "Y[3-6-9-12]")
      )
    )
}

power_df <- add_labels(power_df)
prop_neg_df <- add_labels(prop_neg_df)
bias_df <- add_labels(bias_df)

# -------------------------------------
# Function to create power plot
# -------------------------------------
# -------------------------------------
# Function to create power plot with shapes
# -------------------------------------
plot_power <- function(power_df, setting = "default") {
  
  color_map <- c(
    "(CW)" = "#89CFF0",
    "(ER)" = "#4682B4",          
    "(ER, two-part)" = "#6495ED",
    "(ER, CW)" =  "#191970",
    
    # "(CW)" = "#D8BFD8",
    # "(ER)" = "#DE3163",
    # "(ER, two-part)" = "#E35335",
    # "(ER, CW)" = "#AD002AFF",
    
    "(None)" = "#FFC000" # Unadjusted
  )
  
  shape_map <- c(
    "G-Comp" = 17,     # triangle
    "AIPW" = 16,       # circle
    "Unadjusted" = 15  # square
  )
  
  # extract only the suffix for color mapping
  power_df <- power_df %>%
    mutate(
      color_only = color_label,  # already contains "(ER, CW, ...)" or "(None)"
      shape_only = shape_type
    )
  
  ggplot(power_df, aes(x = n, y = power, color = color_only, 
                       shape = shape_only, linetype = line_type)) +
    geom_line(aes(group = interaction(color_only, line_type)), linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 4, alpha = 0.7) +
    facet_wrap(~ Y_out, scales = "fixed", ncol = 3, labeller = label_parsed) +
    scale_y_continuous(
      limits = c(0, 0.2),
      breaks = seq(0, 0.2, by = 0.05),
      labels = scales::percent_format(accuracy = 1),
      name = "Power"
    ) +
    scale_x_continuous(
      name = "Sample Size",
      labels = scales::comma,
      breaks = scales::pretty_breaks(n = 5)
    ) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_map) +
    scale_linetype_manual(values = c("Naturally Infected" = "solid", "Population" = "dotted")) +
    labs(
      title = paste0("Simulated Shigella vaccine trial power (", setting, " setting)"),
      color = "Assumptions & Modeling Approach",
      shape = "Estimator",
      linetype = "Estimand"
    ) +
    theme_minimal(base_size = 18) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      strip.text = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}


# -------------------------------------------------
# Function to create proportion negative plot
# -------------------------------------------------
plot_prop_neg <- function(prop_neg_df, setting = "default") {
  
  prop_neg_df <- prop_neg_df %>%
    mutate(
      color_label = make_color_label(estimator, estimand, er, cw, two_stage),
      shape_type = ifelse(estimator == "gcomp", "G-Comp",
                          ifelse(estimator == "aipw", "AIPW", "Unadjusted")),
      line_type = ifelse(estimand == "nat_inf", "Naturally Infected", "Population"),
      Y_out = recode(
        Y_out,
        Y_3 = "Y[3]",
        Y_6 = "Y[6]",
        Y_9 = "Y[9]",
        Y_12 = "Y[12]",
        Y_6_12 = "Y[6-12]",
        Y_3_6_9_12 = "Y[3-6-9-12]"
      ),
      Y_out = factor(
        Y_out,
        levels = c("Y[3]", "Y[6]", "Y[9]", "Y[12]", "Y[6-12]", "Y[3-6-9-12]")
      )
    )
  
  color_map <- c(
    "(CW)" = "#89CFF0",
    "(ER)" = "#4682B4",
    "(ER, two-part)" = "#6495ED",
    "(ER, CW)" = "#191970",
    "(None)" = "#FFC000"
  )
  
  shape_map <- c(
    "G-Comp" = 17,
    "AIPW" = 16,
    "Unadjusted" = 15
  )
  
  ggplot(prop_neg_df, aes(x = n, y = prop_neg, color = color_label, 
                          shape = shape_type, linetype = line_type)) +
    geom_line(aes(group = interaction(color_label, line_type)), linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 4, alpha = 0.7) +
    facet_wrap(~ Y_out, scales = "fixed", ncol = 3, labeller = label_parsed) +
    scale_y_continuous(
      name = "Proportion Negative",
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_x_continuous(
      name = "Sample Size",
      labels = scales::comma,
      limits = c(0, 50000),
      breaks = seq(0, 50000, 10000)
    ) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_map) +
    scale_linetype_manual(values = c("Naturally Infected" = "solid", "Population" = "dotted")) +
    labs(
      title = paste0("Proportion Negative (", setting, " setting)"),
      color = "Assumptions & Modeling Approach",
      shape = "Estimator",
      linetype = "Estimand"
    ) +
    theme_minimal(base_size = 18) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

# -------------------------------------------------
# Function to create bias plot
# -------------------------------------------------
plot_bias <- function(bias_df, setting = "default") {
  
  bias_df <- bias_df %>%
    mutate(
      color_label = make_color_label(estimator, estimand, er, cw, two_stage),
      shape_type = ifelse(estimator == "gcomp", "G-Comp",
                          ifelse(estimator == "aipw", "AIPW", "Unadjusted")),
      line_type = ifelse(estimand == "nat_inf", "Naturally Infected", "Population"),
      Y_out = recode(
        Y_out,
        Y_3 = "Y[3]",
        Y_6 = "Y[6]",
        Y_9 = "Y[9]",
        Y_12 = "Y[12]",
        Y_6_12 = "Y[6-12]",
        Y_3_6_9_12 = "Y[3-6-9-12]"
      ),
      Y_out = factor(
        Y_out,
        levels = c("Y[3]", "Y[6]", "Y[9]", "Y[12]", "Y[6-12]", "Y[3-6-9-12]")
      )
    )
  
  color_map <- c(
    "(CW)" = "#89CFF0",
    "(ER)" = "#4682B4",
    "(ER, two-part)" = "#6495ED",
    "(ER, CW)" = "#191970",
    "(None)" = "#FFC000"
  )
  
  shape_map <- c(
    "G-Comp" = 17,
    "AIPW" = 16,
    "Unadjusted" = 15
  )
  
  ggplot(bias_df, aes(x = n, y = bias, color = color_label, 
                      shape = shape_type, linetype = line_type)) +
    geom_line(aes(group = interaction(color_label, line_type)), linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 4, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1, color = "gray40") +
    facet_wrap(~ Y_out, scales = "fixed", ncol = 3, labeller = label_parsed) +
    scale_y_continuous(name = "Bias") +
    scale_x_continuous(
      name = "Sample Size",
      labels = scales::comma,
      limits = c(0, 50000),
      breaks = seq(0, 50000, 10000)
    ) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_map) +
    scale_linetype_manual(values = c("Naturally Infected" = "solid", "Population" = "dotted")) +
    labs(
      title = paste0("Bias (", setting, " setting)"),
      color = "Assumptions & Modeling Approach",
      shape = "Estimator",
      linetype = "Estimand"
    ) +
    theme_minimal(base_size = 18) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}