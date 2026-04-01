# ------------------------------------------------------
# Script to make final results figures for manuscript
# ------------------------------------------------------

# ○ Power plots -- Y_12, all settings
# ○ Proportion negative plots -- Y_12, all settings
# ○ Endpoint comparison -- optimistic setting, all endpoints

# ○ All remaining combinations in the supplement

library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(patchwork)
library(stringr)
library(knitr)
library(purrr)
library(kableExtra)

here::i_am("R/make_results_figures.R")

config_settings <- c("default",
                     "base_12mo",
                     "optimistic_6mo",
                     "optimistic_12mo",
                     "peru_6mo",
                     "peru_12mo")

all_names <- c("General recruitment", "General recruitment",
               "Targeted recruitment", "Targeted recruitment",
               "High early incidnece", "High early incidnece")

config_setting_names <- c("General recruitment: 6 month immunization schedule",
                          "General recruitment: 12 month immunization schedule",
                          "Targeted recruitment: 6 month immunization schedule",
                          "Targeted recruitment: 12 month immunization schedule",
                          "High early incidence: 6 month immunization schedule",
                          "High early incidence: 12 month immunization schedule")

# -------------------------------------
# Add labels for plotting
# -------------------------------------
add_labels <- function(df) {
  df %>%
    mutate(
      estimand = recode(
        estimand,
        nat_inf = "Naturally Infected",
        pop = "Population"
      ),
      line_type = estimand,
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

# -------------------------------------
# Function to create power plot with shapes
# -------------------------------------
plot_power <- function(power_df, setting = "default") {
  
  color_map <- c(
    "Naturally Infected" = "#00468b",
    "Population" = "#ed0000"
  )
  
  linetype_map <- c(
    "Naturally Infected" = "solid",
    "Population" = "dashed"
  )
  
  ggplot(power_df,
         aes(x = n, y = power,
             color = estimand,
             linetype = estimand,
             group = estimand)) +
    
    geom_line(linewidth = 1.2, alpha = 0.9) +
    geom_point(size = 3.5) +
    
    facet_wrap(~ Y_out, ncol = 3, labeller = label_parsed) +
    
    scale_y_continuous(
      limits = c(0,1),
      breaks = seq(0,1,0.2),
      labels = percent
    ) +
    
    scale_x_continuous(
      name = "Sample Size",
      labels = comma
    ) +
    
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    
    labs(
      title = paste0("Simulated Shigella vaccine trial power\n(", setting, " setting)"),
      y = "Power",
      color = "Estimand",
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
# Function to create proportion negative plot
# -------------------------------------------------
plot_prop_neg <- function(prop_neg_df, setting = "default") {
  
  color_map <- c(
    "Naturally Infected" = "#00468b",
    "Population" = "#ed0000"
  )
  
  linetype_map <- c(
    "Naturally Infected" = "solid",
    "Population" = "dashed"
  )
  
  ggplot(prop_neg_df,
         aes(x = n, y = prop_neg,
             color = estimand,
             linetype = estimand,
             group = estimand)) +
    
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    
    facet_wrap(~ Y_out, ncol = 3, labeller = label_parsed) +
    
    scale_y_continuous(
      limits = c(0,1),
      breaks = seq(0,1,0.2),
      labels = percent,
      name = "Proportion Negative"
    ) +
    
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    
    labs(
      title = paste0("Proportion Negative (", setting, " setting)"),
      color = "Estimand",
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
plot_bias <- function(bias_df, add_CI = FALSE, setting = "default") {
  
  color_map <- c(
    "Naturally Infected" = "#00468b",
    "Population" = "#ed0000"
  )
  
  linetype_map <- c(
    "Naturally Infected" = "solid",
    "Population" = "dashed"
  )
  
  p <- ggplot(
    bias_df,
    aes(x = n, y = bias,
        color = estimand,
        linetype = estimand,
        group = estimand)
  )
  
  if (add_CI) {
    p <- p +
      geom_ribbon(
        aes(
          ymin = bias - 1.96 * sd_bias / sqrt(1000),
          ymax = bias + 1.96 * sd_bias / sqrt(1000),
          fill = estimand
        ),
        alpha = 0.15,
        color = NA
      )
  }
  
  p +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    
    facet_wrap(~ Y_out, ncol = 3, labeller = label_parsed) +
    
    scale_color_manual(values = color_map) +
    scale_fill_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    
    labs(
      title = paste0("Bias (", setting, " setting)"),
      y = "Bias",
      color = "Estimand",
      fill = "Estimand",
      linetype = "Estimand"
    ) +
    
    theme_minimal(base_size = 18) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

# -----------------------------------------------------------------------------

# Individual results 

all_truth_df <- dplyr::tibble()

for(i in 1:length(config_settings)){
  
  setting <- config_settings[i]
  setting_name <- config_setting_names[i]
  
  results <- readRDS(here::here(paste0("results/", setting, "_evaluation_results.Rds")))
  
  power_df <- results$power_df 
  prop_neg_df <- results$prop_neg_df 
  bias_df <- results$bias_df 
  
  power_df <- add_labels(power_df)
  bias_df <- add_labels(bias_df)
  prop_neg_df <- add_labels(prop_neg_df)
  
  power_fig <- plot_power(power_df, setting = setting_name)
  bias_fig <- plot_bias(bias_df, setting = setting_name)
  prop_neg_fig <- plot_prop_neg(prop_neg_df, setting = setting_name)
  
  # Extract truth
  truth_df <- results$truth_df |> 
    tibble::enframe(name = "parameter", value = "value") |> 
    dplyr::mutate(
      value = purrr::map_dbl(value, ~ .x),
      setting = setting,
      setting_name = setting_name
    ) |> 
    dplyr::select(setting, setting_name, parameter, value)
  
  all_truth_df <- dplyr::bind_rows(all_truth_df, truth_df)
  
  ggsave(filename = here::here(paste0("results/figures/individual/power_", setting, ".png")), plot = power_fig, width = 14, height = 7)
  ggsave(filename = here::here(paste0("results/figures/individual/bias_", setting, ".png")), plot = bias_fig, width = 14, height = 7)
  ggsave(filename = here::here(paste0("results/figures/individual/prop_neg_", setting, ".png")), plot = prop_neg_fig, width = 14, height = 7)
}

############################################################
# Supplement truth table
############################################################

# Table 1: True effect sizes

truth_long <- all_truth_df %>%
  mutate(
    Estimand = case_when(
      str_detect(parameter, "^nat_inf") ~ "Naturally Infected",
      str_detect(parameter, "^pop") ~ "Population"
    ),
    Outcome = case_when(
      str_detect(parameter, "Y_3_6_9_12")     ~ "Y[3-6-9-12]",
      str_detect(parameter, "Y_6_12")     ~ "Y[6-12]",
      str_detect(parameter, "Y_12")       ~ "Y[12]",
      str_detect(parameter, "Y_9")        ~ "Y[9]",
      str_detect(parameter, "Y_6")        ~ "Y[6]",
      str_detect(parameter, "Y_3")        ~ "Y[3]"
    )
  )


truth_long <- truth_long %>%
  mutate(
    Table_Col = case_when(
      setting == "default"               ~ "general_6mo",
      setting == "base_12mo"              ~ "general_12mo",
      setting == "optimistic_6mo"  ~ "targeted_6mo",
      setting == "optimistic_12mo" ~ "targeted_12mo",
      setting == "peru_6mo"    ~ "early_6mo",
      setting == "peru_12mo"   ~ "early_12mo"
    )
  )

true_effect_df <- truth_long %>%
  select(Estimand, Outcome, Table_Col, value) %>%
  pivot_wider(
    names_from  = Table_Col,
    values_from = value
  ) %>%
  arrange(Estimand)

# Helper fn for rounding
fmt_num <- function(x) {
  ifelse(
    is.na(x),
    NA,
    ifelse(abs(x) < 0.001, "$<$0.001", sprintf("%.3f", x))
  )
}

true_effect_df %>%
  mutate(
    Outcome = case_when(
      Outcome == "Y[3]"           ~ "$Y_{3}$",
      Outcome == "Y[6]"           ~ "$Y_{6}$",
      Outcome == "Y[9]"           ~ "$Y_{9}$",
      Outcome == "Y[12]"          ~ "$Y_{12}$",
      Outcome == "Y[6-12]"        ~ "$Y_{6\\text{--}12}$",
      Outcome == "Y[3-6-9-12]"    ~ "$Y_{3\\text{--}6\\text{--}9\\text{--}12}$"
    ),
    across(
      c(general_6mo, general_12mo,
        targeted_6mo, targeted_12mo,
        early_6mo, early_12mo),
      fmt_num
    )
  ) %>%
  kable(
    format = "latex",
    booktabs = FALSE,
    escape = FALSE,
    align = "c",
    col.names = c(
      "",
      "Outcome",
      "6 mo", "12 mo",
      "6 mo", "12 mo",
      "6 mo", "12 mo"
    )
  ) %>%
  add_header_above(
    c(
      " " = 2,
      "General recruitment" = 2,
      "Targeted recruitment" = 2,
      "High early incidence" = 2
    ),
    bold = TRUE
  ) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    full_width = FALSE,
    position = "center"
  ) %>%
  row_spec(0, bold = TRUE) %>%
  collapse_rows(
    columns = 1,
    valign = "middle",
    latex_hline = "major"
  ) %>%
  column_spec(2:8, width = "2cm")


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

all_settings <- c("default",
                  "base_12mo",
                  "optimistic_6mo",
                  "optimistic_12mo",
                  "peru_6mo",
                  "peru_12mo")

all_schedule <- rep(c("6 month immunization schedule",
                      "12 month immunization schedule"), 3)

all_names <- c("General recruitment", "General recruitment",
               "Targeted recruitment", "Targeted recruitment",
               "High early incidence", "High early incidence")

############################################################
# Main results figure 1: Y_12 power for all settings
############################################################

# grab all Y[12] endpoints for each setting
all_power_df <- pmap_dfr(
  list(
    setting = all_settings,
    setting_name = all_names,
    immunization_schedule = all_schedule
  ),
  function(setting, setting_name, immunization_schedule) {
    results <- readRDS(here::here(paste0("results/", setting, "_evaluation_results.Rds")))
    
    truth_nat_inf <- results$truth_df$nat_inf_Y_12
    truth_pop <- results$truth_df$pop_Y_12
    
    results$power_df %>%
      add_labels() %>%
      filter(Y_out == "Y[12]") %>%
      mutate(
        setting = setting,
        setting_name = setting_name,
        immunization_schedule = immunization_schedule
      ) %>%
      mutate(
        truth = case_when(
          estimand == "Naturally Infected" ~ truth_nat_inf,
          estimand == "Population" ~ truth_pop
        )
      )
  }
)

all_power_df <- all_power_df %>%
  mutate(immunization_schedule = factor(immunization_schedule, levels = c("6 month immunization schedule",
                                                                          "12 month immunization schedule")),
         setting_name = factor(setting_name, levels = c("General recruitment",
                                                        "Targeted recruitment", 
                                                        "High early incidence")))

# add vertical line for power primary endpoint
primary_power_df <- readRDS(here::here("results/primary_endpoint_power.Rds")) %>%
  mutate(
    setting = case_when(
      setting == "default"               ~ "general_6mo",
      setting == "base_12mo"              ~ "general_12mo",
      setting == "optimistic_6mo"  ~ "targeted_6mo",
      setting == "optimistic_12mo" ~ "targeted_12mo",
      setting == "peru_6mo"    ~ "early_6mo",
      setting == "peru_12mo"   ~ "early_12mo"
    ),
    setting_name = case_when(
      grepl("general", setting)  ~ "General recruitment",
      grepl("targeted", setting) ~ "Targeted recruitment",
      grepl("early", setting)    ~ "High early incidence"
    ),
    setting_name = factor(setting_name,
                          levels = c("General recruitment",
                                     "Targeted recruitment",
                                     "High early incidence")),
    immunization_schedule = if_else(grepl("6mo", setting), "6 month immunization schedule", "12 month immunization schedule"),
    immunization_schedule = factor(immunization_schedule,
                                   levels = c("6 month immunization schedule",
                                              "12 month immunization schedule")),
    ages_0_6 = case_when(
      immunization_schedule == "6 month immunization schedule" ~ "6-12 months", 
      immunization_schedule == "12 month immunization schedule" ~ "12-18 months"
    ),
    ages_6_12 = case_when(
      immunization_schedule == "6 month immunization schedule" ~ "12-18 months", 
      immunization_schedule == "12 month immunization schedule" ~ "18-24 months"
    ),
  )


plot_power_combo <- function(df) {
  
  color_map <- c(
    "Naturally Infected" = "#00468b",
    "Population" = "#ed0000"
  )
  
  linetype_map <- c(
    "Naturally Infected" = "solid",
    "Population" = "solid"
  )
  
  ggplot(
    df,
    aes(
      x = n,
      y = power,
      color = estimand,
      linetype = estimand,
      group = estimand
    )
  ) +
    
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    
    # truth annotation
    geom_text(
      aes(
        x = min(n), 
        y = if_else(estimand == "Naturally Infected", 0.55, 0.45),
        label = ifelse(estimand == "Naturally Infected", 
                       paste0("Naturally Infected effect: ", round(truth, 3), " HAZ"),
                       paste0("Population effect: ", round(truth, 3), " HAZ")),
        color = estimand
      ),
      hjust = 0,
      size = 3.8,
      show.legend = FALSE
    ) +
    
    facet_grid(
      setting_name ~ immunization_schedule
    ) +
    
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    
    # vertical line at 90% power sample size (if it's in range)
    # geom_vline(
    #   data = primary_power_df %>% filter(n >= 2500),
    #   aes(xintercept = n),
    #   linetype = "dotted",
    #   color = "black",
    #   alpha = 0.4,
    #   linewidth = 0.8,
    #   inherit.aes = FALSE
    # ) +
    
    # annotation primary endpoint 90% power
    geom_text(
      data = primary_power_df,
      aes(
        x = 2500,
        y = 0.95,
        label = paste0("Expected trial size: ", "n=", scales::comma(n))
      ),
      angle = 0,
      hjust = 0,
      size = 3.8,
      color = "black",
      #alpha = 0.7,
      inherit.aes = FALSE
    ) +
    
    # annotation incidence in placebo arm
    geom_text(
      data = primary_power_df,
      aes(
        x = 2500,
        y = 0.75,
        label = paste0("Placebo incidence (MSD cases / 100 child-years):\n",
                       round(round(inc_sev_0_6,3)*2*100, 1), " (ages ", ages_0_6, ")\n",
                       round(round(inc_sev_6_12,3)*2*100, 1), " (ages ", ages_6_12, ")")
        
        # label = paste0("Placebo incidence: ",
        #        round(inc_0_12*100, 1), " MSD cases / 100 child-years"),
                       
      ),
      angle = 0,
      hjust = 0,
      size = 3.8,
      color = "black",
      #alpha = 0.7,
      inherit.aes = FALSE
    ) +
    
    scale_y_continuous(
      limits = c(0,1),
      breaks = seq(0,1,0.2),
      labels = scales::percent
    ) +
    
    scale_x_log10(
      breaks = c(2500, 5000, 10000, 20000, 40000, 80000),
      labels = scales::comma
    ) +
    
    labs(
      x = "Sample Size",
      y = expression(Power~"to detect vaccine effect on growth at 12 months"),
      color = "Estimand",
      linetype = "Estimand"
    ) +
    
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 14, face = "bold"),
      panel.spacing = unit(14, "pt"),
      panel.grid.minor = element_blank()
    )
}

power_fig <- plot_power_combo(all_power_df)

ggsave(here::here("results/figures/power_figure_Y_12.png"),
       plot = power_fig,
       width = 10,
       height = 10)

############################################################
# Main results figure 2: Y_12 proportion negative for all settings
############################################################

# grab all Y[12] endpoints for each setting
all_prop_neg_df <- pmap_dfr(
  list(
    setting = all_settings,
    setting_name = all_names,
    immunization_schedule = all_schedule
  ),
  function(setting, setting_name, immunization_schedule) {
    results <- readRDS(here::here(paste0("results/", setting, "_evaluation_results.Rds")))
    
    # Get proportion significant -- check if upper bound < 0
    results$results$sig_neg <- ifelse(results$results$estimate + 1.96*results$results$se < 0, 1, 0)
    
    sig_neg_summary <- results$results %>%
      group_by(estimand, n, Y_out) %>%
      summarise(
        prop_sig_neg = mean(sig_neg),
        prop_sig_neg__neg_est = mean(sig_neg[estimate < 0]),
        .groups = "drop"
      )
    
    # stopped here, join with reuslts
    truth_nat_inf <- results$truth_df$nat_inf_Y_12
    truth_pop <- results$truth_df$pop_Y_12
    
    results$prop_neg_df %>%
      left_join(sig_neg_summary, by = c("estimand", "n", "Y_out")) %>%
      add_labels() %>%
      filter(Y_out == "Y[12]") %>%
      mutate(
        setting = setting,
        setting_name = setting_name,
        immunization_schedule = immunization_schedule
      ) %>%
      mutate(
        truth = case_when(
          estimand == "Naturally Infected" ~ truth_nat_inf,
          estimand == "Population" ~ truth_pop
        )
      )
    
  }
)

all_prop_neg_df <- all_prop_neg_df %>%
  mutate(immunization_schedule = factor(immunization_schedule, levels = c("6 month immunization schedule",
                                                                          "12 month immunization schedule")),
         setting_name = factor(setting_name, levels = c("General recruitment",
                                                        "Targeted recruitment", 
                                                        "High early incidence")))


all_prop_neg_df <- all_prop_neg_df %>%
  mutate(
    sig_neg_label = scales::percent(prop_sig_neg__neg_est, accuracy = 0.1)
  )

# For making prop_neg plot all settings
plot_prop_neg_combo <- function(df) {
  
  color_map <- c(
    "Naturally Infected" = "#00468b",
    "Population" = "#ed0000"
  )
  
  linetype_map <- c(
    "Naturally Infected" = "solid",
    "Population" = "solid"
  )
  
  ggplot(
    df,
    aes(
      x = n,
      y = prop_neg,
      color = estimand,
      linetype = estimand,
      group = estimand
    )
  ) +
    
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    
    # geom_point(
    #   aes(y = prop_sig_neg),
    #   shape = 21,
    #   stroke = 1.2,
    #   size = 1,
    #   fill = "white",
    #   position = position_nudge(y = 0.02)
    # ) +
  
    # Add truth annotation in top-left, slightly offset to avoid overlap
    # geom_text(
    #   aes(
    #     x = 40000, 
    #     y = if_else(estimand == "Naturally Infected", 0.45, 0.40),
    #     label = ifelse(estimand == "Naturally Infected", 
    #                    paste0("Naturally Infected: ", round(truth, 3)),
    #                    paste0("Population: ", round(truth, 3))),
    #     color = estimand
    #   ),
    #   hjust = 0,
    #   size = 4,
    #   show.legend = FALSE
    # ) +
    
    facet_grid(
      setting_name ~ immunization_schedule
    ) +
    
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    
    scale_y_continuous(
      limits = c(0,0.50),
      breaks = seq(0,0.5,0.1),
      labels = scales::percent
    ) +
    
    scale_x_log10(
      breaks = c(2500, 5000, 10000, 20000, 40000, 80000),
      labels = scales::comma
    ) +
    
    labs(
      x = "Sample Size",
      y = expression("Proportion negative vaccine effect estimates at 12 months"),
      color = "Estimand",
      linetype = "Estimand"
    ) +
    
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 14, face = "bold"),
      panel.spacing = unit(14, "pt"),
      panel.grid.minor = element_blank()
    )
}

# For making prop_neg plot all settings
plot_prop_sig_neg_combo <- function(df) {
  
  color_map <- c(
    "Naturally Infected" = "#00468b",
    "Population" = "#ed0000"
  )
  
  linetype_map <- c(
    "Naturally Infected" = "solid",
    "Population" = "solid"
  )
  
  ggplot(
    df,
    aes(
      x = n,
      y = prop_sig_neg,
      color = estimand,
      linetype = estimand,
      group = estimand
    )
  ) +
    
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5, fill = "white") +
    
    # geom_point(
    #   aes(y = prop_sig_neg),
    #   shape = 21,
    #   stroke = 1.2,
    #   size = 1,
    #   fill = "white",
    #   position = position_nudge(y = 0.02)
    # ) +
    
    # Add truth annotation in top-left, slightly offset to avoid overlap
    # geom_text(
    #   aes(
    #     x = 40000, 
    #     y = if_else(estimand == "Naturally Infected", 0.45, 0.40),
    #     label = ifelse(estimand == "Naturally Infected", 
    #                    paste0("Naturally Infected: ", round(truth, 3)),
    #                    paste0("Population: ", round(truth, 3))),
    #     color = estimand
    #   ),
    #   hjust = 0,
    #   size = 4,
    #   show.legend = FALSE
    # ) +
    
    facet_grid(
      setting_name ~ immunization_schedule
    ) +
    
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    
    scale_y_continuous(
      limits = c(0,0.04),
      breaks = seq(0,0.04,0.01),
      labels = scales::percent
    ) +
    
    scale_x_log10(
      breaks = c(2500, 5000, 10000, 20000, 40000, 80000),
      labels = scales::comma
    ) +
    
    labs(
      x = "Sample Size",
      y = expression("Proportion significant negative vaccine effect estimates at 12 months"),
      color = "Estimand",
      linetype = "Estimand"
    ) +
    
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 14, face = "bold"),
      panel.spacing = unit(14, "pt"),
      panel.grid.minor = element_blank()
    )
}

prop_neg_fig <- plot_prop_neg_combo(all_prop_neg_df)
prop_neg_sig_fig <- plot_prop_sig_neg_combo(all_prop_neg_df)

ggsave(here::here("results/figures/prop_neg_figure_Y_12.png"),
       plot = prop_neg_fig,
       width = 10,
       height = 10)

ggsave(here::here("results/figures/prop_sig_neg_Y_12.png"),
       plot = prop_neg_sig_fig,
       width = 10,
       height = 10)

############################################################
# Main results figure 3: Comparison of power for all settings
# (and filter differently for supp)
############################################################

# all_settings <- c("optimistic_6mo",
#                   "optimistic_12mo")
# 
# all_schedule <- c("6 month immunization schedule",
#                   "12 month immunization schedule")
# 
# all_names <- c("Targeted recruitment", 
#                "Targeted recruitment")

all_settings <- c("default",
                     "base_12mo",
                     "optimistic_6mo",
                     "optimistic_12mo",
                     "peru_6mo",
                     "peru_12mo")

all_schedule <- rep(c("6 month immunization schedule",
                      "12 month immunization schedule"), 3)

all_names <- c("General recruitment", "General recruitment",
               "Targeted recruitment", "Targeted recruitment",
               "High early incidence", "High early incidence")

all_power_optimistic_df <- pmap_dfr(
  list(
    setting = all_settings,
    setting_name = all_names,
    immunization_schedule = all_schedule
  ),
  function(setting, setting_name, immunization_schedule) {
    results <- readRDS(here::here(paste0("results/", setting, "_evaluation_results.Rds")))
    
    truth_nat_inf_Y_12 <- results$truth_df$nat_inf_Y_12
    truth_pop_Y_12 <- results$truth_df$pop_Y_12
    
    truth_nat_inf_Y_9 <- results$truth_df$nat_inf_Y_9
    truth_pop_Y_9 <- results$truth_df$pop_Y_9
    
    truth_nat_inf_Y_6 <- results$truth_df$nat_inf_Y_6
    truth_pop_Y_6 <- results$truth_df$pop_Y_6
    
    truth_nat_inf_Y_3 <- results$truth_df$nat_inf_Y_3
    truth_pop_Y_3 <- results$truth_df$pop_Y_3
    
    truth_nat_inf_Y_6_12 <- results$truth_df$nat_inf_Y_6_12
    truth_pop_Y_6_12 <- results$truth_df$pop_Y_6_12
    
    truth_nat_inf_Y_3_6_9_12 <- results$truth_df$nat_inf_Y_3_6_9_12
    truth_pop_Y_3_6_9_12 <- results$truth_df$pop_Y_3_6_9_12
    
    results$power_df %>%
      add_labels() %>%
      mutate(
        setting = setting,
        setting_name = setting_name,
        immunization_schedule = immunization_schedule
      ) %>%
      mutate(
        truth = case_when(
          estimand == "Naturally Infected" & Y_out == "Y[12]" ~ truth_nat_inf_Y_12,
          estimand == "Population" & Y_out == "Y[12]" ~ truth_pop_Y_12,
          
          estimand == "Naturally Infected" & Y_out == "Y[9]" ~ truth_nat_inf_Y_9,
          estimand == "Population" & Y_out == "Y[9]" ~ truth_pop_Y_9,
          
          estimand == "Naturally Infected" & Y_out == "Y[6]" ~ truth_nat_inf_Y_6,
          estimand == "Population" & Y_out == "Y[6]" ~ truth_pop_Y_6,
          
          estimand == "Naturally Infected" & Y_out == "Y[3]" ~ truth_nat_inf_Y_3,
          estimand == "Population" & Y_out == "Y[3]" ~ truth_pop_Y_3,
          
          estimand == "Naturally Infected" & Y_out == "Y[6-12]" ~ truth_nat_inf_Y_6_12,
          estimand == "Population" & Y_out == "Y[6-12]" ~ truth_pop_Y_6_12,
          
          estimand == "Naturally Infected" & Y_out == "Y[3-6-9-12]" ~ truth_nat_inf_Y_3_6_9_12,
          estimand == "Population" & Y_out == "Y[3-6-9-12]" ~ truth_pop_Y_3_6_9_12,
          
        )
      )
  }
)

all_power_optimistic_df_text <- all_power_optimistic_df %>%
  mutate(immunization_schedule = factor(immunization_schedule, levels = c("6 month immunization schedule",
                                                                          "12 month immunization schedule")),
         setting_name = factor(setting_name, levels = c("General recruitment",
                                                        "Targeted recruitment", 
                                                        "High early incidence")),
         Y_out = factor(Y_out, levels = c("Y[3]", "Y[6]", "Y[9]", "Y[12]", "Y[6-12]", "Y[3-6-9-12]"), 
                        labels = c("3 month outcome",
                                   "6 month outcome",
                                   "9 month outcome",
                                   "12 month outcome",
                                   "6 and 12 month outcome",
                                   "3,6,9,12 month outcome")))

all_power_optimistic_df <- all_power_optimistic_df %>%
  mutate(immunization_schedule = factor(immunization_schedule, levels = c("6 month immunization schedule",
                                                                          "12 month immunization schedule")),
         setting_name = factor(setting_name, levels = c("General recruitment",
                                                        "Targeted recruitment", 
                                                        "High early incidence")),
         Y_out = factor(Y_out, levels = c("Y[3]", "Y[6]", "Y[9]", "Y[12]", "Y[6-12]", "Y[3-6-9-12]")))

plot_power_all_end <- function(df, text_label = TRUE, plot_title = NULL) {
  
  color_map <- c(
    "Naturally Infected" = "#00468b",
    "Population" = "#ed0000"
  )
  
  linetype_map <- c(
    "Naturally Infected" = "solid",
    "Population" = "solid"
  )
  
  if(text_label){
    ggplot(
      df,
      aes(
        x = n,
        y = power,
        color = estimand,
        linetype = estimand,
        group = estimand
      )
    ) +
      
      geom_line(linewidth = 1.2) +
      geom_point(size = 3.5) +
      
      # Truth annotation
      geom_text(
        aes(
          x = min(n),
          y = if_else(estimand == "Naturally Infected", 0.95, 0.85),
          label = ifelse(
            estimand == "Naturally Infected",
            paste0("Naturally Infected: ", round(truth, 3)),
            paste0("Population: ", round(truth, 3))
          ),
          color = estimand
        ),
        hjust = 0,
        size = 4,
        show.legend = FALSE
      ) +
      
      # primary power 
      # geom_text(
      #   data = primary_power_df,
      #   aes(
      #     x = 2500,
      #     y = 0.75,
      #     label = paste0("90% power primary endpoint: n=", scales::comma(n))
      #   ),
      #   angle = 0,
      #   hjust = 0,
      #   size = 4,
      #   color = "black",
      #   alpha = 0.7,
      #   inherit.aes = FALSE
      # ) +
      
      facet_grid(
        Y_out ~ immunization_schedule#,
        #labeller = labeller(Y_out = label_parsed)
      ) + 
      
      scale_color_manual(values = color_map) +
      scale_linetype_manual(values = linetype_map) +
      
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2),
        labels = scales::percent
      ) +
      
      scale_x_log10(
        breaks = c(2500, 5000, 10000, 20000, 40000, 80000),
        labels = scales::comma
      ) +
      
      labs(
        x = "Sample Size",
        y = "Power to detect vaccine effect",
        color = "Estimand",
        linetype = "Estimand",
        title = plot_title
      ) +
      
      theme_minimal(base_size = 16) +
      theme(
        legend.position = "bottom",
        strip.text = element_text(size = 14, face = "bold"),
        panel.spacing = unit(14, "pt"),
        panel.grid.minor = element_blank()
      )
  } else{
    ggplot(
      df,
      aes(
        x = n,
        y = power,
        color = estimand,
        linetype = estimand,
        group = estimand
      )
    ) +
      
      geom_line(linewidth = 1.2) +
      geom_point(size = 3.5) +
      
      # Truth annotation
      geom_text(
        aes(
          x = min(n),
          y = if_else(estimand == "Naturally Infected", 0.95, 0.80),
          label = ifelse(
            estimand == "Naturally Infected",
            paste0("Naturally Infected: ", round(truth, 3)),
            paste0("Population: ", round(truth, 3))
          ),
          color = estimand
        ),
        hjust = 0,
        size = 4,
        show.legend = FALSE
      ) +
      
      # primary power 
      # geom_text(
      #   data = primary_power_df,
      #   aes(
      #     x = 2500,
      #     y = 0.75,
      #     label = paste0("90% power primary endpoint: n=", scales::comma(n))
      #   ),
      #   angle = 0,
      #   hjust = 0,
      #   size = 4,
      #   color = "black",
      #   alpha = 0.7,
      #   inherit.aes = FALSE
      # ) +
      
      facet_grid(
        Y_out ~ immunization_schedule,
        labeller = labeller(Y_out = label_parsed)
      ) + 
      
      scale_color_manual(values = color_map) +
      scale_linetype_manual(values = linetype_map) +
      
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2),
        labels = scales::percent
      ) +
      
      scale_x_log10(
        breaks = c(2500, 5000, 10000, 20000, 40000, 80000),
        labels = scales::comma
      ) +
      
      labs(
        x = "Sample Size",
        y = "Power to detect vaccine effect",
        color = "Estimand",
        linetype = "Estimand",
        title = plot_title
      ) +
      
      theme_minimal(base_size = 16) +
      theme(
        legend.position = "bottom",
        strip.text = element_text(size = 14, face = "bold"),
        panel.spacing = unit(14, "pt"),
        panel.grid.minor = element_blank()
      )
  }
  
  
}

main_endpoint_power_plot <- plot_power_all_end(all_power_optimistic_df_text %>% 
                                                 filter(Y_out %in%c("6 month outcome", "12 month outcome", "6 and 12 month outcome") &
                                                          setting_name == "Targeted recruitment")) #c("Y[6]",  "Y[12]", "Y[6-12]" )))

ggsave(here::here("results/figures/targeted_all_endpoint.png"),
       plot = main_endpoint_power_plot,
       width = 10,
       height = 10)

# All outcomes for the supplement
targeted_power_plot <- plot_power_all_end(all_power_optimistic_df %>% 
                                                filter(setting_name == "Targeted recruitment"), 
                                              text_label = FALSE,
                                              plot_title = "Targeted recruitment")

general_power_plot <- plot_power_all_end(all_power_optimistic_df %>% 
                                            filter(setting_name == "General recruitment"), 
                                          text_label = FALSE,
                                          plot_title = "General recruitment")

high_power_plot <- plot_power_all_end(all_power_optimistic_df %>% 
                                            filter(setting_name == "High early incidence"), 
                                          text_label = FALSE,
                                          plot_title = "High early incidence")

ggsave(here::here("results/figures/targeted_power_endpoint_supp.png"),
       plot = targeted_power_plot,
       width = 10,
       height = 10)

ggsave(here::here("results/figures/general_power_endpoint_supp.png"),
       plot = general_power_plot,
       width = 10,
       height = 10)

ggsave(here::here("results/figures/high_power_endpoint_supp.png"),
       plot = high_power_plot,
       width = 10,
       height = 10)

############################################################
# Supplement proportion negative figures for all endpoints
############################################################

# grab all Y[12] endpoints for each setting
all_prop_neg_df <- pmap_dfr(
  list(
    setting = all_settings,
    setting_name = all_names,
    immunization_schedule = all_schedule
  ),
  function(setting, setting_name, immunization_schedule) {
    results <- readRDS(here::here(paste0("results/", setting, "_evaluation_results.Rds")))
    
    # Get proportion significant -- check if upper bound < 0
    results$results$sig_neg <- ifelse(results$results$estimate + 1.96*results$results$se < 0, 1, 0)
    
    sig_neg_summary <- results$results %>%
      group_by(estimand, n, Y_out) %>%
      summarise(
        prop_sig_neg = mean(sig_neg),
        prop_sig_neg__neg_est = mean(sig_neg[estimate < 0]),
        .groups = "drop"
      )
    
    # stopped here, join with reuslts
    truth_nat_inf <- results$truth_df$nat_inf_Y_12
    truth_pop <- results$truth_df$pop_Y_12
    
    results$prop_neg_df %>%
      left_join(sig_neg_summary, by = c("estimand", "n", "Y_out")) %>%
      add_labels() %>%
      #filter(Y_out == "Y[12]") %>%
      mutate(
        setting = setting,
        setting_name = setting_name,
        immunization_schedule = immunization_schedule
      ) 
    
  }
)

all_prop_neg_df <- all_prop_neg_df %>%
  mutate(immunization_schedule = factor(immunization_schedule, levels = c("6 month immunization schedule",
                                                                          "12 month immunization schedule")),
         setting_name = factor(setting_name, levels = c("General recruitment",
                                                        "Targeted recruitment", 
                                                        "High early incidence")))


all_prop_neg_df <- all_prop_neg_df %>%
  mutate(
    sig_neg_label = scales::percent(prop_sig_neg__neg_est, accuracy = 0.1)
  )

# For making prop_neg plot all settings
plot_prop_neg_combo <- function(df, plot_title = NULL) {
  
  color_map <- c(
    "Naturally Infected" = "#00468b",
    "Population" = "#ed0000"
  )
  
  linetype_map <- c(
    "Naturally Infected" = "solid",
    "Population" = "solid"
  )
  
  ggplot(
    df,
    aes(
      x = n,
      y = prop_neg,
      color = estimand,
      linetype = estimand,
      group = estimand
    )
  ) +
    
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
  
    facet_grid(
      Y_out ~ immunization_schedule,
      labeller = labeller(Y_out = label_parsed)
    ) + 
    
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    
    scale_y_continuous(
      limits = c(0,0.60),
      breaks = seq(0,0.6,0.2),
      labels = scales::percent
    ) +
    
    scale_x_log10(
      breaks = c(2500, 5000, 10000, 20000, 40000, 80000),
      labels = scales::comma
    ) +
    
    labs(
      x = "Sample Size",
      y = expression("Proportion negative vaccine effect estimates"),
      color = "Estimand",
      linetype = "Estimand",
      title = plot_title
    ) +
    
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 14, face = "bold"),
      panel.spacing = unit(14, "pt"),
      panel.grid.minor = element_blank()
    )
}

# For making prop_neg plot all settings
plot_prop_sig_neg_combo <- function(df, plot_title = NULL) {
  
  color_map <- c(
    "Naturally Infected" = "#00468b",
    "Population" = "#ed0000"
  )
  
  linetype_map <- c(
    "Naturally Infected" = "solid",
    "Population" = "solid"
  )
  
  ggplot(
    df,
    aes(
      x = n,
      y = prop_sig_neg,
      color = estimand,
      linetype = estimand,
      group = estimand
    )
  ) +
    
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5, fill = "white") +
    
    # geom_point(
    #   aes(y = prop_sig_neg),
    #   shape = 21,
    #   stroke = 1.2,
    #   size = 1,
    #   fill = "white",
    #   position = position_nudge(y = 0.02)
    # ) +
    
    # Add truth annotation in top-left, slightly offset to avoid overlap
    # geom_text(
    #   aes(
    #     x = 40000, 
    #     y = if_else(estimand == "Naturally Infected", 0.45, 0.40),
    #     label = ifelse(estimand == "Naturally Infected", 
    #                    paste0("Naturally Infected: ", round(truth, 3)),
    #                    paste0("Population: ", round(truth, 3))),
    #     color = estimand
    #   ),
    #   hjust = 0,
    #   size = 4,
    #   show.legend = FALSE
    # ) +
    
    facet_grid(
      Y_out ~ immunization_schedule,
      labeller = labeller(Y_out = label_parsed)
    ) + 
    
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    
    scale_y_continuous(
      limits = c(0,0.04),
      breaks = seq(0,0.04,0.01),
      labels = scales::percent
    ) +
    
    scale_x_log10(
      breaks = c(2500, 5000, 10000, 20000, 40000, 80000),
      labels = scales::comma
    ) +
    
    labs(
      x = "Sample Size",
      y = expression("Proportion significant negative vaccine effect estimates at 12 months"),
      color = "Estimand",
      linetype = "Estimand",
      title = plot_title
    ) +
    
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 14, face = "bold"),
      panel.spacing = unit(14, "pt"),
      panel.grid.minor = element_blank()
    )
}

prop_neg_fig_general <- plot_prop_neg_combo(all_prop_neg_df %>% filter(setting_name == "General recruitment"), plot_title = "General recruitment")
prop_neg_fig_target <- plot_prop_neg_combo(all_prop_neg_df %>% filter(setting_name == "Targeted recruitment"), plot_title = "Targeted recruitment")
prop_neg_fig_high_early_inc <- plot_prop_neg_combo(all_prop_neg_df %>% filter(setting_name == "High early incidence"), plot_title = "High early incidence")

prop_neg_sig_fig_general <- plot_prop_sig_neg_combo(all_prop_neg_df %>% filter(setting_name == "General recruitment"), plot_title = "General recruitment")
prop_neg_sig_fig_target <- plot_prop_sig_neg_combo(all_prop_neg_df %>% filter(setting_name == "Targeted recruitment"), plot_title = "Targeted recruitment")
prop_neg_sig_fig_high_early_inc <- plot_prop_sig_neg_combo(all_prop_neg_df %>% filter(setting_name == "High early incidence"), plot_title = "High early incidence")

ggsave(here::here("results/figures/prop_neg_figure_general_supp.png"),
       plot = prop_neg_fig_general,
       width = 10,
       height = 10)

ggsave(here::here("results/figures/prop_neg_figure_target_supp.png"),
       plot = prop_neg_fig_target,
       width = 10,
       height = 10)

ggsave(here::here("results/figures/prop_neg_figure_high_early_inc_supp.png"),
       plot = prop_neg_fig_high_early_inc,
       width = 10,
       height = 10)


ggsave(here::here("results/figures/prop_sig_neg_figure_general_supp.png"),
       plot = prop_neg_sig_fig_general,
       width = 10,
       height = 10)

ggsave(here::here("results/figures/prop_sig_neg_figure_target_supp.png"),
       plot = prop_neg_sig_fig_target,
       width = 10,
       height = 10)

ggsave(here::here("results/figures/prop_sig_neg_figure_high_early_inc_supp.png"),
       plot = prop_neg_sig_fig_high_early_inc,
       width = 10,
       height = 10)
