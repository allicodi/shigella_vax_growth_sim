# ------------------------------------------------------------------
# Shigella vaccine trial plots - power, proportion negative
# ------------------------------------------------------------------

library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(patchwork)

library(dplyr)
library(tidyr)
library(stringr)
library(knitr)
library(kableExtra)

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
    # is.na(er) & is.na(cw) & is.na(two_stage), "(None)", # this is the unadjusted version
    is.na(er) & is.na(cw) & is.na(two_stage), "(ER)", # QUESTION unadjusted actually does have ER assumption??
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
    "(CW)" = "#191970",
    "(ER)" = "#6082B6",          
    "(ER, two-part)" = "#AD002AFF",
    "(ER, CW)" =  "#925E9FFF",
    
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
      breaks = seq(0, 0.75, by = 0.25),
      labels = scales::percent_format(accuracy = 1),
      name = "Power"
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_continuous(
      name = "Sample Size",
      labels = scales::comma,
      breaks = scales::pretty_breaks(n = 5)
    ) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_map) +
    scale_linetype_manual(values = c("Naturally Infected" = "solid", "Population" = "dotted")) +
    labs(
      title = paste0("Simulated Shigella vaccine trial power \n(", setting, " setting)"),
      color = "Assumptions & Modeling Approach",
      shape = "Estimator",
      linetype = "Estimand"
    ) +
    theme_minimal(base_size = 18) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.margin = margin(t = 10, b = 10),
      legend.spacing.y = unit(4, "mm"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      strip.text = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10)   
    )+ 
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      shape = guide_legend(nrow = 1),
      linetype = guide_legend(nrow = 1)
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
  
  # color_map <- c(
  #   "(CW)" = "#89CFF0",
  #   "(ER)" = "#4682B4",
  #   "(ER, two-part)" = "#6495ED",
  #   "(ER, CW)" = "#191970",
  #   "(None)" = "#FFC000"
  # )
  
  color_map <- c(
    "(CW)" = "#191970",
    "(ER)" = "#6082B6",          
    "(ER, two-part)" = "#AD002AFF",
    "(ER, CW)" =  "#925E9FFF",
    
    "(None)" = "#FFC000" # Unadjusted
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
plot_bias <- function(bias_df, add_CI = FALSE, setting = "default") {
  
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
    "(CW)" = "#191970",
    "(ER)" = "#6082B6",
    "(ER, two-part)" = "#AD002AFF",
    "(ER, CW)" =  "#925E9FFF",
    "(None)" = "#FFC000"
  )
  
  shape_map <- c(
    "G-Comp" = 17,
    "AIPW" = 16,
    "Unadjusted" = 15
  )
  
  p <- ggplot(bias_df, aes(x = n, y = bias, color = color_label,
                           shape = shape_type, linetype = line_type))
  
  # -------------------------------------------------
  # Add CI ribbon if requested
  # -------------------------------------------------
  if (add_CI) {
    p <- p +
      geom_ribbon(
        aes(
          ymin = bias - 1.96 * sd_bias / sqrt(1000),
          ymax = bias + 1.96 * sd_bias / sqrt(1000),
          fill = color_label,
          group = interaction(color_label, line_type)
        ),
        alpha = 0.15,
        color = NA
      )
  }
  
  # -------------------------------------------------
  # Main plot
  # -------------------------------------------------
  p +
    geom_line(aes(group = interaction(color_label, line_type)),
              linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 4, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed",
               size = 1, color = "gray40") +
    facet_wrap(~ Y_out, scales = "fixed", ncol = 3,
               labeller = label_parsed) +
    scale_y_continuous(name = "Bias") +
    scale_x_continuous(
      name = "Sample Size",
      labels = scales::comma,
      limits = c(0, 50000),
      breaks = seq(0, 50000, 10000)
    ) +
    scale_color_manual(values = color_map) +
    scale_fill_manual(values = color_map) +
    scale_shape_manual(values = shape_map) +
    scale_linetype_manual(values = c("Naturally Infected" = "solid",
                                     "Population" = "dotted")) +
    labs(
      title = paste0("Bias (", setting, " setting)"),
      color = "Assumptions & Modeling Approach",
      fill = "Assumptions & Modeling Approach",
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


# ---------------------------------------------------------------------------

config_settings <- c("default",
                     "base_12mo",
                     # "gambia_6mo",
                     # "gambia_12mo",
                     # "peru_6mo",
                     # "peru_12mo",
                     "gambia_default_VE_6mo",
                     "gambia_default_VE_12mo",
                     "peru_default_VE_6mo",
                     "peru_default_VE_12mo")
                     # "optimistic_6mo",
                     # "optimistic_12mo",
                     # "high_inc_ge_6mo",
                     # "high_inc_ge_12mo")

config_setting_names <- c("6 month trial",
                          "12 month trial",
                          # "6 month trial- The Gambia, optimistic growth effect, higher VE",
                          # "12 month trial- The Gambia, optimistic growth effect, higher VE",
                          # "6 month trial- Peru, optimistic growth effect, higher VE",
                          # "12 month trial- Peru, optimistic growth effect, higher VE",
                          "6 month trial- The Gambia, optimistic growth effect",
                          "12 month trial- The Gambia, optimistic growth effect",
                          "6 month trial- Peru, optimistic growth effect",
                          "12 month trial- Peru, optimistic growth effect")

# config_settings <- c("default",
#                      "vax_12mo",
#                      "vax_6mo_high_inc",
#                      "vax_12mo_high_inc",
#                      "vax_6mo_optimistic_ge",
#                      "vax_12mo_optimistic_ge",
#                      "vax_6mo_optimistic_ge_EM",
#                      "vax_12mo_optimistic_ge_EM",
#                      "vax_6mo_high_inc_optimistic_ge_EM")
# 
# config_setting_names <- c("6 month trial",
#                           "12 month trial",
#                           "6 month trial, high incidence (Peru)",
#                           "12 month trial, high incidence (Peru)",
#                           "6 month trial, optimistic growth effect",
#                           "12 month trial, optimistic growth effect",
#                           "6 month trial, optimistic growth effect,\n effect modification by age",
#                           "12 month trial, optimistic growth effect,\n effect modification by age",
#                           "6 month trial, high incidence (Peru),\n optimistic growth effect, effect modification by age")

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
  
  power_fig <- plot_power(power_df, setting = setting_name)
  bias_fig <- plot_bias(bias_df)
  
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
  
  
  ggsave(filename = here::here(paste0("results/figures/power_", setting, ".png")), plot = power_fig)
  
}

# ------------------------------------------------------------------------
# Figures for paper
# ------------------------------------------------------------------------

# Table 1: True effect sizes

truth_long <- all_truth_df %>%
  mutate(
    Estimand = case_when(
      str_detect(parameter, "^nat_inf") ~ "Naturally Infected",
      str_detect(parameter, "^pop") ~ "Population"
    ),
    Outcome = case_when(
      str_detect(parameter, "Y_3_6_9_12") ~ "Y[3-6-9-12]",
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
      setting == "default"               ~ "Overall_6mo",
      setting == "base_12mo"              ~ "Overall_12mo",
      setting == "gambia_default_VE_6mo"  ~ "The_Gambia_6mo",
      setting == "gambia_default_VE_12mo" ~ "The_Gambia_12mo",
      setting == "peru_default_VE_6mo"    ~ "Peru_6mo",
      setting == "peru_default_VE_12mo"   ~ "Peru_12mo"
    )
  )

true_effect_df <- truth_long %>%
  select(Estimand, Outcome, Table_Col, value) %>%
  pivot_wider(
    names_from  = Table_Col,
    values_from = value
  ) %>%
  arrange(Estimand)



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
      c(Overall_6mo, Overall_12mo,
        The_Gambia_6mo, The_Gambia_12mo,
        Peru_6mo, Peru_12mo),
      ~ sprintf("%.4f", .x)
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
      "Overall" = 2,
      "The Gambia" = 2,
      "Peru" = 2
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
  )


# Figure 1 (2? if babies) - Power for Overall  ------------------------------

# config_settings <- c("default",
#                      "base_12mo")
# 
# config_setting_names <- c("6 month vaccine schedule",
#                           "12 month vaccine schedule")


config_settings <- c("peru_6mo",
                     "peru_12mo")

config_setting_names <- c("Peru 6 month vaccine schedule",
                          "Peru 12 month vaccine schedule")

all_power_df <- purrr::map2_dfr(
  config_settings,
  config_setting_names,
  ~ {
    results <- readRDS(here::here(paste0("results/", .x, "_evaluation_results.Rds")))
    
    power_df <- results$power_df |> add_labels()
    
    power_df |> 
      mutate(
        setting = .x,
        setting_name = .y
      )
  }
)



all_power_df <- all_power_df %>%
  mutate(
    color_group = case_when(
      shape_type == "Unadjusted"  ~ "Unadjusted",
      line_type == "Naturally Infected"                    ~ "Naturally Infected",
      line_type == "Population"                            ~ "Population"
    ),
    setting_name = factor(setting_name, levels = c("Peru 6 month vaccine schedule",
                                                   "Peru 12 month vaccine schedule"))
  )

# FILTER UNADJUSTED
all_power_df <- all_power_df %>%
  filter(shape_type != "Unadjusted")

plot_power_combined <- function(power_df) {
  
  # create a new color_label
  power_df <- power_df %>%
    mutate(
      color_label = case_when(
        line_type == "Naturally Infected" & shape_type == "AIPW"     ~ "Naturally Infected ",
        line_type == "Naturally Infected" & shape_type == "Unadjusted" ~ "Naturally Infected (Unadjusted)",
        line_type == "Population" ~ "Population"
      )
    )
  
  color_map <- c(
    "Naturally Infected " = "#00468b",
    #"Naturally Infected (Unadjusted)"        = "#00468b",
    "Population"                             = "#ed0000"
  )
  
  linetype_map <- c(
    "Naturally Infected " = "solid",
   # "Naturally Infected (Unadjusted)"        = "dotted",
    "Population"                             = "dashed"
  )
  
  ggplot(
    power_df,
    aes(
      x = n,
      y = power,
      color = color_label,
      group = color_label,
      linetype = color_label
    )
  ) +
    geom_line(linewidth = 1.1, alpha = 0.8) +
    geom_point(size = 3.5, alpha = 0.8) +
    
    facet_grid(
      Y_out ~ setting_name,
      scales = "fixed",
      labeller = labeller(
        Y_out = label_parsed,
        setting_name = label_value
      )
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::percent
    ) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 4),
      labels = scales::comma
    ) +
    
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    
    labs(
      x = "Sample size",
      y = "Power",
      color = "Estimand & Estimator",
      linetype = "Estimand & Estimator"
    ) +
    
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      strip.text.y = element_text(size = 14, face = "bold"),
      strip.text.x = element_text(size = 14, face = "bold"),
      panel.spacing = unit(8, "pt")
    ) +
    
    guides(
      color = guide_legend(nrow = 1),
      linetype = guide_legend(nrow = 1)
    )
}


# -----------------------------------------------


# ------------------------------------------------------------------
# Power heatmap at fixed sample size n = 10,000
# ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

# Configuration
config_settings <- c("default",
                     "base_12mo",
                     "gambia_default_VE_6mo",
                     "gambia_default_VE_12mo",
                     "peru_default_VE_6mo",
                     "peru_default_VE_12mo")

config_setting_names <- c("6 month trial - Overall",
                          "12 month trial - Overall",
                          "6 month trial - The Gambia",
                          "12 month trial - The Gambia",
                          "6 month trial - Peru",
                          "12 month trial - Peru")

# Fixed sample size
fixed_n <- 20000

# Collect all power data at fixed n
all_power_heatmap <- purrr::map2_dfr(
  config_settings,
  config_setting_names,
  ~ {
    results <- readRDS(here::here(paste0("results/", .x, "_evaluation_results.Rds")))
    
    power_df <- results$power_df |> 
      add_labels() |>
      filter(n == fixed_n) |>
      mutate(
        setting = .x,
        setting_name = .y
      )
    
    return(power_df)
  }
) 

# Create heatmap data - filter for specific estimator/estimand combinations
heatmap_data <- all_power_heatmap %>%
  # You can filter for specific combinations here
  # For example, just AIPW with naturally infected estimand:
  filter(
    estimator == "aipw",
    estimand == "nat_inf"
  ) %>%
  mutate(
    # # Clean up outcome labels
    # Y_out_clean = case_when(
    #   Y_out == "Y[3]" ~ "Y₃",
    #   Y_out == "Y[6]" ~ "Y₆",
    #   Y_out == "Y[9]" ~ "Y₉",
    #   Y_out == "Y[12]" ~ "Y₁₂",
    #   Y_out == "Y[6-12]" ~ "Y₆₋₁₂",
    #   Y_out == "Y[3-6-9-12]" ~ "Y₃₋₆₋₉₋₁₂"
    # ),
    # # Order outcomes
    # Y_out_clean = factor(
    #   Y_out_clean,
    #   levels = c("Y₃", "Y₆", "Y₉", "Y₁₂", "Y₆₋₁₂", "Y₃₋₆₋₉₋₁₂")
    # ),
    # Order settings
    setting_name = factor(
      setting_name,
      levels = rev(config_setting_names)
    )
  )

# Create heatmap
power_heatmap <- ggplot(
  heatmap_data,
  aes(x = Y_out, y = setting_name, fill = power)
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = sprintf("%.1f%%", power * 100)),
    color = "white",
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = "#d73027",
    mid = "#fee090",
    high = "#1a9850",
    midpoint = 0.3,
    limits = c(0, 0.6),
    breaks = seq(0, 0.6, by = 0.1),
    labels = scales::percent_format(accuracy = 1),
    name = "Power"
  ) +
  labs(
    x = "Outcome",
    y = "Trial Setting",
    title = paste0("Naturally infected AIPW estimator power"),
    subtitle = paste0("n = ", scales::comma(fixed_n))
  ) +
  theme_minimal(base_size = 14) +
  scale_x_discrete(labels = scales::label_parse()) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(face = "bold", size = 13),
    legend.position = "right",
    panel.grid = element_blank()
  )

# Save the plot
ggsave(
  filename = here::here(paste0("results/figures/power_heatmap_n", fixed_n, ".png")),
  plot = power_heatmap,
  width = 10,
  height = 6,
  dpi = 300
)

print(power_heatmap)

# Optional: Create separate heatmaps for different estimator/estimand combinations
create_heatmap_by_group <- function(data, est, estimand_type, n) {
  
  heatmap_subset <- data %>%
    filter(
      estimator == est,
      estimand == estimand_type
    ) %>%
    mutate(
      Y_out_clean = case_when(
        Y_out == "Y[3]" ~ "Y₃",
        Y_out == "Y[6]" ~ "Y₆",
        Y_out == "Y[9]" ~ "Y₉",
        Y_out == "Y[12]" ~ "Y₁₂",
        Y_out == "Y[6-12]" ~ "Y₆₋₁₂",
        Y_out == "Y[3-6-9-12]" ~ "Y₃₋₆₋₉₋₁₂"
      ),
      Y_out_clean = factor(
        Y_out_clean,
        levels = c("Y₃", "Y₆", "Y₉", "Y₁₂", "Y₆₋₁₂", "Y₃₋₆₋₉₋₁₂")
      ),
      setting_name = factor(setting_name, levels = config_setting_names)
    )
  
  estimand_label <- ifelse(estimand_type == "nat_inf", "Naturally Infected", "Population")
  est_label <- ifelse(est == "aipw", "AIPW", ifelse(est == "gcomp", "G-Comp", "Unadjusted"))
  
  p <- ggplot(
    heatmap_subset,
    aes(x = Y_out_clean, y = setting_name, fill = power)
  ) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(
      aes(label = sprintf("%.1f%%", power * 100)),
      color = "white",
      size = 4.5,
      fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = "#d73027",
      mid = "#fee090",
      high = "#1a9850",
      midpoint = 0.4,
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::percent_format(accuracy = 1),
      name = "Power"
    ) +
    labs(
      x = "Outcome",
      y = "Trial Setting",
      title = paste0("Power at n = ", scales::comma(n)),
      subtitle = paste(est_label, "estimator,", estimand_label, "estimand")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(face = "bold", size = 13),
      legend.position = "right",
      panel.grid = element_blank()
    )
  
  ggsave(
    filename = here::here(paste0("results/figures/power_heatmap_", est, "_", estimand_type, "_n", n, ".png")),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

# Create heatmaps for all combinations
create_heatmap_by_group(all_power_heatmap, "aipw", "nat_inf", fixed_n)
create_heatmap_by_group(all_power_heatmap, "aipw", "pop", fixed_n)
create_heatmap_by_group(all_power_heatmap, "gcomp", "nat_inf", fixed_n)
create_heatmap_by_group(all_power_heatmap, "unadj", "nat_inf", fixed_n)

# ----------------------------------------------------------------------

# Make figure of the year version for manuscript

source(here::here("R/parameter_generation_fns.R"))

# this is model for all-ages
shigella_growth_meta_analysis_results <- readRDS(here::here("misc/results/case_control/shigella_growth_effect_data.Rds"))

# Original, unscaled
plot_df <- shigella_growth_meta_analysis_results %>% 
  filter(group != "Any Shigella")

dose_schedule <- "6mo"
scale_growth_effect_0_6 <- -0.25
scale_growth_effect_6_12 <- 0.25
spline_formula <- "y ~ -1 + x + I(pmax(0, x - 4)) + I(pmax(0, x - 8))"

scale_0_6_df <- plot_df %>%
  mutate(
    age_strata = if_else(dose_schedule == "6mo", "6-12 months", "12-18 months"),
    scale_growth_effect_0_6 = scale_growth_effect_0_6
  ) %>%
  mutate(
    pt_est = if_else(scale_growth_effect_0_6 > 0, pt_est + scale_growth_effect_0_6*(upper_ci - pt_est), pt_est + scale_growth_effect_0_6*(pt_est - lower_ci)) 
  ) %>% 
  filter(group != "Any Shigella")

scale_0_6_fits <- scale_0_6_df %>%
  group_by(age_strata, group) %>%
  group_split() %>%
  map( ~ {
    df <- .x
    x <- df$month_num
    y <- df$pt_est
    glm(formula = spline_formula, data = data.frame(x = x, y = y))
  })

names(scale_0_6_fits) <- c("lsd", "msd")

scale_6_12_df <- plot_df %>%
  mutate(
    age_strata = if_else(dose_schedule == "6mo", "12-18 months", "18-24 months"),
    scale_growth_effect_6_12 = scale_growth_effect_6_12
  ) %>%
  mutate(
    pt_est = if_else(scale_growth_effect_6_12 > 0, pt_est + scale_growth_effect_6_12*(upper_ci - pt_est), pt_est + scale_growth_effect_6_12*(pt_est - lower_ci))
  )  %>% 
  filter(group != "Any Shigella")

scale_6_12_fits <- scale_6_12_df %>%
  group_by(age_strata, group) %>%
  group_split() %>%
  map( ~ {
    df <- .x
    x <- df$month_num
    y <- df$pt_est
    glm(formula = spline_formula, data = data.frame(x = x, y = y))
  })

names(scale_6_12_fits) <- c("lsd", "msd")

plot_combined <- bind_rows(
  plot_df %>% mutate(line_type = "12-18 months", size_line = 1.5),
  scale_0_6_df %>% mutate(line_type = "6-12 months", size_line = 0.7),
  scale_6_12_df %>% mutate(line_type = "18-24 months", size_line = 0.7)
)

# Ensure age_strata is consistent for plotting
plot_combined <- plot_combined %>%
  mutate(age_strata = factor(age_strata, levels = c("6-12 months", "12-18 months", "18-24 months")))

# Rename line_type for legend labels
plot_combined <- plot_combined %>%
  mutate(line_type = factor(line_type, levels = c("6-12 months", "12-18 months", "18-24 months")))

facet_plot <- ggplot(plot_combined, aes(x = month_num, y = pt_est, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
  geom_line(data = filter(plot_combined, line_type == "12-18 months"), 
            aes(size = 1.5), linetype = "solid") +
  geom_point(data = filter(plot_combined, line_type == "12-18 months"), size = 4) +
  scale_size_identity() +
  scale_color_manual(values = c(
    "Moderate-to-Severe Shigella" = "#ED0000FF",
    "Less-severe Shigella" = "#42B540FF"
  )) +
  scale_fill_manual(values = c(
    "Moderate-to-Severe Shigella" = "#ED0000FF",
    "Less-severe Shigella" = "#42B540FF"
  )) +
  scale_x_continuous(breaks = 1:12) +
  # Separate splines for each line_type with legend
  stat_smooth(
    data = filter(plot_combined, line_type == "12-18 months"),
    method = "glm",
    formula = spline_formula,
    aes(linetype = "12-18 months"),
    size = 1,
    se = FALSE
  ) +
  stat_smooth(
    data = filter(plot_combined, line_type == "6-12 months"),
    method = "glm",
    formula = spline_formula,
    aes(linetype = "6-12 months"),
    size = 1,
    se = FALSE
  ) +
  stat_smooth(
    data = filter(plot_combined, line_type == "18-24 months"),
    method = "glm",
    formula = spline_formula,
    aes(linetype = "18-24 months"),
    size = 1,
    se = FALSE
  ) +
  scale_linetype_manual(
    name = "Age Group",
    values = c("6-12 months" = "dotted", "12-18 months" = "dashed", "18-24 months" = "dotdash"),
    guide = guide_legend(override.aes = list(size = 0.5, color = "black"))
  ) +
  labs(
    x = "Month",
    y = "Growth Effect",
    # title = "Shigella growth effects in MAL-ED with observed antibiotic use",
    color = "Shigella severity",
    fill = "Shigella severity"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

facet_plot
