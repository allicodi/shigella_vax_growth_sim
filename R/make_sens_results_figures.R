# ----------------------------------------------------------------------------
# Make results figures for sensitivity analyses with unmeasured confounding
# -----------------------------------------------------------------------------

library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(patchwork)
library(stringr)
library(knitr)
library(purrr)
library(kableExtra)

here::i_am("R/make_sens_results_figures.R")

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

sens_settings <- c(
  "S_low_Y_low",
  "S_low_Y_high",
  "S_high_Y_low",
  "S_high_Y_high"
)

sens_setting_names <- c(
  "U decreases risk of infection, \ndecreases HAZ",
  "U decreases risk of infection, \nincreases HAZ",
  "U increases risk of infection, \ndecreases HAZ",
  "U increases risk of infection, \nincreases HAZ"
)

power_settings <- c(
  "optimistic_12mo",
  sens_settings
)

power_setting_names <- c(
  "Original",
  sens_setting_names
)

# -----------------------------------------------------------------------------
# Add labels for plotting
# -----------------------------------------------------------------------------

add_labels <- function(df) {
  df %>%
    mutate(
      estimand = recode(
        estimand,
        nat_inf = "Naturally Infected",
        pop = "Population"
      ),
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
        levels = c(
          "Y[3]", "Y[6]", "Y[9]", "Y[12]",
          "Y[6-12]", "Y[3-6-9-12]"
        )
      ),
      Y_out_text = recode(
        as.character(Y_out),
        "Y[3]" = "3 month outcome",
        "Y[6]" = "6 month outcome",
        "Y[9]" = "9 month outcome",
        "Y[12]" = "12 month outcome",
        "Y[6-12]" = "6 and 12 month outcome",
        "Y[3-6-9-12]" = "3, 6, 9, and 12 month outcome"
      ),
      Y_out_text = factor(
        Y_out_text,
        levels = c(
          "3 month outcome",
          "6 month outcome",
          "9 month outcome",
          "12 month outcome",
          "6 and 12 month outcome",
          "3, 6, 9, and 12 month outcome"
        )
      )
    )
}

# -----------------------------------------------------------------------------
# Individual sensitivity plots
# -----------------------------------------------------------------------------

plot_power <- function(power_df, setting = "default") {
  
  color_map <- c(
    "Naturally Infected" = "#00468b",
    "Population" = "#ed0000"
  )
  
  linetype_map <- c(
    "Naturally Infected" = "solid",
    "Population" = "dashed"
  )
  
  ggplot(
    power_df,
    aes(
      x = n,
      y = power,
      color = estimand,
      linetype = estimand,
      group = estimand
    )
  ) +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    geom_point(size = 3.5) +
    facet_wrap(~ Y_out_text, ncol = 3) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
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

plot_prop_neg <- function(prop_neg_df, setting = "default") {
  
  color_map <- c(
    "Naturally Infected" = "#00468b",
    "Population" = "#ed0000"
  )
  
  linetype_map <- c(
    "Naturally Infected" = "solid",
    "Population" = "dashed"
  )
  
  ggplot(
    prop_neg_df,
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
    facet_wrap(~ Y_out_text, ncol = 3) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
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
    aes(
      x = n,
      y = bias,
      color = estimand,
      linetype = estimand,
      group = estimand
    )
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
    facet_wrap(~ Y_out_text, ncol = 3) +
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
# Run individual sensitivity plots and collect truths
# -----------------------------------------------------------------------------

all_truth_df <- tibble()

for (i in seq_along(sens_settings)) {
  
  setting <- sens_settings[i]
  setting_name <- sens_setting_names[i]
  
  results <- readRDS(here::here(paste0("results/", setting, "_evaluation_results.Rds")))
  
  power_df <- add_labels(results$power_df)
  bias_df <- add_labels(results$bias_df)
  prop_neg_df <- add_labels(results$prop_neg_df)
  
  power_fig <- plot_power(power_df, setting = setting_name)
  bias_fig <- plot_bias(bias_df, setting = setting_name)
  prop_neg_fig <- plot_prop_neg(prop_neg_df, setting = setting_name)
  
  truth_df <- results$truth_df %>%
    tibble::enframe(name = "parameter", value = "value") %>%
    mutate(
      value = purrr::map_dbl(value, ~ .x),
      setting = setting,
      setting_name = setting_name
    ) %>%
    select(setting, setting_name, parameter, value)
  
  all_truth_df <- bind_rows(all_truth_df, truth_df)
  
  ggsave(
    filename = here::here(paste0("results/figures/individual/power_", setting, ".png")),
    plot = power_fig,
    width = 14,
    height = 7
  )
  
  ggsave(
    filename = here::here(paste0("results/figures/individual/bias_", setting, ".png")),
    plot = bias_fig,
    width = 14,
    height = 7
  )
  
  ggsave(
    filename = here::here(paste0("results/figures/individual/prop_neg_", setting, ".png")),
    plot = prop_neg_fig,
    width = 14,
    height = 7
  )
}

# -----------------------------------------------------------------------------
# Truth table for sensitivity settings
# -----------------------------------------------------------------------------

truth_long <- all_truth_df %>%
  mutate(
    Estimand = case_when(
      str_detect(parameter, "^nat_inf") ~ "Naturally Infected",
      str_detect(parameter, "^pop") ~ "Population"
    ),
    Outcome = case_when(
      str_detect(parameter, "Y_3_6_9_12") ~ "Y[3-6-9-12]",
      str_detect(parameter, "Y_6_12") ~ "Y[6-12]",
      str_detect(parameter, "Y_12") ~ "Y[12]",
      str_detect(parameter, "Y_9") ~ "Y[9]",
      str_detect(parameter, "Y_6") ~ "Y[6]",
      str_detect(parameter, "Y_3") ~ "Y[3]"
    ),
    Outcome = factor(
      Outcome,
      levels = c("Y[3]", "Y[6]", "Y[9]", "Y[12]", "Y[6-12]", "Y[3-6-9-12]")
    )
  )

true_effect_df <- truth_long %>%
  select(Estimand, Outcome, setting_name, value) %>%
  pivot_wider(
    names_from = setting_name,
    values_from = value
  ) %>%
  arrange(Estimand, Outcome)

fmt_num <- function(x) {
  ifelse(
    is.na(x),
    NA,
    ifelse(abs(x) < 0.001, "$<$0.001", sprintf("%.3f", x))
  )
}

true_effect_table <- true_effect_df %>%
  mutate(
    Outcome = case_when(
      Outcome == "Y[3]" ~ "$Y_{3}$",
      Outcome == "Y[6]" ~ "$Y_{6}$",
      Outcome == "Y[9]" ~ "$Y_{9}$",
      Outcome == "Y[12]" ~ "$Y_{12}$",
      Outcome == "Y[6-12]" ~ "$Y_{6\\text{--}12}$",
      Outcome == "Y[3-6-9-12]" ~ "$Y_{3\\text{--}6\\text{--}9\\text{--}12}$"
    ),
    across(where(is.numeric), fmt_num)
  ) %>%
  kable(
    format = "latex",
    booktabs = FALSE,
    escape = FALSE,
    align = "c"
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

true_effect_table

# -----------------------------------------------------------------------------
# Main figure: original optimistic_12mo vs sensitivity power curves
# -----------------------------------------------------------------------------

all_power_sens_df <- pmap_dfr(
  list(
    setting = power_settings,
    setting_name = power_setting_names
  ),
  function(setting, setting_name) {
    
    results <- readRDS(here::here(paste0("results/", setting, "_evaluation_results.Rds")))
    
    results$power_df %>%
      add_labels() %>%
      mutate(
        setting = setting,
        setting_name = setting_name
      )
  }
) %>%
  mutate(
    setting_name = factor(setting_name, levels = power_setting_names),
    setting_type = if_else(setting == "optimistic_12mo", "Original", "Sensitivity")
  )

plot_power_original_vs_sens <- function(df,
                                        estimand_keep = "Naturally Infected",
                                        endpoints_keep = NULL,
                                        plot_title = NULL) {
  
  plot_df <- df %>%
    filter(estimand == estimand_keep)
  
  if (!is.null(endpoints_keep)) {
    plot_df <- plot_df %>%
      filter(Y_out %in% endpoints_keep)
  }
  
  ggplot(
    plot_df,
    aes(
      x = n,
      y = power,
      color = setting_name,
      linetype = setting_name,
      group = setting_name
    )
  ) +
    geom_line(linewidth = 1.1, alpha = 0.95) +
    geom_point(size = 3.5) +
    facet_wrap(~ Y_out_text, ncol = 3) +
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
      color = "Setting",
      linetype = "Setting",
      title = plot_title
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      strip.text = element_text(size = 13, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# All endpoints, naturally infected
power_nat_inf_sens_fig <- plot_power_original_vs_sens(
  all_power_sens_df,
  estimand_keep = "Naturally Infected",
  plot_title = "Power curves: original vs unmeasured confounding sensitivity analyses\nNaturally infected estimand"
)

# All endpoints, population
power_pop_sens_fig <- plot_power_original_vs_sens(
  all_power_sens_df,
  estimand_keep = "Population",
  plot_title = "Power curves: original vs unmeasured confounding sensitivity analyses\nPopulation estimand"
)

ggsave(
  here::here("results/figures/power_original_vs_sens_nat_inf_all_endpoints.png"),
  plot = power_nat_inf_sens_fig,
  width = 12,
  height = 8
)

ggsave(
  here::here("results/figures/power_original_vs_sens_pop_all_endpoints.png"),
  plot = power_pop_sens_fig,
  width = 12,
  height = 8
)

# Main endpoints only, naturally infected
power_nat_inf_sens_main_fig <- plot_power_original_vs_sens(
  all_power_sens_df,
  estimand_keep = "Naturally Infected",
  endpoints_keep = c("Y[6]", "Y[12]", "Y[6-12]"),
  plot_title = "Power curves: original vs unmeasured confounding sensitivity analyses\nNaturally infected estimand"
)

# Main endpoints only, population
power_pop_sens_main_fig <- plot_power_original_vs_sens(
  all_power_sens_df,
  estimand_keep = "Population",
  endpoints_keep = c("Y[6]", "Y[12]", "Y[6-12]"),
  plot_title = "Power curves: original vs unmeasured confounding sensitivity analyses\nPopulation estimand"
)

ggsave(
  here::here("results/figures/power_original_vs_sens_nat_inf_main_endpoints.png"),
  plot = power_nat_inf_sens_main_fig,
  width = 12,
  height = 6
)

ggsave(
  here::here("results/figures/power_original_vs_sens_pop_main_endpoints.png"),
  plot = power_pop_sens_main_fig,
  width = 12,
  height = 6
)

# -----------------------------------------------------------------------------
# Bias figure: original optimistic_12mo vs sensitivity bias curves
# -----------------------------------------------------------------------------

all_bias_sens_df <- pmap_dfr(
  list(
    setting = power_settings,
    setting_name = power_setting_names
  ),
  function(setting, setting_name) {
    
    results <- readRDS(here::here(paste0("results/", setting, "_evaluation_results.Rds")))
    
    results$bias_df %>%
      add_labels() %>%
      mutate(
        setting = setting,
        setting_name = setting_name
      )
  }
) %>%
  mutate(
    setting_name = factor(setting_name, levels = power_setting_names),
    setting_type = if_else(setting == "optimistic_12mo", "Original", "Sensitivity")
  )

plot_bias_original_vs_sens <- function(df,
                                       estimand_keep = "Naturally Infected",
                                       endpoints_keep = NULL,
                                       plot_title = NULL) {
  
  plot_df <- df %>%
    filter(estimand == estimand_keep)
  
  if (!is.null(endpoints_keep)) {
    plot_df <- plot_df %>%
      filter(Y_out %in% endpoints_keep)
  }
  
  ggplot(
    plot_df,
    aes(
      x = n,
      y = bias,
      color = setting_name,
      linetype = setting_name,
      group = setting_name
    )
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_line(linewidth = 1.1, alpha = 0.95) +
    geom_point(size = 3.5) +
    facet_wrap(~ Y_out_text, ncol = 3) +
    scale_x_log10(
      breaks = c(2500, 5000, 10000, 20000, 40000, 80000),
      labels = scales::comma
    ) +
    labs(
      x = "Sample Size",
      y = "Bias",
      color = "Setting",
      linetype = "Setting",
      title = plot_title
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      strip.text = element_text(size = 13, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# All endpoints, naturally infected
bias_nat_inf_sens_fig <- plot_bias_original_vs_sens(
  all_bias_sens_df,
  estimand_keep = "Naturally Infected",
  plot_title = "Bias curves: original vs unmeasured confounding sensitivity analyses\nNaturally infected estimand"
)

# All endpoints, population
bias_pop_sens_fig <- plot_bias_original_vs_sens(
  all_bias_sens_df,
  estimand_keep = "Population",
  plot_title = "Bias curves: original vs unmeasured confounding sensitivity analyses\nPopulation estimand"
)

ggsave(
  here::here("results/figures/bias_original_vs_sens_nat_inf_all_endpoints.png"),
  plot = bias_nat_inf_sens_fig,
  width = 12,
  height = 8
)

ggsave(
  here::here("results/figures/bias_original_vs_sens_pop_all_endpoints.png"),
  plot = bias_pop_sens_fig,
  width = 12,
  height = 8
)

# Main endpoints only, naturally infected
bias_nat_inf_sens_main_fig <- plot_bias_original_vs_sens(
  all_bias_sens_df,
  estimand_keep = "Naturally Infected",
  endpoints_keep = c("Y[6]", "Y[12]", "Y[6-12]"),
  plot_title = "Bias curves: original vs unmeasured confounding sensitivity analyses\nNaturally infected estimand"
)

# Main endpoints only, population
bias_pop_sens_main_fig <- plot_bias_original_vs_sens(
  all_bias_sens_df,
  estimand_keep = "Population",
  endpoints_keep = c("Y[6]", "Y[12]", "Y[6-12]"),
  plot_title = "Bias curves: original vs unmeasured confounding sensitivity analyses\nPopulation estimand"
)

ggsave(
  here::here("results/figures/bias_original_vs_sens_nat_inf_main_endpoints.png"),
  plot = bias_nat_inf_sens_main_fig,
  width = 12,
  height = 6
)

ggsave(
  here::here("results/figures/bias_original_vs_sens_pop_main_endpoints.png"),
  plot = bias_pop_sens_main_fig,
  width = 12,
  height = 6
)

# -----------------------------------------------------------------------------
# Combined 2 x 2 figure: power and bias for 12-month endpoint
# -----------------------------------------------------------------------------

plot_power_12mo <- function(df, estimand_keep) {
  
  df %>%
    filter(
      estimand == estimand_keep,
      Y_out == "Y[12]"
    ) %>%
    ggplot(
      aes(
        x = n,
        y = power,
        color = setting_name,
        linetype = setting_name,
        group = setting_name
      )
    ) +
    geom_line(linewidth = 1.2, alpha = 0.95) +
    geom_point(size = 3.5) +
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
      title = NULL,
      x = NULL,
      y = "Power to detect vaccine effect",
      color = "Setting",
      linetype = "Setting"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "bottom",
      # plot.title = element_text(
      #   hjust = 0.5,
      #   face = "bold"
      # ),
      panel.grid.minor = element_blank()
    )
}

# Common symmetric y-axis limits across both 12-month bias panels
bias_12mo_limit <- all_bias_sens_df %>%
  filter(Y_out == "Y[12]") %>%
  summarise(max_abs_bias = max(abs(bias), na.rm = TRUE)) %>%
  pull(max_abs_bias)

bias_12mo_limits <- c(-bias_12mo_limit, bias_12mo_limit)

plot_bias_12mo <- function(df, estimand_keep, y_limits) {
  
  df %>%
    filter(
      estimand == estimand_keep,
      Y_out == "Y[12]"
    ) %>%
    ggplot(
      aes(
        x = n,
        y = bias,
        color = setting_name,
        linetype = setting_name,
        group = setting_name
      )
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "gray40"
    ) +
    geom_line(linewidth = 1.2, alpha = 0.95) +
    geom_point(size = 3.5) +
    scale_y_continuous(
      limits = y_limits
    ) +
    scale_x_log10(
      breaks = c(2500, 5000, 10000, 20000, 40000, 80000),
      labels = scales::comma
    ) +
    labs(
      title = estimand_keep,
      x = "Sample Size",
      y = "Bias",
      color = "Setting",
      linetype = "Setting"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      panel.grid.minor = element_blank()
    )
}


power_nat_inf_12mo <- plot_power_12mo(
  all_power_sens_df,
  estimand_keep = "Naturally Infected"
)

power_pop_12mo <- plot_power_12mo(
  all_power_sens_df,
  estimand_keep = "Population"
)

bias_nat_inf_12mo <- plot_bias_12mo(
  all_bias_sens_df,
  estimand_keep = "Naturally Infected",
  y_limits = bias_12mo_limits
)

bias_pop_12mo <- plot_bias_12mo(
  all_bias_sens_df,
  estimand_keep = "Population",
  y_limits = bias_12mo_limits
)


power_bias_12mo_fig <-
  (
    bias_nat_inf_12mo +
      bias_pop_12mo
  ) /
  (
    power_nat_inf_12mo +
      power_pop_12mo
  ) +
  plot_layout(
    guides = "collect",
    axes = "collect"
  ) +
  plot_annotation(
    title = "Power and bias under unmeasured confounding sensitivity analyses",
    subtitle = "Targeted recruitment, 12 month immunization schedule, 12 month growth endpoint",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        size = 19,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 16,
        hjust = 0.5
      )
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12)
  )


ggsave(
  filename = here::here(
    "results/figures/power_bias_original_vs_sens_12mo.png"
  ),
  plot = power_bias_12mo_fig,
  width = 14,
  height = 10,
  dpi = 300
)
