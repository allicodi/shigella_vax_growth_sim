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

power_df <- results$power_df %>% 
  mutate(estimand = if_else(estimand == "nat_inf", "nat_inf_ER", "pop"))
prop_neg_df <- results$prop_neg_df %>% 
  mutate(estimand = if_else(estimand == "nat_inf", "nat_inf_ER", "pop"))
bias_df <- results$bias_df %>% 
  mutate(estimand = if_else(estimand == "nat_inf", "nat_inf_ER", "pop"))

# Combine other results
no_ER_res <- readRDS(here::here("results/debug_ER_evaluation_results.Rds"))
no_ER_power <- no_ER_res$power_df %>%
  mutate(estimand = "nat_inf_no_ER")
no_ER_prop_neg_df <- no_ER_res$prop_neg_df %>%
  mutate(estimand = "nat_inf_no_ER")
no_ER_bias_df <- no_ER_res$bias_df  %>%
  mutate(estimand = "nat_inf_no_ER")

unadj_res <- readRDS(here::here("results/debug_unadj_evaluation_results.Rds"))

two_stage_res <- readRDS(here::here("results/debug_two_part_evaluation_results.Rds"))
two_stage_power_df <- two_stage_res$power_df %>% 
  mutate(estimand = if_else(estimand == "nat_inf_ER", "nat_inf_ER_2", "pop_2"))
two_stage_prop_neg_df <- two_stage_res$prop_neg_df %>% 
  mutate(estimand = if_else(estimand == "nat_inf_ER", "nat_inf_ER_2", "pop_2"))
two_stage_bias_df <- two_stage_res$bias_df %>% 
  mutate(estimand = if_else(estimand == "nat_inf_ER", "nat_inf_ER_2", "pop_2"))

# combine all
bias_df <- bias_df %>%
  rbind(no_ER_bias_df) %>%
  rbind(unadj_res$bias_df) %>%
  rbind(two_stage_bias_df)

power_df <- power_df %>%
  rbind(no_ER_power) %>%
  rbind(unadj_res$power_df) %>%
  rbind(two_stage_power_df)

prop_neg_df <- prop_neg_df %>%
  rbind(no_ER_prop_neg_df) %>%
  rbind(unadj_res$prop_neg_df) %>%
  rbind(two_stage_prop_neg_df)

# -------------------------------------
# Function to create power plot
# -------------------------------------
plot_power <- function(power_df, setting) {
  
  power_df <- power_df %>%
    mutate(
      estimand = recode(
        estimand,
        nat_inf_ER = "Naturally Infected (ER)",
        nat_inf_no_ER = "Naturally Infected",
        nat_inf_ER_2 = "Naturally Infected (ER; two-part)",
        nat_inf_unadj = "Naturally Infected (unadjusted)",
        pop = "Population",
        pop_2 = "Population (two-part)"
      ),
      estimand_type = if_else(grepl("Population", estimand), "Population", "Naturally Infected"),
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
    "Naturally Infected" = "#2ca02c",
    "Naturally Infected (ER)" = "#98df8a",
    "Naturally Infected (ER; two-part)" = "#1f77b4",
    "Naturally Infected (unadjusted)" = "#aec7e8",
    "Population" = "#003f5c",
    "Population (two-part)" = "#7a5195"
  )
  
  ggplot(power_df, aes(x = n, y = power, color = estimand, group = estimand,
                       linetype = estimand_type)) +
    geom_line(size = 1.5) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.8, linetype = "dashed", size = 1, color = "red") +
    facet_wrap(~ Y_out, scales = "fixed", ncol = 3, labeller = label_parsed) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::percent_format(accuracy = 1),
      name = "Power"
    ) +
    scale_x_continuous(
      name = "Sample Size",
      labels = scales::comma,
      limits = c(0, 50000),
      breaks = seq(0, 50000, 10000)
    ) +
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = c("Naturally Infected" = "solid", "Population" = "dotted")) +
    labs(
      title = paste0("Power (", setting, " setting)"),
      color = "Estimand",
      linetype = "Type"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 18, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}



# -------------------------------------------------
# Function to create proportion negative plot
# -------------------------------------------------
plot_prop_neg <- function(prop_neg_df, setting) {
  
  prop_neg_df <- prop_neg_df %>%
    mutate(
      estimand = recode(
        estimand,
        nat_inf_ER = "Naturally Infected (ER)",
        nat_inf_no_ER = "Naturally Infected",
        nat_inf_ER_2 = "Naturally Infected (ER; two-part)",
        nat_inf_unadj = "Naturally Infected (unadjusted)",
        pop = "Population",
        pop_2 = "Population (two-part)"
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
        levels = c("Y[3]", "Y[6]", "Y[9]", "Y[12]", "Y[6-12]", "Y[3-6-9-12]")
      )
    )
  
  color_map <- c(
    "Naturally Infected" = "#2ca02c",
    "Naturally Infected (ER)" = "#98df8a",
    "Naturally Infected (ER; two-part)" = "#1f77b4",
    "Naturally Infected (unadjusted)" = "#aec7e8",
    "Population" = "#003f5c",
    "Population (two-part)" = "#7a5195"
  )
  
  ggplot(prop_neg_df, aes(x = n, y = prop_neg, color = estimand, group = estimand)) +
    geom_line(size = 1.5) +
    geom_point(size = 3) +
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
    labs(
      title = paste0("Proportion Negative (", setting, " setting)"),
      color = "Estimand"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 18, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}


# -------------------------------------------------
# Function to create bias plot
# -------------------------------------------------
plot_bias <- function(power_df, setting) {
  
  # same transformation
  power_df <- power_df %>%
    mutate(
      estimand = recode(
        estimand,
        nat_inf_ER = "Naturally Infected (ER)",
        nat_inf_no_ER = "Naturally Infected",
        nat_inf_ER_2 = "Naturally Infected (ER; two-part)",
        nat_inf_unadj = "Naturally Infected (unadjusted)",
        pop = "Population",
        pop_2 = "Population (two-part)"
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
        levels = c("Y[3]", "Y[6]", "Y[9]", "Y[12]", "Y[6-12]", "Y[3-6-9-12]")
      )
    )
  
  color_map <- c(
    "Naturally Infected" = "#2ca02c",
    "Naturally Infected (ER)" = "#98df8a",
    "Naturally Infected (ER; two-part)" = "#1f77b4",
    "Naturally Infected (unadjusted)" = "#aec7e8",
    "Population" = "#003f5c",
    "Population (two-part)" = "#7a5195"
  )
  
  ggplot(power_df, aes(x = n, y = bias, color = estimand, group = estimand)) +
    geom_line(size = 1.5) +
    geom_point(size = 3) +
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
    labs(
      title = paste0("Bias (", setting, " setting)"),
      color = "Estimand"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 18, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}


