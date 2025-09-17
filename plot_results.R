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

# -------------------------------------
# Function to create power plot
# -------------------------------------
plot_power <- function(power_df, setting){
  power_results_long<- power_df %>%
    select(n, power_gcomp, power_gcomp_pop_estimand) %>%
    pivot_longer(cols = starts_with("power_"),
                 names_to = "estimand",
                 values_to = "power") %>%
    mutate(estimand = recode(estimand,
                             power_gcomp = "Naturally Infected",
                             power_gcomp_pop_estimand = "Population"))
  
  power_fig <- ggplot(power_results_long, aes(x = n, y = power, color = estimand)) +
    geom_line(size = 3) +
    geom_point(size = 8) +
    geom_hline(yintercept = 0.8, linetype = "dashed", size = 1.5, color = "red") +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = percent_format(accuracy = 1),
      name = "Power"
    ) +
    scale_x_continuous(name = "Sample Size",
                       labels = scales::comma,
                       limits = c(0, 50000),
                       breaks = seq(0,50000,10000)) +
    scale_color_manual(
      values = c(
        "Naturally Infected" = "#2ca02c",  # green
        "Population" = "#003f5c"                  # dark navy blue
      )
    ) +
    labs(
      title = "Default parameterization",
      color = "Estimand"
    ) +
    theme_minimal(base_size = 24) +
    theme(legend.position = "bottom")
  
  return(power_fig)
}

# -------------------------------------------------
# Function to create proportion negative plot
# -------------------------------------------------
plot_prop_neg <- function(prop_neg_df, setting){
  prop_neg_results_long <- prop_neg_df %>%
    select(n, prop_neg_gcomp, prop_neg_gcomp_pop_estimand) %>%
    pivot_longer(cols = starts_with("prop_neg_"),
                 names_to = "estimand",
                 values_to = "prop_neg") %>%
    mutate(estimand = recode(estimand,
                             prop_neg_gcomp = "Naturally Infected",
                             prop_neg_gcomp_pop_estimand = "Population"))
  
  prop_neg_fig <- ggplot(prop_neg_results_long, aes(x = n, y = prop_neg, color = estimand)) +
    geom_line(size = 3) +
    geom_point(size = 8) +
    scale_y_continuous(
      limits = c(0, 0.6),
      breaks = seq(0, 0.6, by = 0.1),
      labels = percent_format(accuracy = 1),
      name = "Negative Point Estimates"
    ) +
    scale_x_continuous(name = "Sample Size",
                       labels = scales::comma,
                       limits = c(0, 50000),
                       breaks = seq(0,50000,10000)) +
    scale_color_manual(
      values = c(
        "Naturally Infected" = "#2ca02c",  # green
        "Population" = "#003f5c"                  # dark navy blue
      )
    ) +
    labs(
      color = "Estimand"
    ) +
    theme_minimal(base_size = 24) +
    theme(legend.position = "bottom")
  
  return(prop_neg_fig)
}

# From SER poster ------------------------------------

# 
# 
# # add_growth_default_bg_results as 'default'
# # optimistic results second panel
# default_results <- readRDS("results/add_growth_default_bg_results.Rds")
# optimistic_results <- readRDS("results/inc_shig_coef_ideal_results.Rds")
# 
# power_df_default <- default_results$power_df
# power_df_optimistic <- optimistic_results$power_df
# 
# prop_neg_default <- default_results$prop_neg_df
# prop_neg_optimistic <- optimistic_results$prop_neg_df
# 
# # Assuming your dataset is called power_df
# # Reshape the data to long format
# power_results_long_default <- power_df_default %>%
#   select(n, power_gcomp, power_gcomp_pop_estimand) %>%
#   pivot_longer(cols = starts_with("power_"),
#                names_to = "estimand",
#                values_to = "power") %>%
#   mutate(estimand = recode(estimand,
#                            power_gcomp = "Naturally Infected",
#                            power_gcomp_pop_estimand = "Population"))
# 
# power_results_long_optimistic <- power_df_optimistic %>%
#   select(n, power_gcomp, power_gcomp_pop_estimand) %>%
#   pivot_longer(cols = starts_with("power_"),
#                names_to = "estimand",
#                values_to = "power") %>%
#   mutate(estimand = recode(estimand,
#                            power_gcomp = "Naturally Infected",
#                            power_gcomp_pop_estimand = "Population"))
# 
# prop_neg_results_long_default <- prop_neg_default %>%
#   select(n, prop_neg_gcomp, prop_neg_gcomp_pop_estimand) %>%
#   pivot_longer(cols = starts_with("prop_neg_"),
#                names_to = "estimand",
#                values_to = "prop_neg") %>%
#   mutate(estimand = recode(estimand,
#                            prop_neg_gcomp = "Naturally Infected",
#                            prop_neg_gcomp_pop_estimand = "Population"))
# 
# prop_neg_results_long_optimistic <- prop_neg_optimistic %>%
#   select(n, prop_neg_gcomp, prop_neg_gcomp_pop_estimand) %>%
#   pivot_longer(cols = starts_with("prop_neg_"),
#                names_to = "estimand",
#                values_to = "prop_neg") %>%
#   mutate(estimand = recode(estimand,
#                            prop_neg_gcomp = "Naturally Infected",
#                            prop_neg_gcomp_pop_estimand = "Population"))
# 
# default_prop_neg_fig <- ggplot(prop_neg_results_long_default, aes(x = n, y = prop_neg, color = estimand)) +
#   geom_line(size = 3) +
#   geom_point(size = 8) +
#   scale_y_continuous(
#     limits = c(0, 0.6),
#     breaks = seq(0, 0.6, by = 0.1),
#     labels = percent_format(accuracy = 1),
#     name = "Negative Point Estimates"
#   ) +
#   scale_x_continuous(name = "Sample Size",
#                      labels = scales::comma,
#                      limits = c(0, 50000),
#                      breaks = seq(0,50000,10000)) +
#   scale_color_manual(
#     values = c(
#       "Naturally Infected" = "#2ca02c",  # green
#       "Population" = "#003f5c"                  # dark navy blue
#     )
#   ) +
#   labs(
#     color = "Estimand"
#   ) +
#   theme_minimal(base_size = 24) +
#   theme(legend.position = "bottom")
# 
# ggsave("results/default_prop_neg.png", default_prop_neg_fig, width = 10, height = 6, dpi = 300, bg = "transparent")
# 
# optimistic_prop_neg_fig <- ggplot(prop_neg_results_long_optimistic, aes(x = n, y = prop_neg, color = estimand)) +
#   geom_line(size = 3) +
#   geom_point(size = 8) +
#   scale_y_continuous(
#     limits = c(0, 0.6),
#     breaks = seq(0, 0.6, by = 0.1),
#     labels = percent_format(accuracy = 1),
#     name = "Negative Point Estimates"
#   ) +
#   scale_x_continuous(name = "Sample Size",
#                      labels = scales::comma,
#                      limits = c(0, 50000),
#                      breaks = seq(0,50000,10000)) +
#   scale_color_manual(
#     values = c(
#       "Naturally Infected" = "#2ca02c",  # green
#       "Population" = "#003f5c"                  # dark navy blue
#     )
#   ) +
#   labs(
#     color = "Estimand"
#   ) +
#   theme_minimal(base_size = 24) +
#   theme(legend.position = "bottom")
# 
# ggsave("results/optimistic_prop_neg.png", optimistic_prop_neg_fig, width = 10, height = 6, dpi = 300, bg = "transparent")
# 
# # Plot
# default_power_fig <- ggplot(power_results_long_default, aes(x = n, y = power, color = estimand)) +
#   geom_line(size = 3) +
#   geom_point(size = 8) +
#   geom_hline(yintercept = 0.8, linetype = "dashed", size = 1.5, color = "red") +
#   scale_y_continuous(
#     limits = c(0, 1),
#     breaks = seq(0, 1, by = 0.2),
#     labels = percent_format(accuracy = 1),
#     name = "Power"
#   ) +
#   scale_x_continuous(name = "Sample Size",
#                      labels = scales::comma,
#                      limits = c(0, 50000),
#                      breaks = seq(0,50000,10000)) +
#   scale_color_manual(
#     values = c(
#       "Naturally Infected" = "#2ca02c",  # green
#       "Population" = "#003f5c"                  # dark navy blue
#     )
#   ) +
#   labs(
#     title = "Default parameterization",
#     color = "Estimand"
#   ) +
#   theme_minimal(base_size = 24) +
#   theme(legend.position = "bottom")
# 
# ggsave("results/default_power_plot.png", default_power_fig, width = 10, height = 6, dpi = 300, bg = "transparent")
# 
# optimistic_power_fig <- ggplot(power_results_long_optimistic, aes(x = n, y = power, color = estimand)) +
#   geom_line(size = 3) +
#   geom_point(size = 8) +
#   geom_hline(yintercept = 0.8, linetype = "dashed", size = 1.5, color = "red") +
#   scale_y_continuous(
#     limits = c(0, 1),
#     breaks = seq(0, 1, by = 0.2),
#     labels = percent_format(accuracy = 1),
#     name = "Power"
#   ) +
#   scale_x_continuous(name = "Sample Size",
#                      labels = scales::comma,
#                      limits = c(0, 50000),
#                      breaks = seq(0,50000,10000)) +
#   scale_color_manual(
#     values = c(
#       "Naturally Infected" = "#2ca02c",  # green
#       "Population" = "#003f5c"                  # dark navy blue
#     )
#   ) +
#   labs(
#     title = "Optimistic parameterization",
#     color = "Estimand"
#   ) +
#   theme_minimal(base_size = 24) +
#   theme(legend.position = "bottom")
# 
# ggsave("results/optimistic_power_plot.png", optimistic_power_fig, width = 10, height = 6, dpi = 300, bg = "transparent")
# 
# combined_plot_prop_neg <- default_prop_neg_fig + optimistic_prop_neg_fig + plot_layout(guides = "collect") & theme(legend.position = "bottom")
# 
# combined_plot <- default_power_fig + default_prop_neg_fig + plot_layout(guides = "collect") & theme(legend.position = "bottom")
# 
# # Save
# ggsave("results/default_power_and_neg_prop_plot.png", combined_plot, width = 20, height = 8, dpi = 300, bg = "transparent")
