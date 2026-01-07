# --------------------------------------------------------------------
# Figure for naive vs naturally infected estimand growth trajectory
# --------------------------------------------------------------------

here::i_am("misc/growth_curve_figure.R")

library(tidyverse)
library(patchwork)

source(here::here("R/parameter_generation_fns.R"))
source(here::here("R/simulate_data.R"))

parameters <- readRDS(here::here("parameters_default.Rds"))

data <- simulate_data(parameters, 
                      n = 5000,
                      VE_mild = 0.4,
                      VE_severe = 0.6, 
                      seed = 123,
                      type = "counterfactual")

data <- data %>%
  # Label strata
  mutate(PS = case_when(
    S_inf_Z0 == 0 & S_inf_Z1 == 0 ~ "Immune",
    S_inf_Z0 == 1 & S_inf_Z1 == 0 ~ "Protected",
    S_inf_Z0 == 1 & S_inf_Z1 == 1 ~ "Doomed"
  )) %>% 
  # Filter to people with baseline HAZ within one SD of the mean (to subset and ideally make curves less extreme)
  filter(X > parameters$mean_X - parameters$sd_X &
         X < parameters$mean_X + parameters$sd_X) %>%
  # cherry pick inf times
  filter(S_inf_time_Z0 < 20 | S_inf_time_Z0 == 52)

# set.seed(5) 
set.seed(3) 
sampled_data <- bind_rows(
  data %>% filter(PS == "Immune")    %>% slice_sample(n = 2),
  data %>% filter(PS == "Protected") %>% slice_sample(n = 2),
  data %>% filter(PS == "Doomed")    %>% slice_sample(n = 1)
) 

# For annotations
pop_outcome_vax_12 <- mean(sampled_data$Y_12_Z1)
pop_outcome_placebo_12 <- mean(sampled_data$Y_12_Z0)

pop_outcome_vax_6 <- mean(sampled_data$Y_6_Z1)
pop_outcome_placebo_6 <- mean(sampled_data$Y_6_Z0)

nat_inf_outcome_vax_12 <- mean(sampled_data$Y_12_Z1[sampled_data$PS == "Doomed" | sampled_data$PS == "Protected"])
nat_inf_outcome_placebo_12 <- mean(sampled_data$Y_12_Z0[sampled_data$PS == "Doomed" | sampled_data$PS == "Protected"])

nat_inf_outcome_vax_6 <- mean(sampled_data$Y_6_Z1[sampled_data$PS == "Doomed" | sampled_data$PS == "Protected"])
nat_inf_outcome_placebo_6 <- mean(sampled_data$Y_6_Z0[sampled_data$PS == "Doomed" | sampled_data$PS == "Protected"])

annotate_df <- data.frame(month = rep(c(6,12),2),
                          Z_cf = c(0,0,1,1),
                          pop = c(pop_outcome_placebo_6, pop_outcome_placebo_12,
                                  pop_outcome_vax_6, pop_outcome_vax_12),
                          nat_inf = c(nat_inf_outcome_placebo_6, nat_inf_outcome_placebo_12,
                                      nat_inf_outcome_vax_6, nat_inf_outcome_vax_12))
                          
long_df <- sampled_data %>%
  pivot_longer(
    cols = matches("^Y_\\d+_Z[01]$"),
    names_to = c("month", "Z_cf"),
    names_pattern = "Y_(\\d+)_Z(\\d)",
    values_to = "Y"
  ) %>%
  mutate(
    month = as.integer(month),
    Z_cf = factor(Z_cf, levels = c("0", "1"),
                  labels = c("Placebo", "Vaccine"))
  ) 

baseline_df <- long_df %>%
  group_by(id, Z_cf) %>%
  slice(1) %>%          # take one row per (id, Z_cf)
  ungroup() %>%
  mutate(
    month = 0,
    Y = X
  )

long_df <- bind_rows(baseline_df, long_df) %>%
  arrange(id, Z_cf, month)


# Infection-time dots: one per id per Z_cf
inf_df <- long_df %>%
  group_by(id, Z_cf) %>%
  summarise(
    inf_month = floor(
      if_else(Z_cf == "Placebo",
              first(S_inf_time_Z0),
              first(S_inf_time_Z1)) / (52 / 12)
    ),
    .groups = "drop"
  ) %>%
  filter(inf_month < 12) %>%
  # get Y value at that month
  left_join(
    long_df %>% select(id, Z_cf, month, Y),
    by = c("id", "Z_cf", "inf_month" = "month")
  ) %>%
  rename(Y_inf = Y) %>%
  distinct(id, .keep_all = TRUE)

# Population dots from annotation df
pop_df <- annotate_df %>%
  mutate(
    Z_cf = factor(Z_cf, levels = c(0,1), labels = c("Placebo","Vaccine")),
    month = month,
    Y = pop
  )

# Population plot
p1 <- ggplot(long_df,
       aes(x = month, y = Y, group = interaction(id, Z_cf), color = Z_cf)) +
  
  # Lines
  geom_line(alpha = 0.5, linewidth = 3) +
  
  # Infection-time black triangles
  geom_point(
    data = inf_df,
    aes(x = inf_month, y = Y_inf, shape = "Infection time"),
    inherit.aes = FALSE,
    color = "black",
    size = 4
  ) +
  
  # Population outcome dots
  geom_point(
    data = pop_df,
    aes(x = month, y = Y, color = Z_cf, shape = "Outcome HAZ"),
    inherit.aes = FALSE,
    size = 4
  ) +
  
  # Shapes & colors
  scale_color_manual(
    values = c("Placebo" = "blue", "Vaccine" = "red"),
    name = "Counterfactual\nTreatment"
  ) +
  scale_shape_manual(
    name = "Point type",
    values = c("Infection time" = 17, "Outcome HAZ" = 16)
  ) +
  
  # Axes
  scale_x_continuous(breaks = 0:12, limits = c(0, 12)) +
  scale_y_continuous(limits = c(-3.25, 1)) +
  
  labs(
    x = "Month",
    y = "HAZ",
    title = "Population Estimand"
  ) +
  theme_minimal()


# Naturally infected plot
nat_inf <- long_df %>% filter(PS %in% c("Doomed", "Protected"))
inf_nat_df <- inf_df %>% filter(id %in% nat_inf$id)

# Population dots from annotation df
nat_inf_annotate_df <- annotate_df %>%
  mutate(
    Z_cf = factor(Z_cf, levels = c(0,1), labels = c("Placebo","Vaccine")),
    month = month,
    Y = nat_inf
  )

p2 <- ggplot(nat_inf,
       aes(x = month, y = Y, group = interaction(id, Z_cf), color = Z_cf)) +
  geom_line(alpha = 0.5, linewidth = 3) +
  geom_point(
    data = inf_nat_df,
    aes(x = inf_month, y = Y_inf, shape = "Infection time"),
    inherit.aes = FALSE,
    color = "black",
    size = 4
  ) +
  geom_point(
    data = nat_inf_annotate_df,
    aes(x = month, y = Y, color = Z_cf, shape = "Outcome HAZ"),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_color_manual(
    values = c("Placebo" = "blue", "Vaccine" = "red"),
    name = "Counterfactual\nTreatment"
  ) +
  scale_shape_manual(
    name = "Point type",
    values = c("Infection time" = 17, "Outcome HAZ" = 16)
  ) +
  # Axes
  scale_x_continuous(breaks = 0:12, limits = c(0, 12)) +
  scale_y_continuous(limits = c(-3.25, 1)) +
  labs(
    x = "Month",
    y = "HAZ",
    title = "Naturally Infected Estimand"
  ) +
  theme_minimal()

# -----------------------------
# Combine plots vertically
# -----------------------------
combined_plot <- p1 / p2 + 
  plot_layout(guides = "collect") +  # share legend
  plot_annotation(
    title = "Counterfactual Growth Trajectories",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

combined_plot
