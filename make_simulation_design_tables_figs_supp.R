# ------------------------------------------------------
# Script to make simulation design table for supplement
# ------------------------------------------------------

here::i_am("make_simulation_design_tables_figs_supp.R")

library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpattern)
library(kableExtra)

source(here::here("R/parameter_generation_fns.R"))

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
               "High early incidnece", "High early incidnece")

params_list <- vector("list", length = length(all_settings))
for(i in 1:length(all_settings)){
  setting <- all_settings[i]  
  params_list[[i]] <- readRDS(here::here(paste0("parameters/parameters_",setting,".Rds")))
}
names(params_list) <- all_settings

## 1. Baseline HAZ (Mean (sd))
haz_df <- data.frame(setting_name = c("General & Targeted recruitment (The Gambia)", "High early incidence (Peru)"),
                     haz_sd_6_12 = c(paste0(round(params_list[["default"]]$mean_X, 2), " (", round(params_list[["default"]]$sd_X, 2), ")"),
                                     #paste0(round(params_list[["optimistic_6mo"]]$mean_X, 2), " (", round(params_list[["optimistic_6mo"]]$sd_X, 2), ")"),
                                     paste0(round(params_list[["peru_6mo"]]$mean_X, 2), " (", round(params_list[["peru_6mo"]]$sd_X, 2), ")")),
                     haz_sd_12_18 = c(paste0(round(params_list[["base_12mo"]]$mean_X, 2), " (", round(params_list[["base_12mo"]]$sd_X, 2), ")"),
                                      #paste0(round(params_list[["optimistic_12mo"]]$mean_X, 2), " (", round(params_list[["optimistic_12mo"]]$sd_X, 2), ")"),
                                      paste0(round(params_list[["peru_12mo"]]$mean_X, 2), " (", round(params_list[["peru_12mo"]]$sd_X, 2), ")")))

haz_df %>%
  rename(
    `Setting` = setting_name,
    `6 month immunization schedule` = haz_sd_6_12,
    `12 month immunization schedule` = haz_sd_12_18
  ) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    align = c("l", "c", "c"),
    caption = "Baseline HAZ, mean (SD)"
  ) %>%
  column_spec(1, bold = TRUE)


## 2. Mean HAZ in absence of infection
# same model for all settings, just plot for all 24 months starting at HAZ 0 

monthly_growth_model <- params_list$default$monthly_growth_model
X <- 0
months <- 6:24

beta_0 <- monthly_growth_model$beta_0[months]
beta_1 <- monthly_growth_model$beta_1[months]

Y_vec <- numeric(length(months))
Y_vec[1] <- X

if(length(months) > 1){
  for(i in 2:length(months)) {
    Y_vec[i] <- beta_0[i] + beta_1[i] * Y_vec[i-1]
  }
}

plot_df <- data.frame(
  month = months,
  Y = Y_vec
)

growth_trajectory_fig <- ggplot(plot_df, aes(x = month, y = Y)) +
  geom_line(
    linewidth = 1.2,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq(min(months),max(months),by=1),
    limits = c(min(months),max(months))
  ) +
  scale_y_continuous(
    breaks = seq(-1.5,0,by=0.5),
    limits = c(-1.5,0)
  ) +
  labs(
    x = "Child age (months)",
    y = "Mean HAZ\n(in absence of infection)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 11),  # increase x-axis font size
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 11)
  )

ggsave( here::here("results/figures/supp_growth_trajectory.png"), growth_trajectory_fig, width = 8, height = 4)


## 3. Incidence
inc_df <- data.frame(setting_name = c("General recruitment", "Targeted recruitment", "High early incidnece"),
                     age_6_12_mad =c(params_list[["default"]]$incidence_S_mean_X_0_6, params_list[["optimistic_6mo"]]$incidence_S_mean_X_0_6, params_list[["peru_6mo"]]$incidence_S_mean_X_0_6),
                     age_6_12_msd =c(params_list[["default"]]$incidence_S_sev_mean_X_0_6, params_list[["optimistic_6mo"]]$incidence_S_sev_mean_X_0_6, params_list[["peru_6mo"]]$incidence_S_sev_mean_X_0_6),
                     age_12_18_mad = c(params_list[["default"]]$incidence_S_mean_X_6_12, params_list[["optimistic_6mo"]]$incidence_S_mean_X_6_12, params_list[["peru_6mo"]]$incidence_S_mean_X_6_12),
                     age_12_18_msd = c(params_list[["default"]]$incidence_S_sev_mean_X_6_12, params_list[["optimistic_6mo"]]$incidence_S_sev_mean_X_6_12, params_list[["peru_6mo"]]$incidence_S_sev_mean_X_6_12),
                     age_18_24_mad = c(params_list[["base_12mo"]]$incidence_S_mean_X_6_12, params_list[["optimistic_12mo"]]$incidence_S_mean_X_6_12, params_list[["peru_12mo"]]$incidence_S_mean_X_6_12),
                     age_18_24_msd = c(params_list[["base_12mo"]]$incidence_S_sev_mean_X_6_12, params_list[["optimistic_12mo"]]$incidence_S_sev_mean_X_6_12, params_list[["peru_12mo"]]$incidence_S_sev_mean_X_6_12))

inc_df <- inc_df %>%
  mutate(across(-setting_name, ~ round(round(.x, 3) * 2 * 100, 1))) # multiply all inc *2*100 and round to 1 decimal 
                     

# Kable latex table where column for setting, spanning headers ages 6-12 months, 12-18 months, 18-24 months with columns for MAD and MSD under each of those. 
# in caption incidence per 100 child years
inc_df %>%
  rename(
    `Setting` = setting_name,
    `MAD` = age_6_12_mad,
    `MSD` = age_6_12_msd,
    `MAD_2` = age_12_18_mad,
    `MSD_2` = age_12_18_msd,
    `MAD_3` = age_18_24_mad,
    `MSD_3` = age_18_24_msd
  ) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    caption = "Incidence per 100 child-years"
  ) %>%
  add_header_above(c(
    " " = 1,
    "6--12 months" = 2,
    "12--18 months" = 2,
    "18--24 months" = 2
  )) %>%
  add_header_above(c(
    " " = 1,
    " " = 2,
    " " = 2,
    " " = 2
  )) %>%
  column_spec(1, bold = TRUE)

# 4. IRR one unit increase

X_grid <- c(-1, -2)

all_IRR_df <- data.frame()

for(i in 1:length(params_list)){
  params <- params_list[[i]]

  # Any shigella diarrhea
  inc_df_0_6 <- expand.grid(
    intercept = params$hazard_S__X_int_0_6,
    haz_coef  = params$hazard_S__X_coef_0_6,
    mean_X    = X_grid
  )

  inc_0_6_neg1 <- cum_inc_by_row(inc_df_0_6[1,])
  inc_0_6_neg2 <- cum_inc_by_row(inc_df_0_6[2,])
  inc_ratio_0_6 <- inc_0_6_neg1 / inc_0_6_neg2

  inc_df_6_12 <- expand.grid(
    intercept = params$hazard_S__X_int_6_12,
    haz_coef  = params$hazard_S__X_coef_6_12,
    mean_X    = X_grid
  )

  inc_ratio_6_12 <- cum_inc_by_row(inc_df_6_12[1,]) / cum_inc_by_row(inc_df_6_12[2,])
  
  inc_6_12_neg1 <- cum_inc_by_row(inc_df_6_12[1,])
  inc_6_12_neg2 <- cum_inc_by_row(inc_df_6_12[2,])
  inc_ratio_6_12 <- inc_6_12_neg1 / inc_6_12_neg2
  
  # severe MSD
  
  # this is conditional on infection/ not from a cox model/ need to get it differently
  # and should probably be checked before it goes in the supp
  
  sev_inc_0_6_neg1 <- inc_0_6_neg1 * plogis(params$hazard_S_sev__X_int_0_6 + params$hazard_S_sev__X_coef_0_6 * X_grid[1])
  sev_inc_0_6_neg2 <- inc_0_6_neg2 * plogis(params$hazard_S_sev__X_int_0_6 + params$hazard_S_sev__X_coef_0_6 * X_grid[2])
  sev_inc_ratio_0_6 <- sev_inc_0_6_neg1 / sev_inc_0_6_neg2
  
  sev_inc_6_12_neg1 <- inc_6_12_neg1 * plogis(params$hazard_S_sev__X_int_6_12 + params$hazard_S_sev__X_coef_6_12 * X_grid[1])
  sev_inc_6_12_neg2 <- inc_6_12_neg2 * plogis(params$hazard_S_sev__X_int_6_12 + params$hazard_S_sev__X_coef_6_12 * X_grid[2])
  sev_inc_ratio_6_12 <- sev_inc_6_12_neg1 / sev_inc_6_12_neg2
  
   if(params$dose_schedule == "6mo"){
    ages <- c("6-12","12-18")
  } else{
    ages <- c("12-18", "18-24")
  }

  df <- data.frame(
    inc_neg1 = round(round(c(inc_0_6_neg1, inc_6_12_neg1),3) * 100 * 2, 1),
    inc_neg2 = round(round(c(inc_0_6_neg2, inc_6_12_neg2),3) * 100 * 2, 1),
    inc_ratio = round(c(inc_ratio_0_6, inc_ratio_6_12), 2),
    sev_inc_neg1 = round(round(c(sev_inc_0_6_neg1, sev_inc_6_12_neg1),3) * 100 * 2, 1),
    sev_inc_neg2 = round(round(c(sev_inc_0_6_neg2, sev_inc_6_12_neg2),3) * 100 * 2, 1),
    sev_inc_ratio = round(c(sev_inc_ratio_0_6, sev_inc_ratio_6_12), 2),
    age_range = ages,
    setting = rep(names(params_list)[i], 2)
  )


  all_IRR_df <- rbind(all_IRR_df, df)

}

all_IRR_df <- all_IRR_df %>%
  mutate(
    recruitment = case_when(
      setting == "default" | setting == "base_12mo" ~ "General recruitment",
      setting == "optimistic_6mo" | setting == "optimistic_12mo" ~ "Targeted recruitment",
      setting == "peru_6mo" | setting == "peru_12mo" ~ "High early incidence"
    ),
    recruitment = factor(recruitment, levels = c("General recruitment", "Targeted recruitment", "High early incidence")),
    age_range = ifelse(setting %in% c("base_12mo", "optimistic_12mo", "peru_12mo") & age_range == "12-18", "12-18 (12mo)", age_range),
    age_range = factor(age_range, levels = c("6-12", "12-18", "12-18 (12mo)", "18-24"), labels = c("6-12 months\n(6mo)", "12-18 months\n(6mo)", "12-18 months\n(12mo)", "18-24 months\n(12mo)"))
  )

# make a barchart grouped by recruitment strategy that has the inc_ratio and sev_inc_ratio for each age range
plot_df <- all_IRR_df %>%
  select(recruitment, age_range, inc_ratio, sev_inc_ratio) %>%
  pivot_longer(
    cols = c(inc_ratio, sev_inc_ratio),
    names_to = "outcome",
    values_to = "irr"
  ) %>%
  mutate(
    outcome = recode(outcome,
                     inc_ratio = "Any diarrhea",
                     sev_inc_ratio = "Moderate-to-severe diarrhea")
  )

irr_plot <- ggplot(
  plot_df,
  aes(
    x = age_range,
    y = irr,
    group = outcome
  )
) +
  geom_col(
    aes(
      fill = outcome
    ),
    position = position_dodge(width = 0.7),
    width = 0.65,
    alpha = 0.9,
    linewidth = 0.6
  ) +
  
  # annotate IRRs above bars
  geom_text(
    aes(
      label = irr
    ),
    position = position_dodge(width = 0.7),
    vjust = -0.3,       # slightly above the bar
    size = 3.5
  ) +
  
  facet_wrap(~ recruitment) +
  
  # match your main color style
  scale_fill_manual(
    values = c(
      "Any diarrhea" = "#00468BFF",
      "Moderate-to-severe diarrhea" = "#ED0000FF"
    )
  ) +
  
  geom_hline(
    yintercept = 1,
    linetype = "dotted",
    color = "gray20"
  ) +
  
  labs(
    x = "Age Group (Immunization schedule)",
    y = "Incidence Rate Ratio\n(HAZ -1 vs HAZ -2)",
    fill = "Severity"
  ) +
  
  coord_cartesian(ylim = c(0.75, 1.25)) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 11),
    legend.position = "bottom"
  )

ggsave(here::here("results/figures/irr_plot_supp.png"), irr_plot, width = 13, height = 6)

