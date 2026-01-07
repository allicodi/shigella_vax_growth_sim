# ------------------------------------------------------
# Script to make simulation design figure
# ------------------------------------------------------

here::i_am("misc/simulation_design_figure.R")

library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpattern)

source(here::here("R/parameter_generation_fns.R"))

# Panel 1: Incidence ------------------------------------------------------

inc_6mo <- get_incidence(dose_schedule = "6mo", enroll_site = "Total")
inc_12mo <- get_incidence(dose_schedule = "12mo", enroll_site = "Total")

# inc_6mo_Peru <- get_incidence(dose_schedule = "6mo", enroll_site = "Peru")
# inc_12mo_Peru <- get_incidence(dose_schedule = "12mo", enroll_site = "Peru")

inc_6mo_Gambia <- get_incidence(dose_schedule = "6mo", enroll_site = "The Gambia")
inc_12mo_Gambia <- get_incidence(dose_schedule = "12mo", enroll_site = "The Gambia")

plot_df <- data.frame(
  age_grp = c(rep(c("6-12 months", "12-18 months", "18-24 months"),2)),
  site = c(rep("Overall",3), rep("The Gambia",3)),
  pt_est_mad = c(inc_6mo$mad_inc_0_6, inc_6mo$mad_inc_6_12, inc_12mo$mad_inc_6_12,
                 inc_6mo_Gambia$mad_inc_0_6, inc_6mo_Gambia$mad_inc_6_12, inc_12mo_Gambia$mad_inc_6_12),
  lower_ci_mad = c(inc_6mo$mad_inc_0_6_lower, inc_6mo$mad_inc_6_12_lower, inc_12mo$mad_inc_6_12_lower,
                   inc_6mo_Gambia$mad_inc_0_6_lower, inc_6mo_Gambia$mad_inc_6_12_lower, inc_12mo_Gambia$mad_inc_6_12_lower),
  upper_ci_mad = c(inc_6mo$mad_inc_0_6_upper, inc_6mo$mad_inc_6_12_upper, inc_12mo$mad_inc_6_12_upper,
                   inc_6mo_Gambia$mad_inc_0_6_upper, inc_6mo_Gambia$mad_inc_6_12_upper, inc_12mo_Gambia$mad_inc_6_12_upper),
  pt_est_msd = c(inc_6mo$msd_inc_0_6, inc_6mo$msd_inc_6_12, inc_12mo$msd_inc_6_12,
                 inc_6mo$msd_inc_0_6_Gambia, inc_6mo$msd_inc_6_12_Gambia, inc_12mo$msd_inc_6_12_Gambia),
  lower_ci_msd = c(inc_6mo$msd_inc_0_6_lower, inc_6mo$msd_inc_6_12_lower, inc_12mo$msd_inc_6_12_lower,
                   inc_6mo_Gambia$msd_inc_0_6_lower, inc_6mo_Gambia$msd_inc_6_12_lower, inc_12mo_Gambia$msd_inc_6_12_lower),
  upper_ci_msd = c(inc_6mo$msd_inc_0_6_upper, inc_6mo$msd_inc_6_12_upper, inc_12mo$msd_inc_6_12_upper,
                   inc_6mo_Gambia$msd_inc_0_6_upper, inc_6mo_Gambia$msd_inc_6_12_upper, inc_12mo_Gambia$msd_inc_6_12_upper)
)

# Ensure age groups are ordered correctly
plot_df <- plot_df %>%
  mutate(
    age_grp = factor(
      age_grp,
      levels = c("6-12 months", "12-18 months", "18-24 months")
    )
  )

plot_long <- plot_df %>%
  mutate(
    age_grp = factor(
      age_grp,
      levels = c("6-12 months", "12-18 months", "18-24 months")
    ),
    site = factor(site, levels = c("Overall", "The Gambia"))
  ) %>%
  pivot_longer(
    cols = starts_with(c("pt_est_", "lower_ci_", "upper_ci_")),
    names_to = c(".value", "severity"),
    names_pattern = "(pt_est|lower_ci|upper_ci)_(mad|msd)"
  ) %>%
  mutate(
    severity = recode(severity,
                      mad = "MAD",
                      msd = "MSD")
  )

inc_fig <- ggplot(
  plot_long,
  aes(
    x = age_grp,
    y = pt_est*2*100,
    color = site,
    linetype = site,
    shape = severity,
    group = interaction(site, severity)
  )
) +
  geom_point(
    position = position_dodge(width = 0.35),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = lower_ci*2*100, ymax = upper_ci*2*100),
    width = 0.15,
    position = position_dodge(width = 0.35)
  ) +
  scale_color_manual(
    values = c(
      "Overall" = "#ED0000FF",
      "The Gambia" = "#00468BFF"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Overall" = "solid",
      "The Gambia" = "solid"
    )
  ) +
  scale_shape_manual(
    values = c(
      "MAD" = 16,  # solid circle
      "MSD" = 17   # triangle
    )
  ) +
  labs(
    x = "Age Group",
    y = "Incidence \n(Shigella episodes\n per 100 child years)",
    color = "Site",
    linetype = "Site",
    shape = "Severity"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

## TODO: Get LSD CIs from Maria

# Panel 2: Baseline HAZ ------------------------------------------------------

# Histogram of baseline HAZ overall & The Gambia, impose Normal curves over top
# there will be a version subset to 6-12

# Histogram of baseline HAZ with age-specific histograms
# Normal curves overlaid for Overall vs The Gambia

efgh_data <- readRDS(here::here("data/efgh/efgh_data.Rds")) %>%
  mutate(
    age_grp = case_when(
      enr_age_months >= 6 & enr_age_months <= 8  ~ "6mo",
      enr_age_months >= 9 & enr_age_months <= 11 ~ "12mo",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(age_grp))

# Histogram fill colors by age only
hist_cols <- c(
  "6mo"  = "grey90",
  "12mo" = "gray65"
)

# Parameter data (site-specific normals)
param_df <- bind_rows(
  data.frame(
    age_grp = "6mo",
    site = "Overall",
    get_baseline_growth(dose_schedule = "6mo", enroll_site = "Overall")
  ),
  data.frame(
    age_grp = "6mo",
    site = "The Gambia",
    get_baseline_growth(dose_schedule = "6mo", enroll_site = "The Gambia")
  ),
  data.frame(
    age_grp = "12mo",
    site = "Overall",
    get_baseline_growth(dose_schedule = "12mo", enroll_site = "Overall")
  ),
  data.frame(
    age_grp = "12mo",
    site = "The Gambia",
    get_baseline_growth(dose_schedule = "12mo", enroll_site = "The Gambia")
  )
)

param_df <- param_df %>%
  mutate(age_grp = factor(age_grp, levels = c("6mo", "12mo")))

# Grid for normal curves
x_grid <- seq(-5, 3, length.out = 500)

norm_df <- param_df %>%
  group_by(age_grp, site) %>%
  do({
    data.frame(
      enr_haz = x_grid,
      density = dnorm(
        x_grid,
        mean = .$mean_enr_haz,
        sd   = .$sd_enr_haz
      )
    )
  }) %>%
  ungroup()

# Corner annotation data
annot_df <- param_df %>%
  arrange(site, age_grp) %>%
  mutate(
    label = paste0(
      site, " (", age_grp, ")\n",
      "\u03BC = ", round(mean_enr_haz, 2),
      ", \u03C3 = ", round(sd_enr_haz, 2)
    ),
    y_pos = c(0.42, 0.32, 0.22, 0.12)
  )

bl_haz_fig <- ggplot() +
  # Age-specific histograms (ALL sites combined)
  # geom_histogram(
  #   data = efgh_data,
  #   aes(
  #     x = enr_haz,
  #     y = after_stat(density),
  #     fill = age_grp
  #   ),
  #   bins = 30,
  #   color = "white",
  #   alpha = 0.6,
  #   position = "identity"
  # ) +
  
  # Normal curves (site-specific)
  geom_line(
    data = norm_df,
    aes(
      x = enr_haz,
      y = density,
      color = site,
      linetype = age_grp,
      group = interaction(site, age_grp)
    ),
    linewidth = 1.2
  ) +
  
  # Mean lines
  geom_vline(
    data = param_df,
    aes(
      xintercept = mean_enr_haz,
      color = site,
      linetype = age_grp
    ),
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  
  # Corner annotations
  geom_text(
    data = annot_df,
    aes(
      x = 2.7,
      y = y_pos,
      label = label,
      color = site
    ),
    hjust = 1,
    vjust = 1,
    size = 2,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(
    values = hist_cols,
    name = "Immunization\nSchedule"
  ) +
  scale_color_manual(
    values = c(
      "Overall" = "#ED0000FF",
      "The Gambia" = "#00468BFF"
    ),
    name = "Site"
  ) +
  scale_linetype_manual(
    values = c(
      "6mo" = "solid",
      "12mo" = "dashed"
    ),
    name = "Immunization\nSchedule",
    guide = guide_legend(
      override.aes = list(
        linewidth = 0.5   # thinner in legend
      )
    )
  ) +
  labs(
    x = "Baseline HAZ",
    y = "Density"
  ) +
  coord_cartesian(xlim = c(-4.5, 3)) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# Panel 3: Effect of BL growth on Shigella infection -----------------------------

parameters_6mo <- readRDS(here::here("parameters/parameters_default.Rds"))
parameters_12mo <- readRDS(here::here("parameters/parameters_base_12mo.Rds"))

parameters_gambia_6mo <- readRDS(here::here("parameters/parameters_gambia_default_VE_6mo.Rds"))
parameters_gambia_12mo <- readRDS(here::here("parameters/parameters_gambia_default_VE_12mo.Rds"))

# Put all parameter sets in a list with names
param_list <- list(
  "Overall (6mo)" = parameters_6mo,
  "Overall (12mo)" = parameters_12mo,
  "The Gambia (6mo)" = parameters_gambia_6mo,
  "The Gambia (12mo)" = parameters_gambia_12mo
)

X_grid <- c(-2, 0)

df_all <- map2_dfr(param_list, names(param_list), function(parameters, param_name) {
  
  inc_df_0_6 <- expand.grid(intercept = parameters$hazard_S__X_int_0_6,
                            haz_coef = parameters$hazard_S__X_coef_0_6,
                            mean_X = X_grid)
  
  inc_ratio_0_6 <- cum_inc_by_row(inc_df_0_6[1,]) / cum_inc_by_row(inc_df_0_6[2,])
  
  inc_df_6_12 <- expand.grid(intercept = parameters$hazard_S__X_int_6_12,
                            haz_coef = parameters$hazard_S__X_coef_6_12,
                            mean_X = X_grid)
  
  inc_ratio_6_12 <- cum_inc_by_row(inc_df_6_12[1,]) / cum_inc_by_row(inc_df_6_12[2,])
  
  df <- data.frame(inc_ratio = c(inc_ratio_0_6, inc_ratio_6_12), 
                   period = c("0–6 months", "6–12 months"),
                   scenario = param_name)
  
  # Recode period into actual age ranges
  df <- df %>%
    mutate(age_range = case_when(
      scenario == "Overall (6mo)" & period == "0–6 months"  ~ "6–12 months",
      scenario == "Overall (6mo)" & period == "6–12 months" ~ "12–18 months",
      scenario == "Overall (12mo)" & period == "0–6 months"  ~ "12–18 months",
      scenario == "Overall (12mo)" & period == "6–12 months" ~ "18–24 months",
      scenario == "The Gambia (6mo)" & period == "0–6 months"    ~ "6–12 months",
      scenario == "The Gambia (6mo)" & period == "6–12 months"   ~ "12–18 months",
      scenario == "The Gambia (12mo)" & period == "0–6 months"   ~ "12–18 months",
      scenario == "The Gambia (12mo)" & period == "6–12 months"  ~ "18–24 months"
    ),
    setting = if_else(scenario == "The Gambia (6mo)" | scenario == "The Gambia (12mo)", "The Gambia", "Overall"))
  
  df
})

df_all <- df_all %>%
  mutate(
    age_range = factor(
      age_range,
      levels = c("6–12 months", "12–18 months", "18–24 months")
    ),
    setting = factor(setting, levels = c("Overall", "The Gambia")),
    scenario = factor(scenario, levels = c("Overall (6mo)",
                                           "The Gambia (6mo)",
                                           "Overall (12mo)",
                                           "The Gambia (12mo)"))
  )

hazard_fig <- ggplot(
  df_all,
  aes(
    x = scenario,
    y = inc_ratio,
    fill = setting,          # color = site
    pattern = age_range      # texture = age group
  )
) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge2(width = 0.8, preserve = "single"),
    width = 0.7,
    color = "gray",
    
    # pattern styling
    pattern_fill = "gray80",
    pattern_colour = "gray90",
    pattern_density = 0.1,
    pattern_spacing = 0.1
  ) +
  scale_fill_manual(
    values = c(
      "Overall" = "#ED0000FF",
      "The Gambia" = "#00468BFF"
    )
  ) +
  scale_pattern_manual(
    values = c(
      "6–12 months"  = "stripe",
      "12–18 months" = "crosshatch",
      "18–24 months" = "circle"
    ), 
    guide = guide_legend(
      override.aes = list(
        pattern_spacing = 0.01,
        pattern_density = 0.05   # thinner in legend
      )
    )
  ) +
  labs(
    x = NULL,
    y = "Incidence rate ratio\n(HAZ −2 vs 0)",
    fill = "Site",
    pattern = "Age group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

hazard_fig


### OLD ###
# X_grid <- seq(-2,1,by=0.01)
# 
# # Put all parameter sets in a list with names
# param_list <- list(
#   "Overall (6mo)" = parameters_6mo,
#   "Overall (12mo)" = parameters_12mo,
#   "The Gambia (6mo)" = parameters_gambia_6mo,
#   "The Gambia (12mo)" = parameters_gambia_12mo
# )
# 
# df_all <- map2_dfr(param_list, names(param_list), function(parameters, param_name) {
#   
#   haz_0_6 <- hazard(
#     intercept = parameters$hazard_S__X_int_0_6,
#     haz_coef = parameters$hazard_S__X_coef_0_6,
#     haz = X_grid
#   )
#   
#   haz_6_12 <- hazard(
#     intercept = parameters$hazard_S__X_int_6_12,
#     haz_coef = parameters$hazard_S__X_coef_6_12,
#     haz = X_grid
#   )
#   
#   df <- tibble(
#     X = rep(X_grid, 2),
#     hazard = c(haz_0_6, haz_6_12),
#     period = rep(c("0–6 months", "6–12 months"), each = length(X_grid)),
#     scenario = param_name
#   )
#   
#   # Recode period into actual age ranges
#   df <- df %>%
#     mutate(age_range = case_when(
#       scenario == "Overall (6mo)" & period == "0–6 months"  ~ "6–12 months",
#       scenario == "Overall (6mo)" & period == "6–12 months" ~ "12–18 months",
#       scenario == "Overall (12mo)" & period == "0–6 months"  ~ "12–18 months",
#       scenario == "Overall (12mo)" & period == "6–12 months" ~ "18–24 months",
#       scenario == "The Gambia (6mo)" & period == "0–6 months"    ~ "6–12 months",
#       scenario == "The Gambia (6mo)" & period == "6–12 months"   ~ "12–18 months",
#       scenario == "The Gambia (12mo)" & period == "0–6 months"   ~ "12–18 months",
#       scenario == "The Gambia (12mo)" & period == "6–12 months"  ~ "18–24 months"
#     ),
#     setting = if_else(scenario == "The Gambia (6mo)" | scenario == "The Gambia (12mo)", "The Gambia", "Overall"))
#   
#   df
# })
# 
# # Plot using age_range as color
# hazard_fig <- ggplot(df_all, aes(X, hazard, color = factor(age_range, levels = c("6–12 months", "12–18 months","18–24 months")), linetype = scenario)) +
#   geom_line(size = 1) +
#   labs(
#     x = "Baseline HAZ",
#     y = "Weekly infection hazard",
#     color = "Age range",
#     linetype = "Trial"
#   ) +
#   theme_minimal()


# Panel 4: Growth trajectory in absence of Shigella ---------------------------

monthly_growth_model_6mo <- get_monthly_growth("6mo")
monthly_growth_model_12mo <- get_monthly_growth("12mo")

# Plot all at X = -1
df <- data.frame(
  trial = factor(c("6mo","12mo"), levels = c("6mo","12mo")),
  X = rep(-1, 2)
)

# df <- data.frame(
#   trial = c("6mo","6mo","12mo","12mo"),
#   site = c("Overall","The Gambia","Overall","The Gambia"),
#   X = c(get_baseline_growth(dose_schedule = "6mo", enroll_site = "Overall")$mean_enr_haz,
#         get_baseline_growth(dose_schedule = "6mo", enroll_site = "The Gambia")$mean_enr_haz,
#         get_baseline_growth(dose_schedule = "12mo", enroll_site = "Overall")$mean_enr_haz,
#         get_baseline_growth(dose_schedule = "12mo", enroll_site = "The Gambia")$mean_enr_haz)
#   )

# Combine all results into a single long dataframe
res_list <- vector("list", nrow(df))

for(j in 1:nrow(df)) {
  X <- df$X[j]
  
  if(df$trial[j] == "6mo") {
    months <- 6:18
    beta_0 <- monthly_growth_model_6mo$beta_0[months]
    beta_1 <- monthly_growth_model_6mo$beta_1[months]
  } else {
    months <- 12:24
    beta_0 <- monthly_growth_model_12mo$beta_0[months]
    beta_1 <- monthly_growth_model_12mo$beta_1[months]
  }
  
  Y_vec <- numeric(length(months))
  Y_vec[1] <- X  # baseline
  if(length(months) > 1){
    for(i in 2:length(months)) {
      Y_vec[i] <- beta_0[i] + beta_1[i] * Y_vec[i-1]
    }
  }
  
  res_list[[j]] <- data.frame(
    trial = df$trial[j],
    #site = df$site[j],
    month = months,
    Y = Y_vec
  )
}

plot_df <- bind_rows(res_list)

# Determine starting points
start_points <- plot_df %>%
  group_by(trial) %>% #, site) %>%
  filter(month == if_else(trial == "6mo", 6L, 12L)) %>%
  ungroup()

# Plot
growth_trajectory_fig <- ggplot(plot_df, aes(x = month, y = Y, 
                                             #color = site, 
                                             linetype = trial)) + #, group = interaction(site, trial))) +
  geom_line(linewidth = 1.2) +
  
  # Bold dot at starting HAZ
  # geom_point(
  #   data = start_points,
  #   aes(x = month, y = Y),
  #   inherit.aes = FALSE,
  #   size = 4,
  #   shape = 16
  # ) +
  
  # Label at starting HAZ
  # geom_text(
  #   data = start_points,
  #   aes(x = month, y = Y, label = round(Y,2)),
  #   vjust = 0.1,hjust = -0.25,
  #   fontface = "bold",
  #   size = 3,
  #   inherit.aes = FALSE
  # ) +
  
  # scale_color_manual(
  #   values = c("Overall" = "#ED0000FF", "The Gambia" = "#00468BFF")
  # ) +
  scale_linetype_manual(
    values = c(
      "6mo" = "solid",
      "12mo" = "dashed"
    ),
    name = "Immunization\nSchedule",
    guide = guide_legend(
      override.aes = list(
        linewidth = 0.5   # thinner in legend
      )
    )
  )+ 
  scale_x_continuous(breaks = seq(6,24,by=1), limits = c(6,24)) +
  
  labs(
    x = "Child age (months)",
    y = "Mean HAZ\n(in absence of infection)",
    #color = "Site",
    linetype = "Immunization\nSchedule"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

# ^^ 
# Plotted at avg HAZ in each setting 
# but that ends up looking weird so maybe just pick something
# and go with same for all settings?

# Panel 5: Shigella growth effect 

# Make figure of the year version for manuscript
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

growth_effect_fig <- ggplot(plot_combined, aes(x = month_num, y = pt_est, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
  geom_line(data = filter(plot_combined, line_type == "12-18 months"), 
            aes(size = 1.5), linetype = "solid") +
  geom_point(data = filter(plot_combined, line_type == "12-18 months"), size = 4) +
  scale_size_identity() +
  scale_color_manual(values = c(
    "Moderate-to-Severe Shigella" = "#925E9FFF",
    "Less-severe Shigella" = "#42B540FF"
  )) +
  scale_fill_manual(values = c(
    "Moderate-to-Severe Shigella" = "#925E9FFF",
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

# Combine all 5 panels using patchwork ----------------------------------------
# Layout: 2 rows of 2 columns, then 1 full-width row at bottom
combined_figure <- (inc_fig + bl_haz_fig) / 
  (hazard_fig + growth_trajectory_fig) / 
  growth_effect_fig +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(plot.tag = element_text(face = 'bold', size = 16))
  )

# Display the combined figure
print(combined_figure)
