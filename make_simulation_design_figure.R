# ------------------------------------------------------
# Script to make simulation design figure with consolidated legends
# ------------------------------------------------------

here::i_am("make_simulation_design_figure.R")

library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpattern)

source(here::here("R/parameter_generation_fns.R"))

cargs <- commandArgs(TRUE)
setting <- cargs[[1]]

params <- readRDS(here::here(paste0("parameters/parameters_",setting,".Rds")))

cfg <- yaml::read_yaml("config.yml")
config <- cfg[[setting]]

severity_colors <- c(
  "MAD" = "#00468BFF",
  "MSD" = "#ED0000FF",
  "LSD" = "#42B540FF"
)

age_linetypes <- c(
  "6-12 months" = "solid",
  "12-18 months" = "longdash",
  "18-24 months" = "twodash"
)

age_patterns <- c(
  "6-12 months" = "none",
  "12-18 months" = "stripe",
  "18-24 months" = "twodash"
)

# Panel 1: Incidence ------------------------------------------------------

inc_res <- get_incidence(dose_schedule = params$dose_schedule, enroll_site = params$site)

if(params$dose_schedule == "6mo"){
  age_grp <- c("6-12 months", "12-18 months")
} else{
  age_grp <-c("12-18 months", "18-24 months")
}

plot_df <- data.frame(
  age_grp = age_grp,
  site = c(rep(params$site,2)),
  pt_est_mad = c(params$incidence_S_mean_X_0_6, params$incidence_S_mean_X_6_12),#, inc_12mo_Gambia$mad_inc_6_12),
  lower_ci_mad = c(inc_res$mad_inc_0_6_lower, inc_res$mad_inc_6_12_lower), #, inc_12mo_Gambia$mad_inc_6_12_lower),
  upper_ci_mad = c(inc_res$mad_inc_0_6_upper, inc_res$mad_inc_6_12_upper), #, inc_12mo_Gambia$mad_inc_6_12_upper),
  pt_est_msd = c(params$incidence_S_sev_mean_X_0_6, params$incidence_S_sev_mean_X_6_12),# inc_12mo_Gambia$msd_inc_6_12),
  lower_ci_msd = c(inc_res$msd_inc_0_6_lower, inc_res$msd_inc_6_12_lower),# inc_12mo_Gambia$msd_inc_6_12_lower),
  upper_ci_msd = c(inc_res$msd_inc_0_6_upper, inc_res$msd_inc_6_12_upper)# inc_12mo_Gambia$msd_inc_6_12_upper)
)

plot_long <- plot_df %>%
  mutate(
    age_grp = factor(
      age_grp,
      levels = age_grp# , "18-24 months")
    )
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
  ) %>%
  mutate(
    severity = factor(
      severity, levels = c("LSD", "MAD", "MSD")
    )
  )

inc_fig <- ggplot(
  plot_long,
  aes(
    x = age_grp,
    y = pt_est*2*100,
    color = severity,
    linetype = age_grp
  )
) +
  geom_point(
    position = position_dodge(width = 0.35),
    size = 3.5
  ) +
  geom_errorbar(
    aes(ymin = lower_ci*2*100, ymax = upper_ci*2*100),
    width = 0.2,
    position = position_dodge(width = 0.35)
  ) +
  scale_color_manual(
    values = severity_colors,
    labels = c(
      "LSD" = "Less-severe diarrhea (LSD)",
      "MAD" = "Any diarrhea",
      "MSD" = "Moderate-to-severe diarrhea (MSD)"
    ),
    name = "Episode severity"
  ) +
  scale_linetype_manual(
    values = age_linetypes
  ) +
  guides(
    linetype = "none"
  ) +
  labs(
    x = "Age Group",
    y = "Incidence \n(Shigella episodes\nper 100 child years)"
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

# Panel 2: Baseline HAZ ------------------------------------------------------

x_grid <- seq(-5, 3, length.out = 500)

norm_df <- data.frame(
  enr_haz = x_grid,
  density = dnorm(
    x_grid,
    mean = params$mean_X,
    sd   = params$sd_X
  )
)

annot_df <- data.frame(
  label = paste0(
    "\u03BC = ", round(params$mean_X, 2),
    ", \u03C3 = ", round(params$sd_X, 2)
  )
)

bl_haz_fig <- ggplot(norm_df, aes(x = enr_haz, y = density)) +
  geom_line(linewidth = 1.2, color = "black") +
  geom_vline(
    xintercept = params$mean_X,
    linewidth = 0.7, 
    linetype = 'dotted'
  ) +
  geom_text(
    data = annot_df,
    aes(
      x = 2.7,
      y = 0.42,
      label = label
    ),
    hjust = 1,
    vjust = 1,
    size = 5
  ) +
  labs(
    x = "Baseline HAZ",
    y = "Density"
  ) +
  coord_cartesian(xlim = c(-4.5, 3)) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 11),  # increase x-axis font size
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 11)
  )

# Panel 3: Effect of BL growth on Shigella infection -----------------------------

X_grid <- c(-1, -2)

inc_df_0_6 <- expand.grid(
  intercept = params$hazard_S__X_int_0_6,
  haz_coef  = params$hazard_S__X_coef_0_6,
  mean_X    = X_grid
)

# Any diarrhea
inc_0_6_neg1 <- cum_inc_by_row(inc_df_0_6[1,])
inc_0_6_neg2 <- cum_inc_by_row(inc_df_0_6[2,])
inc_ratio_0_6 <- inc_0_6_neg1 / inc_0_6_neg2

inc_df_6_12 <- expand.grid(
  intercept = params$hazard_S__X_int_6_12,
  haz_coef  = params$hazard_S__X_coef_6_12,
  mean_X    = X_grid
)

inc_6_12_neg1 <- cum_inc_by_row(inc_df_6_12[1,])
inc_6_12_neg2 <- cum_inc_by_row(inc_df_6_12[2,])
inc_ratio_6_12 <- inc_6_12_neg1 / inc_6_12_neg2

# MSD

sev_inc_0_6_neg1 <- inc_0_6_neg1 * plogis(params$hazard_S_sev__X_int_0_6 + params$hazard_S_sev__X_coef_0_6 * X_grid[1])
sev_inc_0_6_neg2 <- inc_0_6_neg2 * plogis(params$hazard_S_sev__X_int_0_6 + params$hazard_S_sev__X_coef_0_6 * X_grid[2])
sev_inc_ratio_0_6 <- sev_inc_0_6_neg1 / sev_inc_0_6_neg2

sev_inc_6_12_neg1 <- inc_6_12_neg1 * plogis(params$hazard_S_sev__X_int_6_12 + params$hazard_S_sev__X_coef_6_12 * X_grid[1])
sev_inc_6_12_neg2 <- inc_6_12_neg2 * plogis(params$hazard_S_sev__X_int_6_12 + params$hazard_S_sev__X_coef_6_12 * X_grid[2])
sev_inc_ratio_6_12 <- sev_inc_6_12_neg1 / sev_inc_6_12_neg2


df <- data.frame(
  inc_neg1 = round(round(c(inc_0_6_neg1, inc_6_12_neg1),3) * 100 * 2, 1),
  inc_neg2 = round(round(c(inc_0_6_neg2, inc_6_12_neg2),3) * 100 * 2, 1),
  inc_ratio = round(c(inc_ratio_0_6, inc_ratio_6_12), 2),
  sev_inc_neg1 = round(round(c(sev_inc_0_6_neg1, sev_inc_6_12_neg1),3) * 100 * 2, 1),
  sev_inc_neg2 = round(round(c(sev_inc_0_6_neg2, sev_inc_6_12_neg2),3) * 100 * 2, 1),
  sev_inc_ratio = round(c(sev_inc_ratio_0_6, sev_inc_ratio_6_12), 2),
  age_range = age_grp
)

df$age_range <- factor(
  age_grp,
  levels = age_grp
)

plot_df <- df %>%
  select(age_range, inc_ratio, sev_inc_ratio) %>%
  pivot_longer(
    cols = c(inc_ratio, sev_inc_ratio),
    names_to = "severity",
    values_to = "irr"
  ) %>%
  mutate(
    severity = recode(severity,
                     inc_ratio = "Any diarrhea",
                     sev_inc_ratio = "Moderate-to-severe diarrhea")
  )

hazard_fig <- ggplot(
  plot_df,
  aes(
    x = age_range,
    y = irr,
    linetype = age_range,
    group = severity
  )
) +
  geom_col(
    aes(fill = severity),
    color = 'black',
    position = position_dodge(width = 0.7),
    width = 0.65,
    alpha = 0.9,
    linewidth = 0.6
  ) +
  
  geom_text(
    aes(
      label = irr
    ),
    position = position_dodge(width = 0.7),
    vjust = -0.3,       # slightly above the bar
    size = 3.5
  ) +
  
  scale_fill_manual(
    values = c(
      "Any diarrhea" = "#00468BFF",
      "Moderate-to-severe diarrhea" = "#ED0000FF"
    )
  ) + 

  scale_linetype_manual(
    values = age_linetypes,
    guide = "none"     # no legend for linetype
  ) +
  geom_hline(
    yintercept = 1, 
    linetype = 'dotted',
    color = "gray10"
  ) +
  labs(
    x = "Age Group",
    y = "Incidence rate ratio\n(HAZ -1 vs -2)" # one unit increase HAZ -2 to -1
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(0.75, 1.25)) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 11),  # increase x-axis font size
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 11)
  )

# Panel 4: Growth trajectory ---------------------------

monthly_growth_model <- params$monthly_growth_model
X <- params$mean_X

if(params$dose_schedule == "6mo"){
  months <- 6:18
} else{
  months <- 12:24
}

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
    breaks = seq(-1.75,-0.5,by=0.25),
    limits = c(-1.75,-0.5)
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

# Panel 5: Shigella growth effect ---------------------------------------

shigella_growth_meta_analysis_results <- readRDS(here::here("misc/results/case_control/shigella_growth_effect_data.Rds"))

plot_df <- shigella_growth_meta_analysis_results %>% 
  filter(group != "Any Shigella") %>%
  mutate(severity = factor(if_else(group == "Less-severe Shigella", "LSD", "MSD"), levels = c("LSD", "MSD")))

dose_schedule <- params$dose_schedule
scale_growth_effect_0_6 <- config$scale_growth_effect_0_6
scale_growth_effect_6_12 <- config$scale_growth_effect_6_12
spline_formula <- config$effect_shigella_growth_formula

scale_0_6_df <- plot_df %>%
  mutate(
    age_strata = if_else(dose_schedule == "6mo", "6-12 months", "12-18 months"),
    scale_growth_effect_0_6 = scale_growth_effect_0_6
  ) %>%
  mutate(
    pt_est = if_else(scale_growth_effect_0_6 > 0, pt_est + scale_growth_effect_0_6*(upper_ci - pt_est), pt_est + scale_growth_effect_0_6*(pt_est - lower_ci)) 
  ) %>% 
  filter(group != "Any Shigella")

scale_6_12_df <- plot_df %>%
  mutate(
    age_strata = if_else(dose_schedule == "6mo", "12-18 months", "18-24 months"),
    scale_growth_effect_6_12 = scale_growth_effect_6_12
  ) %>%
  mutate(
    pt_est = if_else(scale_growth_effect_6_12 > 0, pt_est + scale_growth_effect_6_12*(upper_ci - pt_est), pt_est + scale_growth_effect_6_12*(pt_est - lower_ci))
  )  %>% 
  filter(group != "Any Shigella")

plot_combined <- bind_rows(
  plot_df %>% mutate(line_type = "All ages", size_line = 1.5),
  scale_0_6_df %>% mutate(line_type = age_grp[1], size_line = 0.7),
  scale_6_12_df %>% mutate(line_type = age_grp[2], size_line = 0.7)
)

plot_combined <- plot_combined %>%
  mutate(
    age_strata = factor(age_strata, levels = c("All ages", age_grp)),
    line_type = factor(line_type, levels = c("All ages", age_grp))
  )

growth_effect_fig <- ggplot(plot_combined, aes(x = month_num, y = pt_est, color = severity, fill = severity)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
  # geom_line(data = filter(plot_combined, line_type == "All ages"), 
  #           aes(size = 1.5), linetype = "solid") +
  geom_point(data = filter(plot_combined, line_type == "All ages"), size = 4) +
  scale_size_identity() +
  scale_color_manual(
    values = severity_colors,
    labels = c(
      "LSD" = "Less-severe diarrhea (LSD)",
      "MSD" = "Moderate-to-severe diarrhea (MSD)"
    ),
    name = "Episode severity"
  ) +
  scale_fill_manual(
    values = severity_colors,
    guide = "none"
  ) +
  scale_x_continuous(breaks = 1:12) +
  stat_smooth(
    data = filter(plot_combined, line_type == age_grp[1]),
    method = "glm",
    formula = spline_formula,
    aes(linetype = age_grp[1]),
    size = 1,
    se = FALSE
  ) +
  stat_smooth(
    data = filter(plot_combined, line_type == age_grp[2]),
    method = "glm",
    formula = spline_formula,
    aes(linetype = age_grp[2]),
    size = 1,
    se = FALSE
  ) +
  scale_linetype_manual(
    name = "Age Group",
    values = age_linetypes,
    guide = guide_legend(
      override.aes = list(
        linewidth = 0.5,   # thinner in legend
        color = 'black'
      )
    )
  ) +
  labs(
    x = "Month",
    y = "HAZ difference (95% CI)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 11),  # increase x-axis font size
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 11)
  )


# Combine all panels with consolidated legend ----------------------------------------

# Dummy legend
legend_df <- expand.grid(
  age_group = factor(age_grp, levels = age_grp),
  severity = factor(c("LSD", "MAD", "MSD"), levels = c("LSD", "MAD", "MSD")),
  x = 1, y = 1
)

dummy_legend <- ggplot(legend_df, aes(x = x, y = y, color = severity, linetype = age_group)) +
  geom_point(size = 4) +
  geom_line() +
  scale_color_manual(
    values = severity_colors,
    labels = c(
      "LSD" = "Less-severe diarrhea (LSD)",
      "MAD" = "Any diarrhea",
      "MSD" = "Moderate-to-severe diarrhea (MSD)"
    ),
    name = "Episode severity"
  ) +
  scale_linetype_manual(
    values = age_linetypes,
    name = "Age group"
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  )


# Create a helper function to extract legend
get_legend <- function(plot) {
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if(length(leg) > 0) {
    legend <- tmp$grobs[[leg]]
  } else {
    legend <- NULL
  }
  return(legend)
}

# Extract legends for consolidation with consistent theme
legend_theme <- theme(
  legend.title = element_text(size = 11, face = "bold"),
  legend.text = element_text(size = 10)
)

legend_dummy <- get_legend(dummy_legend) 

library(gridExtra)
library(grid)

combined_legend <- arrangeGrob(
  legend_dummy, 
  ncol = 1, 
  heights = 1
)

# Combine plots without legends
# combined_figure <- (inc_fig + bl_haz_fig) / 
#   (hazard_fig + growth_trajectory_fig) / 
#   growth_effect_fig +
#   plot_layout(heights = c(1, 1, 3)) +
#   plot_annotation(
#     tag_levels = 'A',
#     theme = theme(plot.tag = element_text(face = 'bold', size = 16))
#   )

combined_figure <- (bl_haz_fig + growth_trajectory_fig) / 
  (inc_fig + hazard_fig) / 
  growth_effect_fig +
  plot_layout(heights = c(1, 1, 3)) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(plot.tag = element_text(face = 'bold', size = 16))
  )

# Combine with legend
final_figure <- wrap_elements(combined_figure) + 
  wrap_elements(combined_legend) + 
  plot_layout(heights = c(5, 0.1))

final_figure

ggsave(here::here(paste0("results/figures/parameterization_", setting, ".png")), 
       plot = final_figure,
       width = 15, height = 10, dpi = 300)

