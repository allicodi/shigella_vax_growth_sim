
library(tidyverse)

here::i_am("misc/plot_growth_tragectories.R")

maled_data <- readRDS(here::here("data/maled/growth_wide_maled.Rds"))

# Pivot to long format
maled_long <- maled_data %>%
  filter(country_id != "PK") %>%
  pivot_longer(
    cols = starts_with("zlen_"),  # adjust if your timepoint columns have a different prefix
    names_to = "month",
    names_prefix = "zlen_",
    names_transform = list(month = as.integer),
    values_to = "zlen"
  )

# Calculate mean by site and month
mean_by_site_month <- maled_long %>%
  group_by(country_id, month) %>%
  summarize(mean_haz = mean(zlen, na.rm = TRUE), .groups = "drop")

# Plot
ggplot(mean_by_site_month, aes(x = month, y = mean_haz, color = country_id)) +
  geom_line() +
  labs(
    title = "Mean HAZ over time by country",
    x = "Month",
    y = "Mean HAZ"
  ) +
  theme_minimal()

# get rid of brazil
# looks comparable trend wise

