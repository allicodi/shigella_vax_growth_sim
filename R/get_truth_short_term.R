# ------------------------------------------------------------------------------
# Function to get true short term growth effect
# ------------------------------------------------------------------------------

here::i_am("R/simulate_data_short_term.R")

source(here::here("R/simulate_data_short_term.R"))

sim_data <- simulate_data(n = 1e5, type = "observed")
data <- sim_data$data
params <- sim_data$params

# mod to be as described in notes
og_data <- data %>%
  mutate(shig = if_else(S_sev == 1, 2, 
                        if_else(S_inf == 1, 1, 0))) %>%
  select(-S_inf, -S_sev)

# make long dataset
long_data <- expand.grid(id = 1:nrow(data), week = 1:52) %>%
  arrange(id) %>%
  left_join(og_data, by = "id") %>%
  filter(!(shig != 0 & week > S_inf_time)) %>%
  mutate(shig = if_else(S_inf_time != 0 & week != S_inf_time, 0, shig)) %>%
  mutate(I_shig = if_else(shig != 0, 1, 0)) %>%
  select(-S_inf_time)

# P(shig in week i | no shig in week i - 1, X) = expit(beta_0 + beta_1*X) 
# == Hazard??
weekly_hazard <- expand.grid(id = 1:nrow(data),
                             wk = 1:52) %>%
  arrange(id) %>%
  left_join(og_data[,c("id", "X", "Z")]) 

weekly_hazard_preds_Z0 <- data.frame(id = 1:nrow(data),
                                  setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = 52)),
                                           paste0("hazard_", 1:52)))

weekly_hazard_preds_Z1 <- data.frame(id = 1:nrow(data),
                                     setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = 52)),
                                              paste0("hazard_", 1:52)))

for(wk in 1:52){
  fit <- glm(I_shig ~ Z + X, 
             data = long_data[long_data$week <= wk,], 
             family = binomial())
  
  pred_data <- weekly_hazard[weekly_hazard$wk == wk,]
  
  pred_hazard_wk_Z1 <- 1 / (1 + exp(-(coef(fit)[1] + coef(fit)[["Z"]]*1 + coef(fit)[["X"]]*pred_data$X)))
  pred_hazard_wk_Z0 <- 1 / (1 + exp(-(coef(fit)[1] + coef(fit)[["Z"]]*0 + coef(fit)[["X"]]*pred_data$X)))
  
  if (wk > 1) {
    # density in week wk = hazard week k * prod(1 - hazard prev week for all prev weeks)
    pred_hazard_wk_Z1 <- pred_hazard_wk_Z1 * apply(as.matrix(weekly_hazard_preds_Z1[, 2:wk]), 1, function(x) prod(1 - x))
    pred_hazard_wk_Z0 <- pred_hazard_wk_Z0 * apply(as.matrix(weekly_hazard_preds_Z0[, 2:wk]), 1, function(x) prod(1 - x))
  }
  
  weekly_hazard_preds_Z1[,wk+1] <- pred_hazard_wk_Z1
  weekly_hazard_preds_Z0[,wk+1] <- pred_hazard_wk_Z0
  
}

# ...now idk what to do
# ...or if any of the above is correct
# also something in notes abt sum? sum hazard to infection time?
# then eventually need to model outcome HAZ itself? & incorporate longitudinal growth measurements?
