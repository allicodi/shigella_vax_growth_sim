# ------------------------------------------------------------------------------
# Function to estimate short term growth effects
# ------------------------------------------------------------------------------

here::i_am("R/estimation_short_term.R")

source(here::here("R/simulate_data_short_term.R"))

# type = counterfactual, observed, both
sim_data <- simulate_data(n = 1e4, type = "observed")
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
                                     X = data$X,
                                     Z = 0,
                                  setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = 52)),
                                           paste0("lambda_", 1:52)))

weekly_hazard_preds_Z1 <- data.frame(id = 1:nrow(data),
                                     X = data$X,
                                     Z = 1,
                                     setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = 52)),
                                              paste0("lambda_", 1:52)))

weekly_hazard_preds <- data.frame(id = 1:nrow(data),
                                   X = data$X,
                                   Z = data$Z,
                                   setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = 52)),
                                            paste0("lambda_", 1:52)))


fit <- glm(I_shig ~ Z + X, 
           data = long_data, 
           family = binomial())

weekly_hazard_preds_Z0$lambda_1 <- predict(fit, newdata = weekly_hazard_preds_Z0, type = 'response')
weekly_hazard_preds_Z1$lambda_1 <- predict(fit, newdata = weekly_hazard_preds_Z1, type = 'response')
weekly_hazard_preds$lambda_1 <- predict(fit, newdata = data, type = 'response')

for(week in 2:52){
  pred_col <- paste0("lambda_", week)
  # constant for every time point
  weekly_hazard_preds_Z0[,pred_col] <- weekly_hazard_preds_Z0$lambda_1 * (1 - weekly_hazard_preds_Z0$lambda_1)^(week - 1)
  weekly_hazard_preds_Z1[,pred_col] <- weekly_hazard_preds_Z1$lambda_1 * (1 - weekly_hazard_preds_Z1$lambda_1)^(week - 1)
  weekly_hazard_preds[,pred_col] <- weekly_hazard_preds$lambda_1 * (1 - weekly_hazard_preds$lambda_1)^(week - 1)
}

# old version too complicated
# for(wk in 1:52){
#   fit <- glm(I_shig ~ Z + X, 
#              data = long_data[long_data$week <= wk,], 
#              family = binomial())
#   
#   pred_data <- weekly_hazard[weekly_hazard$wk == wk,]
#   
#   pred_hazard_wk_Z1 <- 1 / (1 + exp(-(coef(fit)[1] + coef(fit)[["Z"]]*1 + coef(fit)[["X"]]*pred_data$X)))
#   pred_hazard_wk_Z0 <- 1 / (1 + exp(-(coef(fit)[1] + coef(fit)[["Z"]]*0 + coef(fit)[["X"]]*pred_data$X)))
#   
#   if (wk > 1) {
#     # density in week wk = hazard week k * prod(1 - hazard prev week for all prev weeks)
#     pred_hazard_wk_Z1 <- pred_hazard_wk_Z1 * apply(as.matrix(weekly_hazard_preds_Z1[, 2:wk]), 1, function(x) prod(1 - x))
#     pred_hazard_wk_Z0 <- pred_hazard_wk_Z0 * apply(as.matrix(weekly_hazard_preds_Z0[, 2:wk]), 1, function(x) prod(1 - x))
#   }
#   
#   weekly_hazard_preds_Z1[,wk+1] <- pred_hazard_wk_Z1
#   weekly_hazard_preds_Z0[,wk+1] <- pred_hazard_wk_Z0
#   
# }

# for now try:
# For each ID, sum 1:52 weeks, avg of sums?

# P_Sumax_1__Z1_X <- rowSums(weekly_hazard_preds_Z1[,4:ncol(weekly_hazard_preds_Z1)])
# P_Sumax_1__Z1 <- mean(P_Sumax_1__Z1_X )
# 
# P_Sumax_1__Z0_X <- rowSums(weekly_hazard_preds_Z0[,4:ncol(weekly_hazard_preds_Z0)])
# P_Sumax_1__Z0 <- mean(P_Sumax_1__Z0_X )

# ----------------------------------------------------------------------------
# Growth models
# ----------------------------------------------------------------------------

# E[Y_Vu | Z = 1, S_u-1 = 0, X] ----------------------------------------------

# 1. subset to vaccinated
data_Z1 <- data[data$Z == 1,]

# 2. determine which growth measurement we are considering for outcome of the regression 
#.   based on the value of u. ex. if u = 1 we select the three month growth outcome

data_Z1$X_out <- ifelse(data_Z1$S_inf_time <= 4.3, data_Z1$X_3,
                       ifelse(data_Z1$S_inf_time > 4.3 & data_Z1$S_inf_time <= 17.3, data_Z1$X_6,
                              ifelse(data_Z1$S_inf_time > 17.3 & data_Z1$S_inf_time <= 30.4, data_Z1$X_9,
                                     ifelse(data_Z1$S_inf_time > 30.4 & data_Z1$S_inf_time <= 43.5, data_Z1$X_12, 
                                            ifelse(data_Z1$S_inf_time > 43.5 & data_Z1$S_inf == 1, NA, data_Z1$X_12)))))

# QUESTION if infected at >= 43.5, drop??
data_Z1 <- data_Z1[-which(is.na(data_Z1$X_out)),]

# 3. Further subset the data to individuals who have not yet become ill by u (i.e. 
#.   remain illness free through u - 1). This subset of data will potentially include individuals
#.   who become infected either at u or between date and V(u)th visit. That's ok

pred_X_out__Z1_Suminus1_0_X <- data.frame(id = 1:nrow(data),
                                          setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = 52)),
                                                   paste0("pred_", 1:52)))
for(u in 1:52){
  data_Z1_u <- data_Z1[which(data_Z1$S_inf_time >= u), ]
  
  # 4. Regress selected growth outcome on X in this subset of data. This gives your estimate of E[Y_Vu | Z = 1, S_u-1 = 0, X]
  # TODO could be smart and look to see if ppl infected in between
  fit_X_out__Z1_Suminus1_0_X <- glm(X_out ~ X, 
                                    data = data_Z1_u,
                                    family = gaussian())
  
  # QUESTION Now predict on full data??
  pred_X_out__Z1_Suminus1_0_X[,u+1] <- predict(fit_X_out__Z1_Suminus1_0_X, newdata = data)
}

E_X_out__Z1_Suminus1_0_X <- colMeans(pred_X_out__Z1_Suminus1_0_X[,2:ncol(pred_X_out__Z1_Suminus1_0_X)])

# E[Y_Vu | Z = 0, S_u-1 = 1, X] ----------------------------------------------

### If not dense in u (version 1): ###

# unclear if implemented correctly/how to average at the end, think i misinterpreted something in the directions

# 1. subset to unvaccinated
# data_Z0 <- data[data$Z == 0,]
# 
# pred_X_out__Z0_Suminus1_1_X <- vector("list", length = 4)
# names(pred_X_out__Z0_Suminus1_1_X) <- c("v_3", "v_6", "v_9", "v_12")
# 
# # 2. further subset to all individuals who fall ill on any u in U-1(v)
# for(v in c(3,6,9,12)){
#   
#   if(v == 3){
#     umax <- 4.3
#   } else if(v == 6){
#     umax <- 17.3
#   } else if(v == 9){
#     umax <- 30.4
#   } else{
#     umax <- 43.5
#   }
#   
#   data_Z0_Uminus1_v <- data_Z0[which(data_Z0$S_inf_time <= umax),]
#   
#   # 3. create a variable for each individual, say Ti, indicating the time period at which each individual 
#   #    falls ill, ie set Ti equal to the u such that dSu,i = 1
#   
#   # QUESTION is that not already S_inf_time??
#   data_Z0_Uminus1_v$Ti <- data_Z0_Uminus1_v$S_inf_time
#   
#   # QUESTION also assume need to get X_out same way as the other one?
#   data_Z0_Uminus1_v$X_out <- ifelse(data_Z0_Uminus1_v$Ti <= 4.3, data_Z0_Uminus1_v$X_3,
#                           ifelse(data_Z0_Uminus1_v$Ti > 4.3 & data_Z0_Uminus1_v$Ti <= 17.3, data_Z0_Uminus1_v$X_6,
#                                  ifelse(data_Z0_Uminus1_v$Ti > 17.3 & data_Z0_Uminus1_v$Ti <= 30.4, data_Z0_Uminus1_v$X_9,
#                                         ifelse(data_Z0_Uminus1_v$Ti > 30.4 & data_Z0_Uminus1_v$Ti <= 43.5, data_Z0_Uminus1_v$X_12, 
#                                                ifelse(data_Z0_Uminus1_v$Ti > 43.5 & data_Z0_Uminus1_v$S_inf == 1, NA, data_Z0_Uminus1_v$X_12)))))
#   
#   
#   # 4. regress the selected growth outcome on X and T in this subset of data
#   fit_X_out__Z0_Suminus1_1_X_T <- glm(X_out ~ X + Ti,
#                                       data = data_Z0_Uminus1_v,
#                                       family = gaussian())
#   
#   # 5. predicting from this model setting T = u for all individuals provides an estimate of E[Y_Vu | Z = 0, dSu = 1, X]
#   # QUESTION but we have this * 4 V(u)s so unsure how to average at the end??
#   
#   pred_X_out__Z0_Suminus1_1_X_v <- data.frame(id = 1:nrow(data),
#                                               setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = floor(umax))),
#                                                        paste0("pred_", 1:floor(umax))))
#   for(u in 1:floor(umax)){
#     data_u <- data; data$Ti <- u
#     pred_X_out__Z0_Suminus1_1_X_v[,u+1] <- predict(fit_X_out__Z0_Suminus1_1_X_T,
#                                                    newdata = data_u,
#                                                    type = 'response')
#   }
#   
#   pred_X_out__Z0_Suminus1_1_X[[paste0("v_",v)]] <- pred_X_out__Z0_Suminus1_1_X_v
#   
# }

# not sure how to average bc we have set of preds for each v?

### If not dense in u (version 2- pool across visits): ###

# 1. Subset data to unvaccinated individuals
data_Z0 <- data[data$Z == 0,]

# 2. Further subset data to all individuals who fall ill on any u <= umax
# QUESTION is umax == 52? so that's just everyone who is infected?
# or is it 43.5 so eligible to look at end growth?
umax <- 43.5
data_Z0_umax <- data_Z0[which(data_Z0$S_inf == 1 & data_Z0$S_inf_time <= umax),]

# 5. Create a variable, say Ti (as above), indicating the exact time period at which the individual was infected.
data_Z0_umax$Ti <- data_Z0_umax$S_inf_time

# 3. For each individual i, identify the growth outcome that is appropriately proximal to their illness date.
# i.e., let Ti denote the time period at which individual i became ill, then the outcome for this pooled
# model is Yi = YV (Ti). In other words, if you become ill in the first month, we set your outcome to Y3,
# if you become ill in months 2-4, we set your outcome to Y6, etc...

data_Z0_umax$X_out <- ifelse(data_Z0_umax$Ti <= 4.3, data_Z0_umax$X_3,
                                  ifelse(data_Z0_umax$Ti > 4.3 & data_Z0_umax$Ti <= 17.3, data_Z0_umax$X_6,
                                         ifelse(data_Z0_umax$Ti > 17.3 & data_Z0_umax$Ti <= 30.4, data_Z0_umax$X_9,
                                                ifelse(data_Z0_umax$Ti > 30.4 & data_Z0_umax$Ti <= 43.5, data_Z0_umax$X_12, 
                                                       ifelse(data_Z0_umax$Ti > 43.5 & data_Z0_umax$S_inf == 1, NA, data_Z0_umax$X_12)))))

# 4. Create a variable, say Vi, indicating which visit’s outcome was used for each individual.
data_Z0_umax$V <- ifelse(data_Z0_umax$Ti <= 4.3, "v_3",
                                  ifelse(data_Z0_umax$Ti > 4.3 & data_Z0_umax$Ti <= 17.3, "v_6",
                                         ifelse(data_Z0_umax$Ti > 17.3 & data_Z0_umax$Ti <= 30.4, "v_9",
                                                ifelse(data_Z0_umax$Ti > 30.4 & data_Z0_umax$Ti <= 43.5, "v_12", 
                                                       ifelse(data_Z0_umax$Ti > 43.5 & data_Z0_umax$S_inf == 1, NA, "v_12")))))

data_Z0_umax$V <- factor(data_Z0_umax$V, levels = c("v_3", "v_6", "v_9", "v_12"))

# 6. Regress the pooled growth outcome createde on X, V, and T
fit_X_out__Z0_Suminus1_1_X_T_V <- glm(X_out ~ X + Ti + V,
                                    data = data_Z0_umax,
                                    family = gaussian())

# 7. Predicting from this model setting T = u and V = V(u) for all individuals provides 
#.   an estimate of E[Y_V(u) | Z = 0, dS_u = 1, X] 

# QUESTION should this be just combos that apply to each V? so u=1-4 --> V_3, u=5-17 --> V_6, ...

# add Ti and V to full dataset
data$Ti <- data$S_inf_time
data$V <- ifelse(data$Ti <= 4.3, "v_3",
                  ifelse(data$Ti > 4.3 & data$Ti <= 17.3, "v_6",
                         ifelse(data$Ti > 17.3 & data$Ti <= 30.4, "v_9",
                                ifelse(data$Ti > 30.4 & data$Ti <= 43.5, "v_12", 
                                       ifelse(data$Ti > 43.5 & data$S_inf == 1, NA, "v_12")))))

pred_X_out__Z0_Suminus1_1_X <- data.frame(id = 1:nrow(data),
                      setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = floor(umax))),
                               paste0("pred_", 1:floor(umax))))

i <- 2

for(v in c("v_3", "v_6", "v_9", "v_12")){
  if(v == "v_3"){
    for(u in 1:4){
      pred_data <- data
      pred_data$V <- v; pred_data$Ti <- u
      
      pred_X_out__Z0_Suminus1_1_X[,i] <- predict(fit_X_out__Z0_Suminus1_1_X_T_V, 
                             newdata = pred_data, 
                             type = 'response')
      i <- i+1
    }
  } else if(v == "v_6"){
    for(u in 5:17){
      pred_data <- data
      pred_data$V <- v; pred_data$Ti <- u
      
      pred_X_out__Z0_Suminus1_1_X[,i] <- predict(fit_X_out__Z0_Suminus1_1_X_T_V, 
                             newdata = pred_data, 
                             type = 'response')
      i <- i+1
    }
  } else if(v == "v_9"){
    for(u in 18:30){
      pred_data <- data
      pred_data$V <- v; pred_data$Ti <- u
      
      pred_X_out__Z0_Suminus1_1_X[,i] <- predict(fit_X_out__Z0_Suminus1_1_X_T_V, 
                             newdata = pred_data, 
                             type = 'response')
      i <- i+1
    }
  } else{
    for(u in 31:43){
      pred_data <- data
      pred_data$V <- v; pred_data$Ti <- u
      
      pred_X_out__Z0_Suminus1_1_X[,i] <- predict(fit_X_out__Z0_Suminus1_1_X_T_V, 
                             newdata = pred_data, 
                             type = 'response')
      i <- i+1
    }
  }
}

E_X_out__Z0_Suminus1_1_X <- colMeans(pred_X_out__Z0_Suminus1_1_X[,2:ncol(pred_X_out__Z0_Suminus1_1_X)])

# ------------------------------------------------------------------------

# Put the final estimand together

# P(S_umax = 1 | Z = 0) == incidence of infection in the unvaccinated? pre week 43?

# but has to be end of study otherwise 1? so maybe umax shouldn't be 43?
P_S_umax_1__Z_0 <- mean(data$S_inf[data$Z == 0])

# but can only go to 43 here??
estimate <- vector("numeric", length = 43)
for(u in 1:43){
  P_dSu_1__Z_0_X <- weekly_hazard_preds_Z0[,paste0("lambda_",u)]
  E_Y_Vu__Z1_Suminus1_0_X <- pred_X_out__Z1_Suminus1_0_X[,paste0("pred_",u)]
  E_Y_Vu__Z0_Suminus1_1_X <- pred_X_out__Z0_Suminus1_1_X[,paste0("pred_",u)]
  
  estimate[u] <- mean((P_dSu_1__Z_0_X / P_S_umax_1__Z_0) * (E_Y_Vu__Z1_Suminus1_0_X - E_Y_Vu__Z0_Suminus1_1_X))
}

final_growth_effect <- sum(estimate)

# If dense in u (which it won't be): 

# # 1. subset to unvaccinated
# data_Z0 <- data[data$Z == 0,]
# 
# # 2. determine which growth measurement we are considering for outcome of the regression 
# #.   based on the value of u. ex. if u = 1 we select the three month growth outcome
# 
# data_Z0$X_out <- ifelse(data_Z0$S_inf_time <= 4.3, data_Z0$X_3,
#                         ifelse(data_Z0$S_inf_time > 4.3 & data_Z0$S_inf_time <= 17.3, data_Z0$X_6,
#                                ifelse(data_Z0$S_inf_time > 17.3 & data_Z0$S_inf_time <= 30.4, data_Z0$X_9,
#                                       ifelse(data_Z0$S_inf_time > 30.4 & data_Z0$S_inf_time <= 43.5, data_Z0$X_12, 
#                                              ifelse(data_Z0$S_inf_time > 43.5 & data_Z0$S_inf == 1, NA, data_Z0$X_12)))))
# 
# # QUESTION if infected at >= 43.5, drop??
# data_Z0 <- data_Z0[-which(is.na(data_Z0$X_out)),]
# 
# # 3. Further subset data to individuals who fall ill on u
# 
# pred_X_out__Z0_Suminus1_1_X <- data.frame(id = 1:nrow(data),
#                                           setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = 52)),
#                                                    paste0("pred_", 1:52)))
# 
# for(u in 1:52){
#   data_Z0_u <- data_Z0[which(data_Z0$S_inf_time == u), ]
#   
#   # 4. Regress selected growth outcome on X in this subset of data. This gives your estimate of E[Y_Vu | Z = 1, S_u-1 = 0, X]
#   # TODO could be smart and look to see if ppl infected in between
#   fit_X_out__Z0_Suminus1_1_X <- glm(X_out ~ X, 
#                                     data = data_Z0_u,
#                                     family = gaussian())
#   
#   # Now predict on full data??
#   pred_X_out__Z0_Suminus1_1_X[,u+1] <- predict(fit_X_out__Z0_Suminus1_1_X, newdata = data)
# }
