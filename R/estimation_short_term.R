# ------------------------------------------------------------------------------
# Function to estimate short term growth effects
# ------------------------------------------------------------------------------

here::i_am("R/estimation_short_term.R")

library(survival)
library(tidyverse)
source(here::here("R/simulate_parameters.R"))
source(here::here("R/simulate_data_short_term.R"))

# parameters
umax <- 52
pool <- TRUE

# type = counterfactual, observed, both
# params <- simulate_parameters(dose_schedule = "6mo")
# data <- simulate_data(parameters = params, n = 1e5, type = "observed")

params <- simulate_parameters(dose_schedule = "12mo", 
                              incidence_shigella_0_6 = 0.03731, 
                              incidence_shigella_6_12 = 0.05796, 
                              incidence_severe_shigella_0_6 = 0.02025, 
                              incidence_severe_shigella_6_12 = 0.02548)

data <- simulate_data(parameters = params, n = 1e5, type = "observed")

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

# two separate models for first and second half?
long_data_0_6 <- long_data[long_data$week <= 26,]
long_data_6_12 <- long_data[long_data$id %in% unique(long_data$id[long_data$week == 27]),] # remove people who were infected in the first half (they won't have any measurements 27 or later)

# P(shig in week i | no shig in week i - 1, X) = expit(beta_0 + beta_1*X) 
# == Hazard??
weekly_hazard <- expand.grid(id = 1:nrow(data),
                             wk = 1:52) %>%
  arrange(id) %>%
  left_join(og_data[,c("id", "X", "Z")]) 

weekly_density_preds_Z0_0_6 <- data.frame(id = 1:nrow(data),
                                     X = data$X,
                                    setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = 26)),
                                             paste0("lambda_", 1:26)))

weekly_density_preds_Z0_6_12 <- data.frame(id = 1:nrow(data),
                                         X = data$X,
                                         setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = 26)),
                                                  paste0("lambda_", 27:52)))

# weekly_hazard_preds_Z1_0_6 <- data.frame(id = 1:nrow(data),
#                                          X = data$X,
#                                          Z = 1,
#                                          setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = 26)),
#                                                   paste0("lambda_", 1:26)))
# 
# weekly_hazard_preds_Z1_6_12 <- data.frame(id = 1:nrow(data),
#                                           X = data$X,
#                                           Z = 1,
#                                           setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = 26)),
#                                                    paste0("lambda_", 27:52)))
# 
# weekly_hazard_preds_0_6 <- data.frame(id = 1:nrow(data),
#                                    X = data$X,
#                                    Z = data$Z,
#                                    setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = 26)),
#                                             paste0("lambda_", 1:26)))
# 
# weekly_hazard_preds_6_12 <- data.frame(id = 1:nrow(data),
#                                       X = data$X,
#                                       Z = data$Z,
#                                       setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = 26)),
#                                                paste0("lambda_", 27:52)))

fit_0_6_Z0 <- glm(I_shig ~ X, 
               data = long_data_0_6[long_data_0_6$Z == 0, ], 
               family = binomial())

fit_6_12_Z0 <- glm(I_shig ~ X, 
               data = long_data_6_12[long_data_6_12$Z == 0, ], 
               family = binomial())

haz_Z0_0_6 <- predict(fit_0_6_Z0, newdata = weekly_density_preds_Z0_0_6, type = 'response')
density_Z0_0_6 <- do.call(rbind, sapply(haz_Z0_0_6, function(haz){
  haz * (1 - haz)^((1:26) - 1)
}, simplify = FALSE))

haz_Z0_6_12 <- predict(fit_6_12_Z0, newdata = weekly_density_preds_Z0_0_6, type = 'response')
cuminc_Z0_0_6 <- rowSums(density_Z0_0_6)

density_Z0_6_12 <- do.call(rbind, 
  mapply(
    haz_6_12 = haz_Z0_6_12, cum_inc_0_6 = cuminc_Z0_0_6, 
    FUN = function(haz_6_12, cum_inc_0_6){
      haz_6_12 * (1 - haz_6_12)^((1:26) - 1) * (1 - cum_inc_0_6)
    }, SIMPLIFY = FALSE
  )
)

density_Z0_0_12 <- cbind(density_Z0_0_6, density_Z0_6_12)
# cuminc_0_12 <- rowSums(density_Z0_0_12)
# mean(cuminc_0_12)

# weekly_hazard_preds_Z0_6_12$lambda_27 <- predict(fit_6_12_Z0, newdata = weekly_hazard_preds_Z0_6_12, type = 'response')

# weekly_hazard_preds_Z1_0_6$lambda_1 <- predict(fit_0_6, newdata = weekly_hazard_preds_Z1_0_6, type = 'response')
# weekly_hazard_preds_Z1_6_12$lambda_27 <- predict(fit_6_12, newdata = weekly_hazard_preds_Z1_6_12, type = 'response')
# 
# weekly_hazard_preds_0_6$lambda_1 <- predict(fit_0_6, newdata = weekly_hazard_preds_0_6, type = 'response')
# weekly_hazard_preds_6_12$lambda_27 <- predict(fit_6_12, newdata = weekly_hazard_preds_6_12, type = 'response')

for(week in 2:26){
  pred_col <- paste0("lambda_", week)
  # constant for every time point
  weekly_hazard_preds_Z0_0_6[,pred_col] <- weekly_hazard_preds_Z0_0_6$lambda_1 * (1 - weekly_hazard_preds_Z0_0_6$lambda_1)^(week - 1)
  # weekly_hazard_preds_Z1_0_6[,pred_col] <- weekly_hazard_preds_Z1_0_6$lambda_1 * (1 - weekly_hazard_preds_Z1_0_6$lambda_1)^(week - 1)
  # weekly_hazard_preds_0_6[,pred_col] <- weekly_hazard_preds_0_6$lambda_1 * (1 - weekly_hazard_preds_0_6$lambda_1)^(week - 1)
}

for(week in 28:52){
  pred_col <- paste0("lambda_", week)
  # constant for every time point
  weekly_hazard_preds_Z0_6_12[,pred_col] <- weekly_hazard_preds_Z0_6_12$lambda_27 * (1 - weekly_hazard_preds_Z0_6_12$lambda_27)^(week - 1)
  # weekly_hazard_preds_Z1_6_12[,pred_col] <- weekly_hazard_preds_Z1_6_12$lambda_27 * (1 - weekly_hazard_preds_Z1_6_12$lambda_27)^(week - 1)
  # weekly_hazard_preds_6_12[,pred_col] <- weekly_hazard_preds_6_12$lambda_27 * (1 - weekly_hazard_preds_6_12$lambda_27)^(week - 1)
}

# Combine into single dataframe
#weekly_hazard_preds <- cbind(weekly_hazard_preds_0_6, weekly_hazard_preds_6_12[,-which(colnames(weekly_hazard_preds_6_12) %in% c("id","X","Z"))])
weekly_hazard_preds_Z0 <- cbind(weekly_hazard_preds_Z0_0_6, weekly_hazard_preds_Z0_6_12[,-which(colnames(weekly_hazard_preds_Z0_6_12) %in% c("id","X","Z"))])
#weekly_hazard_preds_Z1<- cbind(weekly_hazard_preds_Z1_0_6, weekly_hazard_preds_Z1_6_12[,-which(colnames(weekly_hazard_preds_Z1_6_12) %in% c("id","X","Z"))])

# maybe vectorize??
# for(week in 2:52){
#   pred_col <- paste0("lambda_", week)
#   # constant for every time point
#   weekly_hazard_preds_Z0[,pred_col] <- weekly_hazard_preds_Z0$lambda_1 * (1 - weekly_hazard_preds_Z0$lambda_1)^(week - 1)
#   weekly_hazard_preds_Z1[,pred_col] <- weekly_hazard_preds_Z1$lambda_1 * (1 - weekly_hazard_preds_Z1$lambda_1)^(week - 1)
#   weekly_hazard_preds[,pred_col] <- weekly_hazard_preds$lambda_1 * (1 - weekly_hazard_preds$lambda_1)^(week - 1)
# }

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
# # 
# P_Sumax_1__Z0_X <- rowSums(weekly_hazard_preds_Z0[,4:ncol(weekly_hazard_preds_Z0)])
# P_Sumax_1__Z0 <- mean(P_Sumax_1__Z0_X )

# ----------------------------------------------------------------------------
# Growth models
# ----------------------------------------------------------------------------

# E[Y_Vu | Z = 1, S_u-1 = 0, X] ----------------------------------------------

# 1. subset to vaccinated
data_Z1 <- data[data$Z == 1,]

pred_Y_out__Z1_Suminus1_0_X <- data.frame(id = 1:nrow(data),
                                          setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = umax)),
                                                   paste0("pred_", 1:umax)))
for(u in 1:umax){
  
  # 2. determine which growth measurement we are considering for outcome of the regression 
  #.   based on the value of u. ex. if u = 1 we select the three month growth outcome
  
  data_Z1$u_temp <- u
  data_Z1$Y_out <- ifelse(data_Z1$u_temp <= 4.3, data_Z1$Y_3,
                          ifelse(data_Z1$u_temp > 4.3 & data_Z1$u_temp <= 17.3, data_Z1$Y_6,
                                 ifelse(data_Z1$u_temp > 17.3 & data_Z1$u_temp <= 30.4, data_Z1$Y_9,
                                        ifelse(data_Z1$u_temp > 30.4 & data_Z1$u_temp <= umax, data_Z1$Y_12, NA))))
  
  # 3. Further subset the data to individuals who have not yet become ill by u (i.e. 
  #.   remain illness free through u - 1). This subset of data will potentially include individuals
  #.   who become infected either at u or between date and V(u)th visit. That's ok
  data_Z1_u <- data_Z1[which(data_Z1$S_inf_time >= u), ]
  
  # 4. Regress selected growth outcome on X in this subset of data. This gives your estimate of E[Y_Vu | Z = 1, S_u-1 = 0, X]
  # TODO could be smart and look to see if ppl infected in between
  fit_Y_out__Z1_Suminus1_0_X <- glm(Y_out ~ X, 
                                    data = data_Z1_u,
                                    family = gaussian())
  
  # predict on full data
  pred_Y_out__Z1_Suminus1_0_X[,u+1] <- predict(fit_Y_out__Z1_Suminus1_0_X, newdata = data)
}

# E[Y_Vu | Z = 0, S_u-1 = 1, X] ----------------------------------------------

### If not dense in u (version 1): ###

if(!pool){
  # 1. subset to unvaccinated
  data_Z0 <- data[data$Z == 0,]
  
  pred_Y_out__Z0_Suminus1_1_X <- data.frame(id = 1:nrow(data),
                                            setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = floor(umax))),
                                                     paste0("pred_", 1:floor(umax))))
  
  # # 2. further subset to all individuals who fall ill on any u in U-1(v)
  for(v in c(3,6,9,12)){
    
    # this should be min, max ranges for each v
    if(v == 3){
      v_umin <- 0; v_umax <- 4.3
    } else if(v == 6){
      v_umin <- 4.3; v_umax <- 17.3
    } else if(v == 9){
      v_umin <- 17.3; v_umax <- 30.4
    } else{
      # just make v_umax == umax i guess?? otherwise don't have preds for 44-52
      # v_umin <- 30.4; v_umax <- 43.5
      v_umin <- 30.4; v_umax <- umax
    }
    
    data_Z0_Uminus1_v <- data_Z0[which(data_Z0$S_inf_time > v_umin & data_Z0$S_inf_time <= v_umax),] # bigger than min and less than max
    
    # 3. create a variable for each individual, say Ti, indicating the time period at which each individual
    #    falls ill, ie set Ti equal to the u such that dSu,i = 1
    
    # same as S_inf_time
    data_Z0_Uminus1_v$Ti <- data_Z0_Uminus1_v$S_inf_time
    
    # QUESTION also assume need to get Y_out same way as the other one?
    
    if(v == 3){
      data_Z0_Uminus1_v$Y_out <- data_Z0_Uminus1_v$Y_3
    } else if(v == 6){
      data_Z0_Uminus1_v$Y_out <- data_Z0_Uminus1_v$Y_6
    } else if(v == 9){
      data_Z0_Uminus1_v$Y_out <- data_Z0_Uminus1_v$Y_9
    } else{
      data_Z0_Uminus1_v$Y_out <- data_Z0_Uminus1_v$Y_12
    }
    
    # 4. regress the selected growth outcome on X and T in this subset of data
    fit_Y_out__Z0_Suminus1_1_X_T <- glm(Y_out ~ X + Ti,
                                        data = data_Z0_Uminus1_v,
                                        family = gaussian())
    
    # 5. predicting from this model setting T = u for all individuals provides an estimate of E[Y_Vu | Z = 0, dSu = 1, X]
    # QUESTION but we have this * 4 V(u)s so unsure how to average at the end??
    
    for(u in ceiling(v_umin):floor(v_umax)){
      data_u <- data; data_u$Ti <- u
      pred_Y_out__Z0_Suminus1_1_X[,u+1] <- predict(fit_Y_out__Z0_Suminus1_1_X_T,
                                                   newdata = data_u,
                                                   type = 'response')
    }
  }
}

### If not dense in u (version 2- pool across visits): ###

if(pool){
  # 1. subset to unvaccinated
  data_Z0 <- data[data$Z == 0,]
  
  pred_Y_out__Z0_Suminus1_1_X <- data.frame(id = 1:nrow(data),
                                            setNames(as.data.frame(matrix(NA_real_, nrow = nrow(data), ncol = floor(umax))),
                                                     paste0("pred_", 1:floor(umax))))
  
  pred_data <- data.frame()
  
  # # 2. further subset to all individuals who fall ill on any u in U-1(v)
  for(v in c(3,6,9,12)){
    
    # this should be min, max ranges for each v
    if(v == 3){
      v_umin <- 0; v_umax <- 4.3
    } else if(v == 6){
      v_umin <- 4.3; v_umax <- 17.3
    } else if(v == 9){
      v_umin <- 17.3; v_umax <- 30.4
    } else{
      # just make v_umax == umax i guess?? otherwise don't have preds for 44-52
      # v_umin <- 30.4; v_umax <- 43.5
      v_umin <- 30.4; v_umax <- umax
    }
    
    data_Z0_Uminus1_v <- data_Z0[which(data_Z0$S_inf_time > v_umin & data_Z0$S_inf_time <= v_umax),] # bigger than min and less than max
    
    # 3. create a variable for each individual, say Ti, indicating the time period at which each individual
    #    falls ill, ie set Ti equal to the u such that dSu,i = 1
    
    # same as S_inf_time
    data_Z0_Uminus1_v$Ti <- data_Z0_Uminus1_v$S_inf_time
    
    # QUESTION also assume need to get Y_out same way as the other one?
    
    if(v == 3){
      data_Z0_Uminus1_v$Y_out <- data_Z0_Uminus1_v$Y_3
      data_Z0_Uminus1_v$V <- "v_3"
    } else if(v == 6){
      data_Z0_Uminus1_v$Y_out <- data_Z0_Uminus1_v$Y_6
      data_Z0_Uminus1_v$V <- "v_6"
    } else if(v == 9){
      data_Z0_Uminus1_v$Y_out <- data_Z0_Uminus1_v$Y_9
      data_Z0_Uminus1_v$V <- "v_9"
    } else{
      data_Z0_Uminus1_v$Y_out <- data_Z0_Uminus1_v$Y_12
      data_Z0_Uminus1_v$V <- "v_12"
    }
    
    pred_data <- rbind(pred_data, data_Z0_Uminus1_v)
  }
  
  pred_data$V <- factor(pred_data$V, levels = c("v_3", "v_6", "v_9", "v_12"))
  
  # 6. Regress the pooled growth outcome createde on X, V, and T
  fit_Y_out__Z0_Suminus1_1_X_T_V <- glm(Y_out ~ X + Ti + V ,
                                        data = pred_data,
                                        family = gaussian())
  
  for(u in 1:umax){
    
    temp_data <- data
    temp_data$Ti <- u
    
    # hardcoding u ranges for now??
    if(u %in% 1:4){
      temp_data$V <- "v_3"
    } else if(u %in% 5:17){
      temp_data$V <- "v_6"
    } else if(u %in% 18:30){
      temp_data$V <- "v_9"
    } else{
      temp_data$V <- "v_12"
    }
    
    pred_Y_out__Z0_Suminus1_1_X[,u+1] <- predict(fit_Y_out__Z0_Suminus1_1_X_T_V,
                                                 newdata = temp_data,
                                                 type = 'response')
    
  }
}

# Put the final estimand together ----------------------------------------------

# P(S_umax = 1 | Z = 0) == incidence of infection in the unvaccinated? pre week 43?

#make just preds not 4:end

# cuminc_0_12 <- rowSums(density_Z0_0_12)
# mean(cuminc_0_12)

P_Sumax_1__Z0_X <- rowSums(density_Z0_0_12)
#colnames(P_Sumax_1__Z0_X) <- paste0("lambda_",1:52)

P_Sumax_1__Z0 <- mean(P_Sumax_1__Z0_X )

# P_S_umax_1__Z_0 <- mean(data$S_inf[data$Z == 0 & data$S_inf_time<=umax])

estimate <- vector("numeric", length = umax)
for(u in 1:umax){
  P_dSu_1__Z_0_X <- weekly_hazard_preds_Z0[,paste0("lambda_",u)]
  E_Y_Vu__Z1_Suminus1_0_X <- pred_Y_out__Z1_Suminus1_0_X[,paste0("pred_",u)]
  E_Y_Vu__Z0_Suminus1_1_X <- pred_Y_out__Z0_Suminus1_1_X[,paste0("pred_",u)]
  
  estimate[u] <- mean((P_dSu_1__Z_0_X / P_Sumax_1__Z0) * (E_Y_Vu__Z1_Suminus1_0_X - E_Y_Vu__Z0_Suminus1_1_X))
}

final_growth_effect <- sum(estimate)

# 6 month vax schedule trial ---------

# Short term estimand:
# 0.04572702 with pooling (no interaction)
# 0.03738661 without pooling

# check equivalence with saturated
# 0.03957 with pooling Y_out ~ X * Ti * V
# 0.03957 without pooling Y_out ~ X * Ti

# 12 month vax schedule trial ------

# Short term estimand

# 0.04818 with pooling (no interaction)
# 0.04929 without pooling

# check equivalence with saturated

# 0.04953872 with pooling Y_out ~ X * Ti * V
# 0.04953872 without pooling Y_out ~ X * Ti

# ------------------------------------------------------------------------------

# Just look at 12mo for everyone??