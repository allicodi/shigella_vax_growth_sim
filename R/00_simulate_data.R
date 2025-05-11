# -----------------------------------------------------------------------
# Function to simulate basic data that follows DAG for testing 
# -----------------------------------------------------------------------

# Default numbers from Bangladesh

generate_EFGH_complex <- function(n = 1000,
                                  incidence_shigella = 0.0618,
                                  incidence_severe_shigella = 0.0240,
                                  VE_mild = 0.4,
                                  VE_severe = 0.6,
                                  lfaz_bl_effect_shig = -0.3147,
                                  lfaz_bl_effect_severe_shig = -0.1990,
                                  shig_coef_mild = -0.052,
                                  shig_coef_severe = -0.072,
                                  sd_growth = "baseline", 
                                  zhifei = FALSE,
                                  catch_up = FALSE){

  # Vaccination (RCT 1:1)
  V <- rbinom(n, size=1, prob=0.5)
  
  if(zhifei == TRUE){
    # assign an age group
    age <- sample(c("young", "mid", "old"), n, replace = TRUE)
    
    X <- vector("numeric", length = n)
    mean_X <- vector("numeric", length = n)
    
    # distribution of X varies by age group
    # mean/sd from MAL-ED bangladesh
    # young = 6-11 mo, mid = 12mo-2yrs11mo, old = 3-5years
    mean_X <- ifelse(age == "young", -0.3260, 
                     ifelse(age == "mid", -0.7212,
                            -0.8776))
    sd_X <- ifelse(age == "young", 0.9684, 
                   ifelse(age == "mid", 0.8747,
                          0.8529))
    
    X <- rnorm(n, mean_X, sd_X)
    
  } else {
    # lfazscore at baseline from efgh data - ages 9-12 months
    mean_X <- -0.9395
    sd_X <- 1.1763
    X <- rnorm(n, mean_X, sd_X)
  }

  # Shigella infection (Y_inf) in vaccinated and unvaccinated
  
  # get intercept baseline log-odds 
  bl_log_odds_shigella <- qlogis(incidence_shigella)
  
  prob_Y0is1 <- plogis(bl_log_odds_shigella + lfaz_bl_effect_shig*(X - mean_X))
  Y_inf_V0 <- rbinom(n, 1, prob_Y0is1)
  prob_Y1is1 <- rep(0, n)
  prob_Y1is1[Y_inf_V0 == 1] <- 1 - VE_mild
  Y_inf_V1 <- rbinom(n, 1, prob_Y1is1)
  
  # Shigella severity (Y_sev) in vaccinated and unvaccinated
  bl_log_odds_severe <- qlogis(incidence_severe_shigella/incidence_shigella)
  
  prob_Y_severe_V0 <- ifelse(Y_inf_V0 == 1,
                             plogis(bl_log_odds_severe + lfaz_bl_effect_severe_shig*(X - mean_X)),
                             0)
  Y_sev_V0 <- rbinom(n, 1, prob_Y_severe_V0)
  prob_Y_severe_V1 <- rep(0, n)
  prob_Y_severe_V1[Y_sev_V0 == 1 & Y_inf_V1 == 1] <- (1 - VE_severe) / (1 - VE_mild)
  Y_sev_V1 <- rbinom(n, 1, prob_Y_severe_V1)

  Y_inf <- rep(NA, n)
  Y_inf[V == 1] <- Y_inf_V1[V == 1]
  Y_inf[V == 0] <- Y_inf_V0[V == 0]

  Y_severe <- rep(NA, n)
  Y_severe[V == 1] <- Y_sev_V1[V == 1]
  Y_severe[V == 0] <- Y_sev_V0[V == 0]
  
  # Shigella infection timing - (Y_time)
  Y_time <- sample(seq(1, 12), n, replace = TRUE)
  Y_time[Y_inf == 0] <- 0
  
  infection_adjustment <- shig_coef_mild*I(Y_inf == 1 & Y_severe == 0) + shig_coef_severe*I(Y_inf == 1 & Y_severe == 1)
  
  # Allow or prohibit catch up growth-- smaller effect size if allowing for catch-up
  if(catch_up == FALSE){
    
    # No catch-up growth
    
    # If sd_growth = "baseline", use standard deviation month 24 ~ month 12 MAL-ED Bangladesh (0.4389)
    # Else If sd_growth = "prev-1", use models fit with previous 1 timepoint
    # Else If sd_growth = "prev-2", use models fit with previous 2 timepoints
    # Else assume sd_growth is numeric sd to be used with no intermediate measurements 
    if(sd_growth == "baseline"){
      G <- X + infection_adjustment + rnorm(n, mean = 0, sd = 0.4389)
    } else if (sd_growth == "prev-1") {
      # Growth trajectory from MAL-ED Bangladesh 
      
      # Previous time point
      X_vec <- vector(mode = "list", length = 12) 
      X_vec[[1]]  <- -0.13083351 + 0.9690976*X           + rnorm(n, mean = 0, sd = 0.2120328)
      X_vec[[2]]  <- -0.13359784 + 0.9658743*X_vec[[1]]  + rnorm(n, mean = 0, sd = 0.2096142)
      X_vec[[3]]  <- -0.02840903 + 0.9964731*X_vec[[2]]  + rnorm(n, mean = 0, sd = 0.2143073)
      X_vec[[4]]  <- -0.12824842 + 0.9658842*X_vec[[3]]  + rnorm(n, mean = 0, sd = 0.1941224)
      X_vec[[5]]  <- -0.10182843 + 0.9745330*X_vec[[4]]  + rnorm(n, mean = 0, sd = 0.1732090)
      X_vec[[6]]  <- -0.07077051 + 0.9734801*X_vec[[5]]  + rnorm(n, mean = 0, sd = 0.1712597)
      X_vec[[7]]  <- -0.04753367 + 0.9912921*X_vec[[6]]  + rnorm(n, mean = 0, sd = 0.1616939)
      X_vec[[8]]  <- -0.04551516 + 0.9821661*X_vec[[7]]  + rnorm(n, mean = 0, sd = 0.1704176)
      X_vec[[9]]  <- -0.04299716 + 0.9801243*X_vec[[8]]  + rnorm(n, mean = 0, sd = 0.1597872)
      X_vec[[10]] <- -0.03899147 + 0.9936665*X_vec[[9]]  + rnorm(n, mean = 0, sd = 0.1516835)
      X_vec[[11]] <- -0.04331684 + 0.9741949*X_vec[[10]] + rnorm(n, mean = 0, sd = 0.1652360)
      X_vec[[12]] <- -0.08162686 + 0.9565905*X_vec[[11]] + rnorm(n, mean = 0, sd = 0.1991225)
      
      X_df <- data.frame(X_vec)
      colnames(X_df) <- c("X_1", "X_2", "X_3", "X_4", "X_5", "X_6",
                          "X_7", "X_8", "X_9", "X_10", "X_11", "X_12")
      X_df <- cbind(X_df, Y_time, infection_adjustment)
      
      # Subtract infection adjustment from all time points beyond infection
      adjust_for_inf <- function(row){
        my_Y_time <- row['Y_time']
        if(my_Y_time > 0){
          row[(my_Y_time:12)] <- row[(my_Y_time:12)] + row['infection_adjustment']
        }
        return(row)
      }
      
      X_df <- data.frame(t(apply(X_df, MARGIN = 1, adjust_for_inf)))
      
      G <- X_df$X_12
      
    }else if(sd_growth == "prev-2"){
      X_vec <- vector(mode = "list", length = 12) 
      
      # Previous 2 timepoints
      X_vec[[1]]  <- -0.13083351 + 0.9690976*X          + rnorm(n, mean = 0, sd = 0.2120328)
      X_vec[[2]]  <- -0.13663765 + 0.10708213*X         + 0.8610326*X_vec[[1]]   + rnorm(n, mean = 0, sd = 0.2086736)
      X_vec[[3]]  <- -0.03902150 + 0.01465246*X_vec[[1]]  + 0.9779293*X_vec[[2]]   + rnorm(n, mean = 0, sd = 0.2141136)
      X_vec[[4]]  <- -0.11403457 + 0.20063834*X_vec[[2]]  + 0.7757495*X_vec[[3]]   + rnorm(n, mean = 0, sd = 0.1899936)
      X_vec[[5]]  <- -0.10566722 + 0.04466520*X_vec[[3]]  + 0.9287577*X_vec[[4]]   + rnorm(n, mean = 0, sd = 0.1740689)
      X_vec[[6]]  <- -0.07251206 + 0.04606467*X_vec[[4]]  + 0.9277831*X_vec[[5]]   + rnorm(n, mean = 0, sd = 0.1710707)
      X_vec[[7]]  <- -0.04799386 + 0.06796371*X_vec[[5]]  + 0.9240339*X_vec[[6]]   + rnorm(n, mean = 0, sd = 0.1619633)
      X_vec[[8]]  <- -0.04294224 + 0.13528495*X_vec[[6]]  + 0.8505967*X_vec[[7]]   + rnorm(n, mean = 0, sd = 0.1698145)
      X_vec[[9]]  <- -0.03597587 + 0.14149890*X_vec[[7]]  + 0.8421613*X_vec[[8]]   + rnorm(n, mean = 0, sd = 0.1583651)
      X_vec[[10]] <- -0.03382477 + 0.18068365*X_vec[[8]]  + 0.8170421*X_vec[[9]]   + rnorm(n, mean = 0, sd = 0.1455196)
      X_vec[[11]] <- -0.04020125 + 0.17402791*X_vec[[9]]  + 0.8045709*X_vec[[10]]  + rnorm(n, mean = 0, sd = 0.1630071)
      X_vec[[12]] <- -0.07027040 + 0.24234964*X_vec[[10]] + 0.7192027*X_vec[[11]] + rnorm(n, mean = 0, sd = 0.1905827)
      
      X_df <- data.frame(X_vec)
      colnames(X_df) <- c("X_1", "X_2", "X_3", "X_4", "X_5", "X_6",
                          "X_7", "X_8", "X_9", "X_10", "X_11", "X_12")
      X_df <- cbind(X_df, Y_time, infection_adjustment)
      
      # Subtract infection adjustment from all time points beyond infection
      adjust_for_inf <- function(row){
        my_Y_time <- row['Y_time']
        if(my_Y_time > 0){
          row[(my_Y_time:12)] <- row[(my_Y_time:12)] + row['infection_adjustment']
        }
        return(row)
      }
      
      X_df <- data.frame(t(apply(X_df, MARGIN = 1, adjust_for_inf)))
      
      G <- X_df$X_12
      
    } else{
      G <- X + infection_adjustment + rnorm(n, mean = 0, sd = sd_growth)
    }
  } else{
    # Allow for catch-up growth
    
    if(sd_growth == "baseline"){
      G <- X + shig_coef_mild*I(Y_inf == 1 & Y_severe == 0) + shig_coef_severe*I(Y_inf == 1 & Y_severe == 1) + rnorm(n, mean = 0, sd = 0.4389)
    } else if (sd_growth == "prev-1") {
      # Growth trajectory from MAL-ED Bangladesh 
      
      # Previous time point
      X_1 <- -0.13083351 + 0.9690976*X + I(Y_time == 1)*infection_adjustment + rnorm(n, mean = 0, sd = 0.2120328)
      X_2 <- -0.13359784 + 0.9658743*X_1 + I(Y_time == 2)*infection_adjustment + rnorm(n, mean = 0, sd = 0.2096142)
      X_3 <- -0.02840903 + 0.9964731*X_2 + I(Y_time == 3)*infection_adjustment + rnorm(n, mean = 0, sd = 0.2143073)
      X_4 <- -0.12824842 + 0.9658842*X_3 + I(Y_time == 4)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1941224)
      X_5 <- -0.10182843 + 0.9745330*X_4 + I(Y_time == 5)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1732090)
      X_6 <- -0.07077051 + 0.9734801*X_5 + I(Y_time == 6)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1712597)
      X_7 <- -0.04753367 + 0.9912921*X_6 + I(Y_time == 7)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1616939)
      X_8 <- -0.04551516 + 0.9821661*X_7 + I(Y_time == 8)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1704176)
      X_9 <- -0.04299716 + 0.9801243*X_8 + I(Y_time == 9)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1597872)
      X_10 <- -0.03899147 + 0.9936665*X_9 + I(Y_time == 10)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1516835)
      X_11 <- -0.04331684 + 0.9741949*X_10 + I(Y_time == 11)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1652360)
      
      G <- -0.08162686 + 0.9565905*X_11 + I(Y_time == 12)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1991225)
      
    }else if(sd_growth == "prev-2"){
      # Previous 2 timepoints
      X_1 <- -0.13083351 + 0.9690976*X + I(Y_time == 1)*infection_adjustment + rnorm(n, mean = 0, sd = 0.2120328)
      X_2 <- -0.13663765 + 0.10708213*X + 0.8610326*X_1 + I(Y_time == 2)*infection_adjustment + rnorm(n, mean = 0, sd = 0.2086736)
      X_3 <- -0.03902150 + 0.01465246*X_1 + 0.9779293*X_2 + I(Y_time == 3)*infection_adjustment + rnorm(n, mean = 0, sd = 0.2141136)
      X_4 <- -0.11403457 + 0.20063834*X_2 + 0.7757495*X_3 + I(Y_time == 4)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1899936)
      X_5 <- -0.10566722  + 0.04466520*X_3 + 0.9287577*X_4 + I(Y_time == 5)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1740689)
      X_6 <- -0.07251206 + 0.04606467*X_4 + 0.9277831*X_5 + I(Y_time == 6)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1710707)
      X_7 <- -0.04799386 + 0.06796371*X_5 + 0.9240339*X_6 + I(Y_time == 7)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1619633)
      X_8 <- -0.04294224 + 0.13528495*X_6 + 0.8505967*X_7 + I(Y_time == 8)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1698145)
      X_9 <- -0.03597587 + 0.14149890*X_7 + 0.8421613*X_8 + I(Y_time == 9)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1583651)
      X_10 <- -0.03382477 + 0.18068365*X_8 + 0.8170421*X_9 + I(Y_time == 10)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1455196)
      X_11 <- -0.04020125 + 0.17402791*X_9 + 0.8045709*X_10 + I(Y_time == 11)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1630071)
      
      G <- -0.07027040 + 0.24234964*X_10 + 0.7192027*X_11 + I(Y_time == 12)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1905827)
    } else{
      G <- X + infection_adjustment + rnorm(n, mean = 0, sd = sd_growth)
    }
    
  }
  
  final_df <- data.frame(id = 1:n,
                         V = V,
                         X = X,
                         Y_inf = Y_inf,
                         Y_sev = Y_severe,
                         G = G)
  if(zhifei == TRUE){
    final_df <- cbind(final_df, age)
  }
  
  return(final_df)
}


generate_EFGH_complex_quad <- function(n = 1000,
                                  incidence_shigella = 0.0618,
                                  incidence_severe_shigella = 0.0240,
                                  VE_mild = 0.4,
                                  VE_severe = 0.6,
                                  lfaz_bl_effect_shig = -0.3147,
                                  lfaz_bl_effect_severe_shig = -0.1990,
                                  shig_coef_mild = -0.052,
                                  shig_coef_severe = -0.072,
                                  sd_growth = "baseline", 
                                  zhifei = FALSE,
                                  catch_up = FALSE){
  
  # Vaccination (RCT 1:1)
  V <- rbinom(n, size=1, prob=0.5)
  
  if(zhifei == TRUE){
    # assign an age group
    age <- sample(c("young", "mid", "old"), n, replace = TRUE)
    
    X <- vector("numeric", length = n)
    mean_X <- vector("numeric", length = n)
    
    # distribution of X varies by age group
    # mean/sd from MAL-ED bangladesh
    # young = 6-11 mo, mid = 12mo-2yrs11mo, old = 3-5years
    mean_X <- ifelse(age == "young", -0.3260, 
                     ifelse(age == "mid", -0.7212,
                            -0.8776))
    sd_X <- ifelse(age == "young", 0.9684, 
                   ifelse(age == "mid", 0.8747,
                          0.8529))
    
    X <- rnorm(n, mean_X, sd_X)
    
  } else {
    # lfazscore at baseline from efgh data - ages 9-12 months
    mean_X <- -0.9395
    sd_X <- 1.1763
    X <- rnorm(n, mean_X, sd_X)
  }
  
  # Shigella infection (Y_inf) in vaccinated and unvaccinated
  
  # get intercept baseline log-odds 
  bl_log_odds_shigella <- qlogis(incidence_shigella)
  
  prob_Y0is1 <- plogis(bl_log_odds_shigella + lfaz_bl_effect_shig*(X - mean_X)^2)
  Y_inf_V0 <- rbinom(n, 1, prob_Y0is1)
  prob_Y1is1 <- rep(0, n)
  prob_Y1is1[Y_inf_V0 == 1] <- 1 - VE_mild
  Y_inf_V1 <- rbinom(n, 1, prob_Y1is1)
  
  # Shigella severity (Y_sev) in vaccinated and unvaccinated
  bl_log_odds_severe <- qlogis(incidence_severe_shigella/incidence_shigella)
  
  prob_Y_severe_V0 <- ifelse(Y_inf_V0 == 1,
                             plogis(bl_log_odds_severe + lfaz_bl_effect_severe_shig*(X - mean_X)^2),
                             0)
  Y_sev_V0 <- rbinom(n, 1, prob_Y_severe_V0)
  prob_Y_severe_V1 <- rep(0, n)
  prob_Y_severe_V1[Y_sev_V0 == 1 & Y_inf_V1 == 1] <- (1 - VE_severe) / (1 - VE_mild)
  Y_sev_V1 <- rbinom(n, 1, prob_Y_severe_V1)
  
  Y_inf <- rep(NA, n)
  Y_inf[V == 1] <- Y_inf_V1[V == 1]
  Y_inf[V == 0] <- Y_inf_V0[V == 0]
  
  Y_severe <- rep(NA, n)
  Y_severe[V == 1] <- Y_sev_V1[V == 1]
  Y_severe[V == 0] <- Y_sev_V0[V == 0]
  
  # Shigella infection timing - (Y_time)
  Y_time <- sample(seq(1, 12), n, replace = TRUE)
  Y_time[Y_inf == 0] <- 0
  
  infection_adjustment <- shig_coef_mild*I(Y_inf == 1 & Y_severe == 0) + shig_coef_severe*I(Y_inf == 1 & Y_severe == 1)
  
  # Allow or prohibit catch up growth-- smaller effect size if allowing for catch-up
  if(catch_up == FALSE){
    
    # No catch-up growth
    
    # If sd_growth = "baseline", use standard deviation month 24 ~ month 12 MAL-ED Bangladesh (0.4389)
    # Else If sd_growth = "prev-1", use models fit with previous 1 timepoint
    # Else If sd_growth = "prev-2", use models fit with previous 2 timepoints
    # Else assume sd_growth is numeric sd to be used with no intermediate measurements 
    if(sd_growth == "baseline"){
      G <- X + infection_adjustment + rnorm(n, mean = 0, sd = 0.4389)
    } else if (sd_growth == "prev-1") {
      # Growth trajectory from MAL-ED Bangladesh 
      
      # Previous time point
      X_vec <- vector(mode = "list", length = 12) 
      X_vec[[1]]  <- -0.13083351 + 0.9690976*X           + rnorm(n, mean = 0, sd = 0.2120328)
      X_vec[[2]]  <- -0.13359784 + 0.9658743*X_vec[[1]]  + rnorm(n, mean = 0, sd = 0.2096142)
      X_vec[[3]]  <- -0.02840903 + 0.9964731*X_vec[[2]]  + rnorm(n, mean = 0, sd = 0.2143073)
      X_vec[[4]]  <- -0.12824842 + 0.9658842*X_vec[[3]]  + rnorm(n, mean = 0, sd = 0.1941224)
      X_vec[[5]]  <- -0.10182843 + 0.9745330*X_vec[[4]]  + rnorm(n, mean = 0, sd = 0.1732090)
      X_vec[[6]]  <- -0.07077051 + 0.9734801*X_vec[[5]]  + rnorm(n, mean = 0, sd = 0.1712597)
      X_vec[[7]]  <- -0.04753367 + 0.9912921*X_vec[[6]]  + rnorm(n, mean = 0, sd = 0.1616939)
      X_vec[[8]]  <- -0.04551516 + 0.9821661*X_vec[[7]]  + rnorm(n, mean = 0, sd = 0.1704176)
      X_vec[[9]]  <- -0.04299716 + 0.9801243*X_vec[[8]]  + rnorm(n, mean = 0, sd = 0.1597872)
      X_vec[[10]] <- -0.03899147 + 0.9936665*X_vec[[9]]  + rnorm(n, mean = 0, sd = 0.1516835)
      X_vec[[11]] <- -0.04331684 + 0.9741949*X_vec[[10]] + rnorm(n, mean = 0, sd = 0.1652360)
      X_vec[[12]] <- -0.08162686 + 0.9565905*X_vec[[11]] + rnorm(n, mean = 0, sd = 0.1991225)
      
      X_df <- data.frame(X_vec)
      colnames(X_df) <- c("X_1", "X_2", "X_3", "X_4", "X_5", "X_6",
                          "X_7", "X_8", "X_9", "X_10", "X_11", "X_12")
      X_df <- cbind(X_df, Y_time, infection_adjustment)
      
      # Subtract infection adjustment from all time points beyond infection
      adjust_for_inf <- function(row){
        my_Y_time <- row['Y_time']
        if(my_Y_time > 0){
          row[(my_Y_time:12)] <- row[(my_Y_time:12)] + row['infection_adjustment']
        }
        return(row)
      }
      
      X_df <- data.frame(t(apply(X_df, MARGIN = 1, adjust_for_inf)))
      
      G <- X_df$X_12
      
    }else if(sd_growth == "prev-2"){
      X_vec <- vector(mode = "list", length = 12) 
      
      # Previous 2 timepoints
      X_vec[[1]]  <- -0.13083351 + 0.9690976*X          + rnorm(n, mean = 0, sd = 0.2120328)
      X_vec[[2]]  <- -0.13663765 + 0.10708213*X         + 0.8610326*X_vec[[1]]   + rnorm(n, mean = 0, sd = 0.2086736)
      X_vec[[3]]  <- -0.03902150 + 0.01465246*X_vec[[1]]  + 0.9779293*X_vec[[2]]   + rnorm(n, mean = 0, sd = 0.2141136)
      X_vec[[4]]  <- -0.11403457 + 0.20063834*X_vec[[2]]  + 0.7757495*X_vec[[3]]   + rnorm(n, mean = 0, sd = 0.1899936)
      X_vec[[5]]  <- -0.10566722 + 0.04466520*X_vec[[3]]  + 0.9287577*X_vec[[4]]   + rnorm(n, mean = 0, sd = 0.1740689)
      X_vec[[6]]  <- -0.07251206 + 0.04606467*X_vec[[4]]  + 0.9277831*X_vec[[5]]   + rnorm(n, mean = 0, sd = 0.1710707)
      X_vec[[7]]  <- -0.04799386 + 0.06796371*X_vec[[5]]  + 0.9240339*X_vec[[6]]   + rnorm(n, mean = 0, sd = 0.1619633)
      X_vec[[8]]  <- -0.04294224 + 0.13528495*X_vec[[6]]  + 0.8505967*X_vec[[7]]   + rnorm(n, mean = 0, sd = 0.1698145)
      X_vec[[9]]  <- -0.03597587 + 0.14149890*X_vec[[7]]  + 0.8421613*X_vec[[8]]   + rnorm(n, mean = 0, sd = 0.1583651)
      X_vec[[10]] <- -0.03382477 + 0.18068365*X_vec[[8]]  + 0.8170421*X_vec[[9]]   + rnorm(n, mean = 0, sd = 0.1455196)
      X_vec[[11]] <- -0.04020125 + 0.17402791*X_vec[[9]]  + 0.8045709*X_vec[[10]]  + rnorm(n, mean = 0, sd = 0.1630071)
      X_vec[[12]] <- -0.07027040 + 0.24234964*X_vec[[10]] + 0.7192027*X_vec[[11]] + rnorm(n, mean = 0, sd = 0.1905827)
      
      X_df <- data.frame(X_vec)
      colnames(X_df) <- c("X_1", "X_2", "X_3", "X_4", "X_5", "X_6",
                          "X_7", "X_8", "X_9", "X_10", "X_11", "X_12")
      X_df <- cbind(X_df, Y_time, infection_adjustment)
      
      # Subtract infection adjustment from all time points beyond infection
      adjust_for_inf <- function(row){
        my_Y_time <- row['Y_time']
        if(my_Y_time > 0){
          row[(my_Y_time:12)] <- row[(my_Y_time:12)] + row['infection_adjustment']
        }
        return(row)
      }
      
      X_df <- data.frame(t(apply(X_df, MARGIN = 1, adjust_for_inf)))
      
      G <- X_df$X_12
      
    } else{
      G <- X + infection_adjustment + rnorm(n, mean = 0, sd = sd_growth)
    }
  } else{
    # Allow for catch-up growth
    
    if(sd_growth == "baseline"){
      G <- X + shig_coef_mild*I(Y_inf == 1 & Y_severe == 0) + shig_coef_severe*I(Y_inf == 1 & Y_severe == 1) + rnorm(n, mean = 0, sd = 0.4389)
    } else if (sd_growth == "prev-1") {
      # Growth trajectory from MAL-ED Bangladesh 
      
      # Previous time point
      X_1 <- -0.13083351 + 0.9690976*X + I(Y_time == 1)*infection_adjustment + rnorm(n, mean = 0, sd = 0.2120328)
      X_2 <- -0.13359784 + 0.9658743*X_1 + I(Y_time == 2)*infection_adjustment + rnorm(n, mean = 0, sd = 0.2096142)
      X_3 <- -0.02840903 + 0.9964731*X_2 + I(Y_time == 3)*infection_adjustment + rnorm(n, mean = 0, sd = 0.2143073)
      X_4 <- -0.12824842 + 0.9658842*X_3 + I(Y_time == 4)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1941224)
      X_5 <- -0.10182843 + 0.9745330*X_4 + I(Y_time == 5)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1732090)
      X_6 <- -0.07077051 + 0.9734801*X_5 + I(Y_time == 6)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1712597)
      X_7 <- -0.04753367 + 0.9912921*X_6 + I(Y_time == 7)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1616939)
      X_8 <- -0.04551516 + 0.9821661*X_7 + I(Y_time == 8)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1704176)
      X_9 <- -0.04299716 + 0.9801243*X_8 + I(Y_time == 9)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1597872)
      X_10 <- -0.03899147 + 0.9936665*X_9 + I(Y_time == 10)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1516835)
      X_11 <- -0.04331684 + 0.9741949*X_10 + I(Y_time == 11)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1652360)
      
      G <- -0.08162686 + 0.9565905*X_11 + I(Y_time == 12)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1991225)
      
    }else if(sd_growth == "prev-2"){
      # Previous 2 timepoints
      X_1 <- -0.13083351 + 0.9690976*X + I(Y_time == 1)*infection_adjustment + rnorm(n, mean = 0, sd = 0.2120328)
      X_2 <- -0.13663765 + 0.10708213*X + 0.8610326*X_1 + I(Y_time == 2)*infection_adjustment + rnorm(n, mean = 0, sd = 0.2086736)
      X_3 <- -0.03902150 + 0.01465246*X_1 + 0.9779293*X_2 + I(Y_time == 3)*infection_adjustment + rnorm(n, mean = 0, sd = 0.2141136)
      X_4 <- -0.11403457 + 0.20063834*X_2 + 0.7757495*X_3 + I(Y_time == 4)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1899936)
      X_5 <- -0.10566722  + 0.04466520*X_3 + 0.9287577*X_4 + I(Y_time == 5)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1740689)
      X_6 <- -0.07251206 + 0.04606467*X_4 + 0.9277831*X_5 + I(Y_time == 6)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1710707)
      X_7 <- -0.04799386 + 0.06796371*X_5 + 0.9240339*X_6 + I(Y_time == 7)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1619633)
      X_8 <- -0.04294224 + 0.13528495*X_6 + 0.8505967*X_7 + I(Y_time == 8)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1698145)
      X_9 <- -0.03597587 + 0.14149890*X_7 + 0.8421613*X_8 + I(Y_time == 9)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1583651)
      X_10 <- -0.03382477 + 0.18068365*X_8 + 0.8170421*X_9 + I(Y_time == 10)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1455196)
      X_11 <- -0.04020125 + 0.17402791*X_9 + 0.8045709*X_10 + I(Y_time == 11)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1630071)
      
      G <- -0.07027040 + 0.24234964*X_10 + 0.7192027*X_11 + I(Y_time == 12)*infection_adjustment + rnorm(n, mean = 0, sd = 0.1905827)
    } else{
      G <- X + infection_adjustment + rnorm(n, mean = 0, sd = sd_growth)
    }
    
  }
  
  final_df <- data.frame(id = 1:n,
                         V = V,
                         X = X,
                         Y_inf = Y_inf,
                         Y_sev = Y_severe,
                         G = G)
  if(zhifei == TRUE){
    final_df <- cbind(final_df, age)
  }
  
  return(final_df)
}