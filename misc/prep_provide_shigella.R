
here::i_am("misc/prep_provide_shigella.R")

library(tidyverse)

load("misc/provide_data/PROVIDE diarrhea.RData")
load("misc/provide_data/provide_f10_r1-3bk7-.RData")

# Get cases
shigella_diarrhea <- diarrhea4 %>%
  filter(shigella_eiec_afe > 0.5) %>%
  mutate(case = 1) %>%
  select(sid, 
         case,
         dedt, 
         daysbirth,
         adenovirus_40_41,
         aeromonas, 
         astrovirus, 
         c_jejuni_coli, 
         cryptosporidium, 
         cyclospora, 
         e_histolytica, 
         isospora, 
         norovirus_gii, 
         rotavirus, 
         salmonella, 
         sapovirus, 
         shigella_eiec, 
         ST_ETEC,
         TEPEC, 
         v_cholerae, 
         ETECLT, 
         EAEC)

# Get severity information for cases

# GEMS MSD definition = 
# sunken eyes  --- in dehydration assessment (dehyd)
# loss of skin tugor (slow return) --- in dehydration assessment
# IV rehydration admin --- rehyd
# dysentery - no 
# hospitalized with diarrhea or dysentery -- we have SAEs??

# MAL-ED variables
# "dysentery" -- no
# "fever", -- yes
# "fever_days", -- no, but we have epidays
# "dehyd", -- yes
# "lsstools", -- episevrt (number of stool / 24hr)
# "daysvomit" -- yes

severity_df <- mgmt_bv_ddep_f10 %>%
  filter(sid %in% unique(shigella_diarrhea$sid)) %>%
  mutate(MSD = if_else(dehyd > 1, 1, 0)) %>%    # dehyd >= 2 means skin tenting, sunken eyes, and/or dry mucus membrane
  mutate(MSD = if_else(rehyd == 3, 1, MSD)) %>% # rehyd == 3 means IV rehydration 
  mutate(high_fever = if_else(feversc >= 2, 1, 0),
         dehyd_bin = if_else(dehyd == -9, NA,
                             if_else(dehyd > 1, 1, 0)),
         epidate = as.Date(epidate)) %>%
  select(sid,epidate, daysbirth, MSD, high_fever, dehyd_bin, vomtdays, episevrt, epidays) %>%
  rename('lsstools' = episevrt,
         'daysvomit' = vomtdays, 
         'epi_duration' = epidays)

shigella_diarrhea <- left_join(shigella_diarrhea, severity_df, by = c("sid","daysbirth"))

# df of all hospitalizations for case SIDs
hosp_df <- mgmt_rts_rotatrial_sae_f10 %>%
  filter(saevnt == 3) %>% # hospitalization
  filter(sid %in% unique(shigella_diarrhea$sid))

# if there is a hospitalization within a week of the shigella episode, attribute to that
hosp_match <- lapply(1:nrow(shigella_diarrhea), function(i){
  row <- shigella_diarrhea[i,,drop = FALSE]
  
  date_min <- as.Date(row$epidate) - days(7)
  date_max <- as.Date(row$epidate) + days(7)
  shig_sid <- row$sid
  
  hosp_sid <- hosp_df %>%
    filter(sid == shig_sid,
           saespdt >= date_min & saespdt <= date_max)
  
  if(nrow(hosp_sid > 0)) return(1)
  else return(0)
}) 

# incorporate hosp into MSD definition
shigella_diarrhea$hosp <- as.numeric(hosp_match)
shigella_diarrhea$MSD <- ifelse(shigella_diarrhea$hosp, 1, shigella_diarrhea$MSD)

# make smaller to cases only
haz_long <- mgmt_bv_danth_f10 %>%
  filter(sid %in% unique(shigella_diarrhea$sid)) %>%
  mutate(haz = if_else(haz == -9, NA, haz),
         ageday = if_else(ageday == -9, NA, ageday),
         agemonth = ageday / 30.44,
         vstwk_num = str_remove(vstwk, "^Week\\s+"))

targets_months <- c(3, 6, 9, 12)

shigella_diarrhea <- lapply(1:nrow(shigella_diarrhea), function(i) {
  row <- shigella_diarrhea[i, , drop = FALSE]
  agedays_inf <- row$daysbirth
  
  # Baseline measurement (closest pre-episode)
  bl_row <- haz_long %>%
    filter(sid == row$sid,
           ageday != -9,
           ageday <= agedays_inf) %>%
    slice_max(ageday, n = 1, with_ties = FALSE) %>%
    select(sid, dob, ageday, agemonth, gender, haz) %>%
    rename(bl_haz = haz,
           ageday_bl = ageday,
           agemonth_bl = agemonth)
  
  # Follow-up measurements
  followup_row <- lapply(targets_months, function(m) {
    target_day <- agedays_inf + m * 30.44
    
    match_row <- haz_long %>%
      filter(sid == row$sid,
             ageday != -9,
             ageday > agedays_inf,
             abs(ageday - target_day) <= 30) %>%
      slice_min(abs(ageday - target_day), n = 1, with_ties = FALSE)
    
    if (nrow(match_row) == 0) {
      tibble(sid = row$sid,
             haz = NA_real_,
             ageday = NA_real_,
             agemonth = NA_real_,
             target_month = m)
    } else {
      match_row %>%
        mutate(target_month = m) %>%
        select(sid, haz, ageday, agemonth, target_month)
    }
  }) %>%
    bind_rows() %>%
    pivot_wider(
      id_cols = sid,
      names_from = target_month,
      values_from = c(haz, ageday, agemonth),
      names_glue = "{.value}_fu_{target_month}m"
    ) %>%
    select(-sid)  # drop sid to avoid duplicates
  
  # Episode info
  epi_row <- row %>%
    rename(agedays_episode = daysbirth) %>%
    mutate(agemonths_episode = agedays_episode / 30.44) %>%
    select(-sid)  # drop sid to avoid duplicates
  
  # Combine
  cbind(bl_row, epi_row, followup_row)
}) %>% bind_rows()

shigella_diarrhea$case_id <- 1:nrow(shigella_diarrhea)

control_data <- lapply(1:nrow(shigella_diarrhea), function(i) {
  
  row <- shigella_diarrhea[i, , drop = FALSE]
  
  # Matching window for DOB ±15 days
  date_min_dob <- as.Date(row$dob) - days(15)
  date_max_dob <- as.Date(row$dob) + days(15)
  
  # Find potential controls matching gender, enrollment visit, DOB window, exclude case itself
  potential_controls <- mgmt_bv_danth_f10 %>%
    filter(gender == row$gender,
           vstwk == "Enrollment") %>%
    mutate(dob = as.Date(dob)) %>%
    filter(dob >= date_min_dob & dob <= date_max_dob,
           sid != row$sid)
  
  # Remove controls with diarrhea in week before case episode
  date_min_epi <- as.Date(row$epidate) - days(7)
  diarrhea_wk_before <- diarrhea4 %>%
    filter(sid %in% potential_controls$sid,
           dedt >= date_min_epi & dedt <= row$dedt)
  
  controls <- potential_controls %>%
    filter(!(sid %in% diarrhea_wk_before$sid))
  
  if (nrow(controls) == 0) return(NULL)
  
  # For each control, get baseline and follow-up measurements aligned to case time points
  controls_long <- lapply(1:nrow(controls), function(j) {
    ctrl <- controls[j, , drop = FALSE]
    
    ctrl_haz <- mgmt_bv_danth_f10 %>%
      filter(sid == ctrl$sid) %>%
      mutate(haz = if_else(haz == -9, NA_real_, haz),
             ageday = if_else(ageday == -9, NA_real_, ageday),
             agemonth = ageday / 30.44)
    
    # Baseline closest to case baseline HAZ age ±30 days (no before restriction)
    bl_row <- ctrl_haz %>%
      filter(!is.na(ageday),
             abs(ageday - row$ageday_bl) <= 30) %>%
      slice_min(abs(ageday - row$ageday_bl), n = 1, with_ties = FALSE) %>%
      mutate(dob = as.Date(dob)) %>%  
      select(sid, dob, ageday, agemonth, gender, haz) %>%
      rename(bl_haz = haz,
             ageday_bl = ageday,
             agemonth_bl = agemonth)
    
    # If no baseline found, fill with NA
    if (nrow(bl_row) == 0) {
      bl_row <- tibble(
        sid = ctrl$sid,
        dob = as.Date(NA),
        ageday_bl = NA_real_,
        agemonth_bl = NA_real_,
        gender = ctrl$gender,
        bl_haz = NA_real_
      )
    }
    
    # Follow-up target ages from case
    target_months <- c(3, 6, 9, 12)
    target_days <- bl_row$ageday_bl + target_months * 30.44
    names(target_days) <- c("3m", "6m", "9m", "12m")
    
    followup_row <- lapply(seq_along(target_days), function(idx) {
      target_day <- target_days[idx]
      target_name <- names(target_days)[idx]
      
      if (is.na(target_day)) {
        tibble(
          sid = ctrl$sid,
          haz = NA_real_,
          ageday = NA_real_,
          agemonth = NA_real_,
          target_month = target_name
        )
      } else {
        match_row <- ctrl_haz %>%
          filter(!is.na(ageday),
                 abs(ageday - target_day) <= 30) %>%
          slice_min(abs(ageday - target_day), n = 1, with_ties = FALSE)
        
        if (nrow(match_row) == 0) {
          tibble(
            sid = ctrl$sid,
            haz = NA_real_,
            ageday = NA_real_,
            agemonth = NA_real_,
            target_month = target_name
          )
        } else {
          match_row %>%
            mutate(target_month = target_name) %>%
            select(sid, haz, ageday, agemonth, target_month)
        }
      }
    }) %>%
      bind_rows() %>%
      pivot_wider(
        id_cols = sid,
        names_from = target_month,
        values_from = c(haz, ageday, agemonth),
        names_glue = "{.value}_fu_{target_month}"
      ) %>%
      select(-sid)
    
    # Episode info from case plus control_sid for linking
    epi_row <- row %>%
      mutate(control_sid = ctrl$sid) %>%
      select(case_id, agedays_episode, agemonths_episode)
    
    cbind(control_sid = ctrl$sid, bl_row, epi_row, followup_row)
  }) %>% bind_rows()
  
  controls_long
}) %>% bind_rows()

# NEED TO COMBINE CASE & CONTROL, THEN GET BASELINE COVARIATES (ses, edu, etc)

# maled covariates --> provide covariates
# "site" -- NA?
# "sex" -- sex
# "agemonths" -- age at infection? / control matched
# "baseline_haz" -- get HAZ at 6months
# "mated_bin" -- yes
# "drinkimp" -- yes
# "sanitimp" -- yes
# wami_quintile -- make SES quintile

# "adenovirus_40_41_new",      - yes            
# "aeromonas_new",             - yes         
# "astrovirus_new",            - yes       
# "campylobacter_pan_new",     - c jejuni
# "cryptosporidium_new",       - yes           
# "cyclospora_new",            - yes  
# "e_histolytica_new",         - yes          
# "isospora_new",              - yes    
# "norovirus_new",             - yes         
# "rotavirus_new",             - duh       
# "salmonella_new",            - yes        
# "sapovirus_new",             - yes
# "shigella_new",              - yes
# "st_etec_new",               - yes         
# "tEPEC_new",                 - yes           
# "v_cholerae_new",            - yes       
# "ETEC_new",                  - yes       
# "e_bieneusi_new",            - no       
# "eaec_new",                  - yes

# "dysentery",
# "fever", -- yes
# "fever_days", 
# "dehyd", -- yes
# "lsstools", -- episevrt (number of stool / 24hr)
# "daysvomit" -- yes

# GEMS MSD definition = 
# sunken eyes  --- in dehydration assessment
# loss of skin tugor (slow return) --- in dehydration assessment
# IV rehydration admin --- rehyd
# dysentery - no 
# hospitalized with diarrhea or dysentery -- we have SAEs?? 

# get cases
# make controls matching on same criteria as MAL-ED
# Get growth outcomes X time later 

# growth measured at 

# approx every three months with some flexibility 

# enrollment
# week 6 
# week 10 
# week 12 
# week 14 
# week 17 
# week 18 
# week 24 
# week 39 
# week 40 
# week 52 
# week 65
# week 78
# week 91
# week 104


# MAL-ED Control matching

# matched_controls <- lapply(1:nrow(case_data), function(i, case_data, all_controls, tac_data){
#   row <- case_data[i,]
#   
#   # get age range to match depending on case age
#   # age range = 0-11 months (1-364 days) --> +- 2mo (30.44*2 mo= 61 days)
#   # age range = 12+ months (365 days +) --> +- 4mo (30.44*4 mo = 122 days)
#   if(row$agedays < 365){
#     min_age <- max(0, row$agedays - 61)
#     max_age <- min(row$agedays + 61, 364)
#   } else{
#     min_age <- max(365, row$agedays - 122)
#     max_age <- row$agedays + 122
#   }
#   
#   # match sex, site, time, age, not their own control
#   matching_controls <- all_controls %>%
#     filter(cafsex == row$cafsex) %>%
#     filter(country_id == row$country_id) %>%
#     filter(date < (row$date + days(15)) & date > row$date - days(15)) %>%
#     filter(agedays >= min_age & agedays <= max_age) %>%
#     filter(pid != row$pid) %>%
#     group_by(pid) %>%
#     slice_max(order_by = date, n = 1) %>% # Keep the latest sample per individual
#     ungroup()
#   
#   # for each matching control, make sure no diarrhea 7 days prior
#   control_eligible <- rep(TRUE, nrow(matching_controls))
#   for(j in 1:nrow(matching_controls)){
#     control_row <- matching_controls[j,]
#     tac_data_match <- tac_data %>%
#       filter(pid == control_row$pid) %>%
#       filter(date <= control_row$date & date > control_row$date - days(7))
#     
#     if(any(tac_data_match$stooltype == "D1")){
#       control_eligible[j] <- FALSE
#     } 
#   }
#   
#   # eliminate ineligible controls
#   matching_controls <- matching_controls[control_eligible,]
#   
#   matching_controls$case_pid <- row$pid
#   matching_controls$case_sid <- row$sid
#   
#   return(matching_controls)
#   
# }, case_data = case_data, all_controls = all_controls, tac_data = tac_data)

# get age range to match depending on case age
# age range = 0-11 months (1-364 days) --> +- 2mo (30.44*2 mo= 61 days)
# age range = 12+ months (365 days +) --> +- 4mo (30.44*4 mo = 122 days)