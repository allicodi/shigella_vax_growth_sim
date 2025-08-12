
here::i_am("misc/prep_provide_shigella.R")

library(tidyverse)

load("misc/provide_data/PROVIDE diarrhea.RData")
load("misc/provide_data/provide_f10_r1-3bk7-.RData")

# Get cases
shigella_diarrhea <- diarrhea4 %>%
  filter(shigella_eiec_afe > 0.5) %>%
  mutate(case = 1) %>%
  mutate(ruuska_sev = if_else(ruuska >= 11, 1, 0)) %>%
  select(sid, 
         case,
         ruuska_sev,
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
  mutate(GEMS_MSD = if_else(dehyd > 1, 1, 0)) %>%    # dehyd >= 2 means skin tenting, sunken eyes, and/or dry mucus membrane
  mutate(GEMS_MSD = if_else(rehyd == 3, 1, GEMS_MSD)) %>% # rehyd == 3 means IV rehydration 
  mutate(high_fever = if_else(feversc >= 2, 1, 0),
         dehyd_bin = if_else(dehyd == -9, NA,
                             if_else(dehyd > 1, 1, 0)),
         epidate = as.Date(epidate)) %>%
  select(sid,epidate, daysbirth, GEMS_MSD, high_fever, dehyd_bin, vomtdays, episevrt, epidays) %>%
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

# incorporate hosp into GEMS_MSD definition
shigella_diarrhea$hosp <- as.numeric(hosp_match)
shigella_diarrhea$GEMS_MSD <- ifelse(shigella_diarrhea$hosp, 1, shigella_diarrhea$GEMS_MSD)

# Get antibiotics info for cases

# Antibiotic use

abx_key <- c(
  'azithromycin',
  'erythromycin',
  'metronidazole',
  'ciprofloxacin',
  'nalidixic acid',
  'pivmecillinam',
  'cephradine',
  'co-trimoxazole',
  'amoxycillin',
  'cefixime',
  'ceftazidime',
  'cefuroxime',
  'ceftriaxone',
  'fluclocaxillin',
  'other'
)

abx_who <- c('azithromycin','ciprofloxacin','ceftriaxone','pivmecillinam')
abx_maybe <- c('erythromycin','nalidixic acid','amoxycillin','cefixime','ceftazidime', 'cefuroxime','co-trimoxazole')
abx_no <- c('metronidazole','cephradine','fluclocaxillin','other')

bv_dep_diarrheal_episode_f10$dedt <- as.Date(bv_dep_diarrheal_episode_f10$dedt)

abx_df <- shigella_diarrhea %>%
  left_join(bv_dep_diarrheal_episode_f10, by = c("sid", "dedt")) %>%
  mutate(antb = factor(antb, levels = 1:3, labels = c("None", "Pre-visit", "At-visit")),
         antb1 = if_else(antb1 == 99, NA, antb1),
         antb2 = if_else(antb2 == 99, NA, antb2), 
         antb1 = factor(antb1, levels = 1:15, labels = abx_key),
         antb2 = factor(antb2, levels = 1:15, labels = abx_key)) %>%
  group_by(sid, dedt) %>%
  summarise(all_abx = if_else(any(antb1 %in% abx_who) | any(antb2 %in% abx_who), "Guideline recommended",
                               if_else(any(antb1 %in% abx_maybe) | any(antb2 %in% abx_maybe), "Possibly effective", "No or ineffective")))


shigella_diarrhea <- left_join(shigella_diarrhea, abx_df, by = c("sid", "dedt") )
shigella_diarrhea$all_abx <- factor(shigella_diarrhea$all_abx)

# make smaller to cases only
haz_long <- mgmt_bv_danth_f10 %>%
  filter(sid %in% unique(shigella_diarrhea$sid)) %>%
  mutate(haz = if_else(haz == -9, NA, haz),
         ageday = if_else(ageday == -9, NA, ageday),
         agemonth = ageday / 30.44,
         vstwk_num = str_remove(vstwk, "^Week\\s+"))

targets_months <- c(1, 3, 6, 9, 12)

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
             abs(ageday - target_day) <= 45) %>%
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
      names_glue = "{.value}_{target_month}m"
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

# Get rid of controls that only have 2 or less observations (early dropout?)
mgmt_bv_danth_f10_filtered <- mgmt_bv_danth_f10 %>%
  filter(haz != -9) %>%
  group_by(sid) %>%
  filter(n() >= 3 ) %>%
  ungroup()

control_data <- lapply(1:nrow(shigella_diarrhea), function(i) {
  
  row <- shigella_diarrhea[i, , drop = FALSE]
  
  # Matching window for DOB ±15 days
  date_min_dob <- as.Date(row$dob) - days(15)
  date_max_dob <- as.Date(row$dob) + days(15)
  
  # Find potential controls matching gender, enrollment visit, DOB window, exclude case itself
  potential_controls <- mgmt_bv_danth_f10_filtered %>%
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
    
    ctrl_haz <- mgmt_bv_danth_f10_filtered %>%
      filter(sid == ctrl$sid) %>%
      mutate(haz = if_else(haz == -9, NA_real_, haz),
             ageday = if_else(ageday == -9, NA_real_, ageday),
             agemonth = ageday / 30.44)
    
    # Baseline closest to case baseline HAZ age ±45 days (no before restriction)
    bl_row <- ctrl_haz %>%
      filter(!is.na(ageday),
             abs(ageday - row$ageday_bl) <= 45) %>%
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
    target_months <- c(1, 3, 6, 9, 12)
    target_days <- bl_row$ageday_bl + target_months * 30.44
    names(target_days) <- c("1m", "3m", "6m", "9m", "12m")
    
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
                 abs(ageday - target_day) <= 45) %>%
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
        names_glue = "{.value}_{target_month}"
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
# Prepare cases dataframe
# Prepare cases: rename sid to case_sid
cases_combined <- shigella_diarrhea %>%
  mutate(case = 1) %>%
  select(case_id, sid, case, everything()) %>%
  mutate(dob = as.Date(dob))

controls_combined <- control_data %>%
  mutate(case = 0) %>%
  select(case_id, sid, case, everything()) %>%
  mutate(dob = as.Date(dob))

combined_df <- bind_rows(
  cases_combined,
  controls_combined
) %>%
  arrange(case_id, desc(case))



# Baseline SES related covariates

# Number of siblings <5years
# Number of people sleeping in household
# Mother education
# Father education
# Total monthly income
# Electricity
# Cooking gas
# TV
# Drinking water source
# Toilet
# Food availability

cols_for_pca <- c("elec",  "gas",   "phon" , "almr"  ,"tabl",  "chair", "bench" ,"clock" ,"bed"   ,"radio","tv"  ,  "bcycl" ,"mcycl" ,"sewm"  ,"fan")

ses_pca <- prcomp(bv_ses_water_f10[complete.cases(bv_ses_water_f10[,cols_for_pca]), cols_for_pca], 
                  center = TRUE, scale = TRUE)

pca_score <- as.matrix(bv_ses_water_f10[, cols_for_pca]) %*% ses_pca$rotation[, 1]

bv_ses_water_f10$ses_score <- pca_score[, 1]  
bv_ses_water_f10 <- bv_ses_water_f10 %>%
  mutate(ses_quintile = ntile(ses_score, 5)) 

sub_ses <- bv_ses_water_f10 %>%
  select(sid,
         sb5y,
         pepl,
         medu,
         ses_quintile,
         watr,
         toil) %>% 
  rename('num_hh_lt_5' = sb5y,
         'num_hh_sleep' = pepl)

data <- left_join(combined_df, sub_ses, by = "sid")

# Transform factors 
data$gender <- factor(data$gender)
data$medu <- factor(data$medu, levels = 1:18, labels = c("No formal education",
                                                         "1st year",
                                                         "2nd year",
                                                         "3rd year",
                                                         "4th year",
                                                         "5th year",
                                                         "6th year",
                                                         "7th year",
                                                         "8th year",
                                                         "9th year",
                                                         "SSC_Dakhil passed",
                                                         "HSC_Fazil passed",
                                                         "Vocational_diploma_homeopathy_LMF_etc",
                                                         "Degree_Alim passed",
                                                         "Hons passed_3 or 4yr hons",
                                                         "Master_Kamil passed",
                                                         "MBBS_MD_FCPS_FRCP",
                                                         "BSC Engineer_MSc_PhD_etc"))
# Make edu vars with fewer categories
data <- data %>%
  mutate(
    medu_cat = case_when(
      medu %in% c("No formal education") ~ 1,
      medu %in% c("1st year", "2nd year", "3rd year", 
                  "4th year", "5th year") ~ 2,
      medu %in% c("6th year", "7th year") ~ 3,
      medu %in% c("8th year", "9th year", "SSC_Dakhil passed", "Degree_Alim passed", "Vocational_diploma_homeopathy_LMF_etc", 
                  "HSC_Fazil passed", "Hons passed_3 or 4yr hons", 
                  "Master_Kamil passed", "MBBS_MD_FCPS_FRCP", "BSC Engineer_MSc_PhD_etc") ~ 4,
      TRUE ~ NA_real_
    ), medu_cat = factor(medu_cat, levels = 1:4, labels = c("No formal education",
                                                           "Primary",
                                                           "Secondary",
                                                           "Higher")))

data$watr <- factor(data$watr, levels = 1:3, labels = c("Municipality supply_piped",
                                                        "Own arrangement by pump",
                                                        "Tube well"))

data$toil <- factor(data$toil, levels = 1:5, labels = c("Septic tank or toilet",
                                                        "Water sealed or slap latrine",
                                                        "Pit latrine",
                                                        "Open latrine",
                                                        "Hanging latrine"))

saveRDS(data, here::here("misc/provide_data/shigella_provide_case_control.Rds"))
