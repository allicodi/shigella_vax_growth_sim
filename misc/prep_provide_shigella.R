
here::i_am("misc/prep_provide_shigella.R")

library(tidyverse)

load("misc/provide_data/PROVIDE diarrhea.RData")
load("misc/provide_data/provide_f10_r1-3bk7-.RData")

shigella_diarrhea <- diarrhea4 %>%
  filter(shigella_eiec_afe > 0.5) 

anthro_data <- mgmt_bv_danth_f10 

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


# get age range to match depending on case age
# age range = 0-11 months (1-364 days) --> +- 2mo (30.44*2 mo= 61 days)
# age range = 12+ months (365 days +) --> +- 4mo (30.44*4 mo = 122 days)