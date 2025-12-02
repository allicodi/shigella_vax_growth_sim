here::i_am("misc/high_inc_site.R")

source("R/parameter_generation_fns.R")

trials <- c("6mo", "12mo")
sites <- c("Bangladesh", "Kenya", "Malawi", "Mali",
           "Pakistan", "Peru", "The Gambia")

combos <- expand.grid(trial = trials,
                      site = sites)

results <- data.frame()

for(i in 1:nrow(combos)){
  results <- rbind(results, get_incidence(dose_schedule = combos$trial[i],enroll_site = combos$site[i]))
}

results <- cbind(combos, results)

# Peru seems to be the highest for both 6mo and 12mo trials

