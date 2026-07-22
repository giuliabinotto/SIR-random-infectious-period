###########################################################
#### Stochastic epidemics models with random duration #####
######### Influenza: Hong Kong University study ###########
####### Simulation of the HKU study using our model #######
###########################################################
################# X. Bardina, G. Binotto ##################
###########################################################


### Libraries ###
library(dplyr)
library(here)
  here::i_am("influenza/Influ_p.R")


### Sources ###
source(here::here("R", "Functions.R"))
source(here::here("influenza", "Influ_q.R"))
source(here::here("influenza", "Influ_SAR.R"))
source(here::here("influenza", "Influ_p.R"))

  
##### Simulations with our model #####

secondary_infections <- function(n.trials) {
    
    # Secondary infections by household size
    results_list <- lapply(hh_sizes_influ, function(n_hh) {
        
        invisible(capture.output(
            results <- simulations(
                n.trials = n.trials, n = n_hh, Y0 = 1, p = 0.03, r = 1, 
                d = 9, q.distr = "HKstudy", q.param = NULL, lambda = 1, 
                Time = 9, N.distr = "constant", N.param1 = n_hh-1, 
                N.param2 = NULL
            )
        ))
        
        mean_infected <- mean(results$Ytot - 1)
        
        data.frame(
            hh_size = n_hh,
            mean_infected = mean_infected,
            proportion_infected = mean_infected / (n_hh-1)
        )
    })
    
    df <- do.call(rbind, results_list)
    
    # Simulated SAR
    probs <- freq_hh_size_noindex / sum(freq_hh_size_noindex)
    
    SAR_simul <- sum(df$proportion_infected * probs)
    
    list(
        secondary_infections = df,
        SAR_simul = SAR_simul
    )
}

set.seed(123)
secondary_infections(n.trials = 1000)


###################### EXTRA ###########################
##### Best seed #####

seeds <- 1:100

SAR_results <- data.frame(
    seed = seeds,
    SAR = sapply(seeds, function(s) {
        set.seed(s)
        secondary_infections(n.trials = 1000)$SAR_simul
    })
)

summary(SAR_results$SAR)
hist(SAR_results$SAR)

# Seed with minimum SAR
SAR_results[which.min(SAR_results$SAR), ]   # Seed 49

SAR_results %>%
    filter(SAR < 0.142)
