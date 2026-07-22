###########################################################
#### Stochastic epidemics models with random duration #####
######### Influenza: Hong Kong University study ###########
######## Estimation of the contagion probability p ########
###########################################################
################# X. Bardina, G. Binotto ##################
###########################################################

# This script estimates the influenza transmission probability p using the infectious 
# period parameters (r, d and q) and the secondary attack rate (SAR) estimated 
# from the combined HKU 2008 and 2009 datasets.
#
# Source this file to make the p available, e.g.:
#   source("influenza/Influ_p.R")
#
# This script depends on (sourced below):
#   Influ_q.R   -> r_influ, d_influ, q_influ  (parameters of the infectious period)
#   Influ_SAR.R -> SAR_influ                  (secondary attack rate)


### Libraries ###
library(here)
  here::i_am("influenza/Influ_p.R")


### Sources ###
source(here::here("influenza", "Influ_q.R"))
source(here::here("influenza", "Influ_SAR.R"))

  
### Parameters ###
# Infectious period
r <- r_influ
d <- d_influ
q <- q_influ
D.prob <- q-c(q[2:length(q)],0)   # Probability function of D
T.values <- r:(r+d)   # Support (values) of T=r+D

# SAR
SAR <- SAR_influ

# R0 (from literature)
R0 <- 1.28


### Function: estimate the contagion probability p ###

#' @param D.prop      Probability function of D
#' @param T.values    Support of T=r+D
#' @param SAR         Estimated SAR for influenza
#' @param R0          R0 for influenza (from literature)
#'
#' @return A list with:
#'   p    - contagion probability
#'   N    - (constant) number of daily contacts

estimate_p <- function(D.prob, T.values, SAR, R0){
  
  ### Expected values of T and T^2 ###
  E <- sum(D.prob*T.values)
  E2 <- sum(D.prob*T.values^2)
  
  ### Estimation of p ###
  # Taylor's expansion of second order: (1-p)^T=1-pT+p^2*T*(T-1)/2+o(p^3)
  # SAR=p*E(T)-p^1*(E(T^2)-E(T))/2 (aproximation)
  p <- Re(polyroot(c(SAR, -E, (E2-E)/2))[1])
  
  ### Estimation of N (constant number of daily contacts) ###
  N <- R0/SAR
  
  ### Output ###
  result <- data.frame(p = round(p, 4), N = round(N, 4))
  
  return(result)
}

res <- estimate_p(D.prob, T.values, SAR, R0)


### Output: probability p ###
# This is the object other scripts depend on (influ_HKStudy.R).
p_influ <- res$p
