###########################################################
#### Stochastic epidemics models with random duration #####
######### Influenza: Hong Kong University study ###########
## Estimation of the duration probabilities q_i=P(D>=i) ###
###########################################################
################# X. Bardina, G. Binotto ##################
###########################################################

# This script defines the values of r and d and the vector q (q_i = P(D >= i) 
# with i = 0, ..., d), estimated from the combined HKU 2008 and 2009 datasets.
#
# Source this file to make the r, d and q available, e.g.:
#   source("influenza/Influ_q.R")
#
# Required data (NOT included in this repo, see README.md "Data"):
#   data/influenza/HKUstudy2008/
#   data/influenza/HKUstudy2009/
#     - home_pcr.csv / qPCR.csv : hhID, member, visit, qPCR
#     - hchar_h.csv             : hhID, clinic_day, v1_day, v2_day, v3_day
#     - antiviral_m.csv         : hhID, member, av, dayfrom, dayto
#     - incomplete_m.csv        : hhID, member


### Libraries ###
library(dplyr)
library(ggplot2)
library(here)
  here::i_am("influenza/Influ_q.R")


### Variables 
qPCR_threshold <- 900  # Detection threshold used in the HKU studies (copies/mL)
r_tolerance    <- 0.94 # Minimum value for duration probabilities to be treated as "certainly infectious" when estimating r


### Function: estimate r, q_i and d from empirical duration probabilities ###

#' Estimation of the parameters of the infectious period from empirical duration probabilities
#'
#' @param durationprobs_vector  Numeric vector of the empirical durations probabilities P(duration>=i)
#' @param tolerance Minimum value of duration probabilities to be treated as "certainly infectious"
#'
#' @return A list with:
#'   r    - estimated fixed minimum infectious period
#'   q    - duration probabilities q_i = P(D >= i), i = 0, ..., d
#'   d    - stochastic duration of the infection

# reduce_to_model_q     
estimate_rqd <- function(durationprobs_vector, tolerance = r_tolerance) {
  n <- length(durationprobs_vector)
  
  ### Estimation of r ###
  below_tol <- which(durationprobs_vector[-1] < tolerance)
  r <- if (length(below_tol) == 0) n - 1 else below_tol[1] - 1
  
  ### Estimation of q_i ###
  if (r + 1 >= n) {
    q <- 1  # every duration probabilities is within tolerance: no stochastic period
  } else {
    q <- c(1, durationprobs_vector[(r + 2):n])
  }
  
  ### Estimation of d ###
  d <- length(q) - 1
  
  ### Output ###
  list(r = r, d = d, q = round(q, 4))
}


### Function: estimate the parameters of the infectious duration for a single study ###

#' Estimation of the parameters r, d and q from HKU data
#'
#' @param pcr_path        Path to the PCR results csv (hhID, member, visit, qPCR)
#' @param hchar_path      Path to the household visit-day csv
#' @param antiviral_path  Path to the antiviral-use csv
#' @param incomplete_path Path to the incomplete-follow-up csv
#'
#' @return A list with:
#'   r          - estimated fixed minimum infectious period
#'   d          - estimated stochastic infectious period
#'   q          - data.frame(i, qi) with the estimated q_i
#'   durations  - data.frame with the observed infectious duration per individual
#'   n          - number of confirmed secondary infections used in the estimation

estimate_durationparam <- function(pcr_path, hchar_path, antiviral_path, incomplete_path) {
  
  ### Read data ###
  pcr        <- read.csv(pcr_path,        stringsAsFactors = FALSE)
  hchar      <- read.csv(hchar_path,      stringsAsFactors = FALSE)
  antivirals <- read.csv(antiviral_path,  stringsAsFactors = FALSE)
  incomplete <- read.csv(incomplete_path, stringsAsFactors = FALSE)
  
  ### Assign the real sample day to each PCR measurement ###
  # visit 0 = clinic visit (clinic_day); visits 1-3 = home visits (v1-v3_day).
  # All days are measured from the index case's symptom onset.
  visit_days <- hchar %>%
    select(hhID, clinic_day, v1_day, v2_day, v3_day)
  
  pcr_with_day <- pcr %>%
    left_join(visit_days, by = "hhID") %>%
    mutate(
      sample_day = case_when(
        visit == 0 ~ clinic_day,
        visit == 1 ~ v1_day,
        visit == 2 ~ v2_day,
        visit == 3 ~ v3_day,
        TRUE       ~ NA_real_
      )
    ) %>%
    filter(!is.na(sample_day), !is.na(qPCR))
  
  ### Keep secondary contacts only (exclude index cases, member == 0) ###
  pcr_secondary <- pcr_with_day %>%
    filter(member != 0)
  
  ### Exclude individuals with incomplete follow-up ###
  pcr_secondary <- pcr_secondary %>%
    anti_join(incomplete, by = c("hhID", "member"))
  
  ### Exclude individuals who received antivirals ###
  antiviral_individuals <- antivirals %>%
    distinct(hhID, member)
  
  pcr_secondary <- pcr_secondary %>%
    anti_join(antiviral_individuals, by = c("hhID", "member"))
  
  ### Identify PCR-confirmed infections ###
  # Criterion: at least one qPCR reading above qPCR_threshold.
  infected <- pcr_secondary %>%
    group_by(hhID, member) %>%
    summarise(
      any_positive = any(qPCR > qPCR_threshold, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(any_positive)
  
  pcr_infected <- pcr_secondary %>%
    semi_join(infected, by = c("hhID", "member"))
  
  n_individuals <- nrow(infected)
  
  ### Compute the observed infectious duration per individual ###
  # Duration = day of the last PCR-positive sample (> qPCR_threshold).
  durations <- pcr_infected %>%
    group_by(hhID, member) %>%
    summarise(
      duration = max(sample_day[qPCR > qPCR_threshold], na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    filter(is.finite(duration))
  
  ### Estimate q_i ###
  i_max  <- max(durations$duration)
  i_vals <- 0:i_max
  
  qi_vals <- sapply(i_vals, function(i) mean(durations$duration >= i, na.rm = TRUE))
  
  qi_est <- round(estimate_rqd(qi_vals)$q, 4)
  i_est <- 0:(length(qi_est)-1)
  
  r_est  <- estimate_rqd(qi_vals)$r
  d_est <- estimate_rqd(qi_vals)$d
  q_est <- data.frame(i = i_est, qi = qi_est)
  
  return(list(r = r_est, d = d_est, q = q_est, durations = durations, n = n_individuals))
}


### Function: estimate the parameters of the infectious duration for the combinations of the HKU studies ###

#' Estimation of the parameters r, d and q from HKU data (2008 and 2009)
#'
#' @param study2008    Application of the function estimate_durationparam() to the data of the HKU study of 2008
#' @param study2009    Application of the function estimate_durationparam() to the data of the HKU study of 2009
#'
#' @return A list with:
#'   r    - estimated fixed minimum infectious period
#'   d    - estimated stochastic infectious period
#'   q    - data.frame(i, qi) with the estimated q_i
#'   n    - number of confirmed secondary infections used in the estimation

combined_durationparam <- function(study2008, study2009){
  
  ### Combine the duration probabilities of the two studies ###
  combined_data <- bind_rows(
    res_2008$durations %>% mutate(study = "2008"),
    res_2009$durations %>% mutate(study = "2009")
  )
  
  ### Define the global duration ###
  i_max  <- max(combined_data$duration)
  i <- 0:i_max
  
  ### Estimate the global duration probabilities ###
  combined_durationprobs <- round(sapply(i, function(k) mean(combined_data$duration >= k, na.rm = TRUE)), 4)
  
  ### Estimate the parameters of the infectious duration ###
  result <- estimate_rqd(durationprobs_vector = combined_durationprobs, tolerance = r_tolerance)
  
  return(list(r = result$r, d = result$d, q = result$q, n = study2008$n + study2009$n))
  
}


### Check that the data has been downloaded ###
# The raw data is NOT included in this repository (see README.md 
# for download links and licensing).

influenza_data_dir <- here("data", "influenza")

if (!file.exists(here(influenza_data_dir, "HKUstudy2008", "home_pcr.csv")) ||
    !file.exists(here(influenza_data_dir, "HKUstudy2009", "qPCR.csv"))) {
  stop(
    "Could not find the HKU 2008/2009 study data. ",
    "Download it following the instructions in the 'Data' section of ",
    "README.md and place it in data/influenza/HKUstudy2008/ and ",
    "data/influenza/HKUstudy2009/ respectively."
  )
}


### Application to both HKU studies ###

res_2008 <- estimate_durationparam(
  pcr_path        = here(influenza_data_dir, "HKUstudy2008", "home_pcr.csv"),
  hchar_path      = here(influenza_data_dir, "HKUstudy2008", "hchar_h.csv"),
  antiviral_path  = here(influenza_data_dir, "HKUstudy2008", "antiviral_m.csv"),
  incomplete_path = here(influenza_data_dir, "HKUstudy2008", "incomplete_m.csv")
)

res_2009 <- estimate_durationparam(
  pcr_path        = here(influenza_data_dir, "HKUstudy2009", "qPCR.csv"),
  hchar_path      = here(influenza_data_dir, "HKUstudy2009", "hchar_h.csv"),
  antiviral_path  = here(influenza_data_dir, "HKUstudy2009", "antiviral_CREAT.csv"),
  incomplete_path = here(influenza_data_dir, "HKUstudy2009", "incomplete_CREAT.csv")
)


### Combination of both studies ###

durationparam <- combined_durationparam(res_2008, res_2009)


### Output: r, d and q ###
# These are the objects other scripts depend on (BB_SEM24_Functions.R,
# Influ_p.R, influ_HKUStudy.R).
r_influ <- durationparam$r
q_influ <- durationparam$q
d_influ <- durationparam$d


### Plot (only when run interactively, not when sourced by other scripts) ###
if (interactive()) {
  
  color1 <- "#00798C"
  color2 <- "#FE7F2D"
  color3 <- "#FCCA46"
  colors <- c("Combined" = color1, "2008" = color2, "2009" = color3)
  
  x_max <- max(res_2008$q$i, res_2009$q$i, d_influ)
    
  ggplot() +
    geom_line(aes(x = res_2008$q$i, y = res_2008$q$qi, color = "2008"), linewidth = 1.5, linetype = "dashed") +
    geom_line(aes(x = res_2009$q$i, y = res_2009$q$qi, color = "2009"), linewidth = 1.5, linetype = "dashed") +
    geom_line(aes(x = 0:d_influ, y = q_influ, color = "Combined"), linewidth = 1.5) +
    geom_point(aes(x = 0:d_influ, y = q_influ, color = "Combined"), size = 4) +
    scale_x_continuous(breaks = seq(0, x_max, 2)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    theme_grey(base_size = 25) +
    labs(x = expression(i), y = expression(q[i]), color = "Legend") +
    scale_color_manual(values = colors, breaks = names(colors)[c(1, 2, 3)]) +
    guides(color = guide_legend(override.aes = list(shape = NA))) +
    theme(legend.position = c(0.9, 0.83))
}
