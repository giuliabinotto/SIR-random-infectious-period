###########################################################
#### Stochastic epidemics models with random duration #####
######### Influenza: Hong Kong University study ###########
###### Estimation of the Secondary Attack Rate (SAR) ######
###########################################################
################# X. Bardina, G. Binotto ##################
###########################################################


# This script estimate the secondary attack rate from the combined HKU 2008 and 2009 datasets.
#
# Source this file to make the SAR available, e.g.:
#   source("influenza/Influ_SAR.R")
#
# Required data (NOT included in this repo, see README.md "Data"):
#   data/influenza/HKUstudy2008/
#   data/influenza/HKUstudy2009/
#     - home_pcr.csv / qPCR.csv : hhID, member, qPCR


### Libraries ###
library(dplyr)
library(here)
  here::i_am("influenza/Influ_SAR.R")


### Sources ###
#source("R/BB_SEM24_Functions.R")


### Variables 
qPCR_threshold <- 900  # Detection threshold used in the HKU studies (copies/mL)


### Function: estimate secondary attack rate from the HKU 2008 and 2009 datasets ###

#' @param study2008_path    Path to the PCR results csv (hhID, member, qPCR) from HKU 2008 study
#' @param study2009_path    Path to the PCR results csv (hhID, member, qPCR) from HKU 2009 study
#'
#' @return A list with:
#'   SAR          - estimated global SAR
#'   SAR_hh_size  - estimated SAR by household size
#'   hh_sizes     - household sizes

estimate_SAR <- function(study2008_path, study2009_path) {
  
  ### Read data ###
  pcr_2008 <- read.csv(study2008_path, stringsAsFactors = FALSE)
  pcr_2009 <- read.csv(study2009_path, stringsAsFactors = FALSE)
  
  ### Create merged dataset (columns: hhID, member, qPCR) ### 
  cols <- c("hhID", "member", "qPCR")
  pcr <- bind_rows(
    select(pcr_2008, all_of(cols)),
    select(pcr_2009, all_of(cols))
  )
  
  ### Modify database pcr to add missing index cases ###
  # Check for each household if index case (member==0) is present
  no_index_household <- pcr %>%
    group_by(hhID) %>%
    filter(!any(member == 0)) %>%
    distinct(hhID) %>%
    pull()
  
  # Add index case where missing
  pcr <- pcr %>%
    bind_rows(
      tibble(
        hhID = no_index_household,
        member = 0,
        qPCR = NA
      )
    ) %>%
    arrange(hhID, member)
  
  ### Estimate global SAR (Secondary Attack Rate) ###
  # Consider only household contacts (without index cases)
  pcr_contacts <- pcr %>%
    filter(member != 0)
  
  # Total number of unique contacts that were followed up
  total_contacts <- pcr_contacts %>%
    distinct(hhID, member) %>%
    nrow()
  
  # Secondary infections (viral loads > qPCR_threshold)
  sec_infections <- pcr_contacts %>%
    filter(qPCR > qPCR_threshold) %>%
    distinct(hhID, member) %>%
    nrow()
  
  # Global SAR
  SAR <- (sec_infections / total_contacts)
  
  ### Estimate SAR by household size ###
  # Frequency table of household size
  hh_size <- aggregate(member ~ hhID, data = pcr, FUN = function(x) length(unique(x)))
  freq_hh_size <- table(hh_size$member)   # Household size (including index case)
  freq_hh_size_noindex <- table(hh_size$member-1)   # Household size (excluding index case)
  
  # SAR by household size
  SAR_hh_size <- pcr_contacts %>%
    group_by(hhID) %>% 
    summarise(
      hh_size = n_distinct(member) + 1,   # Household size: unique contacts + 1 (index case)
      sec_infected = n_distinct(member[qPCR > qPCR_threshold & !is.na(qPCR)]),   # Unique members with at least one PCR > qPCR_threshold
      .groups = "drop"
    ) %>%
    group_by(hh_size) %>%
    summarise(
      SAR_hh = sum(sec_infected) / sum(hh_size - 1),
      .groups = "drop"
    ) %>%
    arrange(hh_size)
  
  ### Output: SAR (global and by household size) and  frequency table of household size ###
  return(list(SAR = SAR, SAR_hh_size = SAR_hh_size, hh_sizes = sort(unique(hh_size$member))))
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


### SAR estimation ###
res_SAR <- estimate_SAR(
  study2008_path = here(influenza_data_dir, "HKUstudy2008", "home_pcr.csv"),
  study2009_path = here(influenza_data_dir, "HKUstudy2009", "qPCR.csv")
)


### Output: SAR ###
# These are the objects other scripts depend on (Influ_p.R, influ_HKUStudy.R).
SAR_influ <- res_SAR$SAR
hh_sizes_influ <- res_SAR$hh_sizes


