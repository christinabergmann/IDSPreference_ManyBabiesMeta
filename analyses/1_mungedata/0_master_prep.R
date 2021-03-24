

# This file runs all data prep scripts in order (#1-4) and handles the need to run 
#  some files multiple times to create different versions of datasets (for age-matched analyses
# and analyses that change the inclusion criteria).

library(here)
setwd( here("analyses/1_mungedata") )

# only needs to be run once regardless of which datasets are to be created
source("1_get_tidy_MA_data.R")

# number of trials that must be passed
# There are 8 pairs in total, so 4 pairs = 50%, 6 pairs = 75%
# main analysis is "1"
criteria.vec = c(1, 4, 6)


for ( .c %in% criteria.vec ) {
  
  n_trial_pairs_criterion = .c
  
  # for main analysis, run as both age-matched and not matched
  if ( .c == 1 ) {
    
    for ( .a in c(TRUE, FALSE) ) {
      # for replications, should we look at the much smaller, age-matched dataset?
      age.matched = .a
      source("2_get_tidy_MB_data.R")
    }
    
  } else {
    age.matched == FALSE
  }
}
