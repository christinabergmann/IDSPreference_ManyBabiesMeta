

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

for ( .c in criteria.vec ) {
  
  n_trial_pairs_criterion = .c
  
  # for main analysis, run as both age-matched and not matched
  if ( .c == 1 ) {
    for ( .a in c(FALSE, TRUE) ) {
      # for replications, should we look at the much smaller, age-matched dataset?
      age.matched = .a
      setwd( here("analyses/1_mungedata") )
      suppressMessages( source("2_get_tidy_MB_data.R") )
    }
  }
  
  # for sensitivity analyses with other inclusion criteria, 
  #  don't also do age-matching
  if ( .c != 1 ) {
    age.matched = FALSE
    setwd( here("analyses/1_mungedata") )
    suppressMessages( source("2_get_tidy_MB_data.R") )
  }
}


#bm: now need to run 3_merge_MA_MB.R, which had the following args:s

# for replications, should we make the sensitivity-analysis dataset with more
#  stringent inclusion (set to TRUE), or the main-analysis dataset?
ic.dataset = FALSE

# for replications, should we look at the much smaller, age-matched dataset instead?
age.matched = FALSE
  
  

