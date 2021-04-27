

# This file runs all data prep scripts in order (#1-4) and handles the need to run 
#  some files multiple times to create different versions of datasets (for age-matched analyses
# and analyses that change the inclusion criteria).

library(here)
library(testthat)
setwd( here("analyses/1_mungedata") )

# STEP 1: Get Tidy MA Data --------------------------------------------------------
# only needs to be run once regardless of which datasets are to be created

# 2 versions: using original or corrected Dunst data
for ( .u in c(FALSE, TRUE) ) {
  use.corrected.Dunst = .u
  source("1_get_tidy_MA_data.R")
}


# STEP 2: Get Tidy MB Data --------------------------------------------------------

# this step happens for each possible inclusion criterion
# and within the main-analysis inclusion criterion, makes both
#  age-matched and non-matched datasets

# number of trials that must be passed
# There are 8 pairs in total, so 4 pairs = 50%, 6 pairs = 75%
# main analysis is 1 = 0.125*8

# only .125 and .75 are used in analyses
criteria.vec = c(.125, .25, .5, .75)*8


# even though this part of script preps the MB data,
#  depends on use.corrected.dunst because of age-matching
for ( .u in c(FALSE, TRUE) ) {
  
  use.corrected.dunst = .u
  
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
    #  avoid age-matching
    if ( .c != 1 ) {
      age.matched = FALSE
      setwd( here("analyses/1_mungedata") )
      suppressMessages( source("2_get_tidy_MB_data.R") )
    }
  }
  
}

# STEP 3-4: Merge MA and MB Data and Prep for Analysis --------------------------------------------------------

# .ic: for replications, should we make the sensitivity-analysis dataset with more
#  stringent inclusion (set to TRUE), or the main-analysis dataset?

for ( .u in c(FALSE, TRUE) ) {
  
  use.corrected.dunst = .u
  
  for ( .ic in c(FALSE, TRUE) ) {
    
    ic.dataset = .ic
    
    # for main analysis, run as both age-matched and not matched
    if ( .ic == FALSE ) {
      for ( .a in c(FALSE, TRUE) ) {
        # for replications, should we look at the much smaller, age-matched dataset?
        age.matched = .a
        setwd( here("analyses/1_mungedata") )
        suppressMessages( source("3_merge_MA_MB_data.R") )
        setwd( here("analyses/1_mungedata") )
        suppressMessages( source("4_prep_for_analysis.R") )
      }
    }
    
    # for sensitivity analyses with other inclusion criteria, 
    #  don't also do age-matching
    if ( .ic == TRUE ) {
      age.matched = FALSE
      setwd( here("analyses/1_mungedata") )
      suppressMessages( source("3_merge_MA_MB_data.R") )
      setwd( here("analyses/1_mungedata") )
      suppressMessages( source("4_prep_for_analysis.R") )
    }
  }
  
}
