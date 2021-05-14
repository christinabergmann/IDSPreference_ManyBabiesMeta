

library(here)
library(testthat)
setwd( here("analyses/2_analyze") )

# 2 versions: using original or corrected Dunst data
for ( .u in c(FALSE, TRUE) ) {
  use.corrected.dunst = .u
  source("1_analyze.R")
}
