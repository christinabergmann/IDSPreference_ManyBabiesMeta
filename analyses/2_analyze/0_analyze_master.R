
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here,
               testthat)

setwd( here("analyses/2_analyze") )

# 2 versions: using original or corrected Dunst data
for ( .u in c(FALSE, TRUE) ) {
  use.corrected.dunst = .u
  # ~ Code-Running Parameters ------------------------------------------------------------------
  # should we remove existing results file instead of overwriting individual entries? 
  start.res.from.scratch = TRUE
  # should we use the grateful package to scan and cite packages?
  cite.packages.anew = TRUE
  # should we bootstrap from scratch or read in old resamples?
  boot.from.scratch = TRUE
  # make plots from scratch?
  redo.plots = TRUE
  
  source("1_analyze.R")
}
