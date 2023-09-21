library(here)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here,
               testthat)

set.seed(4711)

# 4 versions: using original, corrected Dunst, extended community-augmented meta-analysis data
ma.versions <- c("Dunst_original","Dunst_corrected","augmented_ma_extended")
for (ma.version in ma.versions) {
    # ~ Code-Running Parameters ------------------------------------------------------------------
    # should we remove existing results file instead of overwriting individual entries? 
    start.res.from.scratch = TRUE
    # should we use the grateful package to scan and cite packages?
    cite.packages.anew = TRUE
    # should we bootstrap from scratch or read in old resamples?
    boot.from.scratch = TRUE
    # make plots from scratch?
    redo.plots = TRUE
    
    setwd( here("analyses","2_analyze") )
    source("1_analyze.R")
}

