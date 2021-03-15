

# This file checks key analysis results by refitting models "manually"
#   i.e., without calling the fancy helper fns


# PRELIMINARIES ------------------------------------------------------------------

library(tidyverse) 
library(knitr)
library(here)
library(tableone)
library(corrr)
library(robumeta)
library(fastDummies)
library(weightr)
library(PublicationBias)
library(xtable)
library(boot)
library(testthat)
library(ggplot2)
library(metafor)
library(MatchIt)