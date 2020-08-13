
############################## PRELIMINARIES ############################## 

library(tidyverse) 
library(knitr)
library(here)
library(tableone)

data.dir = here("data")
# where to save results
results.dir = here("results_from_R")
# results.dir = "~/Dropbox/Personal computer/Independent studies/2020/Christina's ManyBabiesMeta (MB-Meta)/IDSPreference_ManyBabiesMeta/results_from_R"
overleaf.dir = "~/Dropbox/Apps/Overleaf/MB-Meta/R_objects"
code.dir = here("analyses/1_mungedata")

setwd(code.dir)
source("prep_helper.R")

# should we use the grateful package to scan and cite packages?
cite.packages.anew = FALSE

# helper code
# source( here("analyses/2_analyze/analyze_helper.R") )

# read in dataset
# scrambling was done by the "outline" Rmd file
setwd(data.dir)
# codebook: https://github.com/langcog/metalab2/blob/master/metadata/spec.yaml
d = read_csv("mb_ma_combined_scrambled.csv")


############################## RECODE MODERATORS ############################## 

# list of moderators
mods = c( "mean_agec",
          "test_lang",  # whether stimuli were in native language; almost constant in meta
          "method",
          
          # constant in RRR:
          "speech_type",
          "own_mother",
          "presentation",
          "dependent_measure",
          "main_question_ids_preference",
          
          # varied in RRR:
          #"stimulus_set", # ~~~ not in the dataset# ~~~ not in the datasete
          "trial_control" )
#"human_coded", # ~~~ not in the dataset

# which are continuous?
contMods = "mean_agec"


# fix inconsistent capitalization
d = d %>% mutate_at( .vars = c("method", "speech_type", "speaker", "presentation"),
                 .funs = tolower )

# collapse categories of speaker
d$own_mother = (d$speaker == "child’s mother")

# dummies for being meta-analysis and being replication
d$isMeta = (d$study_type == "MA")
d$isRep = (d$study_type == "MB")


# center continuous moderators
# age in months
d$mean_agec = d$mean_age/12 - mean( d$mean_age[ d$isMeta == TRUE ]/12, na.rm = TRUE)


# distribution of moderators in RRR and MA
# mean of mean_agec should be 0 in the MA but nonzero in MB
CreateTableOne(vars = mods, 
               strata = "study_type",
               data = d)


# recode all moderators so reference levels are mode in meta-analysis
#  (recode via alphabetization)
d = d %>% mutate_at( .vars = mods[ !mods == contMods ],
                     .funs = function(.v) code_mode_as_ref(vec = .v,
                                                           isMeta = isMeta) )


# same table again, after recoding
CreateTableOne(vars = mods, 
               strata = "study_type",
               data = d)



############################## MAKE OTHER VARIABLES ############################## 

# rename for simplicity
d$yi = d$d_calc
d$vi = d$d_var_calc
d$sei = sqrt(d$vi)

# p-values and affirmative status for publication bias analyses
d$pval = 2 * ( 1 - pnorm( abs(d$yi/d$sei) ) )
d$affirm = ( d$yi > 0 ) & ( d$pval < 0.05 )


############################## SAVE DATASET ############################## 

setwd(data.dir)
write.csv(d, "mb_ma_combined_scrambled_prepped.csv")




