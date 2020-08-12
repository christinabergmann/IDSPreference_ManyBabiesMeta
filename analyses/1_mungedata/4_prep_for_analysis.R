
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


# should we use the grateful package to scan and cite packages?
cite.packages.anew = FALSE

# helper code
# source( here("analyses/2_analyze/analyze_helper.R") )

# read in dataset
# scrambling was done by the "outline" Rmd file
setwd(data.dir)
# codebook: https://github.com/langcog/metalab2/blob/master/metadata/spec.yaml
d = read_csv("mb_ma_combined_scrambled.csv")


############################## RECODE VARIABLES ############################## 


# look at moderators
mods = c( "mean_agec",
          "test_lang",  # whether stimuli were in native language; almost constant in meta
          "method",
          
          # constant in RRR:
          "speech_type",
          "speaker",
          "presentation",
          "dependent_measure",
          "main_question_ids_preference",
          
          # varied in RRR:
          #"stimulus_set", # ~~~ not in the dataset# ~~~ not in the datasete
          "trial_control" )
#"human_coded", # ~~~ not in the dataset

# fix inconsistent capitalization
d = d %>% mutate_at( .vars = c("method", "speech_type", "speaker", "presentation"),
                 .funs = tolower )

# center continuous moderators
d$mean_agec = d$mean_age - mean(d$mean_age, na.rm = TRUE)

# distribution of moderators in RRR and MA
CreateTableOne(vars = mods, 
               strata = "study_type",
               data = d)



# reference levels to use:
# mean_age = 0 (mean)
# test_lang = 0 because constant in meta
# method = "cf" (most common in meta)
# 


# rename for simplicity
d$yi = d$d_calc
d$vi = d$d_var_calc
d$sei = sqrt(d$vi)


d$isMeta = (d$study_type == "MA")
d$isRep = (d$study_type == "MB")




############################## SAVE DATASET ############################## 

setwd(data.dir)
write.csv(d, "mb_ma_combined_scrambled_prepped.csv")




