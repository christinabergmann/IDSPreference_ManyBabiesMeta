
  
# PRELIMINARIES --------------------------------------------------------

library(tidyverse) 
library(knitr)
library(here)
library(tableone)

# expects global vars set by master prep script: ic.dataset, age.matched

# to reproduce main analyses, set to FALSE
# If working on the preregistration, set this to T, otherwise F to use the veridical dataset without scrambling.
prereg = FALSE


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
setwd(data.dir)
# codebook: https://github.com/langcog/metalab2/blob/master/metadata/spec.yaml


# For preregistration, we used a scrambled dataset
if(prereg){
  
  COMBINED_DATASET <- here("data/mb_ma_combined.csv")
  OUTFILE <- here("data/mb_ma_combined_scrambled.csv")
  
  full_dataset <- read_csv(COMBINED_DATASET)
  
  set.seed(1819)
  full_dataset_shuffled <- full_dataset %>%
    mutate(d_calc = sample(d_calc, replace = F), 
           d_var_calc = sample(d_var_calc, replace = F))
  
  write_csv(full_dataset_shuffled, OUTFILE)
  
  d = full_dataset_shuffled
  
} else {
  if ( ic.dataset == FALSE & age.matched == FALSE ) d = read_csv("mb_ma_combined_0.125.csv")
  if ( ic.dataset == TRUE & age.matched == FALSE ) d = read_csv("mb_ma_combined_0.75.csv")
  if ( ic.dataset == FALSE & age.matched == TRUE ) d = read_csv("mb_ma_combined_0.125_age_matched.csv")
  if ( ic.dataset == TRUE & age.matched == TRUE ) stop("Case not handled")
}



# RECODE MODERATORS --------------------------------------------------------

# list of moderators (not yet created)
mods = c( "mean_agec_mos",
          "test_lang",  # whether stimuli were in native language
          "method",
          
          # constant in RRR:
          "speech_type",
          "own_mother",
          "presentation",
          "dependent_measure",
          "main_question_ids_preference"#,
          
          # varied in RRR:
          #"stimulus_set", # ~~~ not in the dataset 
          #"trial_control" # ~~~ not in the dataset
          ) 
#"human_coded", # ~~~ not in the dataset

# which are continuous?
contMods = "mean_agec_mos"


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


# make mean_agec_mos in MONTHS
# mean_age is coded in days
daysPerMonth = 30.44
d$mean_agec_mos = d$mean_age/daysPerMonth - mean( d$mean_age[ d$isMeta == TRUE ]/daysPerMonth, na.rm = TRUE)

# sanity check: mean age by source
d %>% group_by(study_type) %>%
  summarise(mean(mean_agec_mos))

# Now fix method names
d = d  %>% mutate(method = ifelse(method %in% c("singlescreen", "eyetracking"), "cf", 
                                  ifelse(method %in% c("fc", "cht"), "other", method)))

# And collapse dependent variable 
d = d %>% mutate(dependent_measure = ifelse(dependent_measure == "facial_expression", "affect", "preference"))


# distribution of moderators in RRR and MA
# mean of mean_agec_mos should be 0 in the MA but nonzero in MB
CreateTableOne(vars = mods, 
               strata = "study_type",
               data = d)

# recode all moderators so reference levels are mode in meta-analysis
#  (recode via alphabetization)
d = d %>% mutate_at( .vars = mods[ !mods == contMods ],
                     .funs = function(.v) code_mode_as_ref(vec = .v,
                                                           isMeta = d$isMeta) ) #Needed to hack this


# also recode as dummy variables for meta-regression joy
dummyMods = d %>% select(mods) %>%
  dummy_cols() %>%
  # remove pre-recoding moderators to avoid duplicated names below
  select(-mods)
d = bind_cols( d, dummyMods)


# same table again, after recoding
CreateTableOne(vars = mods, 
               strata = "study_type",
               data = d)



# MAKE OTHER VARIABLES --------------------------------------------------------

# unique ID
# don't use expt_num here because it's only defined for MA studies
d = d %>% group_by(study_id) %>%
  mutate( unique = paste( study_id, row_number(), sep = ", #" ) )


# rename for simplicity
d$yi = d$d_calc
d$vi = d$d_var_calc
d$sei = sqrt(d$vi)

# CI limits (for plotting)
d$lo = d$yi - qnorm(0.975) * sqrt(d$vi)
d$hi = d$yi + qnorm(0.975) * sqrt(d$vi)

# recode for plotting joy
d$studyTypePretty = NA
d$studyTypePretty[ d$isMeta == TRUE ] = "Meta-analysis"
d$studyTypePretty[ d$isMeta == FALSE ] = "Replications"

d$sourcePretty = "Meta-analysis"
d$sourcePretty[ d$isMeta == FALSE ] = "Replications"

# p-values and affirmative status for publication bias analyses
d$pval = 2 * ( 1 - pnorm( abs(d$yi/d$sei) ) )
d$affirm = ( d$yi > 0 ) & ( d$pval < 0.05 )


# for additional publication bias analysis:
# add indicator for whether effect was *reported* to be significant
if ( ic.dataset == FALSE & age.matched == FALSE ) {
  d$pvalSignif = as.numeric(d$pval < 0.05)
  
  d$reportedSignif = NA
  d$reportedSignif[ d$effect_significance_reported == "significant" ] = 1 
  d$reportedSignif[ d$effect_significance_reported == "non-significant" ] = 0
}


# SAVE DATASET --------------------------------------------------------

setwd(data.dir)

if ( prereg == TRUE ) write.csv(d, "mb_ma_combined_scrambled_prepped.csv")

if ( prereg == FALSE ) {
  if ( ic.dataset == FALSE & age.matched == FALSE ) write.csv(d, "mb_ma_combined_prepped_0.125.csv")
  if ( ic.dataset == TRUE & age.matched == FALSE ) write.csv(d, "mb_ma_combined_prepped_0.75.csv")
  if ( ic.dataset == FALSE & age.matched == TRUE ) write.csv(d, "mb_ma_combined_prepped_0.125_age_matched.csv")
  
}




