if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               here,
               janitor)

source(here("analyses","1_mungedata","compute_es_IDS.R"))

MONTH_IN_DAYS <- 365.25/12

if ( use.corrected.dunst == FALSE ) {
  MA1_PATH <- here("data","from_data_team","Dunst_original.csv")
  MA_OUT_PATH <- here("data","prepped_with_original_dunst","ma_data_tidy.csv")
}

if ( use.corrected.dunst == TRUE ) {
  MA1_PATH <- here("data","from_data_team","Dunst_corrected.csv")
  MA_OUT_PATH <- here("data","prepped_with_corrected_dunst","ma_data_tidy.csv")
}

RECOMPUTE_D <- use.corrected.dunst 
#if we use the original (uncorrected) meta-analysis, d is *never* recomputed
#if we use the corrected meta-analysis, d is recomputed where possible

TARGET_VARS <- c("study_id", "short_cite", "expt_num", "original_ma", "main_question_ids_preference", 
                 "response_mode", "exposure_phase", "method", "dependent_measure", "participant_design",
                 "native_lang", "test_lang", "infant_type", "group_name_1", "group_name_2",
                 "t", "n_1", "n_2", "mean_age_1", "same_infant", "mean_age_2", "x_1", "x_2", "sd_1", "sd_2", "d", "d_var",
                 "r", "corr", "gender_1", "num_trials", "speaker_fam", "speaker_experience", "speaker_female",
                 "speaker", "presentation", "setting", "speech_type", "effect_significance_reported")


ma_data_raw <- read_csv(file = MA1_PATH) %>%
  clean_names()

ma_data_tidy <- ma_data_raw %>%
  select(all_of(TARGET_VARS)) %>%
  mutate(id = 1:n()) 

ma_data_tidy_with_es <-  ma_data_tidy %>%
  group_by(id) %>%
  nest() %>%
  mutate(es_data = map(data, compute_es,recompute_d=RECOMPUTE_D)) %>%
  unnest(cols = c(data, es_data)) %>%
  ungroup()

ma_data_tidy <- ma_data_tidy_with_es %>%
  filter(!is.na(d_calc)) %>% # we don't have all effect sizes
  filter(!is.na(d_var_calc)) %>% # remove if we don't have the variance
  select(-c(d, d_var)) %>%
  rename(d = d_calc, d_var = d_var_calc) # If d was recomputed, we replace the reported d column with the recalculated ones. 
# When we couldn't recalculate d (or for the uncorrected meta-analysis analysis), we use the one reported in the meta-analyses. 
# We always compute d_var from scratch (not reported in the main meta-analysis).
# See "es_method" in the intermediate dataset.

remove(ma_data_tidy_with_es)


# get study characteristics
study_moderators <- ma_data_tidy %>%
        group_by(id) %>%
        mutate(n = sum(c(n_1, n_2), na.rm = TRUE),
              mean_age = mean(c(mean_age_1, mean_age_2), na.rm = TRUE),
              month = round(mean_age/MONTH_IN_DAYS,2),
              age_group =  case_when(month >= 0 & month <= 3.5 ~ "0-3 mo",
                                     month > 3 & month <= 6.5 ~ "3-6 mo",
                                     month > 6.5 & month <= 9.5 ~ "6-9 mo",
                                     month > 9.5 & month <= 12.5 ~ "9-12 mo",
                                     month > 12.5 & month <= 15.5 ~ "12-15 mo",
                                     month > 15.5 & month <= 18.5 ~ "15-18 mo")) %>%
  rename(gender = gender_1) %>%
  select(-contains("group_name"), -contains("_1"), -contains("_2"),
         -month, -r,-t, -corr, -original_ma)


ma_data <- full_join(ma_data_tidy, study_moderators) %>%
  select(-id) %>%
  rename(prop_female = gender) %>% 
  mutate(speech_type = case_when(
    speech_type %in% c("Filtered", "Synthesized") ~ "filtered/synthesized", 
    TRUE ~ speech_type))

write_csv(ma_data, MA_OUT_PATH)

# # 2022-1-5: MM doing a sanity check given unexpected changes to results of corrected
# #  Dunst analysis
# # fit a simple subset model to ma_data
# robu( d ~ 1, 
#       data = ma_data, 
#       studynum = as.factor(study_id),
#       var.eff.size = d_var,
#       modelweights = "HIER",
#       small = TRUE)
# # yes, this still yields 0.448, which is much less than the naive 0.70










