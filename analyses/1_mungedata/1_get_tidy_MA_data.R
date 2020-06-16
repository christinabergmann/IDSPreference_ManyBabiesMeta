# The final ManyBabies1 dataset from the public github repository, plus the processing pipeline from the paper
library(tidyverse)
library(here)
library(janitor)
source(here("analyses/1_mungedata/compute_es_IDS.R"))

MONTH_IN_DAYS <- 365.25/12
MA1_PATH <- here("data/IDS Preference Meta-Analytic Dataset - Sheet1.csv")
MA_OUT_PATH <- here("data/ma_data_tidy.csv")

TARGET_VARS <- c("study_id", "short_cite", "main_question_ids_preference", "trial_control",
                 "response_mode", "exposure_phase", "method", "dependent_measure", "participant_design",
                 "native_lang", "test_lang", "infant_type", "group_name_1", "group_name_2",
                 "t", "n_1", "n_2", "mean_age_1", "same_infant", "mean_age_2", "x_1", "x_2", "sd_1", "sd_2", "d", "d_var",
                 "r", "corr", "gender_1", "num_trials", "speaker_fam", "speaker_experience", "speaker_female",
                 "speaker", "presentation", "setting", "speech_type")

ma_data_raw <- read_csv(file = MA1_PATH) %>%
  clean_names()

ma_data_tidy <- ma_data_raw %>%
  select(all_of(TARGET_VARS)) %>%
  mutate(id = 1:n())

# calculate effect sizes
ma_data_tidy_with_es <-  ma_data_tidy %>%
  group_by(id) %>%
  nest() %>%
  mutate(es_data = map(data, compute_es)) %>%
  unnest() %>%
  select(id, d_calc, d_var_calc, es_method) %>%
  ungroup()

# get study characteristics
study_moderators <- ma_data_tidy %>%
        group_by(id) %>%
        mutate(n = mean(c(n_1, n_2), na.rm = TRUE),
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
         -month, -r,-t, -d, -d_var, -corr)

ma_data <- full_join(ma_data_tidy_with_es, study_moderators) %>%
  select(-id) %>%
  rename(prop_female = gender)

write_csv(ma_data, MA_OUT_PATH)
# rename variables so they match up
