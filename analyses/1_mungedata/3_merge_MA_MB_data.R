# Merge study level mb and ma data

library(tidyverse)
library(here)

MB_DATA_PATH <- here("data/mb_data_tidy.csv")
MA_DATA_PATH <- here("data/ma_data_tidy.csv")
OUTFILE <- here("data/mb_ma_combined.csv")

mb_data_raw <- read_csv(MB_DATA_PATH)
ma_data_raw <- read_csv(MA_DATA_PATH)


mb_data_tidy <- mb_data_raw %>%
  rename(study_id = lab,
         d_calc = d_z,
         d_var_calc = d_z_var,
         native_lang = modal_lang1) %>%
  mutate(study_type = "MB",
         short_cite = "ManyBabies Consortium (2020)", 
         main_question_ids_preference = "yes",
         dependent_measure = "looking_time",
         exposure_phase = "test_only",
         test_lang = ifelse(prop_nae>.9, "native", "nonnative"),
         infant_type = case_when(prop_preterm > .25 | prop_curr_earinfection > .25 |
                                   prop_hearing_vision > .25 | prop_cognitive_developmental > .25 |
                                   prop_monolingual < .75 | mean_lang1_exposure < .75 ~ "non-typical",
                                 TRUE ~ "typical")) %>% # define typicality based on demographic variables
  select(-prop_preterm, -prop_curr_earinfection,
         -prop_hearing_vision, -prop_cognitive_developmental,
         -prop_monolingual, -mean_lang1_exposure) %>%
  group_by(study_id) %>%
  mutate(same_infant = 1:n())


ma_data_tidy <- ma_data_raw %>%
  rename(d_calc = d,
         d_var_calc = d_var) %>%
  mutate(study_type = "MA",
         prop_nae = case_when(native_lang %in% c("American English", "Canadian English") ~ 1,
                              TRUE ~ 0),
         native_lang =case_when(str_detect(native_lang, "English") ~ "english",
                                TRUE ~ tolower(native_lang)),
         infant_type = case_when(infant_type != "typical" ~ "non-typical",
                                 TRUE ~ infant_type),
         prop_caregiver_seat = NA,
         parent_education = NA)

tidy_df <- bind_rows(mb_data_tidy, ma_data_tidy) %>%
  select(study_type, study_id, short_cite,  mean_age, age_group,
         n, d_calc, d_var_calc, native_lang, prop_female, prop_nae, infant_type,
         main_question_ids_preference, same_infant, method, dependent_measure,  presentation, response_mode,
         everything())


write_csv(tidy_df, OUTFILE)
