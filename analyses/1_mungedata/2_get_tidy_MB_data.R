# The final ManyBabies1 dataset from the public github repository, plus the processing pipeline from the paper
library(tidyverse)
library(here)

MB1_PATH <- "https://raw.githubusercontent.com/manybabies/mb1-analysis-public/master/processed_data/03_data_diff_main.csv"
MB_OUT_PATH <- here("data/mb_data_tidy.csv")

TARGET_VARS <- c("lab", "subid_unique", "trial_num", "method", "age_days", "age_group",
                 "lang_group", "lang1", "lang1_exposure", "parenta_education", "parentb_education",
                 "caregiver_seat", "preterm",
                 "curr_earinfection", "hearing_vision", "cognitive_developmental",
                  "gender", "beard", "nae", "ADS", "IDS", "diff")

mb_data_raw <- read_csv(MB1_PATH)

mb_data_tidy <- mb_data_raw %>%
  filter(!is.na(diff)) %>%
  select(all_of(TARGET_VARS))

# tidy factors
mb_data_tidy_fct <- mb_data_tidy %>%
    mutate_at(vars(lang1, caregiver_seat, preterm, curr_earinfection, hearing_vision, cognitive_developmental,
                   beard), tolower) %>%
    mutate(lang_group = case_when(
            lang_group %in%  c("monolingual", "monilingual", "monolingual, not english", "monolinugal") ~ "monolingual",
            TRUE ~ lang_group),
         lang1 = case_when(
            str_detect(lang1, "english") ~ "english",
            TRUE ~ lang1),
         caregiver_seat = case_when(
           str_detect(caregiver_seat, "caregiver") ~ "caregiver",
           TRUE ~ "seat"),
         preterm = case_when(
           preterm %in% c("preterm", "y") ~ "preterm",
           preterm %in% c("full", "n", "0") ~ "full"),
         curr_earinfection = case_when(
           curr_earinfection == "na" ~ NA_character_,
           str_detect(curr_earinfection, "n") | str_detect(curr_earinfection, "0")  ~ "no",
           str_detect(curr_earinfection, "y") | str_detect(curr_earinfection, "1.0")  ~ "yes" ,
           TRUE ~ NA_character_),
         hearing_vision = case_when(
           hearing_vision == "na" ~ NA_character_,
           str_detect(hearing_vision, "n")  ~ "no",
           str_detect(hearing_vision, "y")   ~ "yes" ,
           TRUE ~ NA_character_),
         cognitive_developmental = case_when(
           str_detect(cognitive_developmental, "n")  ~ "no",
           str_detect(cognitive_developmental, "y")   ~ "yes" ,
           TRUE ~ NA_character_),
         beard  = case_when(
           str_detect(beard, "n")  ~ "no",
           str_detect(beard, "beard") ~ "yes",
           is.na(beard) ~ NA_character_)
         ) %>%
  mutate_if(is.character, as.factor)


# calculate es by "study"
d_var_calc <- function(n, d) {
  (2/n) + (d ^ 2 / (4 * n))
}

es_by_study <- mb_data_tidy_fct %>%
  group_by(lab, age_group, subid_unique) %>%
  summarise(d = mean(diff, na.rm = TRUE)) %>%
  group_by(lab, age_group) %>%
  summarise(d_z = mean(d)/ sd(d),
            n = n(),
            d_z_var = d_var_calc(n, d_z))

# get study characteristics
study_moderators <-  mb_data_tidy_fct %>%
  group_by(lab, age_group) %>%
  nest() %>%
  mutate(method = map(data, ~unique(.$method)),
         mean_age =  map(data, ~mean(.$age_days, na.rm = T)),
         prop_monolingual = map(data, ~sum(.$lang_group == "monolingual")/nrow(.[!is.na("lab_group"),])),
         modal_lang1 = map(data, ~count(., lang1) %>% arrange(-n) %>% slice(1) %>% pull(lang1)),
         mean_lang1_exposure = map(data, ~mean(.$lang1_exposure, na.rm = T)),
         parent_education = map(data, ~mean(c(mean(.$parenta_education, na.rm = T), mean(.$parentb_education, na.rm = T)))),
         prop_caregiver_seat = map(data, ~sum(.$caregiver_seat == "caregiver", na.rm = T)/nrow(.[!is.na("caregiver"),])),
         prop_preterm = map(data, ~sum(.$preterm == "preterm")/nrow(.[!is.na("preterm"),])),
         prop_curr_earinfection = map(data, ~sum(.$curr_earinfection == "yes", na.rm = T)/nrow(.[!is.na("curr_earinfection"),])),
         prop_hearing_vision = map(data, ~sum(.$hearing_vision == "yes", na.rm = T)/nrow(.[!is.na("hearing_vision"),])),
         prop_cognitive_developmental = map(data, ~sum(.$cognitive_developmental == "yes",na.rm = T)/nrow(.[!is.na("cognitive_developmental"),])),
         prop_female = map(data, ~sum(.$gender == "F", na.rm = T)/nrow(.[!is.na("gender"),])),
         prop_beard = map(data, ~sum(.$beard == "yes", na.rm = T)/nrow(.[!is.na("beard"),])),
         prop_nae = map(data, ~sum(.$nae, na.rm = T)/nrow(.[!is.na("nae"),]))) %>%
  select(-data) %>%
  unnest()

mb_data <- full_join(es_by_study, study_moderators)

#add methodological variables
methodological_vars <- read_csv(here("data/mb_methodological_variables.csv"))

mb_data <- full_join(mb_data, methodological_vars)

write_csv(mb_data, MB_OUT_PATH)
