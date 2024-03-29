# The final ManyBabies1 dataset from the public github repository, plus the processing pipeline from the paper
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               here)

# script expects these global vars from master: n_trial_pairs_criterion, age.matched

MB1_PATH <- "https://raw.githubusercontent.com/manybabies/mb1-analysis-public/master/processed_data/03_data_diff_main.csv"
ma_file_outname <- paste("mb_data_tidy_", n_trial_pairs_criterion/8, ".csv", sep = "")

if (ma_version == "Dunst_original") {
  data.dir = here("data","prepped_with_original_dunst")
} else if (ma_version == "Dunst_corrected") {
  data.dir = here("data","prepped_with_corrected_dunst")
} else {
  data.dir = here("data",paste0("prepped_with_",ma_version))
}

MB_OUT_PATH = here(data.dir,ma_file_outname)

TARGET_VARS <- c("lab", "subid_unique", "trial_num", "method", "age_days", "age_group",
                 "lang_group", "lang1", "lang1_exposure",
                 "preterm", "curr_earinfection", "hearing_vision", "cognitive_developmental",
                 "nae", "ADS", "IDS", "diff")

mb_data_raw <- read_csv(MB1_PATH)

mb_data_tidy <- mb_data_raw %>%
  #filter(!is.na(diff)) %>% #Removing this line because it is not in the MB scripts
  select(all_of(TARGET_VARS))

# tidy factors
mb_data_tidy_fct <- mb_data_tidy %>%
  mutate_at(vars(lang1, 
                 preterm, curr_earinfection, hearing_vision, cognitive_developmental
  ), tolower) %>%
  mutate(lang_group = case_when(
    lang_group %in%  c("monolingual", "monilingual", "monolingual, not english", "monolinugal") ~ "monolingual",
    TRUE ~ lang_group),
    lang1 = case_when(
      str_detect(lang1, "english") ~ "english", 
      TRUE ~ lang1),
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
      TRUE ~ NA_character_)
  ) %>%
  mutate_if(is.character, as.factor)


# calculate es by "study"
d_var_calc <- function(n, d) {
  (2/n) + (d ^ 2 / (4 * n))
}

# ManyBabies github code version for reference
# Source: https://github.com/manybabies/mb1-analysis-public/blob/master/paper/mb1-paper.Rmd
# l 537f
 
#es_by_study <- mb_data_tidy_fct %>%
#  group_by(lab, age_group, method, nae, subid_unique) %>%
#  summarise(d = mean(diff, na.rm = TRUE)) %>%
#  group_by(lab, age_group, method, nae) %>%
#  summarise(d_z = mean(d, na.rm = TRUE) / sd(d, na.rm = TRUE), 
#            n = length(unique(subid_unique)), 
#            d_z_var = d_var_calc(n, d_z)) %>%
#  filter(n >= 10) %>%
#  # left_join(ages) %>%
#  filter(!is.na(d_z)) 

# V 1 by Molly et al - with additional inclusion criterion update
es_by_participant <- mb_data_tidy_fct %>%
  group_by(lab, age_group, method,subid_unique) %>%
  filter(!is.na(diff)) %>% # MZ/CC: first filter diff NAs to avoid multiple trial pair counting issue
  summarise(d = mean(diff, na.rm = TRUE), n_trial_pairs = n()) %>%
  filter(n_trial_pairs >= n_trial_pairs_criterion)

es_by_study <- es_by_participant %>%
  group_by(lab, age_group, method) %>%
  summarise(d_z = mean(d, na.rm=TRUE)/ sd(d, na.rm=TRUE),
            n = n(),
            d_z_var = d_var_calc(n, d_z)) %>%
   filter(n>9) %>% # Match ManyBabies1 dataset by adding this inclusion criterion
   filter(!is.na(d_z)) # MZ: I think this shouldn't be needed but keeping to be safe

# get study characteristics
study_moderators <-  mb_data_tidy_fct %>%
  group_by(lab, age_group) %>%
  nest() %>%
  mutate(method = map(data, ~unique(.$method)),
         mean_age =  map(data, ~mean(.$age_days, na.rm = T)),
         prop_monolingual = map(data, ~sum(.$lang_group == "monolingual")/nrow(.[!is.na("lang_group"),])),
         modal_lang1 = map(data, ~count(., lang1) %>% arrange(-n) %>% slice(1) %>% pull(lang1)),
         mean_lang1_exposure = map(data, ~mean(.$lang1_exposure, na.rm = T)),
         prop_preterm = map(data, ~sum(.$preterm == "preterm")/nrow(.[!is.na("preterm"),])),
         prop_curr_earinfection = map(data, ~sum(.$curr_earinfection == "yes", na.rm = T)/nrow(.[!is.na("curr_earinfection"),])),
         prop_hearing_vision = map(data, ~sum(.$hearing_vision == "yes", na.rm = T)/nrow(.[!is.na("hearing_vision"),])),
         prop_cognitive_developmental = map(data, ~sum(.$cognitive_developmental == "yes",na.rm = T)/nrow(.[!is.na("cognitive_developmental"),])),
         prop_nae = map(data, ~sum(.$nae, na.rm = T)/nrow(.[!is.na("nae"),]))) %>%
  select(-data) %>%
  unnest(cols = c(method, mean_age, prop_monolingual, modal_lang1, mean_lang1_exposure, 
                  prop_preterm, prop_curr_earinfection, prop_hearing_vision, 
                  prop_cognitive_developmental, prop_nae))

mb_data <- full_join(es_by_study, study_moderators)

#add methodological variables
methodological_vars <- read_csv(here("data/mb_methodological_variables.csv"))

mb_data <- full_join(mb_data, methodological_vars)

if (age.matched == TRUE) {
  
  mb_data <- mb_data %>%
    unite(unique_study_id,lab,method,age_group,sep="_",remove=FALSE)
  
  # in main analysis, MB subjects were not well-matched with MA subjects
  
  # check mean age in MA: 248 days (median 186 days)
  setwd(data.dir)
  dma_original <- read_csv("ma_data_tidy.csv")
  summary(dma_original$mean_age)
  
  #check mean age of subset MA data: 216 days (median 186 days)
  dma_mb_ages <- read_csv("ma_data_tidy_mb_ages.csv")
  summary(dma_mb_ages$mean_age)
  #select subset MA restricted to MB age groups (3-15 months)
  dma <- dma_mb_ages
  
  # without further restriction, mean age in MB is ~280-290 days
  summary(mb_data$mean_age)
  
  # CC: proposal to carry out age matching; sample proportions from each age group.
  #calculate proportion of subjects in MA:
  proportion_subjects_in_age_groups_MA <- dma %>%
    group_by(age_group) %>%
    dplyr::summarise(n = sum(n), nprop = sum(n)/sum(dma$n))
  
  proportion_subjects_in_age_groups_MB <- mb_data %>%
    ungroup() %>%
    filter(!is.na(d_z)) %>%
    mutate(total_n = sum(n,na.rm=TRUE)) %>%
    group_by(age_group,total_n) %>%
    dplyr::summarise(n = sum(n,na.rm=TRUE)) %>%
    mutate(nprop = n/total_n)
  
  #sample MB dataset to include all subjects between 3 and 6 months and a proportional number
  #of participants (to the MA) from infants between 6 and 9 months and 9 to 12 months.
  
  sample_studyids_3_6 <- mb_data %>%
    filter(!is.na(d_z)) %>%
    filter(age_group=="3-6 mo") %>%
    ungroup() %>%
    distinct(unique_study_id)
  
  set.seed(40)
  
  sample_studyids_6_9 <- mb_data %>%
    ungroup() %>%
    filter(!is.na(d_z)) %>%
    filter(age_group=="6-9 mo") %>%
    slice_sample(
      n=round(filter(proportion_subjects_in_age_groups_MA,age_group=="6-9 mo")$nprop/filter(proportion_subjects_in_age_groups_MA,age_group=="3-6 mo")$nprop*nrow(sample_studyids_3_6))) %>%
    distinct(unique_study_id)
  
  #current studyids to use for subsample
  sample_studyids <- bind_rows(sample_studyids_3_6,sample_studyids_6_9)

  
  # check whether the MA has a 9-12 month age group
  # if yes sample 9-12 month olds proportionally
  if ("9-12 mo" %in% proportion_subjects_in_age_groups_MA$age_group) {
    sample_studyids_9_12 <- mb_data %>%
      ungroup() %>%
      filter(!is.na(d_z)) %>%
      filter(age_group=="9-12 mo") %>%
      slice_sample(
        n=round(filter(proportion_subjects_in_age_groups_MA,age_group=="9-12 mo")$nprop/filter(proportion_subjects_in_age_groups_MA,age_group=="3-6 mo")$nprop*nrow(sample_studyids_3_6))) %>%
      distinct(unique_study_id)
    #subids to use for subsample
    sample_studyids <- bind_rows(sample_studyids,sample_studyids_9_12)
  }
  
  # check whether the MA has a 12-15 month age group
  # if yes sample 12-15 month olds proportionally
  if ("12-15 mo" %in% proportion_subjects_in_age_groups_MA$age_group) {
    sample_studyids_12_15 <- mb_data %>%
      ungroup() %>%
      filter(!is.na(d_z)) %>%
      filter(age_group=="12-15 mo") %>%
      slice_sample(
        n=round(filter(proportion_subjects_in_age_groups_MA,age_group=="12-15 mo")$nprop/filter(proportion_subjects_in_age_groups_MA,age_group=="3-6 mo")$nprop*nrow(sample_studyids_3_6))) %>%
      distinct(unique_study_id)
    #subids to use for subsample
    sample_studyids <- bind_rows(sample_studyids,sample_studyids_12_15)
  }
  
  #select final age-matched sample
  sample_ages_full <- mb_data %>%
    filter(unique_study_id %in% sample_studyids$unique_study_id)
  
  summary(sample_ages_full$mean_age)
  summary(dma$mean_age)
  
  #show overlapping age plot - much better match than the "cutoff method" previously used
  #the MB looks a bit more "peak-y" mainly because it does not have any participants in the 0-3 mo age range (which affects the peak of the density)
  age_matching_plot <- ggplot() +
    geom_density(data = dma, aes(x = mean_age,weight=n/sum(n), fill = "#FC4E07"), show.legend = T, 
                 alpha = 0.8, color = "black") +
    geom_density(data = sample_ages_full, aes(x = mean_age,weight=n/sum(n), fill = "steelblue"), 
                 alpha = 0.8, color = "black") +
    xlim(c(-50, 500)) +
    xlab('Mean Age in Days') +
    ylab('Density') +
    scale_fill_manual(name = "", labels = c("Meta-analysis", "ManyBabies"), 
                      values = c("#FC4E07", "steelblue")) +
    ggtitle('Age-Matched Samples across MA and MB') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size=15),
          legend.position = "right",
          axis.text.x = element_text(size = 13),
          axis.title.x = element_text(size = 13),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 13))
  
  #save plot for supplementary materials:
  ggsave(plot = age_matching_plot, file = here('misc_figures',paste0("age_matching_plot_","prepped_with_",ma_version,'.png')), height = 5, width = 8)
  
  mb_data <- sample_ages_full
  
  # retitle the dataset
  MB_OUT_PATH <- paste(data.dir, "/mb_data_tidy_", n_trial_pairs_criterion/8, "_age_matched.csv", sep = "")
}

write_csv(mb_data, MB_OUT_PATH)




