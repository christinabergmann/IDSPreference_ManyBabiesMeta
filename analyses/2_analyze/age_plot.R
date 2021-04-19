library(tidyverse)
library(here)
library(robumeta)
library(clubSandwich)

d <- read_csv(here("data/mb_ma_combined_prepped_0.125.csv"))

moderators = c( "isMeta",  # code this way since we expect meta to have larger effect sizes
           "mean_agec_mos",
           "test_lang",  # whether stimuli were in native language; almost constant in meta
           "method",
           
           # constant in RRR:
           "speech_type",
           "own_mother",
           "presentation",
           "dependent_measure",  # causes singularity
           "main_question_ids_preference" )

naive_model <- robu(yi ~ isMeta,
                        data = d, 
                        studynum = as.factor(study_id),
                        var.eff.size = vi,
                        modelweights = "HIER",
                        small = TRUE)


moderated_model <- robu(yi ~ isMeta + mean_agec_mos + test_lang + method,
      data = d, 
      studynum = as.factor(study_id),
      var.eff.size = vi,
      modelweights = "HIER",
      small = TRUE)

pred <- predict(moderated_model, 
                pred.vector = c(1, 0, 0, 0, 0, 0, 0), level = .95)

ggplot(d, 
       aes(x = mean_age, y = d_calc, col = study_type)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, lty = 2) + 
  geom_hline(yintercept = 0, lty = 2, col = "black") + 
  langcog::theme_mikabr() + 
  langcog::scale_color_solarized(name = "Study Type") + 
  theme(legend.position = "bottom") + 
  ylab("Effect size (d)") +
  xlab("Age (days)")
