# calculate ES for IDS MA

library(tidyverse)
library(assertthat)

complete <- function(...) {
  args = list(...)
  !any(unlist(map(args, ~(is.null(.x) || is.na(.x)))))
}

compute_es <- function(ma_df) {

  participant_design <- ma_df$participant_design
  x_1 <- ma_df$x_1
  x_2 <- ma_df$x_2
  SD_1 <- ma_df$sd_1
  SD_2 <- ma_df$sd_2
  n_1 <- ma_df$n_1
  n_2 <- ma_df$n_2
  d <- ma_df$d
  d_var <- ma_df$d_var

  assert_that(participant_design %in% c("between", "within_two"))

  d_calc <- NA
  d_var_calc <- NA
  es_method <- "missing"

  #effect size calculation
 if (participant_design == "between" & complete(x_1, x_2, n_1, n_2, SD_1, SD_2)) {
    pooled_SD <- sqrt((((n_1 - 1) * SD_1 ^ 2) + ((n_2 - 1) * SD_2 ^ 2)) / (n_1 + n_2 - 2)) # "classic" Cohen's d with marginal SD
    d_calc <- (x_1 - x_2) / pooled_SD
    d_var_calc <- ((n_1 + n_2) / (n_1 * n_2)) + (d_calc ^ 2 / (2 * (n_1 + n_2)))
    es_method  <- "classic_cohen_d_between"

  } else if (complete(d, d_var)) { # if raw data are not available do we want to use d already reported (using who knows what method)? I'd vote for no.
    d_calc <- d
    d_var_calc <- d_var
    es_method  <- "reported_d"
  } else if (participant_design == "within_two" & complete(x_1, x_2, n_1, SD_1, SD_2)) {
    pooled_SD <- sqrt((SD_1 ^ 2) + (SD_2 ^ 2) / 2) 
    d_calc <- (x_1 - x_2) / pooled_SD
    d_var_calc <- ((2 / n_1)  + (d_calc ^ 2 / (2 * n_1)))
    es_method  <- "classic_cohen_d_within"
    
  } 
  

  data_frame("d_calc" = d_calc,
             "d_var_calc" = d_var_calc,
             "es_method" = es_method)

}
