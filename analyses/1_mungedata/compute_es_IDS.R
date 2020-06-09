#This script calculates effect sizes for entries in metalab
#metalab.stanford.edu
#for questions please contact team members on metalab website

library(dplyr)
library(purrr)
library(assertthat)

complete <- function(...) {
  args = list(...)
  !any(unlist(map(args, ~(is.null(.x) || is.na(.x)))))
}


compute_es <- function(ma_df) {

  participant_design <- ma_df$participant_design
  x_1 <- ma_df$x_1
  x_2 <- ma_df$x_2
  x_dif <- ma_df$x_dif
  SD_1 = ma_df$sd_1
  SD_2 = ma_df$sd_2
  SD_dif = ma_df$SD_dif
  n_1 =ma_df$n_1
  n_2 = ma_df$n_2
  t = ma_df$t
  f = ma_df$f
  d = ma_df$d
  d_var = ma_df$d_var
  corr = ma_df$corr
  r = ma_df$xr2
  r_var = ma_df$r_var


    assert_that(participant_design %in% c("between", "within_two"))

    #we introduce variables calles d_calc and d_var_calc to distiguish them from the fields d and d_var, which are fields where effect sizes were already available from the source of the data
    d_calc <- NA
    d_var_calc <- NA
    es_method <- "missing"

    if (participant_design == "between") {

      #effect size calculation
      if (complete(x_1, x_2, SD_1, SD_2)) {
        pooled_SD <- sqrt(((n_1 - 1) * SD_1 ^ 2 + (n_2 - 1) * SD_2 ^ 2) / (n_1 + n_2 - 2)) # Lipsey & Wilson, 3.14
        d_calc <- (x_1 - x_2) / pooled_SD # Lipsey & Wilson (2001)
      } else if (complete(t)) {
        d_calc <- t * sqrt((n_1 + n_2) / (n_1 * n_2)) # Lipsey & Wilson, (2001)
      } else if (complete(f)) {
        d_calc <- sqrt(f * (n_1 + n_2) / (n_1 * n_2)) # Lipsey & Wilson, (2001)
      } else if (complete(r)) {
        d_calc <- 2 * r / sqrt(1 - r ^ 2)
      }
      if (complete(n_1, n_2, d_calc)) {
        #now that effect size are calculated, effect size variance is calculated
        d_var_calc <- ((n_1 + n_2) / (n_1 * n_2)) + (d_calc ^ 2 / (2 * (n_1 + n_2)))
      } else if (complete(r, r_var)) {
        #if r instead of d is reported, transform for standardization
        d_var_calc <- 4 * r_var / ((1 - r ^ 2) ^ 3)
      } else if (complete(d, d_var)) {
        #if d and d_var were already reported, use those values
        d_calc <- d
        d_var_calc <- d_var
      }

      es_method  <- "between"

    } else if (participant_design == "within_two") {
      if (is.na(corr) & complete(x_1, x_2, SD_1, SD_2, t)) {
        # Use raw means, SD, and t-values to calculate correlations
        corr <- (SD_1^2 + SD_2^2 - (n_1 * (x_1 - x_2)^2 / t^2)) / (2 * SD_1 * SD_2)
      }
      if (is.na(corr) | corr > .99 | corr < .01){
        #if correlation between two measures is not reported, use an imputed correlation value
        #we also account for the observation that some re-calculated values are impossible and replace those
        corr <- .5
      }

      #effect size calculation
      if (complete(x_1, x_2, SD_1, SD_2)) {
        pooled_SD <- sqrt((SD_1 ^ 2 + SD_2 ^ 2) / 2) # Lipsey & Wilson (2001)
        d_calc <- (x_1 - x_2) / pooled_SD # Lipsey & Wilson (2001)
        es_method  <- "group_means_two"
      } else if (complete(x_1, x_2, SD_dif)) {
        within_SD <- SD_dif / sqrt(2 * (1 - corr)) # Lipsey & Wilson (2001); Morris & DeShon (2002)
        d_calc <- (x_1 - x_2) / within_SD # Lipsey & Wilson (2001)
        es_method  <- "group_means_two"
      } else if (complete(t)) {
        wc <- sqrt(2 * (1 - corr))
        d_calc <- (t / sqrt(n_1)) * wc #Dunlap et al., 1996, p.171
        es_method  <- "t_two"
      }

      if (complete(n_1, d_calc)) {
        #now that effect size are calculated, effect size variance is calculated
        #d_var_calc <- ((1 / n_1) + (d_calc ^ 2 / (2 * n_1))) * 2 * (1 - corr) #we used this until 4/7/17
      #  d_var_calc <- (2 * (1 - corr)/ n_1) + (d_calc ^ 2 / (2 * n_1)) # Lipsey & Wilson (2001)
        d_var_calc <-   (2/n_1) + (d_calc ^ 2 / (4 * n_1)) # This is what is used in the MB data
      } else if (complete(d, d_var)) {
        #if d and d_var were already reported, use those values
        d_calc <- d
        d_var_calc <- d_var
        es_method  <- "d_two"
      }

    }

  data_frame("d_calc" = d_calc,
             "d_var_calc" = d_var_calc,
             "es_method" = es_method)

}
