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
          es_method  <- "group_means_between"

        } else if (complete(d)) {
          d_calc <- d
          es_method  <- "d_between"
        }

        #effect size var calculation
        if (complete(n_1, n_2, d_calc)) {
          d_var_calc <- ((n_1 + n_2) / (n_1 * n_2)) + (d_calc ^ 2 / (2 * (n_1 + n_2)))
        } else if (complete(d_var)) {
          d_var_calc <- d_var
        }

    } else if (participant_design == "within_two") {

        #effect size calculation
       if (complete(x_1, x_2, SD_1, SD_2)) {
            pooled_SD <- sqrt((SD_1 ^ 2 + SD_2 ^ 2) / 2) # Lipsey & Wilson (2001)
            d_calc <- (x_1 - x_2) / pooled_SD # Lipsey & Wilson (2001)
            es_method  <- "group_means_two"

       } else if (complete(d)) {
         d_calc <- d
         es_method  <- "d_two"
       }

      #effect size var calculation
      if (complete(n_1,  d_calc)) {
        d_var_calc <-   (2/n_1) + (d_calc ^ 2 / (4 * n_1))

      } else if (complete(d_var)){
        d_var_calc <- d_var
      }

    }

  data_frame("d_calc" = d_calc,
             "d_var_calc" = d_var_calc,
             "es_method" = es_method)

}
