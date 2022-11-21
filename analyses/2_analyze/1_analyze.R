
# META-NOTES ------------------------------------------------------------------ 


# ~ Usage notes  ------------------------------------------------------------------

# - Key points about the assumptions and interpretation of various fn outputs are marked with "**" in this file and in analyze_helper.R. If you make changes to fns themselves or to how variables are coded, make sure you read those to make sure things won't break.

# - Code writes graphics, tables, and stats_for_paper.csv straight to relative directories (results.dir) and also to a fixed directory on MBM's machine
#  (overleaf.dir) that pipes into the LaTeX manuscript.

# - Code also writes numerical results to stats_for_paper.csv in a format that can be piped into LaTeX. To see the results as you run the code, use vr(). You can wipe that file using wr().

# - Names of important model objects:
# - naive.MA.only and naive.reps.only: meta-analyses within subsets; no moderators
# - cond.MA.only and cond.reps.only: meta-analyses within subsets; 
#   conditions moderators to averages in meta-analysis

# - Remember not to use asterisks in file names bc they cause syncing trouble for CB.

# - Acronyms: "MB" or "MLR" refer to the replications. "MA" refers to meta-analysis.




# 0. PRELIMINARIES ------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")

# This script uses renv to preserve the R environment specs (e.g., package versions.)
pacman::p_load(renv)
# run this if you want to reproduce results using the R environment we had:
# library(here); setwd(here)
# renv::restore()
pacman::p_load(
  tidyverse,
  knitr,
  data.table,
  here,
  tableone,
  corrr,
  robumeta,
  fastDummies,
  weightr,
  PublicationBias,
  xtable,
  boot,
  testthat,
  ggplot2,
  metafor,
  MatchIt,
  table1)

# run this only if you want to update the R environment specs
# library(here); setwd(here())
# renv::snapshot()

if ( use.corrected.dunst == FALSE ) {
  data.dir = here("data/prepped_with_original_dunst")
  # where to save results
  results.dir = here("results_from_R/results_with_original_dunst")
  overleaf.dir = "~/Dropbox/Apps/Overleaf/MB-Meta/R_objects"
}

if ( use.corrected.dunst == TRUE ) {
  data.dir = here("data/prepped_with_corrected_dunst")
  # where to save results
  results.dir = here("results_from_R/results_with_corrected_dunst")
  overleaf.dir = "~/Dropbox/Apps/Overleaf/MB-Meta/R_objects/corrected_dunst"
}

if(!dir.exists(results.dir)){
  dir.create(results.dir)
} 


code.dir = here("analyses/2_analyze")

# helper fns
setwd(code.dir)
source("analyze_helper.R")

# package update not yet on CRAN, but we need the cluster-bootstrap functionality
source("MetaUtility development functions.R")


# ~ Code-Running Parameters ------------------------------------------------------------------
# should we remove existing results file instead of overwriting individual entries? 
start.res.from.scratch = TRUE
# should we use the grateful package to scan and cite packages?
cite.packages.anew = FALSE
# should we bootstrap from scratch or read in old resamples?
boot.from.scratch = TRUE
# make plots from scratch?
redo.plots = TRUE

# if (redo.mod.selection == FALSE) {
#   # read in the surviving moderators
#   setwd(results.dir)
#   modsS = read.csv("surviving_mods.csv")[,2]
# }

# wipe results csvs if needed
if ( start.res.from.scratch == TRUE ) wr()

# ~ Constants of Universe ------------------------------------------------------------------
digits = 2
pval.cutoff = 10^-4  # threshold for using "<"
boot.reps = 1000 # for all bootstrapped inference 
excelMax = 8
# max digits in Excel; only used for running sanity checks against the results that are auto-written to a csv file

# plot colors
colors = c("darkgray", "red2")  # replications, originals

# stop sci notation
options(scipen=999)

# moderator names
# as used in table
mods = c( "study_type",
          "mean_agec_mos",
          "test_lang",  # whether stimuli were in native language; almost constant in meta
          "native_lang", #Won't be used in main analyses
          "method",
          
          # constant in RRR:
          "speech_type",
          "own_mother",
          "presentation",
          "dependent_measure",
          "main_question_ids_preference" )

# and as used in meta-regression (prior to selection to get model to converge)
# same list as mods except has isMeta instead of study_type
#  and does not have native_lang
#  also used in matching
mods2 = c( "isMeta",  # code this way since we expect meta to have larger effect sizes
           "mean_agec_mos",
           "test_lang",  # whether stimuli were in native language; almost constant in meta
           "method",
           
           # constant in RRR:
           "speech_type",
           "own_mother",
           "presentation",
           "dependent_measure",  # causes singularity
           "main_question_ids_preference" )

# ~ Read Datasets ------------------------------------------------------------------
setwd(data.dir)
d = suppressMessages( suppressWarnings( read_csv("mb_ma_combined_prepped_0.125.csv") ) ) %>%
  filter(!is.na(d_calc))

# dataset with just the meta-analysis
dma = d %>% filter(isMeta == TRUE)

# dataset with just the replications
dr = d %>% filter(isMeta == FALSE)

# dataset with more stringent inclusion criteria in replications
dic = suppressMessages( suppressWarnings( read_csv("mb_ma_combined_prepped_0.75.csv") ) ) %>%
  filter(!is.na(d_calc))
dric = dic %>% filter(isMeta == FALSE)%>%
  filter(!is.na(d_calc))

# dataset with IPD age-matching in MB
dage = suppressMessages( suppressWarnings( read_csv("mb_ma_combined_prepped_0.125_age_matched.csv") ) ) %>%
  filter(!is.na(d_calc))
drage = dage %>% filter(isMeta == FALSE)%>%
  filter(!is.na(d_calc))



# 1. CHARACTERISTICS OF INCLUDED STUDIES ------------------------------------------------------------------           


# for updating result csv
section = 1

# ~ Basics ------------------------------------------------------------------

t = d %>% group_by(study_type) %>%
  summarise( m = n(),
             k = length(unique(study_id) ),
             nSubjMed = median(n),
             nSubjMin = min(n),
             nSubjMax = max(n), )
t


update_result_csv( name = paste( "m ests", t$study_type, sep = " "),
                   value = t$m )

update_result_csv( name = "m ests total",
                   value = sum(t$m) )

update_result_csv( name = paste( "k studies", t$study_type, sep = " "),
                   value = t$k )

update_result_csv( name = paste( "median n", t$study_type, sep = " "),
                   value = t$nSubjMed )

update_result_csv( name = paste( "min n", t$study_type, sep = " "),
                   value = t$nSubjMin )

update_result_csv( name = paste( "max n", t$study_type, sep = " "),
                   value = t$nSubjMax )


# ~ Moderators ------------------------------------------------------------------

# distribution of moderators in RRR and MA
# t = CreateTableOne(vars = mods, 
#                    strata = "study_type",
#                    data = d)

#xtable( print(t, noSpaces = TRUE, printToggle = FALSE) )

# Helpful for table https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html

# clean up variable names and values for nicer presentation in Table 1
t <- d %>%
  select(study_type, mods) %>%
  mutate(study_type = as.factor(study_type)) %>%
  mutate(`centered age (months)` = mean_agec_mos) %>%
  mutate(`test language` = factor(test_lang, labels = c("native", "non-native", "artificial"))) %>%
  mutate(`native language` = factor(native_lang)) %>%
  mutate(method = factor(method, labels = c("central fixation", "headturn preference procedure", "other"))) %>%
  mutate(`speech type` = factor(speech_type, labels = c("simulated", "naturalistic", "filtered/synthesized"))) %>%
  mutate(`own mother` = factor( own_mother, labels = c("no", "yes"))) %>%
  mutate(presentation = factor(presentation, labels = c("tape recording", "video recording"))) %>%
  mutate(`main question was about IDS prerence` = factor(main_question_ids_preference, labels = c("yes", "no"))) %>%
  select(study_type, `centered age (months)`, `test language`, `native language`, method, `speech type`, `own mother`, presentation, `main question was about IDS prerence`)


my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(round(x,2)), digits=3), c("",
                                                                    "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}

table1(~ `centered age (months)` + `test language` + `native language` + method + `speech type` + `own mother` + presentation + `main question was about IDS prerence`| study_type, data = t, overall = F, render.continuous=my.render.cont)



# ~~ Moderator Correlation Matrix ------------------------------------------------------------------

# temporarily recode all categorical moderators as dummies
# dataset already has dummies, but this is easier than subsetting the right names
cat.mods = mods[ !mods == "mean_agec_mos" ]
temp = dummy_cols( d[,cat.mods],
                   select_columns = cat.mods,
                   remove_selected_columns = TRUE,
                   ignore_na =  TRUE )  
temp$mean_agec_mos = d$mean_agec_mos 
names(temp)

# make correlation matrix
corrs = temp %>%
  correlate( use = "pairwise.complete.obs" ) %>%
  stretch() %>%
  arrange(desc(r)) %>%
  group_by(r) %>%
  filter(row_number()==1)

# save it

setwd(results.dir)
write.csv(corrs, "moderator_cormat.csv")
# too long to put in paper




# 2. PRIMARY MODERATOR ANALYSES (PREREG'D) ------------------------------------------------------------------

section = 2

# ~ Subset meta-analyses for MA an MB ------------------------------------------------------------------
# this fn also writes stats to results csv
( naive.MA.only = fit_subset_meta( .dat = dma,
                                   .mods = "1",
                                   .label = "MA subset naive" ) )


( naive.reps.only = fit_subset_meta( .dat = dr,
                                     .mods = "1",
                                     .label = "Reps subset naive" ) )



# sanity check: refit one of these models manually
# resCSV is assigned as a global var by fns that call update_result_csv
if ( exists("resCSV") ) {
  
  temp = robu( yi ~ 1, 
               data = dma, 
               studynum = as.factor(study_id),
               var.eff.size = vi,
               modelweights = "HIER",
               small = TRUE)
  
  expect_equal( resCSV$value[ resCSV$name == "MA subset naive est X.Intercept." ],
                as.character( round(temp$b.r, 2) ) )
  
  expect_equal( resCSV$value[ resCSV$name == "MA subset naive lo X.Intercept." ],
                as.character( round(temp$reg_table$CI.L, 2) ) )
  
  expect_equal( resCSV$value[ resCSV$name == "MA subset naive hi X.Intercept." ],
                as.character( round(temp$reg_table$CI.U, 2) ) )
  
  expect_equal( resCSV$value[ resCSV$name == "MA subset naive pval X.Intercept." ],
                as.character( round(temp$reg_table$prob, excelMax) ) )
  
}




# ~ Both sources: Naive and moderated regressions ------------------------------------------------------------------ 

# cannot skip this section because needs naiveRes and modRes later

# order of moderator importance given in prereg:
# age, test_lang, method, speaker, speech_type, own_mother, presentation, DV, main question

# per prereg, if model isn't estimable, the moderators are to be removed in the opposite order
#  of this importance list

# list of moderators to attempt for each of 2 models
mod.sets = list( c("isMeta"),  # naive model
                 mods2 )  # moderated model

labels = c("naive",
           "moderated")

# fit the naive model
# fit_mr automatically writes results to the results csv file and table
# caps indicates it's a primary model
naiveRes = fit_mr( .dat = d,
                   .label = "NAIVE",
                   .mods = mod.sets[[1]],
                   .write.to.csv = TRUE,
                   .write.table = TRUE,
                   .simple.return = FALSE )

# sanity check: should be similar to naive meta-regression model
#  but not necessarily equal
expect_equal( as.numeric(naive.MA.only$b.r), round(naiveRes$est.ma, digits), tol = 0.03 )
expect_equal( as.numeric(naive.MA.only$reg_table$CI.L), round(naiveRes$est.ma.lo, digits), tol = 0.03 )
expect_equal( as.numeric(naive.MA.only$reg_table$CI.U), round(naiveRes$est.ma.hi, digits), tol = 0.04 )
#Increased tolerance from 0.03

# fit the meta-regression with all covariates 
#  and remove them in prespecified order if needed until 
#  model becomes estimables

section = 1  # loop will break otherwise

gotError = TRUE  # initialize so the while-loop is entered


# this will print an xtable with the moderator estimates for use in the paper
while ( gotError == TRUE ) {
  

  tryCatch({
    mod1Res = fit_mr( .dat = d,
                      .label = "MOD1",
                      .mods = mod.sets[[2]],
                      .write.to.csv = TRUE,
                      .write.table = TRUE,
                      .simple.return = FALSE )
    gotError = FALSE
    message("MODEL WAS ESTIMABLE. OH YASSSSS <3")
    
  }, error = function(err){
    gotError <<- TRUE
    
    # remove one moderator from the end of the list
    message( paste( "\n Removing ", 
                    mod.sets[[2]][ length(mod.sets[[2]]) ],
                    " from moderators and trying again" ) )
    mod.sets[[2]] <<- mod.sets[[2]][ 1 : ( length(mod.sets[[2]]) - 1 ) ]
    
    
  })
  
}

# above call will automatically print an xtable for the final mod1Res
# that table is in paper

# look at the surviving moderators and save for later use
( modsS = mod.sets[[2]] )

##### Sanity checks
# sanity check: refit one of these models manually
# resCSV is assigned as a global var by fns that call update_result_csv
if ( exists("resCSV") ) {
  
  # refit naive model
  temp = robu( yi ~ isMeta, 
               data = d, 
               studynum = as.factor(study_id),
               var.eff.size = vi,
               modelweights = "HIER",
               small = TRUE)
  
  expect_equal( resCSV$value[ resCSV$name == "NAIVE k" ],
                as.character( round( nrow(d), 0) ) )
  
  expect_equal( resCSV$value[ resCSV$name == "NAIVE n subj" ],
                as.character( round( sum(d$n), 0) ) )
  
  expect_equal( resCSV$value[ resCSV$name == "NAIVE tau" ],
                as.character( round( sqrt(temp$mod_info$tau.sq), 2) ) )
  
  # intercept estimate and inference
  expect_equal( resCSV$value[ resCSV$name == "NAIVE est X.Intercept." ],
                as.character( round(temp$b.r[1], 2) ) )
  
  expect_equal( resCSV$value[ resCSV$name == "NAIVE lo X.Intercept." ],
                as.character( round(temp$reg_table$CI.L[1], 2) ) )
  
  expect_equal( resCSV$value[ resCSV$name == "NAIVE hi X.Intercept." ],
                as.character( round(temp$reg_table$CI.U[1], 2) ) )
  
  expect_equal( resCSV$value[ resCSV$name == "NAIVE pval X.Intercept." ],
                format.pval(temp$reg_table$prob[1], eps = pval.cutoff) )              
  # moderator estimate and inference
  expect_equal( resCSV$value[ resCSV$name == "NAIVE est isMetaTRUE" ],
                as.character( round(temp$b.r[2], 2) ) )
  
  expect_equal( resCSV$value[ resCSV$name == "NAIVE lo isMetaTRUE" ],
                as.character( round(temp$reg_table$CI.L[2], 2) ) )
  
  expect_equal( resCSV$value[ resCSV$name == "NAIVE hi isMetaTRUE" ],
                as.character( round(temp$reg_table$CI.U[2], 2) ) )
  
  expect_equal( resCSV$value[ resCSV$name == "NAIVE pval isMetaTRUE" ],
                format.pval(temp$reg_table$prob[2], eps = pval.cutoff) )                   
  
  # average for meta-analysis obtained by recoding isMeta
  temp = robu( yi ~ isRep, 
               data = d, 
               studynum = as.factor(study_id),
               var.eff.size = vi,
               modelweights = "HIER",
               small = TRUE)
  
  expect_equal( resCSV$value[ resCSV$name == "NAIVE avg for meta" ],
                as.character( round(temp$b.r[1], 2) ) )
  
  expect_equal( resCSV$value[ resCSV$name == "NAIVE avg lo for meta" ],
                as.character( round(temp$reg_table$CI.L[1], 2) ) )
  
  # this one is in scientific notation; what a PITA
  format.pval(temp$reg_table$prob[1], eps = pval.cutoff)
  resCSV$value[ resCSV$name == "NAIVE avg pval for meta" ]
  
  # moderated model
  # for this one, maybe just spot-check some of the moderators?
  temp = robu( yi ~ isMeta + mean_agec_mos + test_lang + method, 
               data = d, 
               studynum = as.factor(study_id),
               var.eff.size = vi,
               modelweights = "HIER",
               small = TRUE)
  
  # for checking all the coeff estimates as a vector
  suffix = c("X.Intercept.", "isMetaTRUE", "mean_agec_mos", "test_langb.nonnative", "test_langc.artificial", "methodb.hpp", "methodc.other")
  
  expect_equal( resCSV$value[ resCSV$name %in% paste( "MOD1 est", suffix, sep = " ") ],
                as.character( round(temp$b.r, 2) ) )
  
  expect_equal( resCSV$value[ resCSV$name %in% paste( "MOD1 lo", suffix, sep = " ") ],
                as.character( round(temp$reg_table$CI.L, 2) ) )
  
  expect_equal( resCSV$value[ resCSV$name %in% paste( "MOD1 hi", suffix, sep = " ") ],
                as.character( round(temp$reg_table$CI.U, 2) ) )
  
  expect_equal( resCSV$value[ resCSV$name %in% paste( "MOD1 pval", suffix, sep = " ") ],
                format.pval(temp$reg_table$prob, eps = pval.cutoff) )  
  # oh yeahhhhhh
}


# 2. DENSITY PLOT OF META-ANALYSIS VS. REPLICATION CALIBRATED ESTIMATES ------------------------------------------------------------------


# ~ Get Marginal Calibrated Estimates for Plot ------------------------------------------------------------------

calibR = calib_ests(yi = dr$yi, sei = sqrt(dr$vi) )
mean(calibR>0.2)  # quick check
mean(calibR)

calibMA = calib_ests(yi = dma$yi, sei = sqrt(dma$vi) )
mean(calibMA>0.2)  # quick check
mean(calibMA)

d$calibNaive = NA
d$calibNaive[ d$isMeta == FALSE ] = calibR
d$calibNaive[ d$isMeta == TRUE ] = calibMA


# ~ Get Conditional Calibrated Estimates for Plot ------------------------------------------------------------------
# fit the final, moderated models to each subset
# to see what the heterogeneity looks like in each case
( cond.MA.only = fit_subset_meta( .dat = dma,
                                  .mods = modsS[ !modsS == "isMeta" ],
                                  .label = "MA subset mods" ) )

( cond.reps.only = fit_subset_meta( .dat = dr,
                                    .mods = modsS[ !modsS == "isMeta" ],
                                    .label = "Reps subset mods" ) )


# get conditional calibrated estimates
# this fn already has unit tests in analyze_helper.Rs
d$calibCond[ d$isMeta == FALSE ] = conditional_calib_ests(cond.reps.only)$calib.shift
d$calibCond[ d$isMeta == TRUE ] = conditional_calib_ests(cond.MA.only)$calib.shift

# ***important: MB heterogeneity estimate is 0 after conditioning on mods
#  but MA heterogeneity estimate actually slightly increases


if ( redo.plots == TRUE ) {
  
  ##### Marginal Calibrated Estimates
  # choose axis scaling
  summary(d$yi)
  xmin = -1
  xmax = 3.5
  tickJump = 0.5  # space between tick marks
  
  ggplot( data = d,
          aes( x = calibNaive,
               fill = studyTypePretty,
               color = studyTypePretty ) ) +
    
    # mean estimates from subset models
    geom_vline( xintercept = naive.MA.only$b.r,
                color = colors[2],
                lty = 2 ) +
    
    geom_vline( xintercept = naive.reps.only$b.r,
                color = colors[1],
                lty = 2) +
    
    # null
    geom_vline(xintercept = 0,
               color = "gray",
               lty = 1) +
    
    geom_density(alpha = 0.3) +
    
    theme_bw() +
    
    xlab("Marginal population SMDs") +
    scale_x_continuous( limits = c(xmin, xmax), breaks = seq(xmin, xmax, tickJump)) +
    
    ylab("Estimated density") +
    
    scale_color_manual( values = rev(colors), name = "Source" ) +
    scale_fill_manual( values = rev(colors), name = "Source" ) +
    
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  
  my_ggsave( name = "naive_calibrated_plot.pdf",
             width = 8,
             height = 5 )
  
  
  ##### Conditional Calibrated Estimates
  ggplot( data = d,
          aes( x = calibCond,
               fill = studyTypePretty,
               color = studyTypePretty ) ) +
    
    # conditional mean estimates from subset models
    geom_vline( xintercept = cond.MA.only$b.r[1],
                color = colors[2],
                lty = 2 ) +
    
    geom_vline( xintercept = cond.reps.only$b.r[1],
                color = colors[1],
                lty = 2) +
    
    # null
    geom_vline(xintercept = 0,
               color = "gray",
               lty = 1) +
    
    # ensemble estimates shifted to Z=0
    geom_density(alpha = 0.3) +
    
    theme_bw() +
    
    xlab("Conditional population SMDs") +
    scale_x_continuous( limits = c(xmin, xmax), breaks = seq(xmin, xmax, tickJump)) +
    
    ylab("Estimated density") +
    
    scale_color_manual( values = rev(colors), name = "Source" ) +
    scale_fill_manual( values = rev(colors), name = "Source" ) +
    
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  
  my_ggsave( name = "cond_calibrated_plot.pdf",
             width = 8,
             height = 5 )
  
  
}


# 3. NAIVE AND MODERATED MODEL INFERENCE------------------------------------------------------------------

# MM audited this section 2021-3-22

# get inference for conditional Phats and their difference for the moderated model

# structure of the code below:
# - for each resample, fits both naive model and moderated model
# - via fit_mr, gets 9 stats of interest for each model
# - returns results of naive model, moderated model, AND their difference
#    in a (3*9)-length vector

# here is the order of fit_mr's 9 returned stats:
# c( est.ma,
#    est.rep,
#    est.ma - est.rep,
#    
#    Phat0.ma,
#    Phat0.rep,
#    Phat0.diff,
#    
#    Phat0.2.ma,
#    Phat0.2.rep,
#    Phat0.2.diff )

# ~ Create the Resamples ------------------------------------------------------------------

# with boot.reps - 1,000, takes about 10 min
if ( boot.from.scratch == TRUE ) {
  boot.res = boot( data = d, 
                   parallel = "multicore",
                   R = boot.reps, 
                   statistic = function(original, indices) {
                     # draw resample with replacement
                     # ignore the indices passed by boot because we're
                     #  doing clustered bootstrapping, oh yeah
                     
                     b = cluster_bt(.dat = original,
                                    .clustervar = "study_id")
                     
                     # multi-argument returns need to be via c(), not list or df or whatever
                     tryCatch({
                       # **key part: fit the two models
                       # get their stats and also the cross-model differences
                       .naiveRes = fit_mr( .dat = b, .mods = "isMeta" )
                       .modRes = fit_mr( .dat = b, .mods = modsS )
                       .diff = .naiveRes - .modRes
                       
                       c(.naiveRes, .modRes, .diff)
                       
                     }, error = function(err){
                       # number of NA's here needs to match the match .simple.return = TRUE structure of fit_mr
                       # x 3 separate calls to fit_mr
                       return( rep(NA, 9 * 3) )
                     })
                     
                   } )
  
  # save bootstraps
  setwd(results.dir)
  base::save(boot.res, file = "saved_bootstraps.RData")
}

if ( boot.from.scratch == FALSE ){
  # read in existing bootstraps
  setwd(results.dir)
  base::load("saved_bootstraps.RData")
}



# ~ Process the Resamples ------------------------------------------------------------------

# number of non-failed bootstrap reps
( boot.reps.successful = sum( !is.na( boot.res$t[,1] ) ) )


# prep for BCa interval construction: get stats on original data in same format as returned by boot()
# "2" suffix is because we already have naiveRes and modRes from the model selection
#   part, but those objects use the full return structure rather than the boot-friendly one
naiveRes2 = fit_mr( .dat = d, .mods = "isMeta" )
modRes2 = fit_mr( .dat = d, .mods = modsS )
diff = naiveRes2 - modRes2
( t0 = c(naiveRes2, modRes2, diff) )

# **NB: because of the way we hacked the "statistic" argument of boot()
#  above st it actually creates a resample rather than only calculating a statistics,
#  the below will NOT match the actual estimates from original dataset
# the bootstrap "bias" estimates will also be wrong for this reason
# boot.res$t0

# sanity check: agreement between fit_mr's two
#  return structures
expect_equal( naiveRes$est.ma, naiveRes2[1] )
expect_equal( naiveRes$est.rep, naiveRes2[2] )
expect_equal( naiveRes$avgDiff, naiveRes2[3] )
expect_equal( naiveRes$Phat0.ma, naiveRes2[4] )
expect_equal( naiveRes$Phat0.rep, naiveRes2[5] )
expect_equal( naiveRes$Phat0.ma - naiveRes$Phat0.rep, naiveRes2[6] )
expect_equal( naiveRes$Phat0.2.ma, naiveRes2[7] )
expect_equal( naiveRes$Phat0.2.rep, naiveRes2[8] )
expect_equal( naiveRes$Phat0.2.ma - naiveRes$Phat0.2.rep, naiveRes2[9] )


# **critical step: replace boot()'s incorrect t0 with the correct one
#  to allow boot.ci to work correctly
boot.res$t0 = t0


# sanity check: bootstrap bias estimate
# should be reasonably close to 0
round( colMeans(boot.res$t, na.rm = TRUE) - t0, 2 )

# get the bootstrapped CIs for all 3 statistics respectively using BCa method
( CIs = get_boot_CIs(boot.res, n.ests = ncol(boot.res$t) ) )
# convert to df
res = as.data.frame( do.call( what = rbind, args = CIs ) )
names(res) = c("lo", "hi")

# add in sample point estimates and names of variables
res = res %>% add_column( est = t0, .before = 1 )
# from return structure of fit_mr:
# **must match return structure of fit_mr
varNames = c( "AvgM",
              "AvgR",
              "AvgDiff",
              
              "Phat0M",
              "Phat0R",
              "Phat0Diff",
              
              "Phat0.2M",
              "Phat0.2R",
              "Phat0.2Diff" )

# return structure of fit_mr for reference:
# c( est.ma,
#    est.rep,
#    est.ma - est.rep,
#    
#    Phat0.ma,
#    Phat0.rep,
#    Phat0.diff,
#    
#    Phat0.2.ma,
#    Phat0.2.rep,
#    Phat0.2.diff )


# 3 for number of calls to fit_mr
res = res %>% add_column( stat = rep( varNames, 3 ), .before = 1 )
res = res %>% add_column( model = rep( c("naive", "mod", "modelDiff" ), each = length(varNames) ),
                          .before = 1 )


# ~~ Clean Up CIs ------------------------------------------------------------------
# sanity check: look for point estimates that aren't in their CIs
res %>% filter( est > hi | est < lo )
# there is only one, and it makes sense b/c estimate was at ceiling

# set such CIs to NA
res$lo[ res$est < res$lo | res$est > res$hi ] = NA
res$hi[ res$est < res$lo | res$est > res$hi ] = NA

# if only one CI limit was estimable, set both to NA
res$lo[ is.na(res$hi) ] = NA
res$hi[ is.na(res$lo) ] = NA

# **all of the CIs in this df are bootstrapped

# # sanity check: compare model-based CIs (for coeff estimates) to boot ones
# # model-based ones are consistently a lot wider, actually
# res[ res$model == "naive" & res$stat == "AvgM", c("lo", "hi") ]; c( naiveRes$est.ma.lo, naiveRes$est.ma.hi )
# 
# res[ res$model == "naive" & res$stat == "AvgR", c("lo", "hi") ]; c( naiveRes$est.rep.lo, naiveRes$est.ma.hi )
# 
# res[ res$model == "naive" & res$stat == "AvgDiff", c("lo", "hi") ]; c( naiveRes$avgDiffLo, naiveRes$avgDiffHi )
# 
# # why are boot CIs so different from model-based ones?
# boot.res
# boot.ci(boot.res, index=1)  # all types of boot CIs are pretty similar
# c( naiveRes$est.ma.lo, naiveRes$est.ma.hi )  # model-based much wider
# # I went through the browser of fit_br and confirmed that above model-based ones
# #  were extracted correctly
# 
# # look at bootstraps themselves
# x = boot.res$t[,1]
# hist(x)
# mean(x, na.rm = TRUE); t0[1]
# # mean of bootstraps is definitely lower than the sample one, which could indicate bias
# # in MRM, we did see relative bias of 5% on average for conditional expectation
# # so this amount of bias is actually plausible


# **for coefficient estimates (i.e., not Phats), use the model-based CIs instead 
# for consistency with earlier table showing meta-regression estimates
res[ res$model == "naive" & res$stat == "AvgM", c("lo", "hi") ] = c( naiveRes$est.ma.lo, naiveRes$est.ma.hi )

res[ res$model == "naive" & res$stat == "AvgR", c("lo", "hi") ] = c( naiveRes$est.rep.lo, naiveRes$est.ma.hi )

res[ res$model == "naive" & res$stat == "AvgDiff", c("lo", "hi") ] = c( naiveRes$avgDiffLo, naiveRes$avgDiffHi )

res[ res$model == "mod" & res$stat == "AvgM", c("lo", "hi") ] = c( mod1Res$est.ma.lo, mod1Res$est.ma.hi )

res[ res$model == "mod" & res$stat == "AvgR", c("lo", "hi") ] = c( mod1Res$est.rep.lo, mod1Res$est.ma.hi )

res[ res$model == "mod" & res$stat == "AvgDiff", c("lo", "hi") ] = c( mod1Res$avgDiffLo, mod1Res$avgDiffHi )



# ~~ Prettify and Save Results Tables ------------------------------------------------------------------

# ~~~ Unrounded numeric table (just for reproducibility) ------------------------------------------------------------------
setwd(results.dir)
setwd("table_model_diffs")
fwrite( res, "table_model_diffs_unrounded.csv" )

# ~~~ Rounded numeric table (for piping numbers into manuscript) ------------------------------------------------------------------

res2 = res

# list of numeric variables
numVars = c("est", "lo", "hi")
res2 = res2 %>% mutate_if(is.numeric, function(x) round(x,2))

# turn Phats into percentages
inds = grepl( pattern = "Phat", x = res2$stat )
res2[ inds, numVars ] = 100*res2[ inds, numVars ]
res2$stat[ inds ] = str_replace( string = res2$stat[ inds ],
                                 pattern = "Phat",
                                 replacement = "Perc" )

# make single key column for use with rlgetrum in TeX
res2 = res2 %>% add_column( unique = paste( res2$model, res2$stat, sep = "_") )

fwrite( res2, "table_model_diffs_rounded.csv" )

# also save to Overleaf
if ( dir.exists(overleaf.dir) ) {
  setwd(overleaf.dir)
  fwrite( res2, "table_model_diffs_rounded.csv" )
}


# ~~~ Rounded character table (for manuscript table) ------------------------------------------------------------------
res3 = res2
res3 = res3 %>% mutate_at( numVars, as.character )

# combine estimates and CIs into single string
res3$valueString = paste( res3$est,
                          " [",
                          res3$lo, 
                          ", ", 
                          res3$hi, 
                          "]",
                          sep = "" ) 
# if CI is [NA, NA], delete it from string
res3$valueString = str_replace( string = res3$valueString,
                                pattern = " \\[NA, NA\\]",
                                replacement = "" )

res3 = res3 %>% select(-numVars)
res3 = res3 %>% select(-unique)


# reshape table for readability
res3 = res3 %>% pivot_wider( names_from = model,
                             values_from = valueString )

setwd(results.dir)
setwd("table_model_diffs")
fwrite( res3, "table_model_diffs_rounded_pretty.csv" )

# TeX for paper
print( xtable(res3), include.rownames = FALSE)


# ~~ Sanity checks on conditional Phats ------------------------------------------------------------------

# don't move this section
# depends on having object res2 above
# spot-check the conditional Phats

# check Perc0.2M in naive model
temp = robu( yi ~ 1, 
             data = dma, 
             studynum = as.factor(study_id),
             var.eff.size = vi,
             modelweights = "HIER",
             small = TRUE)
expect_equal( round( 100 * mean(calibMA > 0.2) ),
              res2$est[ res2$model == "naive" & res2$stat == "Perc0.2M"] )


# check Perc0.2R in naive model
temp = robu( yi ~ 1, 
             data = dr, 
             studynum = as.factor(study_id),
             var.eff.size = vi,
             modelweights = "HIER",
             small = TRUE)
expect_equal( round( 100 * mean(calibR > 0.2) ),
              res2$est[ res2$model == "naive" & res2$stat == "Perc0.2R"] )


# check Perc0.2M in moderated model
temp = robu( yi ~ mean_agec_mos + test_lang + method, 
             data = dma, 
             studynum = as.factor(study_id),
             var.eff.size = vi,
             modelweights = "HIER",
             small = TRUE)

# conditional_calib_ests fn already has unit tests in analyze_helper.R :)
calib = conditional_calib_ests(.model = temp)$calib.shift
expect_equal( round( 100 * mean(calib > 0.2) ),
              res2$est[ res2$model == "mod" & res2$stat == "Perc0.2M"] )


# check Perc0.2R in moderated model
temp = robu( yi ~ mean_agec_mos + test_lang + method, 
             data = dr, 
             studynum = as.factor(study_id),
             var.eff.size = vi,
             modelweights = "HIER",
             small = TRUE)

# conditional_calib_ests fn already has unit tests in analyze_helper.R :)
calib = conditional_calib_ests(.model = temp)$calib.shift
expect_equal( round( 100 * mean(calib > 0.2) ),
              res2$est[ res2$model == "mod" & res2$stat == "Perc0.2R"] )




# 4. PUBLICATION BIAS (SOME PREREG'D AND SOME POST HOC) ------------------------------------------------------------------

section = 4


# ~ Affirmative and Nonaffirmative Counts ------------------------------------------------------------------

t = d %>% group_by(isMeta) %>%
  summarise( k.affirm = sum(affirm),
             k.nonaffirm = sum(affirm == FALSE),
             Paffirm = mean(affirm) )

update_result_csv( name = "MA k nonaffirmative",
                   value = t$k.nonaffirm[ t$isMeta == TRUE ] )

update_result_csv( name = "MA k affirmative",
                   value = t$k.affirm[ t$isMeta == TRUE ] )

# ~ Hedges Selection Model ------------------------------------------------------------------
# be careful about inference due to correlated point estimates
# can't fit model with 3 cutoffs because there are no significant negative studies
( m1 = weightfunct( effect = dma$yi,
                    v = dma$vi,
                    steps = c(0.025, 1),
                    table = TRUE ) )
# actually makes the estimate larger
# get SEs via the Hessian
H = m1[[2]]$hessian
ses = sqrt( diag( solve(H) ) )


statCI_result_csv( "weightr mu",
                   c(m1[[2]]$par[2],
                     m1[[2]]$par[2] - qnorm(.975) * ses[2],
                     m1[[2]]$par[2] + qnorm(.975) * ses[2]) )


update_result_csv( name = "weightr mu pval",
                   value = format_stat( 2 * ( 1 - pnorm( abs(m1[[2]]$par[2]) / ses[2] ) ), cutoffs = c(0.10, pval.cutoff) ) )


# ~ Worst-Case Meta-Analysis------------------------------------------------------------------

# meta-analyze only the nonaffirmatives
# 2-sided pval
( meta.worst = robu( yi ~ 1, 
                     data = dma[ dma$affirm == FALSE, ], 
                     studynum = study_id,
                     var.eff.size = vi,
                     modelweights = "HIER",
                     small = TRUE) )
mu.worst = meta.worst$b.r
t2.worst = meta.worst$mod_info$tau.sq
mu.lo.worst = meta.worst$reg_table$CI.L
mu.hi.worst = meta.worst$reg_table$CI.U
mu.se.worst = meta.worst$reg_table$SE
pval.worst = meta.worst$reg_table$prob


# sanity check
corrected_meta(yi = dma$yi,
               vi = dma$vi,
               eta = 1000,
               clustervar = dma$study_id,
               favor.positive = TRUE,
               model = "robust")


statCI_result_csv( "Worst mu",
                   c(meta.worst$b.r,
                     meta.worst$reg_table$CI.L,
                     meta.worst$reg_table$CI.U) )


update_result_csv( name = "Worst mu pval",
                   value = round(pval.worst, 3) )


# ~ Meta-Analysis with Eta = 4.70 (Post Hoc) ------------------------------------------------------------------
# code modified from PublicationBias::significance_funnel innards
eta = 4.70
pvals = dma$pval
A = dma$affirm
yi = dma$yi
vi = dma$vi
clustervar = dma$study_id
dat = dma
small = TRUE

metaCorr = corrected_meta( yi = yi,
                           vi = vi, 
                           eta = eta,
                           clustervar = clustervar,
                           model = "robust",
                           favor.positive = TRUE )

# **wow! quite close to replication mean


update_result_csv( name = "Corr MA est",
                   value = round( metaCorr$est, 2) )

update_result_csv( name = "Corr MA lo",
                   value = round( metaCorr$lo, 2) )

update_result_csv( name = "Corr MA hi",
                   value = round( metaCorr$hi, 2) )

update_result_csv( name = "Corr MA pval",
                   value = format.pval( metaCorr$pval, eps = pval.cutoff) )


# # N.P. for both point estimate and CI
# update_result_csv( name = "sval est to 1",
#                    section = 2,
#                    value = res$sval.est,
#                    print = FALSE )
# update_result_csv( name = "sval CI to 1",
#                    section = 2,
#                    value = res$sval.ci,
#                    print = FALSE )

#@not possible even though worst-case estimate has CI crossing null
# must be the weight thing again



# ~ Diagnostics for assumption that publication bias is one-tailed  ------------------------------------------------------------------

# p-value plot
if ( redo.plots == TRUE ) {
  
  pval_plot( yi = dma$yi,
             vi = dma$vi )
  
  my_ggsave("meta_pval_plot.png",
            width = 5,
            height = 4)
}

# percent of one-tailed p-values below 0.025 and above 0.975

# one-tailed p-vals
pvalOne = 1 - pnorm( dma$yi / sqrt(dma$vi) )

update_result_csv( name = "Perc pvals <0.025",
                   value = round( 100*mean(pvalOne < 0.025) ) )

update_result_csv( name = "Perc pvals >0.975",
                   value = round( 100*mean(pvalOne > 0.975) ) )



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                       3. SUBSET MODELS AND PLOTS (POST HOC)            
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

section = 3

# look at direction of moderator associations in MB vs. MA
# *a key difference: HPP (vs. CF) is associated with LARGER ES in replications, but SMALLER ES in meta-analysis

# fit pruned model within MA and MB separately
#  to look heuristically at interactions
# these xtables are in paper
fit_mr( .dat = dma,
        .label = "mod1_MA_only",
        .mods = c("mean_agec_mos", "test_lang", "method"),
        #.write.to.csv = TRUE,
        .write.table = TRUE )

fit_mr( .dat = dr,
        .label = "mod1_MLR_only",
        .mods = c("mean_agec_mos", "test_lang", "method"),
        #.write.to.csv = TRUE,
        .write.table = TRUE )


# sanity check: refit models manually
robu( yi ~ isMeta + mean_agec_mos + test_lang + (method=="b.hpp"), 
      data = d, 
      studynum = as.factor(study_id),
      var.eff.size = vi,
      modelweights = "HIER",
      small = TRUE)

# adding interactions doesn't improve fit
robu( yi ~ isMeta*mean_agec_mos + test_lang + isMeta*(method=="b.hpp"), 
      data = d, 
      studynum = as.factor(study_id),
      var.eff.size = vi,
      modelweights = "HIER",
      small = TRUE)


# 5. TABLE AND FOREST PLOT OF SUBSET ESTIMATES ------------------------------------------------------------------

# for each categorical moderator, fit subset meta-analysis


# ~ Fit subset model to each level of each categorical model ------------------------------------------------------------------
modsCat = mods2[ !mods2 %in% c("isMeta", "mean_agec_mos" ) ]

for ( i in 1:length(modsCat) ) {
  
  .mod = modsCat[[i]] 
  .levels = unique( d[[.mod]] )
  
  for ( l in .levels ) {
    
    # subset to chosen level of moderator and fit marginal meta-analysis
    .metaMA = fit_subset_meta( .dat = dma[ dma[,.mod] == l, ],
                               .mods = "1",
                               .label = NA,
                               .simple.return = TRUE )
    
    .metaR = fit_subset_meta( .dat = dr[ dr[,.mod] == l, ],
                              .mods = "1",
                              .label = NA,
                              .simple.return = TRUE ) 
    
    new.rows = data.frame( mod = .mod, 
                           level = l,
                           source = c( "Meta-analysis", "Replications"),
                           k = c( .metaMA$k, .metaR$k ),
                           est = c( .metaMA$intercept, .metaR$intercept ),
                           lo = c( .metaMA$intercept.lo, .metaR$intercept.lo ),
                           hi = c( .metaMA$intercept.hi, .metaR$intercept.hi ) )
    
    
    if ( i == 1 & l == .levels[1] ) res = new.rows else res = rbind( res, new.rows )
    
  }
  
}

setwd( results.dir )
setwd( "tables_to_prettify" )
write.csv( res, "subsets_by_moderator.csv")


### Sanity check: check one subset
temp = robu( yi ~ 1, 
             data = dr[ dr$method == "b.hpp", ], 
             studynum = as.factor(study_id),
             var.eff.size = vi,
             modelweights = "HIER",
             small = TRUE)

# results from calling fit_subset_meta
resRow = res %>% filter( mod == "method" & level == "b.hpp" & source == "Replications" )

expect_equal( resRow$est,
              temp$b.r[1] )

expect_equal( resRow$lo,
              temp$reg_table$CI.L )

expect_equal( resRow$hi,
              temp$reg_table$CI.U )

expect_equal( resRow$k,
              sum( dr$method == "b.hpp" ) )


# ~ Plot it like it's hot ------------------------------------------------------------------

if ( redo.plots == TRUE ) {
  # choose reasonable axis limits
  summary(res$est, na.rm = TRUE)
  min(res$lo, na.rm = TRUE)
  max(res$hi, na.rm = TRUE)
  breaks = seq(0, 1.6, 0.2)
  
  subsetplot <- res
  
  
  #prettyLabels = c( main_question_ids_preference = "Primary analysis")
  
  subsetplot$mod <- factor(subsetplot$mod, levels = 
                             c("test_lang", "method", "speech_type", "own_mother", "presentation", 
                               "dependent_measure", "main_question_ids_preference"), 
                           labels = 
                             c("TestLang", "Method", "SpeechType", "Mother", "Presentation", 
                               "DV", "MainQ:IDS"))
  
  subsetplot$level <- factor(subsetplot$level, levels = 
                               c("b.nonnative", "a.native", "c.artificial", "a.cf", "b.hpp", 
                                 "c.other", "b.naturalistic", "a.simulated", "c.filtered/synthesized", 
                                 "a.no", "b.yes", "a.tape recording", "b.video recording", "a.preference", 
                                 "b.affect", "a.yes", "b.no"), 
                             labels = 
                               c("non-native", "native", "artificial", "Central Fixation", "HPP", 
                                 "Other", "naturalistic", "simulated", "filtered/synthesized", 
                                 "no", "yes", "Audio", "Video", "Preference", "Affect", "yes", "no"))
  
  
  
  
  #@to do in prep code: should make pretty versions of moderator variables (var names and levels) in prep code
  
  ggplot( data = subsetplot, 
          aes( x = est,
               y = level,
               color = source ) ) + 
    geom_vline( xintercept = naive.MA.only$b.r,
                lty = 2,
                color = colors[2]) + 
    geom_vline( xintercept = naive.reps.only$b.r,
                lty = 2,
                color = colors[1]) + 
    
    geom_point( size = 3,
                position=ggstance::position_dodgev(height=-0.4) ) + 
    geom_errorbar( aes(xmin = lo,
                       xmax = hi),
                   width = 0,
                   position=ggstance::position_dodgev(height=-0.4) ) +
    xlab( "Subset estimate (SMD)") +
    ylab( "Moderator subset") + 
    scale_x_continuous( breaks = breaks ) +
    scale_color_manual( name = "Source", values = rev(colors) ) + 
    theme_classic() +
    facet_grid( rows = vars(mod),
                scales = "free_y",
                switch = "both") +
    # https://michaelbach.de/2012/07/22/R-ggplot2-and-axis-limits-unexpected-behaviour-solution.html
    coord_cartesian( xlim = c(breaks[1], breaks[length(breaks)] ) ) 
  
  my_ggsave( name = "subset_forest.pdf",
             width = 10,
             height = 7 )
  
}

# 6. MATCHING AND IPW (POST HOC) ------------------------------------------------------------------

section = 6

# ~ Matching  ------------------------------------------------------------------

# ~~ Look at age distribution in each source ------------------------------------------------------------------
# this is the only continuous covariate to be matched
# and it's one with little overlap
age_densities(d)


# ~~ Make matches - coarsened exact matching (CEM) ------------------------------------------------------------------

# start from all candidate moderators, not just the survivors from multivariable meta-regression
mods3 = mods2[ !mods2 %in% c("isMeta") ]

# # quintiles
# ageCuts = quantile( d$mean_agec_mos,
#                     probs = seq(0, 1, 1/4) )
# # how big are the differences between bins (months)?
# diff(ageCuts)

# 6-month age bins
( ageCuts = seq( min(d$mean_agec_mos), max(d$mean_agec_mos), 6 ) )

# exact matching on all categorical variables and coarsened exact matching, based on ageCuts above,
# for age
# MatchIt::matchit docs: "if a categorical variable does not appear in grouping, it will not be coarsened, so exact matching will take place on it"
# regarding subclasses: "setting method = "cem" performs coarsened exact matching. With coarsened exact matching, covariates are coarsened into bins, and a complete cross of the coarsened covariates is used to form subclasses defined by each combination of the coarsened covariate levels. Any subclass that doesn't contain both treated and control units is discarded, leaving only subclasses containing treatment and control units that are exactly equal on the coarsened covariates."
string = paste( "isMeta ~ ",
                paste( mods3, collapse=" + "),
                collapse = "")

x = matchit( formula = eval( parse( text = string ) ),
             data = d,
             method = "cem",
             cutpoints = list( mean_agec_mos =  ageCuts) )
x

# 11 if using 6-month bins (8 MLR and 3 MA)
# 4 if using up to 3-month bins
summary(x)



# # nearest neighbor based on PS score: 102 matches
# x = matchit( formula = eval( parse( text = string ) ),
#              data = d,
#              method = "nearest" )
# x
# 
# # optimal pair matching based on PS score: also 102 matches
# library(optmatch)
# x = matchit( formula = eval( parse( text = string ) ),
#              data = d,
#              method = "optimal" )
# x

### How close are the matches on age?
# obviously the other variables are exactly matched

# look at match diagnostics (only interesting for age)
#  and number of matches in MA and MLR
summary(x)


# matched pairs only
dmt = match.data(x)
age_densities(dmt)  # definitely looks better


# how many subclasses?
# NAs are for unmatched observations
summary(x$subclass)

# look at characteristics of each subclass
CreateTableOne(vars = mods2, 
               strata = "study_type",
               data = dmt[ dmt$subclass == 1, ] ) 

CreateTableOne(vars = mods2, 
               strata = "study_type",
               data = dmt[ dmt$subclass == 2, ] ) 
# subclass 1 is younger and method = a.cf
# subclass 2 is older and method = b.hpp
# otherwise the other covariates are the same between strata




# ~ Analyze the CEM-matched dataset ------------------------------------------------------------------

# simple comparison of sources' average point estimates after matching
( t1 = dmt %>% group_by(isMeta) %>%
    summarise( k = n(),
               mean(yi),
               mAgec = mean(mean_agec_mos) ) )
# wow - matching seems to exacerbate the between-source discrepancy

# look within subclasses
( t2 = dmt %>% group_by(subclass, isMeta) %>%
    summarise( k = n(),
               mean(yi),
               mAgec = mean(mean_agec_mos) ) )


# the real analysis
matchRes = fit_mr( .dat = dmt,
                   .label = "matched",
                   .mods = "isMeta",
                   .write.to.csv = TRUE,
                   .write.table = FALSE,
                   .simple.return = FALSE )


update_result_csv( name = "matched k total",
                   value = nrow(dmt) )

update_result_csv( name = "matched k from MA",
                   value = sum( t2$k[ t2$isMeta == TRUE ] ) )

update_result_csv( name = "matched k from R",
                   value = sum( t2$k[ t2$isMeta == FALSE ] ) )

update_result_csv( name = "matched age diff mos",
                   value = diff( t1$mAgec ) )

update_result_csv( name = "pre-matching age diff mos",
                   value = mean( dr$mean_agec_mos ) - mean( dma$mean_agec_mos ) )


update_result_csv( name = paste( "matched subclass", t2$subclass, "isMeta", t2$isMeta, "k", sep = " " ),
                   value = t2$k )

update_result_csv( name = "pre-matching mean_age mos",
                   value = mean( d$mean_age / 30.44 ) )

update_result_csv( name = "matched mean_age mos subclass 1",
                   value = mean( dmt$mean_age[ dmt$subclass == 1 ] / 30.44 ) )

update_result_csv( name = "matched mean_age mos subclass 2",
                   value = mean( dmt$mean_age[ dmt$subclass == 2 ] / 30.44 ) )



# ~ IPW ------------------------------------------------------------------


# fit PS model: logit{ P(isMeta = 1) } = XB

string = paste( "isMeta ~ ",
                paste( mods3, collapse=" + "),
                collapse = "")

PSmod = glm( eval( parse( text = string ) ),
             family = binomial(link = "logit"),
             data = d )

# as expected, younger mean age is strongest predictor of being in MA vs. MLR
# can heuristically sanity-check these against Table 1 breakdown of moderators by source
summary(PSmod)

# add propensity scores to dataset
d$propScore = predict(PSmod, type = "response")
d$PSweight = 1/d$propScore
summary(d$propScore)


# ~~ Look at covariate overlap  ------------------------------------------------------------------
# not surprisingly, there is very poor overlap, which is why we get so few matches
ggplot( data = d[ d$isMeta == TRUE, ],
        aes(x = propScore,
            fill = isMeta ) ) +
  geom_histogram(aes(y = - ..density..)) + # negative sign for mirroring
  
  geom_histogram(data = d[ d$isMeta == FALSE, ],
                 aes(x = propScore,
                     y = ..density..,
                     fill = factor(isMeta)))


# ~~ Look at weights  ------------------------------------------------------------------

#bm
summary(d$PSweight)
# very extreme weights, as expected given lack of overlap

# # trimming weights doesn't really work because so many weights are equal to the max:
# table(d$PSweight == max(d$PSweight))
# q5 = quantile(d$PSweight, 0.1)
# q95 = quantile(d$PSweight, 0.9)
# d$PSweightTrim = pmin( d$PSweight, q95 )
# d$PSweightTrim = pmax( d$PSweightTrim, q5 )
# summary(d$PSweightTrim)

# fit weighted robust model
# just like SAPB, theoretically
# use the naive t2 from MA as the naive estimate
IPW.robu = robu( yi ~ isMeta,
                 studynum = as.factor(study_id),
                 data = d,
                 userweights = PSweight / (vi + c(naive.MA.only$mod_info$tau.sq) ),
                 var.eff.size = vi,
                 modelweights = "HIER",
                 small = TRUE )

est = as.numeric(IPW.robu$b.r)
lo = IPW.robu$reg_table$CI.L
hi = IPW.robu$reg_table$CI.U
pval.est = IPW.robu$reg_table$prob


update_result_csv( name = paste( "IPW robu est", IPW.robu$labels ),
                   value = est )

update_result_csv( name = paste( "IPW robu lo", IPW.robu$labels ),
                   value = round(lo, digits) )

update_result_csv( name = paste( "IPW robu hi", IPW.robu$labels ),
                   value = round(hi, digits) )

update_result_csv( name = paste( "IPW robu pval", IPW.robu$labels ),
                   value = pval.est )



# 7. FOREST PLOT ------------------------------------------------------------------


# don't move this section to be earlier! 
# must be after matching 
# relies on having the matched dataset for coloring those points

if ( redo.plots == TRUE ) {
  
  # should we color-code studies in the matched set?
  color.subclasses = TRUE
  
  # ~ Make plotting dataframe ------------------------------------------------------------------
  
  # relative weight of each study in meta-analysis (within group)
  d = d %>% group_by(isMeta) %>%
    mutate( rel.wt = 100 * (1/vi) / sum(1/vi) )
  
  # plotting df
  dp = d
  
  
  # strings with information about each meta-analysis
  dp$info.strings = paste( "n = ", d$n,
                           ", ", round( d$yi, 2 ), " ",
                           format_CI( d$lo, d$hi, 2 ),
                           sep = "" )
  
  
  # labels for each meta-analysis, including stats
  dp$label = paste(dp$unique, " (", dp$info.strings, ")", sep = "")
  
  
  # sort by point estimate
  dp = dp[ order( dp$sourcePretty,
                  dp$yi,
                  decreasing = FALSE), ]
  
  
  # add pooled point estimates after each group
  # these have arbitrary relative weights
  
  for ( l in unique(dp$sourcePretty) ) {
    # row ID after which to insert the pooled row
    this.group.rows = which(dp$sourcePretty == l)
    
    pooled.label = paste( "POOLED - ", toupper(l), sep = "" )
    
    dp = add_row( as.data.frame(dp),
                  .after = this.group.rows[ length(this.group.rows) ],
                  yi = ifelse( l == "Meta-analysis", naive.MA.only$b.r, naive.reps.only$b.r ),
                  lo = ifelse( l == "Meta-analysis", naive.MA.only$reg_table$CI.L, naive.reps.only$reg_table$CI.L ),
                  hi = ifelse( l == "Meta-analysis", naive.MA.only$reg_table$CI.U, naive.reps.only$reg_table$CI.U ),
                  label = pooled.label,
                  sourcePretty = l,
                  rel.wt = 5 )
  }
  
  
  # for aes on plot
  dp$is.pooled = grepl("POOL", dp$label)
  if ( color.subclasses == TRUE ) dp$is.pooled[ dp$unique %in% dmt$unique ] = 2
  
  # for pooled estimates, don't include number of studies in string
  ind = grepl("POOL", dp$label)
  dp$info.strings[ind] = paste( round( dp$yi[ind], 2 ), " ",
                                format_CI( dp$lo[ind], dp$hi[ind], 2 ),
                                sep = "" )
  
  # to force y-axis ordering, turn into a factor whose levels match the 
  #  actual order in the dataset
  dp$label = factor( dp$label, levels = rev(dp$label) )
  levels(dp$label)
  
  
  # regular circle: 19
  # 2 is open triangle
  shapes = c(19,17)
  if ( color.subclasses == TRUE ) shapes = c(19,17,19)
  
  # ~ Make the Plot ------------------------------------------------------------------
  
  # now color-coding by whether it's the pooled estimate or not
  colors2 = c("black", "red2")
  if( color.subclasses == TRUE ) colors2 = c("black", "red2", "green3" )
  
  # choose good axis breaks
  min(d$lo)
  max(d$hi)
  xBreaks = seq( -1.5, 4.5, 0.5)
  
  p = ggplot( data = dp, aes( x = yi, 
                              y = label, 
                              size = rel.wt,
                              shape = as.factor(is.pooled),
                              color = as.factor(is.pooled) ) ) +
    
    geom_errorbarh( aes(xmin = lo,
                        xmax = hi),
                    lwd = .5,
                    height = .001 ) +
    
    geom_point() +
    
    # geom_point( data = dp, aes( x = yi,
    #                                y = label ),
    #             size = 3,
    #             shape = 124,
    #             color = "black") +  # "4" for X
    
    xlab( "Point estimate (SMD)" ) +
    ylab("") +
    
    geom_vline(xintercept = 0, lty = 2) +
    
    guides(size = "none") +
    
    scale_shape_manual(values = shapes,
                       name = "",
                       guide="none") +
    
    scale_color_manual(values = colors2,
                       name = "",
                       guide="none") +
    
    scale_x_continuous(limits = c(min(xBreaks), max(xBreaks)),
                       breaks = xBreaks) +
    
    
    facet_grid(sourcePretty ~ .,
               scales = "free",
               space = "free_y") +  # allows y-axes to drop levels from other groups
    
    theme_classic() + 
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 20),
          strip.text.y = element_text(size = 20))
  
  p
  
  
  
  my_ggsave( name = "basic_forest.pdf",
             width = 14,
             height = 28 )
  
  
}



# 8. WITH MODIFIED SUBJECT INCLUSION CRITERIA ------------------------------------------------------------------


section = 8
setwd(data.dir)

quick_sens_analysis( .dat = dic,
                     .suffix = "IC" )


# sanity check: check just the dataset dimension counts
# other than that, quick_sens_analysis is just a straightforward set of calls to fit_mr
# resCSV is assigned as a global var by fns that call update_result_csv
if ( exists("resCSV") ) {
  
  expect_equal( resCSV$value[ resCSV$name == "Reps subset naiveIC k" ],
                as.character( nrow( dic[ dic$isMeta == FALSE, ] ) ) )
  
  expect_equal( resCSV$value[ resCSV$name == "Reps subset naiveIC n subj" ],
                as.character( round( sum( dic[ dic$isMeta == FALSE, "n" ], 0 ) ) ) )
}

# 9. WITH IPD AGE-MATCHING IN MLR ------------------------------------------------------------------

section = 9
setwd(data.dir)

# some basic stats about age-matched data
# this should be very close to 0 now
update_result_csv( name = "mean mean_agec_mos drage",
                   value = mean(drage$mean_agec_mos) )

# raw age (transformed to months)
update_result_csv( name = "mean mean_age mos drage",
                   value = round( mean(drage$mean_age / 30.44), 1 ) )

# raw age in MA for comparison (months)
update_result_csv( name = "mean mean_age mos dma",
                   value = round( mean(dma$mean_age / 30.44), 1 ) )


quick_sens_analysis( .dat = dage,
                     .suffix = "ageMatch" )


# 10. OTHER SENSITIVITY ANALYSES (SUPPLEMENT)  ------------------------------------------------------------------

# ~ Exclude between-subjects designs  ------------------------------------------------------------------

# look at design types
table(d$participant_design,
      d$sourcePretty,
      useNA = "ifany")

update_result_csv( name = "k between-S studies",
                   value = sum(d$participant_design == "between") )

#bm
quick_sens_analysis( .dat = d[ d$participant_design != "between", ],
                     .suffix = "WithinS" )



# ~ Look at age linearity assumption  ------------------------------------------------------------------

# not in manuscript
# omit SEs because they'll be wrong
# use calibrated estimates to somewhat account for studies' precisions
ggplot(d, aes(mean_age, calibNaive, color=sourcePretty)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  scale_color_manual(values=colors) +
  theme_classic()

# apparent nonlinearity in replications seems largely driven by the outlying replication
#  study with largest mean age
# try excluding that one
temp = d %>% filter( !(isMeta == TRUE & mean_age > 300 ) )
expect_equal(nrow(temp), nrow(d) - 1)

ggplot(temp, aes(mean_age, calibNaive, color=sourcePretty)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  scale_color_manual(values=colors) +
  theme_classic()




# ~ Publication bias with reported (not calculated) significance  ------------------------------------------------------------------

# keep only MA studies that reported whether results were significant
dma2 = dma %>% filter( !is.na(reportedSignif) )
nrow(dma2)

# proportions of estimates reported significant,
#  stratified by whether calculated p-value was < 0.05
( t = dma2 %>%
    group_by(pvalSignif) %>%
    summarise( n(),
               mean(reportedSignif) ) )

update_result_csv( name = paste( "Perc reportedSignif when pvalSignif was", t$pvalSignif ),
                   value = round( 100*t$`mean(reportedSignif)`, 0) )

# huh! not a very strong correlation between 
#  calculated and reported significance
cor.test(dma2$reportedSignif, dma2$pvalSignif)

# new affirmative indicator
dma2$affirm2 = (dma2$reportedSignif == 1) & (dma2$yi > 0)

# this worst-case estimate is a bit smaller than the one in main analysis
( meta.worst2 = robu( yi ~ 1, 
                      data = dma2[ dma2$affirm2 == FALSE, ], 
                      studynum = study_id,
                      var.eff.size = vi,
                      modelweights = "HIER",
                      small = TRUE) )

mu.worst = meta.worst2$b.r
t2.worst = meta.worst2$mod_info$tau.sq
mu.lo.worst = meta.worst2$reg_table$CI.L
mu.hi.worst = meta.worst2$reg_table$CI.U
mu.se.worst = meta.worst2$reg_table$SE
pval.worst = meta.worst2$reg_table$prob


# sanity check
corrected_meta_2(yi = dma2$yi,
                 vi = dma2$vi,
                 reportedSignif = dma2$reportedSignif,
                 eta = 1000,
                 clustervar = dma2$study_id,
                 favor.positive = TRUE,
                 model = "robust")


statCI_result_csv( "Worst reportedSignif mu",
                   c(meta.worst2$b.r,
                     meta.worst2$reg_table$CI.L,
                     meta.worst2$reg_table$CI.U) )


update_result_csv( name = "Worst reportedSignif mu pval",
                   value = round(pval.worst, 3) )



# ~ S-values -----
# turns out similarly to main results
# s-values to reduce to null
( Sval0 = svalue( yi = dma2$yi,
                  vi = dma2$vi,
                  q = 0,
                  clustervar = dma2$study_id,
                  favor.positive = TRUE,
                  model = "robust" ) )

# s-values to reduce effect size to match replications
( SvalR = svalue( yi = dma$yi,
                  vi = dma$vi,
                  q = naive.reps.only$b.r,
                  clustervar = dma$study_id,
                  favor.positive = TRUE,
                  model = "robust" ) )

# N.P. for estimate
update_result_csv( name = "sval est to reps",
                   value = round( SvalR$sval.est, 2 ),
                   print = FALSE )
update_result_csv( name = "sval CI to reps",
                   value = round( SvalR$sval.ci, 2 ),
                   print = FALSE )


# ~ Meta-Analysis with Eta = 4.70 (Post Hoc) ------------------------------------------------------------------
# turns out very similarly to main results
eta = 4.70


metaCorr2 = corrected_meta_2( yi = dma2$yi,
                              vi = dma2$vi,
                              reportedSignif = dma2$reportedSignif,
                              eta = eta,
                              clustervar = dma2$study_id,
                              model = "robust",
                              favor.positive = TRUE )

# again quite close to replication mean
# not reported in manuscript except to say that results were similar
update_result_csv( name = "Corr MA reportedSignif est",
                   value = round( metaCorr2$est, 2) )

update_result_csv( name = "Corr MA reportedSignif lo",
                   value = round( metaCorr2$lo, 2) )

update_result_csv( name = "Corr MA reportedSignif hi",
                   value = round( metaCorr2$hi, 2) )

update_result_csv( name = "Corr MA reportedSignif pval",
                   value = format.pval( metaCorr2$pval, eps = pval.cutoff) )

