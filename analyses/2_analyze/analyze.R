
# META-NOTES ------------------------------------------------------------------ 

# ~ To do  ------------------------------------------------------------------


# - Think about doing something IPD for age (look at raw data to see if this makes sense)

# - Put analysis that doesn't match on age but matches on the others into manuscript. 

# stop using asterisks in file names bc they cause syncing trouble (for now, I just manually removed them)

# section variable needs updating after code structure is done

# ~ Ask CB, et al ------------------------------------------------------------------

# in the 0.75 dataset with the more stringent inclusion criterion, am I right in thinking that there are fewer effect sizes because some age groups are dropped completely? (but the mean age in MB doesn't change much at all)

# the 0.125 dataset corresponds with main analysis, right?


# ~ Other notes  ------------------------------------------------------------------

# names of important model objects:
# - naive.MA.only and naive.reps.only: meta-analyses within subsets; no moderators


# PRELIMINARIES ------------------------------------------------------------------

library(tidyverse) 
library(knitr)
library(here)
library(tableone)
library(corrr)
library(robumeta)
library(fastDummies)
library(weightr)
library(PublicationBias)
library(xtable)
library(boot)
library(testthat)
library(ggplot2)
library(metafor)
library(MatchIt)

data.dir = here("data")
# where to save results
results.dir = here("results_from_R")
# results.dir = "~/Dropbox/Personal computer/Independent studies/2020/Christina's ManyBabiesMeta (MB-Meta)/IDSPreference_ManyBabiesMeta/results_from_R"
overleaf.dir = "~/Dropbox/Apps/Overleaf/MB-Meta/R_objects"
code.dir = here("analyses/2_analyze")

# helper fns
setwd(code.dir)
source("analyze_helper.R")

# package update not yet on CRAN, but we need the cluster-bootstrap functionality
source("MetaUtility development functions.R")

# ~ Code-Running Parameters ------------------------------------------------------------------
# should we remove existing results file instead of overwriting individual entries? 
start.res.from.scratch = FALSE
# should we use the grateful package to scan and cite packages?
cite.packages.anew = FALSE
# should we bootstrap from scratch or read in old resamples?
boot.from.scratch = FALSE
# make plots from scratch?
redo.plots = FALSE
# redo iterative selection of moderators to get model to converge?
redo.mod.selection = FALSE

if (redo.mod.selection == FALSE) {
  # read in the surviving moderators
  setwd(results.dir)
  modsS = read.csv("surviving_mods.csv")[,2]
}

# wipe results csvs if needed
if ( start.res.from.scratch == TRUE ) wr()

#~ Constants of Universe ------------------------------------------------------------------
digits = 2
pval.cutoff = 10^-4  # threshold for using "<"
boot.reps = 1000 # for all bootstrapped inference 
# plot colors
colors = c("darkgray", "red")  # replications, originals

# stop sci notation
options(scipen=999)

# moderator names
mods = c( "study_type",
          "mean_agec",
          "test_lang",  # whether stimuli were in native language; almost constant in meta
          "native_lang", #Won't be used in main analyses
          "method",
          
          # constant in RRR:
          "speech_type",
          "own_mother",
          "presentation",
          "dependent_measure",
          "main_question_ids_preference" )

# ~ Read Datasets ------------------------------------------------------------------
setwd(data.dir)
d = suppressMessages( suppressWarnings( read_csv("mb_ma_combined_prepped.csv") ) )

# dataset with just the meta-analysis
dma = d %>% filter(isMeta == TRUE)

# dataset with just the replications
dr = d %>% filter(isMeta == FALSE)

# dataset with more stringent inclusion criteria in replications
dic = suppressMessages( suppressWarnings( read_csv("mb_ma_combined_prepped_0.75.csv") ) )
dric = dic %>% filter(isMeta == FALSE)

# dataset with IPD age-matching in MLR
dage = suppressMessages( suppressWarnings( read_csv("mb_ma_combined_prepped_0.125_age_matched.csv") ) )
drage = dage %>% filter(isMeta == FALSE)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 0. CHARACTERISTICS OF INCLUDED STUDIES ------------------------------------------------------------------           


# for updating result csv
section = 0

# ~ Basics ------------------------------------------------------------------

t = d %>% group_by(study_type) %>%
  summarise( m = n(),
             k = length(unique(study_id) ),
             nSubjMed = median(n) )
t
# seems like way fewer than what we expected?


update_result_csv( name = paste( "m ests", t$study_type, sep = " "),
                   value = t$m )

update_result_csv( name = "m ests total",
                   value = sum(t$m) )

update_result_csv( name = paste( "k studies", t$study_type, sep = " "),
                   value = t$k )

update_result_csv( name = paste( "median n", t$study_type, sep = " "),
                   value = t$nSubjMed )


# ~ Moderators ------------------------------------------------------------------

# distribution of moderators in RRR and MA
t = CreateTableOne(vars = mods, 
                   strata = "study_type",
                   data = d)

xtable( print(t, noSpaces = TRUE, printToggle = FALSE) )


# ~~ Moderator Correlation Matrix ------------------------------------------------------------------

# temporarily recode all categorical moderators as dummies
cat.mods = mods[ !mods == "mean_agec" ]
temp = dummy_cols( d[,cat.mods],
                   select_columns = cat.mods,
                   remove_selected_columns = TRUE,
                   ignore_na =  TRUE )  
temp$mean_agec = d$mean_agec 
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




# 1. PRIMARY MODERATOR ANALYSES (PREREG'D) ------------------------------------------------------------------

section = 1


# create MA subset and MLR subset
# this fn also writes stats to results csv
( naive.MA.only = fit_subset_meta( .dat = dma,
                                   .mods = "1",
                                   .label = "MA subset naive" ) )


( naive.reps.only = fit_subset_meta( .dat = dr,
                                     .mods = "1",
                                     .label = "Reps subset naive" ) )


# ~ Fit naive and moderated regressions ------------------------------------------------------------------ 


if ( redo.mod.selection == TRUE ) {
  
  # order of importance given in prereg:
  # age, test_lang, method, speaker, speech_type, own_mother, presentation, DV, main question
  
  # per prereg, if model isn't estimable, the moderators are to be removed in the opposite order
  #  of this importance list
  
  # try to fit meta-regression
  # same list as mods except has isMeta instead of study_type
  #  and does not have native_lang
  mods2 = c( "isMeta",  # code this way since we expect meta to have larger effect sizes
             "mean_agec",
             "test_lang",  # whether stimuli were in native language; almost constant in meta
             "method",
             
             # constant in RRR:
             "speech_type",
             "own_mother",
             "presentation",
             "dependent_measure",  # causes singularity
             "main_question_ids_preference" )
  
  
  
  mod.sets = list( c("isMeta"),  # naive model
                   mods2 )
  
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
      
    }, error = function(err){
      gotError <<- TRUE
      
      # remove one moderator from the end of the list
      message( paste( "\n Removing ", 
                      mod.sets[[2]][ length(mod.sets[[2]]) ],
                      " from moderators and trying again" ) )
      mod.sets[[2]] <<- mod.sets[[2]][ 1 : ( length(mod.sets[[2]]) - 1 ) ]
      
      
    })
    
  }
  
  # look at the surviving moderators
  mod.sets[[2]]
  
  # write the list of moderators  so we don't have to do this again
  setwd(results.dir)
  write.csv(mod.sets[[2]], file = "surviving_mods.csv")
  
  # for later use
  modsS = mod.sets[[2]]
  
  ##### Sanity checks
  # sanity check: refit the naive and pruned models manually
  robu( yi ~ isMeta, 
        data = d, 
        studynum = as.factor(study_id),
        var.eff.size = vi,
        modelweights = "HIER",
        small = TRUE)
  
  robu( yi ~ isMeta + mean_agec + test_lang + method, 
        data = d, 
        studynum = as.factor(study_id),
        var.eff.size = vi,
        modelweights = "HIER",
        small = TRUE)
  
} # end if ( redo.mod.selection == TRUE )


# 2. DENSITY PLOT OF META-ANALYSIS VS. REPLICATION CALIBRATED ESTIMATES ------------------------------------------------------------------


# ~ Get Marginal Calibrated Estimates for Plot ------------------------------------------------------------------
# SAVE calibR, calibM because needed for plot
calibR = calib_ests(yi = dr$yi, sei = sqrt(dr$vi) )
mean(calibR>0.2)
mean(calibR)

calibMA = calib_ests(yi = dma$yi, sei = sqrt(dma$vi) )
mean(calibMA>0.2)
mean(calibMA)

d$calibNaive = NA
d$calibNaive[ d$isMeta == FALSE ] = calibR
d$calibNaive[ d$isMeta == TRUE ] = calibMA


# ~ Get Conditional Calibrated Estimates for Plot ------------------------------------------------------------------
# fit the final, moderated models to each subset
# to see what the heterogeneity looks like in each case
# @move this?
( cond.MA.only = fit_subset_meta( .dat = dma,
                                  .mods = modsS[ !modsS == "isMeta" ],
                                  .label = "MA subset mods" ) )

( cond.reps.only = fit_subset_meta( .dat = dr,
                                    .mods = modsS[ !modsS == "isMeta" ],
                                    .label = "Reps subset mods" ) )


# get conditional calibrated estimates
d$calibCond[ d$isMeta == FALSE ] = conditional_calib_ests(cond.reps.only)$calib.shift
d$calibCond[ d$isMeta == TRUE ] = conditional_calib_ests(cond.MA.only)$calib.shift

# ***important: MLR heterogeneity estimate is 0 after conditioning on mods
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

# get inference for conditional Phats and their difference FOR the moderated model

# for each resample, fits both naive model and moderated model
# via fit_mr, gets 7 stats of interest for each model
# returns results of naive model, moderated model, AND their difference
# in a 21-length vector

# order of fit_mr's 9 returned stats:
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
                     
                     # # this works fine
                     # b = original[indices,]
                     # mean(rnorm(20))
                     
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
                       # increase number of NA's if needed to match .simple.return = TRUE structure of fit_mr
                       # x 3 separate calls to fit_mr
                       return( rep(NA, 9 * 3) )
                     })
                     
                   } )
  
  # save bootstraps
  setwd(results.dir)
  save(boot.res, file = "saved_bootstraps.RData")
}

if ( boot.from.scratch == FALSE ){
  # read in existing bootstraps
  setwd(results.dir)
  load("saved_bootstraps.RData")
}



# ~ Process the Resamples ------------------------------------------------------------------

# number of non-failed bootstrap reps
( boot.reps.successful = sum( !is.na( boot.res$t[,1] ) ) )


# get stats on original data
# "2" suffix is because we already have naiveRes and modRes from the model selection
#   part, but those objects use the full return structure rather than the boot-friendly one
naiveRes2 = fit_mr( .dat = d, .mods = "isMeta" )
modRes2 = fit_mr( .dat = d, .mods = modsS )
diff = naiveRes2 - modRes2
( t0 = c(naiveRes2, modRes2, diff) )

# NB: because of the way we hacked the "statistic" argument of boot()
#  above st it actually creates a resample, the below will NOT match
#  the actual estimates from original dataset
# the bootstrap "bias" estimates will also be wrong for this reason
# boot.res$t0

# sanity check: agreement between fit_mr's two
#  return structures
#@make sure make notes about reason for needing both return structures
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

# return structure of fit_mr:
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
# this one makes sense b/c estimate was at ceiling

# set such CIs to NA
res$lo[ res$est < res$lo | res$est > res$hi ] = NA
res$hi[ res$est < res$lo | res$est > res$hi ] = NA

# if only one CI limit was estimable, set both to NA
res$lo[ is.na(res$hi) ] = NA
res$hi[ is.na(res$lo) ] = NA

# **all of the CIs in this df are bootstrapped

# ### COMPARE BOOT VS. MODEL-BASED CIs FOR COEFF ESTIMATES:
# # model-based ones are consistently a lot wider
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
# 
# # *for this reason, perhaps it actually makes more sense to use the boot CIs here
# # also for comparability within the table
# ### END CI COMPARISON



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
library(data.table)
fwrite( res, "table_model_diffs_unrounded.csv" )

# ~~~ Rounded numeric table (for piping numbers into manuscript) ------------------------------------------------------------------
res2 = res
numVars = c("est", "lo", "hi")
res2[ , numVars ] = round( res2[ , numVars ], 2 )
# also turn Phats into percentages
inds = grepl( pattern = "Phat", x = res2$stat )
res2[ inds, numVars ] = 100*res2[ inds, numVars ]
res2$stat[ inds ] = str_replace( string = res2$stat[ inds ],
                                 pattern = "Phat",
                                 replacement = "Perc" )

# make single key column for use with rlgetrum in TeX
res2 = res2 %>% add_column( unique = paste( res2$model, res2$stat, sep = "_") )

fwrite( res2, "table_model_diffs_rounded.csv" )

# also save to Overleaf
setwd(overleaf.dir)
fwrite( res2, "table_model_diffs_rounded.csv" )


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


# ~~ Sanity checks ------------------------------------------------------------------

#@also sanity-check the simple models using calls to quick_phat above

# get_and_write_phat( .dat = dma,
#                     .q = 0,
#                     label = "Phat0MA")
# 
# get_and_write_phat( .dat = dma,
#                     .q = 0.2,
#                     label = "Phat0.2MA")
# 
# get_and_write_phat( .dat = dr,
#                     .q = 0,
#                     label = "Phat0R")
# 
# get_and_write_phat( .dat = dr,
#                     .q = 0.2,
#                     label = "Phat0.2R")



# 4. PUBLICATION BIAS (SOME PREREG'D AND SOME POST HOC) ------------------------------------------------------------------

section = 4

# @MOVE THIS
# sanity check: should be similar to naive meta-regression model
#  but not necessarily equal
expect_equal( as.numeric(naive.MA.only$b.r), round(naiveRes$est.ma, digits), tol = 0.03 )
expect_equal( as.numeric(naive.MA.only$reg_table$CI.L), round(naiveRes$est.ma.lo, digits), tol = 0.03 )
expect_equal( as.numeric(naive.MA.only$reg_table$CI.U), round(naiveRes$est.ma.hi, digits), tol = 0.03 )

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


# ~ S-values ------------------------------------------------------------------
# s-values to reduce to null
( Sval0 = svalue( yi = dma$yi,
                  vi = dma$vi,
                  q = 0,
                  clustervar = dma$study_id,
                  favor.positive = TRUE,
                  model = "robust" ) )



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


# ~ Significance Funnel ------------------------------------------------------------------

if ( redo.plots == TRUE ) {
  # modified from package
  yi = d$yi
  vi = d$vi
  xmin = min(yi)
  xmax = max(yi)
  ymin = 0  # so that pooled points are shown
  ymax = max( sqrt(vi) )
  xlab = "Point estimate (SMD)"
  ylab = "Estimated standard error"
  favor.positive = TRUE
  alpha.select = 0.05
  plot.pooled = TRUE
  
  # pooled points to plot
  est.MA = naive.MA.only$b.r  # naive est in MA
  est.R = naive.reps.only$b.r  # naive est in reps
  est.SAPBE = metaCorr$est  # not actually worst-case, but rather SAPB-E estimate
  est.W = mu.worst
  
  
  d$sei = sqrt(vi)
  
  # calculate p-values
  d$pval = 2 * ( 1 - pnorm( abs(yi) / sqrt(vi) ) )
  
  # which direction of effects are favored?
  # if we have the pooled point estimate, but not the favored direction,
  #  assume favored direction matches sign of pooled estimate (but issue warning)
  if ( !is.na(est.all) & is.na(favor.positive) ) {
    favor.positive = (est.all > 0)
    warning("favor.positive not provided, so assuming publication bias favors estimates whose sign matches est.all")
  }
  if ( is.na(est.all) & is.na(favor.positive) ) {
    stop("Need to specify favor.positive")
  }
  
  # affirmative vs. nonaffirmative indicator
  d$affirm = rep(NA, nrow(d))
  
  if ( favor.positive == TRUE ) {
    d$affirm[ (d$yi > 0) & (d$pval < alpha.select) ] = "Affirmative"
    d$affirm[ (d$yi < 0) | (d$pval >= alpha.select) ] = "Non-affirmative"
  }
  if ( favor.positive == FALSE ) {
    d$affirm[ (d$yi < 0) & (d$pval < alpha.select) ] = "Affirmative"
    d$affirm[ (d$yi > 0) | (d$pval >= alpha.select) ] = "Non-affirmative"
  }
  
  # reorder levels for plotting joy
  d$affirm = factor( d$affirm, c("Non-affirmative", "Affirmative") )
  
  # stop if no studies in either group
  if ( sum( d$affirm == "Non-affirmative" ) == 0 ) {
    stop("There are no non-affirmative studies. The plot would look silly.")
  }
  
  if ( sum( d$affirm == "Affirmative" ) == 0 ) {
    stop("There are no affirmative studies. The plot would look silly.")
  }
  
  
  # set up pooled estimates for plotting
  pooled.pts = data.frame( yi = c(est.MA, est.R, est.SAPBE, est.W),
                           sei = c(0,0,0,0) )
  
  # for a given SE (y-value), return the "just significant" point estimate value (x-value)
  just_signif_est = function( .sei ) .sei * qnorm(1 - alpha.select/2)
  
  # calculate slope and intercept of the "just affirmative" line
  # i.e., 1.96 = (just affirmative estimate) / se
  if (favor.positive == TRUE) sl = 1/qnorm(1 - alpha.select/2)
  if (favor.positive == FALSE) sl = -1/qnorm(1 - alpha.select/2)
  int = 0
  # # sanity check: should be exactly alpha.select
  # 2 * ( 1 - pnorm( abs(1) / sl ) )
  
  
  # make the plot
  p.funnel = ggplot( data = d, aes( x = yi,
                                    y = sei,
                                    color = isMeta ) )
  
  if ( plot.pooled == TRUE ) {
    
    # plot the pooled points
    # outer part of diamonds
    p.funnel = p.funnel + geom_point(
      data = pooled.pts,
      aes( x = yi, y = sei ),
      size = 4,
      shape = 5,
      fill = NA,
      color = c("red", "darkgray", NA, NA)
    ) +
      
      # inner part of diamonds
      geom_point(
        data = pooled.pts,
        aes( x = yi, y = sei ),
        size = 4,
        shape = 18,
        color =  c("red", "darkgray", "red", "orange"),
        alpha = 1
      ) +
      
      # just for visual separation of pooled ests
      geom_hline( yintercept = 0 ) +
      
      # diagonal "just significant" line
      geom_abline(slope=sl,intercept = int, color = "gray")
  }
  
  p.funnel = p.funnel +
    
    # # semi-transparent points with solid circles around them
    # geom_point( size = 3, alpha=.4) +
    # geom_point( size = 3, shape = 1) +
    
    geom_point(size = 3, shape = 1) +
    
    scale_color_manual(values = colors) +
    #scale_shape_manual(values = c(1,2)) +
    
    xlab(xlab) +
    ylab(ylab) +
    
    scale_x_continuous( limits = c(xmin, xmax), breaks = seq(-0.5, 3, .5) ) +
    scale_y_continuous( limits = c(ymin, ymax), breaks = seq(0, .6, .1) ) +
    
    theme_classic() +
    #theme(legend.title=element_blank())
    theme(legend.position = "none")
  
  plot(p.funnel)
  
  
  my_ggsave("mb_signif_funnel.png",
            width = 5,
            height = 4)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                       3. SUBSET MODELS AND PLOTS (POST HOC)            
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

section = 3

robu( yi ~ mean_agec + test_lang + method,
      data = dma,
      studynum = as.factor(study_id),
      var.eff.size = vi,
      modelweights = "HIER",
      small = TRUE)

robu( yi ~ mean_agec + test_lang + method,
      data = dr,
      studynum = as.factor(study_id),
      var.eff.size = vi,
      modelweights = "HIER",
      small = TRUE)



# *a key difference: HPP (vs. CF) is associated with LARGER ES in replications, but SMALLER ES in meta-analysis

# sanity check: fit pruned model within MA and MLR separately
#  to look heuristically at interactions
fit_mr( .dat = dma,
        .label = "mod1_MA_only",
        .mods = c("mean_agec", "test_lang", "method"),
        #.write.to.csv = TRUE,
        .write.table = TRUE )

fit_mr( .dat = dr,
        .label = "mod1_MLR_only",
        .mods = c("mean_agec", "test_lang", "method"),
        #.write.to.csv = TRUE,
        .write.table = TRUE )




robu( yi ~ isMeta + mean_agec + test_lang + (method=="b.hpp"), 
      data = d, 
      studynum = as.factor(study_id),
      var.eff.size = vi,
      modelweights = "HIER",
      small = TRUE)

# adding interactions doesn't improve fit
robu( yi ~ isMeta*mean_agec + test_lang + isMeta*(method=="b.hpp"), 
      data = d, 
      studynum = as.factor(study_id),
      var.eff.size = vi,
      modelweights = "HIER",
      small = TRUE)

# # save:
# 5. TABLE AND FOREST PLOT OF SUBSET ESTIMATES ------------------------------------------------------------------

# for each categorical moderator, fit subset meta-analysis


# ~ Fit subset model to each level of each categorical model ------------------------------------------------------------------
modsCat = mods2[ !mods2 %in% c("isMeta", "mean_agec" ) ]

for ( i in 1:length(modsCat) ) {
  
  .mod = modsCat[[i]] 
  .levels = unique( d[[.mod]] )
  
  for ( l in .levels ) {
    
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



# ~ Plot it like it's hot ------------------------------------------------------------------

if ( redo.plots == TRUE ) {
  # choose reasonable axis limits
  summary(res$est, na.rm = TRUE)
  min(res$lo, na.rm = TRUE)
  max(res$hi, na.rm = TRUE)
  breaks = seq(0, 1.6, 0.2)
  
  
  prettyLabels = c( main_question_ids_preference = "Primary analysis")
  
  
  #@to do in prep code: should make pretty versions of moderator variables (var names and levels) in prep code
  
  ggplot( data = res, 
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



# ~~ Make matches - CEM ------------------------------------------------------------------
mods3 = mods2[ !mods2 %in% c("isMeta") ]

string = paste( "isMeta ~ ",
                paste( mods3, collapse=" + "),
                collapse = "")


# # quintiles
# ageCuts = quantile( d$mean_agec,
#                     probs = seq(0, 1, 1/4) )
# # how big are the differences between bins (months)?
# diff(ageCuts)

# 12-month bins
( ageCuts = seq( min(d$mean_agec), max(d$mean_agec), 12 ) )

# exact matching on all categorical variables and coarsened exact matching, based on ageCuts above,
# for age
# MatchIt::matchit docs: "if a categorical variable does not appear in grouping, it will not be coarsened, so exact matching will take place on it"
# regarding subclasses: "setting method = "cem" performs coarsened exact matching. With coarsened exact matching, covariates are coarsened into bins, and a complete cross of the coarsened covariates is used to form subclasses defined by each combination of the coarsened covariate levels. Any subclass that doesn't contain both treated and control units is discarded, leaving only subclasses containing treatment and control units that are exactly equal on the coarsened covariates."
x = matchit( formula = eval( parse( text = string ) ),
             data = d,
             method = "cem",
             cutpoints = list( mean_agec =  ageCuts) )
x

# 49 matches if not including mean_agec (but still only 3 from MA)
# 14 if using mean_agec with exact quintile matching
# 19 if using age tertiles
# 4 if using 6-8 month bins
# 17 if using 12-month bins
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

# how close are the matches on age?
# obviously the other variables are exactly matched


# look at match diagnostics (only interesting for age)
#  and number of matches in MA and MLR
summary(x)


# matched pairs only
dmt = match.data(x)
age_densities(dmt)  # definitely looks better


# how many subclasses?
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

# simple
t1 = dmt %>% group_by(isMeta) %>%
  summarise( k = n(),
             mean(yi),
             mAgec = mean(mean_agec) )
# wow - really exacerbates the difference

# look within subclasses
t2 = dmt %>% group_by(subclass, isMeta) %>%
  summarise( k = n(),
             mean(yi),
             mAgec = mean(mean_agec) )



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

update_result_csv( name = "matched age diff",
                   value = diff( t1$mAgec ) )

update_result_csv( name = "pre-matching age diff",
                   value = mean( dr$mean_agec ) - mean( dma$mean_agec ) )


update_result_csv( name = paste( "matched subclass", t2$subclass, "isMeta", t2$isMeta, "k", sep = " " ),
                   value = t2$k )

update_result_csv( name = "pre-matching mean_age",
                   value = mean( d$mean_age ) )

update_result_csv( name = "matched mean_age subclass 1",
                   value = mean( dmt$mean_age[ dmt$subclass == 1 ] ) )

update_result_csv( name = "matched mean_age subclass 2",
                   value = mean( dmt$mean_age[ dmt$subclass == 2 ] ) )



# ~ IPW ------------------------------------------------------------------


# fit PS model
# logit{ P(isMeta = 1) } = Xb + eps

string = paste( "isMeta ~ ",
                paste( mods3, collapse=" + "),
                collapse = "")

PSmod = glm( eval( parse( text = string ) ),
             family = binomial(link = "logit"),
             data = d )

# add propensity scores to dataset
d$propScore = predict(PSmod, type = "response")
d$PSweight = 1/d$propScore



# ~~ Look at covariate overlap  ------------------------------------------------------------------
# very poor overlap, which is why we get so few matches
ggplot( data = d[ d$isMeta == TRUE, ],
        aes(x = propScore,
            fill = isMeta ) ) +
  geom_histogram(aes(y = - ..density..)) + # negative sign for mirroring
  
  geom_histogram(data = d[ d$isMeta == FALSE, ],
                 aes(x = propScore,
                     y = ..density..,
                     fill = factor(isMeta)))


# ~~ Look at weights  ------------------------------------------------------------------
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
                 small = small )

est = as.numeric(IPW.robu$b.r)[2]  # isMeta coeff only
lo = IPW.robu$reg_table$CI.L[2]
hi = IPW.robu$reg_table$CI.U[2]
pval.est = IPW.robu$reg_table$prob[2]


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
# relies on having the matched dataset for coloring those points

if ( redo.plots == TRUE ) {
  
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
  # have breaks approximately evenly spaced on log scale
  #  but round for prettiness
  #breaks = unique( round( 10^seq(-3, 3, .2), 1 ) )
  
  
  # 
  # # re-level group.pretty so legend order matches display order
  # #  (i.e., from high to low pooled estimate)
  # agg2p$group.pretty = factor( agg2p$group.pretty,
  #                              levels =  rev(c( "Metalab",
  #                                               "Top psychology",
  #                                               "Top medical",
  #                                               "PLOS One" )  ))
  
  
  # ~ Make the Plot ------------------------------------------------------------------
  
  # now color-coding by whether it's the pooled estimate or not
  colors2 = c("black", "red")
  if( color.subclasses == TRUE ) colors2 = c("black", "red", "green" )
  
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
    
    geom_vline(xintercept = 1, lty = 2) +
    
    guides(size = FALSE ) +
    
    scale_shape_manual(values = shapes,
                       name = "",
                       guide=FALSE) +
    
    scale_color_manual(values = colors2,
                       name = "",
                       guide=FALSE) +
    
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

# replication dataset that excludes 
# see n_trial_pairs_criterion variable in "2_get_tidy_MB_data.R" for how this is created

# note that there are now fewer effect sizes, presumably because some age groups get lost completely upon excluding badly behaved subjects
dim(dric)

# ~ Subset model in MLR ------------------------------------------------------------------

# with more stringent inclusion criteria, estimate in MLR increases from 0.35 to 0.43
# still less than in meta-analysis
# this fn also writes stats to results csv
( naiveIC.reps.only = fit_subset_meta( .dat = dric,
                                     .mods = "1",
                                     .label = "Reps subset naiveIC" ) )

# ~ Naive model with both sources ------------------------------------------------------------------
naiveResIC = fit_mr( .dat = dic,
                   .label = "naiveIC",
                   .mods = "isMeta",
                   .write.to.csv = TRUE,
                   .write.table = TRUE,
                   .simple.return = FALSE )

# ~ Moderated model with both sources ------------------------------------------------------------------
# again does not improve, and in fact somewhat worsens, the discrepancy
modResIC = fit_mr( .dat = dic,
                     .label = "modIC",
                     .mods = modsS,
                     .write.to.csv = TRUE,
                     .write.table = TRUE,
                     .simple.return = FALSE )


# 9. WITH IPD AGE-MATCHING IN MLR ------------------------------------------------------------------

section = 9
setwd(data.dir)

# only 14 effect sizes now
nrow(drage)

# sanity check: should be very close to 0 now
mean(drage$mean_agec)

# ~ Subset model in MLR ------------------------------------------------------------------

# estimate in MLR further DECREASES from 0.35 to 0.162
# still less than in meta-analysis
# this fn also writes stats to results csv
( naiveAge.reps.only = fit_subset_meta( .dat = drage,
                                       .mods = "1",
                                       .label = "Reps subset naiveAge" ) )

# ~ Naive model with both sources ------------------------------------------------------------------
naiveResAge = fit_mr( .dat = dage,
                     .label = "naiveAge",
                     .mods = "isMeta",
                     .write.to.csv = TRUE,
                     .write.table = TRUE,
                     .simple.return = FALSE )

# ~ Moderated model with both sources ------------------------------------------------------------------
# discrepancies might decrease a little bit this time, but still 0.41
modResAge = fit_mr( .dat = dage,
                   .label = "modAge",
                   .mods = modsS,
                   .write.to.csv = TRUE,
                   .write.table = TRUE,
                   .simple.return = FALSE )

#bm: this analysis seems useful; should clean up the 2_XXX prep file
#  and merge this branch into master


# (Obsolete?) SUBSET MODELS ------------------------------------------------------------------ 
# 
# # ~ Sanity check: Subsets that resemble each other  ------------------------------------------------------------------
# 
# # MCF: "Could you do something like a "modal subsample" where you pull native language HPP (or CF?) results and plot those across MA/MB1, possibly with age included?"
# fit_subset_meta( .dat = dma %>% filter( method == "b.hpp" & test_lang == "a.native" ),
#                  .mods = "1",
#                  .label = NA )
# 
# fit_subset_meta( .dat = dr %>% filter( method == "b.hpp" & test_lang == "a.native" ),
#                  .mods = "1",
#                  .label = NA )
# 
# # interesting!!! these agree pretty well (replications actually a little bigger, 0.58 vs 0.51)
# # bm
# 
# 
# # **complements of above - do not agree well at all (difference 0.43, more like in meta-regression below)
# # 0.29 replications vs. 0.72 MA
# fit_subset_meta( .dat = dma %>% filter( method != "b.hpp" | test_lang != "a.native" ),
#                  .mods = "1",
#                  .label = NA )
# 
# fit_subset_meta( .dat = dr %>% filter( method != "b.hpp" | test_lang != "a.native" ),
#                  .mods = "1",
#                  .label = NA )
# 
# 
# # continues to agree very well if we look just at HPP subset
# fit_subset_meta( .dat = dma %>% filter( method == "b.hpp" ),
#                  .mods = "1",
#                  .label = NA )
# 
# fit_subset_meta( .dat = dr %>% filter( method == "b.hpp" ),
#                  .mods = "1",
#                  .label = NA )
# 
# 
# # surviving mods were mean age, language, and method
# # compare to meta-regression
# 
# # first one is estimated difference of 0.24 
# fit_mr( .dat = d,
#         .mods = c("isMeta", "method", "test_lang"),
#         .simple.return = FALSE)
# 
# # recode method and language as binary
# d$method.hpp = d$method == "b.hpp"
# d$test_lang.native = d$test_lang == "a.native"
# 
# fit_mr( .dat = d,
#         .mods = c("isMeta", "method.hpp", "test_lang.native"),
#         .simple.return = FALSE)
# 

