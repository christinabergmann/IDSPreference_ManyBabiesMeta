

 

############################## PRELIMINARIES ############################## 

library(tidyverse) 
library(knitr)
library(here)
library(tableone)
library(corrr)
library(robumeta)
library(MetaUtility)
library(fastDummies)
library(weightr)
library(PublicationBias)
library(xtable)
library(boot)
library(testthat)

data.dir = here("data")
# where to save results
results.dir = here("results_from_R")
# results.dir = "~/Dropbox/Personal computer/Independent studies/2020/Christina's ManyBabiesMeta (MB-Meta)/IDSPreference_ManyBabiesMeta/results_from_R"
overleaf.dir = "~/Dropbox/Apps/Overleaf/MB-Meta/R_objects"
code.dir = here("analyses/2_analyze")

# helper fns
setwd(code.dir)
source("analyze_helper.R")
# # source internal fns from boot package
# # it's a long story
# # see my_boot() in helper fns if you like long stories
# source( here("analyses/2_analyze/bootfuns.R") )

# should we remove existing results file instead of overwriting individual entries? 
start.res.from.scratch = FALSE
# should we use the grateful package to scan and cite packages?
cite.packages.anew = FALSE
# should we bootstrap from scratch or read in old resamples?
boot.from.scratch = FALSE

# wipe results csvs if needed
if ( start.res.from.scratch == TRUE ) wr()

# constants of universe
digits = 2
pval.cutoff = 10^-4  # threshold for using "<"
boot.reps = 2000 # for cross-model comparisons 

# read in dataset
setwd(data.dir)
d = suppressMessages( read_csv("mb_ma_combined_scrambled_prepped.csv") )

# dataset with just the meta-analysis
dma = d %>% filter(isMeta == TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                       0. CHARACTERISTICS OF INCLUDED STUDIES            
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# for updating result csv
section = 0

############################## BASICS ############################## 

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


############################## MODERATORS ##############################

# look at moderators
mods = c( "study_type",
          "mean_agec",
          "test_lang",  # whether stimuli were in native language; almost constant in meta
          "native_lang",
          #"prop_nae", # ~~~ not motivated, was basis for test_lang
          "method",
          
          # constant in RRR:
          "speech_type",
          "own_mother",
          "presentation",
          "dependent_measure",
          "main_question_ids_preference",
         #"stimulus_set", # ~~~ not in the dataset 
          
          # varied in RRR and MA:          
          #"trial_control", # ~~~ not in the dataset
         #"human_coded", # ~~~ not in the dataset
)
# distribution of moderators in RRR and MA
t = CreateTableOne(vars = mods, 
                   strata = "study_type",
                   data = d)

xtable( print(t, noSpaces = TRUE, printToggle = FALSE) )


##### Moderator Correlation Matrix #####

# temporarily recode all categorical moderators as dummies
cat.mods = mods[ !mods == "mean_agec" ]
temp = dummy_cols( d[,cat.mods],
                   select_columns = cat.mods,
                   remove_selected_columns = TRUE,
                   ignore_na =  TRUE )  # removes the "parent" column
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                       1. NAIVE AND MODERATED META-REGRESSIONS           
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

section = 1

############################## FIT NAIVE AND MODERATED META-REGRESSIONS ############################## 

# try to fit meta-regression
mods2 = c( "isMeta",  # code this way since we expect meta to have larger effect sizes
           "mean_agec",
           "test_lang",  # whether stimuli were in native language; almost constant in meta
           "method",
           
           # constant in RRR:
           "speech_type",
           "own_mother",
           "presentation",
           #"dependent_measure",  # causes singularity
           "main_question_ids_preference" )

# varied in RRR:
#"trial_control" )  # causes singularity

mod.sets = list( c("isMeta"),  # naive model
                 mods2 )

labels = c("naive",
           "mod1")

# fit the naive and moderated models
# these write their results to the results csv file and table
naiveRes = fit_mr( .dat = d,
                   .label = "naive",
                   .mods = mod.sets[[1]],
                   .write.to.csv = TRUE,
                   .write.table = TRUE,
                   .simple.return = FALSE )

mod1Res = fit_mr( .dat = d,
                  .label = "mod1",
                  .mods = mod.sets[[2]],
                  .write.to.csv = TRUE,
                  .write.table = TRUE,
                  .simple.return = FALSE )


############################## CROSS-MODEL INFERENCE ##############################

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
                     
                     # WORKS with boot
                     # multi-argument returns need to be via c(), not list or df or whatever
                     tryCatch({
                       fit_mr( .dat = b,
                               .mods = mod.sets[[2]] )
                     }, error = function(err){
                       return( c(NA, NA, NA) )
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


# number of non-failed bootstrap reps
( boot.reps.successful = sum( !is.na( boot.res$t[,1] ) ) )

# get the bootstrapped CIs for all 3 statistics respectively using BCa method
( bootCIs = get_boot_CIs(boot.res, n.ests = 3) )

# order of statistics in the above: 
# between-source difference in estimated means
# between-source difference in Phat0
# between-source difference in Phat0.2

# # naive model metrics of source discrepancy
# # CIS WRONG HERE
# statCI_result_csv( "naive avgDiff", c(naiveRes$avgDiff, bootCIs[[1]][1], bootCIs[[1]][2]) )
# 
# statCI_result_csv( "naive Phat0Diff", c(naiveRes$Phat0Diff, bootCIs[[2]][1], bootCIs[[2]][2]) )
# 
# statCI_result_csv( "naive Phat0.2Diff", c(naiveRes$Phat0.2Diff, bootCIs[[3]][1], bootCIs[[3]][2]) )
# 
# 
# # moderated model metrics of metrics of source discrepancy
# statCI_result_csv( "mod1 avgDiff", c(mod1Res$avgDiff, bootCIs[[1]][1], bootCIs[[1]][2]) )
# 
# statCI_result_csv( "mod1 Phat0Diff", c(mod1Res$Phat0Diff, bootCIs[[2]][1], bootCIs[[2]][2]) )
# 
# statCI_result_csv( "mod1 Phat0.2Diff", c(mod1Res$Phat0.2Diff, bootCIs[[3]][1], bootCIs[[3]][2]) )


# between-model difference
statCI_result_csv( "moderator reduction avgDiff", c(naiveRes$avgDiff - mod1Res$avgDiff, bootCIs[[1]][1], bootCIs[[1]][2]) )

statCI_result_csv( "moderator reduction Phat0Diff", c(naiveRes$Phat0Diff - mod1Res$Phat0Diff, bootCIs[[2]][1], bootCIs[[2]][2]) )

statCI_result_csv( "moderator reduction Phat0.2Diff", c(naiveRes$Phat0.2Diff - mod1Res$Phat0.2Diff, bootCIs[[3]][1], bootCIs[[3]][2]) )


# ***p-values would require resampling under H0...

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                           2. PUBLICATION BIAS           
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

section = 2

##### Sanity Check: Estimate in Meta-Analysis Alone #####
( meta = robu( yi ~ 1, 
               data = dma, 
               studynum = as.factor(study_id),
               var.eff.size = vi,
               modelweights = "HIER",
               small = TRUE) )

# sanity check: should be similar to naive meta-regression model
expect_equal( as.numeric(meta$b.r), round(naiveRes$est.ma, digits), tol = 0.01 )
expect_equal( as.numeric(meta$reg_table$CI.L), round(naiveRes$est.ma.lo, digits), tol = 0.01 )
expect_equal( as.numeric(meta$reg_table$CI.U), round(naiveRes$est.ma.hi, digits), tol = 0.01 )

##### Affirmative and Nonaffirmative Counts #####

t = d %>% group_by(isMeta) %>%
  summarise( k.affirm = sum(affirm),
             k.nonaffirm = sum(affirm == FALSE),
             Paffirm = mean(affirm) )

update_result_csv( name = "meta k nonaffirmative",
                   value = t$k.nonaffirm[ t$isMeta == TRUE ] )

update_result_csv( name = "meta k affirmative",
                   value = t$k.affirm[ t$isMeta == TRUE ] )

##### Hedges Selection Model #####
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



##### Worst-Case Meta-Analysis ######

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



statCI_result_csv( "Worst mu",
                   c(meta.worst$b.r,
                     meta.worst$reg_table$CI.L,
                     meta.worst$reg_table$CI.U) )


update_result_csv( name = "Worst mu pval",
                   value = round(pval.worst, 3) )


# ##### S-values ######
# omitted because not very informative given the worst-case results above
# # s-values to reduce to null
# ( res = svalue( yi = d$logRR,
#                 vi = d$varlogRR,
#                 q = log(1), 
#                 clustervar = d$authoryear,
#                 model = "robust" ) )
# # N.P. for both point estimate and CI
# update_result_csv( name = "sval est to 1",
#                    section = 2,
#                    value = res$sval.est,
#                    print = FALSE )
# update_result_csv( name = "sval CI to 1",
#                    section = 2,
#                    value = res$sval.ci,
#                    print = FALSE )
# 
# 
# # s-values to reduce effect size to RR=1.1
# ( res = svalue( yi = d$logRR,
#                 vi = d$varlogRR,
#                 q = log(1.1), 
#                 clustervar = d$authoryear,
#                 model = "robust" ) )
# # N.P. for estimate
# update_result_csv( name = "sval est to 1.1",
#                    section = 2,
#                    value = res$sval.est,
#                    print = FALSE )
# update_result_csv( name = "sval CI to 1.1",
#                    section = 2,
#                    value = round( res$sval.ci, 2 ),
#                    print = FALSE )


##### Significance Funnel ######
significance_funnel(yi = dma$yi,
                    vi = dma$vi,
                    est.N = mu.worst,
                    est.all = ests,
                    favor.positive = TRUE)

setwd(results.dir)
ggsave( "funnel.pdf",
        width = 8,
        height = 6 )

setwd(overleaf.dir)
ggsave( "funnel.pdf",
        width = 8,
        height = 6 )


# *** If worst-case meta-analysis estimate was less than replication estimate, report about amount of publication bias required to shift meta-analysis to match replications.



