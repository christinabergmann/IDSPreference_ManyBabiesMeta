

 
############################## PRELIMINARIES ############################## 

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


# should we remove existing results file instead of overwriting individual entries? 
start.res.from.scratch = FALSE
# should we use the grateful package to scan and cite packages?
cite.packages.anew = FALSE
# should we bootstrap from scratch or read in old resamples?
boot.from.scratch = TRUE

# wipe results csvs if needed
if ( start.res.from.scratch == TRUE ) wr()

# constants of universe
digits = 2
pval.cutoff = 10^-4  # threshold for using "<"
boot.reps = 2000 # for cross-model comparisons 

# stop sci notation
options(scipen=999)

# read in dataset
setwd(data.dir)
d = suppressMessages( read_csv("mb_ma_combined_prepped.csv") )

# dataset with just the meta-analysis
dma = d %>% filter(isMeta == TRUE)

# dataset with just the replications
dr = d %>% filter(isMeta == FALSE)

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
          "native_lang", #Won't be used in main analyses
          "method",
          
          # constant in RRR:
          "speech_type",
          "own_mother",
          "presentation",
          "dependent_measure",
          "main_question_ids_preference" )

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




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                       1. NAIVE AND MODERATED META-REGRESSIONS           
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

section = 1

############################## FIT NAIVE AND MODERATED META-REGRESSIONS ############################## 
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
naiveRes = fit_mr( .dat = d,
                   .label = "naive",
                   .mods = mod.sets[[1]],
                   .write.to.csv = TRUE,
                   .write.table = TRUE,
                   .simple.return = FALSE )



# fit the meta-regression with all covariates 
#  and remove them in prespecified order if needed until 
#  model becomes estimables

gotError = TRUE  # initialize so the while-loop is entered


# this will print an xtable with the moderator estimates for use in the paper
while ( gotError == TRUE ) {
  
  tryCatch({
    mod1Res = fit_mr( .dat = d,
                      .label = "mod1",
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




############################## SUBSET MODELS ##############################

# subset to MA and replications separately

( naive.MA.only = fit_subset_meta( .dat = dma,
                               .mods = "1",
                               .label = "MA subset naive" ) )

( naive.reps.only = fit_subset_meta( .dat = dr,
                                 .mods = "1",
                                 .label = "Reps subset naive" ) )



############################## PROPORTION MEANINGFULLY STRONG EFFECTS ##############################

get_and_write_phat( .dat = dma,
                     .q = 0,
                     label = "Phat0MA")

get_and_write_phat( .dat = dma,
                    .q = 0.2,
                    label = "Phat0.2MA")

get_and_write_phat( .dat = dr,
                    .q = 0,
                    label = "Phat0R")

get_and_write_phat( .dat = dr,
                    .q = 0.2,
                    label = "Phat0.2R")




calibR = calib_ests(yi = dr$yi, sei = sqrt(dr$vi) )
mean(calibR>0.2)
mean(calibR)

calibMA = calib_ests(yi = dma$yi, sei = sqrt(dma$vi) )
mean(calibMA>0.2)
mean(calibMA)

d$calib = NA
d$calib[ d$isMeta == FALSE ] = calibR
d$calib[ d$isMeta == TRUE ] = calibMA


############################## DENSITY PLOT OF META-ANALYSIS VS. REPLICATION CALIBRATED ESTIMATES ##############################

colors = c("black", "orange")

#@move to prep
d$studyTypePretty = NA
d$studyTypePretty[ d$isMeta == TRUE ] = "Meta-analysis"
d$studyTypePretty[ d$isMeta == FALSE ] = "Replications"

# choose axis scaling
summary(d$yi)
xmin = -1
xmax = 3.5
tickJump = 0.5  # space between tick marks

ggplot( data = d,
        aes( x = calib,
             fill = studyTypePretty,
             color = studyTypePretty ) ) +
  
  # mean estimates from subset models
  geom_vline( xintercept = naive.MA.only$b.r,
              color = colors[1],
              lty = 2 ) +
  
  geom_vline( xintercept = naive.reps.only$b.r,
              color = colors[2],
              lty = 2) +
  
  # shifted threshold for Z=z
  geom_vline(xintercept = 0,
             color = "gray",
             lty = 1) +
  
  # ensemble estimates shifted to Z=0
  geom_density(alpha = 0.3) +
  
  theme_bw() +
  
  xlab("Calibrated study estimate (SMD)") +
  scale_x_continuous( limits = c(xmin, xmax), breaks = seq(xmin, xmax, tickJump)) +
  
  ylab("Density") +
  
  scale_color_manual( values = colors, name = "Source" ) +
  scale_fill_manual( values = colors, name = "Source" ) +
  
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


my_ggsave( name = "calibrated_plot.pdf",
           width = 8,
           height = 5 )



############################## CROSS-MODEL INFERENCE ##############################


# refit mod1Res for each boot iterate, each time estimating avgDiff and two currently fake stats
#  for Phat

#bm
# @temp only
# @increase later
boot.reps = 500

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
# between-source difference in Phat0 (**FAKE)
# between-source difference in Phat0.2 (**FAKE)

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
# framed as the INCREASED discrepancy after controlling for moderators
statCI_result_csv( "moderator increase avgDiff", c(mod1Res$avgDiff - naiveRes$avgDiff, bootCIs[[1]][1], bootCIs[[1]][2]) )

# SAVE
# currently these are fake
# statCI_result_csv( "moderator reduction Phat0Diff", c(naiveRes$Phat0Diff - mod1Res$Phat0Diff, bootCIs[[2]][1], bootCIs[[2]][2]) )
# 
# statCI_result_csv( "moderator reduction Phat0.2Diff", c(naiveRes$Phat0.2Diff - mod1Res$Phat0.2Diff, bootCIs[[3]][1], bootCIs[[3]][2]) )


# ***p-values would require resampling under H0...

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                           2. PUBLICATION BIAS           
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

section = 2

# @MOVE THIS
# sanity check: should be similar to naive meta-regression model
#  but not necessarily equal
expect_equal( as.numeric(naive.MA.only$b.r), round(naiveRes$est.ma, digits), tol = 0.03 )
expect_equal( as.numeric(naive.MA.only$reg_table$CI.L), round(naiveRes$est.ma.lo, digits), tol = 0.03 )
expect_equal( as.numeric(naive.MA.only$reg_table$CI.U), round(naiveRes$est.ma.hi, digits), tol = 0.03 )

##### Affirmative and Nonaffirmative Counts #####

t = d %>% group_by(isMeta) %>%
  summarise( k.affirm = sum(affirm),
             k.nonaffirm = sum(affirm == FALSE),
             Paffirm = mean(affirm) )

update_result_csv( name = "MA k nonaffirmative",
                   value = t$k.nonaffirm[ t$isMeta == TRUE ] )

update_result_csv( name = "MA k affirmative",
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


##### S-values ######
# s-values to reduce to null
( Sval0 = svalue( yi = dma$yi,
                vi = dma$vi,
                q = 0,
                clustervar = dma$study_id,
                favor.positive = TRUE,
                model = "robust" ) )



##### Meta-Analysis with Eta = 4.70 (Post Hoc) ######
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


##### Significance Funnel ######
# significance_funnel(yi = dma$yi,
#                     vi = dma$vi,
#                     est.N = mu.worst,
#                     est.all = as.numeric(naive.MA.only$b.r),
#                     favor.positive = TRUE)
# 
# setwd(results.dir)
# ggsave( "funnel.pdf",
#         width = 8,
#         height = 6 )
# 
# setwd(overleaf.dir)
# ggsave( "funnel.pdf",
#         width = 8,
#         height = 6 )


# modified version
yi = d$yi
vi = d$vi
xmin = min(yi)
xmax = max(yi)
ymin = 0  # so that pooled points are shown
ymax = max( sqrt(vi) )
xlab = "Point estimate"
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


##### Make the Plot #####
# replications, originals
colors = c("darkgray", "red")

p.funnel = ggplot( data = d, aes( x = d$yi,
                                  y = d$sei,
                                  color = d$isMeta,
                                  shape = d$affirm ) )

if ( plot.pooled == TRUE ) {
  
  # plot the pooled points
  p.funnel = p.funnel + geom_point(
    data = pooled.pts,
    aes( x = yi, y = sei ),
    size = 4,
    shape = 5,
    fill = NA,
    color = c("red", "black", "orange", "darkgray")
  ) +
    
    geom_point(
      data = pooled.pts,
      aes( x = yi, y = sei ),
      size = 4,
      shape = 18,
      color =  c("red", "black", "orange", "darkgray"),
      alpha = 1
    ) +
    
    # just for visual separation of pooled ests
    geom_hline( yintercept = 0 ) +
    
    # diagonal "just significant" line
    geom_abline(slope=sl,intercept = int, color = "gray")
}

p.funnel = p.funnel +
  
  # semi-transparent points with solid circles around them
  geom_point( size = 3, alpha=.4) +
  geom_point( size = 3, shape = 1) +
  
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(1,19)) +
  
  xlab(xlab) +
  ylab(ylab) +
  
  scale_x_continuous( limits = c(xmin, xmax), breaks = seq(-0.5, 3, .5) ) +
  scale_y_continuous( limits = c(ymin, ymax), breaks = seq(0, .6, .1) ) +
  
  theme_classic() +
  #theme(legend.title=element_blank())
  theme(legend.position = "none")

plot(p.funnel)

setwd(overleaf.dir)
ggsave(plot = p.funnel,
       filename = "mb_signif_funnel.png",
       width = 5,
       height = 4)




