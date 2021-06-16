
################################ MISCELLANEOUS ################################

# save ggplot to results.dir and overleaf.dir for piping joy
my_ggsave = function(name,
                     width,
                     height,
                     .results.dir = results.dir,
                     .overleaf.dir = overleaf.dir) {
  
  if(!dir.exists(.results.dir)){
    dir.create(.results.dir)
  } 
  setwd(.results.dir)
  ggsave( name,
          width = width, 
          height = height)
  
  if(dir.exists(.overleaf.dir)){
    setwd(.overleaf.dir)
    ggsave( name,
            width = width, 
            height = height)
  }
}

# wrapper for update_result_csv to easily write stat and confidence interval
statCI_result_csv = function(name,
                             vec){
  
  update_result_csv( name = paste( name, c("", "lo", "hi"), sep = " " ),
                     value = round( vec, digits ) )
}



# for reproducible manuscript-writing
# adds a row to the file "stats_for_paper" with a new statistic or value for the manuscript
# optionally, "section" describes the section of code producing a given result
# expects global vars: overleaf.dir, results.dir
update_result_csv = function( name,
                              .section = if( exists("section") ) section else "",  # defaults to global var
                              value = NA,
                              print = FALSE,
                              .results.dir = results.dir,
                              .overleaf.dir = overleaf.dir ) {
  if(!dir.exists(.results.dir)){
    dir.create(.results.dir)
  } 
  setwd(.results.dir)
  
  new.rows = data.frame( name,
                         value = as.character(value),
                         section = as.character(.section) )
  
  # to avoid issues with variable types when overwriting
  new.rows$name = as.character(new.rows$name)
  new.rows$value = as.character(new.rows$value)
  new.rows$section = as.character(new.rows$section)
  
  
  if ( "stats_for_paper.csv" %in% list.files() ) {
    res = read.csv( "stats_for_paper.csv",
                    stringsAsFactors = FALSE,
                    colClasses = rep("character", 3 ) )
    
    # if this entry is already in the results file, overwrite the
    #  old one
    if ( all(name %in% res$name) ) res[ res$name %in% name, ] = new.rows
    else res = rbind(res, new.rows)
  }
  
  if ( !"stats_for_paper.csv" %in% list.files() ) {
    res = new.rows
  }
  
  write.csv( res, 
             "stats_for_paper.csv",
             row.names = FALSE,
             quote = FALSE )
  
  # also write to Overleaf
  if ( dir.exists(.overleaf.dir) ) {
    setwd(.overleaf.dir)
    write.csv( res, 
               "stats_for_paper.csv",
               row.names = FALSE,
               quote = FALSE )
  }
  
  if ( print == TRUE ) {
    View(res)
  }
  
  # also save as global object 
  #  useful for running sanity checks
  resCSV <<- res
}


# stands for "wipe results"
# removes file stats_for_paper.csv
wr = function(){
  setwd(results.dir)
  if( "stats_for_paper.csv" %in% list.files() ) system("rm stats_for_paper.csv")
  if(dir.exists(.overleaf.dir)){
    setwd(overleaf.dir)
    if( "stats_for_paper.csv" %in% list.files() ) system("rm stats_for_paper.csv")
  }
}

# stands for "view results"
# displays contents of stats_for_paper.csv
vr = function(){
  setwd(results.dir)
  View( read.csv("stats_for_paper.csv") )
}



################################ FIT META-REGRESSION AND SAVE RESULTS ################################

# For a given set of moderators (.mods):
#  - fits meta-regression
#  - gets estimate for meta-analysis and for replications
#  - estimates the 3 differences of interest (source coef, difference in Phats)
#  and optionally return them
#  - writes a table with all the meta-regression estimates (optionally)
#  - writes certain desired stats directly to results csv (optionally)

# Arguments:
# - .dat: dataset (used as an argument only to allow for bootstrapping,
#    because fn already  internally makes MA and MLR subsets)

# - .label: name of analysis for labeling files and results that will be written

# - .mods: moderators to adjust (use only "isMeta" for naive) 

# - .write.table: do you want a csv file of regression coeffs to be written to results.dir?

# - .write.to.scv: do you want estimates written to stats_for_paper.csv?

# - .simple.return: should fn return only the stats to be bootstrapped
#   (as an unlabeled numeric vector) or a more informative dataframe that's better for human use?

# **Phats for MA and reps are based on SUBSET models
# **In returned stats, differences are always meta-analysis - replications

# Sanity checks for this fn are in analyze.R

fit_mr = function( .dat,
                   .label = NA,  
                   .mods, 
                   .write.table = FALSE,
                   .write.to.csv = FALSE,
                   .simple.return = TRUE ) {
  
  # TEST ONLY:
  # .dat = d
  # .mods = modsS
  # .label = "naive"
  # .write.table = FALSE
  # .write.to.csv = FALSE
  # .simple.return = FALSE
  
  # find out if dataset has both MLR and MA
  if ( length( unique(.dat$isMeta) ) != 2 ) {
    message("FYI, dataset has only one source (MLR or MA), so not looking at isMeta discrepancy")
    hasBoth = FALSE
  } else {
    hasBoth = TRUE
  }
  
  ##### 1. Fit Meta-Regression (All Data) #####
  linpred.string = paste( .mods, collapse=" + ")
  string = paste( "yi ~ ", linpred.string, collapse = "")
  
  ( meta = robu( eval( parse( text = string ) ), 
                 data = .dat, 
                 studynum = as.factor(study_id),
                 var.eff.size = vi,
                 modelweights = "HIER",
                 small = TRUE) )
  
  est = meta$b.r
  t2 = meta$mod_info$tau.sq
  mu.lo = meta$reg_table$CI.L
  mu.hi = meta$reg_table$CI.U
  mu.se = meta$reg_table$SE
  pval = meta$reg_table$prob
  V = meta$VR.r  # variance-covariance matrix
  
  # save this one separately since it will ultimately be returned
  if (hasBoth == TRUE) est.rep = est[ meta$labels == "X.Intercept."]
  
  # rounded and formatted estimates for text
  # expects pval.cutoff to be a global var
  ests = round( est, 2 )
  pvals2 = format.pval(pval, eps = pval.cutoff)
  
  # save results to csv file
  if ( .write.to.csv == TRUE ){
    update_result_csv( name = paste(.label, "tau"),
                       value = round( sqrt(t2), 2 ) )
    
    update_result_csv( name = paste(.label, "k"),
                       value = nrow(meta$data.full) )
    
    #@CHECK THIS because I'm not sure n is what I think it is
    # temporarily rounding because it's fractional in study_id Kaplan1995a
    update_result_csv( name = paste(.label, "n subj"),
                       value = round( sum(.dat$n) ) )
    
    update_result_csv( name = paste( .label, "est", meta$labels ),
                       value = ests )
    
    update_result_csv( name = paste( .label, "lo", meta$labels ),
                       value = round(mu.lo, digits) )
    
    update_result_csv( name = paste( .label, "hi", meta$labels ),
                       value = round(mu.hi, digits) )
    
    update_result_csv( name = paste( .label, "pval", meta$labels ),
                       value = pvals2 )
  }
  
  # also save selected results to a df to be returned
  if ( hasBoth == TRUE & .simple.return == FALSE ) .res = data.frame( est.rep = est.rep,
                                                                      est.rep.lo = mu.lo[ meta$labels == "X.Intercept." ], 
                                                                      est.rep.hi = mu.hi[ meta$labels == "X.Intercept." ])
  
  
  
  ##### 3. Write Meta-Regression Table #####
  if ( .write.table == TRUE ){
    CIs = format_CI( mu.lo,
                     mu.hi )
    temp = data.frame( Moderator = meta$labels, 
                       EstCI = paste( ests, CIs, sep = " " ),
                       Pval = pvals2 )
    
    # save results
    setwd(results.dir)
    if(!dir.exists("tables_to_prettify")){
      dir.create("tables_to_prettify")
    } 
    setwd("tables_to_prettify")
    
    write.csv( temp,
               paste("model_", .label, "_table.csv", sep = "" ),
               row.names = FALSE)
    
    # also print it as an xtable
    #browser()
    print( xtable(temp), include.rownames = FALSE)
    
  }
  
  
  if ( hasBoth == TRUE ) {
    
    ##### 3. Get Estimated Mean for Meta-Analysis #####
    # could use linear combo of coefficients and vars, but df are complicated
    #  for robumeta
    # so easiest way is to just refit the model, reversing coding of study type variable
    
    .mods2 = .mods
    .mods2[ .mods == "isMeta" ] = "isRep"  # now intercept will be for meta-analysis
    linpred.string2 = paste( .mods2, collapse=" + ")
    string2 = paste( "yi ~ ", linpred.string2, collapse = "")
    ( meta2 = robu( eval( parse( text = string2 ) ), 
                    data = .dat, 
                    studynum = as.factor(study_id),
                    var.eff.size = vi,
                    modelweights = "HIER",
                    small = TRUE) )
    
    
    pval = meta2$reg_table$prob[meta$labels == "X.Intercept."]
    pval2 = format_stat(pval)
    
    est.ma = meta2$b.r[meta2$labels == "X.Intercept."]
    
    # save results to csv
    if ( .write.to.csv == TRUE ){
      update_result_csv( name = paste( .label, "avg for meta" ),
                         value = round( est.ma, digits ) )
      
      update_result_csv( name = paste( .label, "avg lo for meta" ),
                         value = round( meta2$reg_table$CI.L[meta$labels == "X.Intercept."], digits ) )
      
      update_result_csv( name = paste( .label, "avg hi for meta" ),
                         value = round( meta2$reg_table$CI.U[meta$labels == "X.Intercept."], digits ) )
      
      update_result_csv( name = paste( .label, "avg pval for meta" ),
                         value = pval2 )
      
      
    }
    
    # also save selected results to a df to be returned
    if ( .simple.return == FALSE ){
      .res$est.ma = est.ma
      .res$est.ma.lo = meta2$reg_table$CI.L[meta$labels == "X.Intercept."]
      .res$est.ma.hi = meta2$reg_table$CI.U[meta$labels == "X.Intercept."]
    }
    
    # sanity check
    expect_equal( est.ma - est.rep, meta$b.r[ meta$labels == "isMetaTRUE"] )
    
    
    ##### 4. Get Phat for Meta, Replications, and Difference #####
    # from subset models
    
    # remove "isMeta" from moderators
    if ( length(.mods) > 1 ) {
      .mods3 = .mods[ !.mods == "isMeta" ]
      linpred.string3 = paste( .mods3, collapse=" + ")
      string3 = paste( "yi ~ ", linpred.string3, collapse = "")
    } else {
      # handle naive model
      string3 = "yi ~ 1"
    }
    
    # replications only
    # if there are moderators, get conditional calibrated estimates:
    if ( length(.mods) > 1 ) {
      mR = robu( eval( parse( text = string3 ) ), 
                 data = .dat %>% filter(isMeta == FALSE), 
                 studynum = as.factor(study_id),
                 var.eff.size = vi,
                 modelweights = "HIER",
                 small = TRUE)
      
      calibR = conditional_calib_ests(mR)$calib.shift
      
      # if no moderators, get marginal calibrated estimates
    } else {
      calibR = calib_ests( yi = .dat$yi[ .dat$isMeta == FALSE ],
                           sei = .dat$sei[ .dat$isMeta == FALSE ] ) 
    }
    
    Phat0.rep = mean( calibR > 0 )
    Phat0.2.rep = mean( calibR > 0.2 )
    
    
    # meta-analysis only
    # if there are moderators, get conditional calibrated estimates:
    if ( length(.mods) > 1 ) {
      mM = robu( eval( parse( text = string3 ) ), 
                 data = .dat %>% filter(isMeta == TRUE), 
                 studynum = as.factor(study_id),
                 var.eff.size = vi,
                 modelweights = "HIER",
                 small = TRUE )
      
      calibM = conditional_calib_ests(mM)$calib.shift
      
      # if no moderators, get marginal calibrated estimates
    } else {
      calibM = calib_ests( yi = .dat$yi[ .dat$isMeta == TRUE ],
                           sei = .dat$sei[ .dat$isMeta == TRUE ] ) 
    }
    
    Phat0.ma = mean( calibM > 0 )
    Phat0.2.ma = mean( calibM > 0.2 )
    
    Phat0.diff = Phat0.ma - Phat0.rep
    Phat0.2.diff = Phat0.2.ma - Phat0.2.rep
    
    if ( .simple.return == FALSE ) {
      .res$avgDiff = meta$b.r[ meta$labels == "isMetaTRUE"]
      .res$avgDiffLo = meta$reg_table$CI.L[ meta$labels == "isMetaTRUE"]
      .res$avgDiffHi = meta$reg_table$CI.U[ meta$labels == "isMetaTRUE"]
      
      .res$Phat0.rep = Phat0.rep
      .res$Phat0.ma = Phat0.ma
      .res$Phat0.2.rep = Phat0.2.rep
      .res$Phat0.2.ma = Phat0.2.ma
      .res$Phat0Diff = Phat0.diff
      .res$Phat0.2Diff = Phat0.2.diff
      
    }
    
    # be careful about changing this return structure
    #  various parts of the bootstrapping in analyze.R expect this exact structure
    if ( .simple.return == TRUE ){
      # return as a numeric vector for compatibility with boot()
      return( c( est.ma,
                 est.rep,
                 est.ma - est.rep,
                 
                 Phat0.ma,
                 Phat0.rep,
                 Phat0.diff,
                 
                 Phat0.2.ma,
                 Phat0.2.rep,
                 Phat0.2.diff ) ) 
      
    } else return(.res)
  }
  
}


# fit meta-analysis to arbitrary subset, .dat
# see fit_mr for args
# if you want no moderators at all, use .mods = "1"
# this fn has sanity checks in analyze.R
fit_subset_meta = function( .dat,
                            .label = NA,  
                            .mods,
                            .simple.return = FALSE) {
  
  # catch case where subset is empty
  if ( nrow(.dat) == 0 ) {
    message("Subset is empty")
    
    # return something whose structure is compatible with .simple.return = TRUE
    return( data.frame( intercept = NA,
                        intercept.lo = NA,
                        intercept.hi = NA,
                        k = NA ) )
  }
  
  
  ##### 1. Fit Meta-Regression #####
  linpred.string = paste( .mods, collapse=" + ")
  string = paste( "yi ~ ", linpred.string, collapse = "")
  
  ( meta = robu( eval( parse( text = string ) ), 
                 data = .dat, 
                 studynum = as.factor(study_id),
                 var.eff.size = vi,
                 modelweights = "HIER",
                 small = TRUE) )
  
  est = meta$b.r
  t2 = meta$mod_info$tau.sq
  mu.lo = meta$reg_table$CI.L
  mu.hi = meta$reg_table$CI.U
  mu.se = meta$reg_table$SE
  pval = meta$reg_table$prob
  V = meta$VR.r  # variance-covariance matrix
  
  # rounded and formatted estimates for text
  # expects pval.cutoff to be a global var
  ests = round( est, 2 )
  pvals2 = format.pval(pval, eps = pval.cutoff)
  
  # save results to csv file
  if ( !is.na(.label) ) {
    update_result_csv( name = paste(.label, "tau"),
                       value = round( sqrt(t2), 2 ) )
    
    update_result_csv( name = paste(.label, "k"),
                       value = nrow(meta$data.full) )
    
    #@CHECK THIS because I'm not sure n is what I think it is
    # temporarily rounding because it's fractional in study_id Kaplan1995a
    update_result_csv( name = paste(.label, "n subj"),
                       value = round( sum(.dat$n) ) )
    
    update_result_csv( name = paste( .label, "est", meta$labels ),
                       value = ests )
    
    update_result_csv( name = paste( .label, "lo", meta$labels ),
                       value = round(mu.lo, digits) )
    
    update_result_csv( name = paste( .label, "hi", meta$labels ),
                       value = round(mu.hi, digits) )
    
    update_result_csv( name = paste( .label, "pval", meta$labels ),
                       value = pvals2 )
  }
  
  if ( .simple.return == FALSE ) return(meta)
  
  if ( .simple.return == TRUE ) {
    
    # if you update this return structure, also update the placeholder at top
    #  for when subset is empty
    return( data.frame( intercept = est[1],
                        intercept.lo = mu.lo[1],
                        intercept.hi = mu.hi[1],
                        k = nrow(.dat) ) )
  }
  
  
}


# **IMPORTANT: this fn hard-codes the moderators to be conditioned
#   and their level names and orders

# give conditional calibrated estimates when setting all moderators
#  to their means/modes in MA
# sanity checks below
conditional_calib_ests = function(.model){
  # get data from robu object
  dat = .model$data
  
  bhat = .model$b.r  # vector of coefficient estimates 
  t2 = .model$mod_info$tau.sq  # heterogeneity estimate
  
  # calculate the linear predictor (minus intercept) for each study
  # would need to be modified for other datasets
  otherBhat = bhat[-1]  # non-intercept terms
  otherBhatVars = .model$labels[-1]
  
  # for each bhat, check if it's actually in the model
  # important because some subsets (e.g., boot resamples) could be 
  #  homo on any given moderator
  linpred = bhat[1]
  
  for ( i in 1:length(otherBhatVars) ) {
    
    b = otherBhatVars[i]
    if ( b == "mean_agec_mos" ) {
      linpred = linpred + dat$mean_agec_mos * otherBhat[i]
      
    } else if (b == "test_langb.nonnative" ) {
      linpred = linpred + ( dat$test_lang == "b.nonnative" ) * otherBhat[i]
      
    } else if (b == "test_langc.artificial" ) {
      linpred = linpred + ( dat$test_lang == "c.artificial" ) * otherBhat[i]
      
    } else if (b == "methodb.hpp" ) {
      linpred = linpred + ( dat$method == "b.hpp" ) * otherBhat[i]
      
    } else if (b == "methodc.other" ) {
      linpred = linpred + ( dat$method == "c.other" ) * otherBhat[i]
      
    } else {
      browser()
      stop( "Variable name was not found in dataset" )
    }
  }
  
  
  ##### Calculate Calibrated Estimates #####
  # calibrated estimate, shifted to set effect modifiers to 0
  calib.shift = c(bhat[1]) + sqrt( c(t2) / ( c(t2) + dat$vi) ) * ( dat$yi - linpred )
  
  # sanity check
  #var(calib.shift); t2
  
  return( data.frame( calib.shift, linpred ) )
}


# ##### sanity check #1:
# x = conditional_calib_ests(cond.reps.only)
# mod = cond.reps.only  # see analyze.R
# dat = dr
# t2 = mod$mod_info$tau.sq
# 
# myLinpred = mod$b.r[1] + mod$b.r[2] * dat$mean_agec_mos +
#   mod$b.r[3] * ( dat$test_lang == "b.nonnative" ) +
#   mod$b.r[4] * ( dat$method == "b.hpp" )
# 
# # check linear predictor calculation
# expect_equal( myLinpred, x$linpred)
# 
# # calculate conditional calibrated estimates manually
# myCalib = mod$b.r[1] + sqrt( c(t2) / ( c(t2) + dat$vi) ) * ( dat$yi - myLinpred ) 
# expect_equal(myCalib, x$calib.shift)
# # yesss
# 
# ##### sanity check #2: different dataset that is homo
# #  on different covariates
# x = conditional_calib_ests(cond.MA.only)
# mod = cond.MA.only  # see analyze.R
# dat = dma
# t2 = mod$mod_info$tau.sq
# 
# myLinpred = mod$b.r[1] + mod$b.r[2] * dat$mean_agec_mos +
#   mod$b.r[3] * ( dat$test_lang == "b.nonnative" ) +
#   mod$b.r[4] * ( dat$test_lang == "c.artificial" ) +
#   mod$b.r[5] * ( dat$method == "b.hpp" ) + 
#   mod$b.r[6] * ( dat$method == "c.other" )
# 
# # check linear predictor calculation
# expect_equal( myLinpred, x$linpred)
# 
# # calculate conditional calibrated estimates manually
# myCalib = mod$b.r[1] + sqrt( c(t2) / ( c(t2) + dat$vi) ) * ( dat$yi - myLinpred ) 
# expect_equal(myCalib, x$calib.shift)
# # yesss


# do a quick sensitivity analysis and write the following results to csv:
#  - average effect size in replication subset
#  - average effect size in MA subset
#  - fit_mr stats from naive model
#  - fit_mr stats from moderated model
# expects a few variables to be defined globally: section, modsS
# has minimal sanity check in analyze.R
quick_sens_analysis = function( .dat,
                                .suffix ) {
  # # test only
  # .dat = dic
  # .suffix = "IC"
  
  ### Naive, replication subset
  .naive.reps.only = fit_subset_meta( .dat = .dat[ .dat$isMeta == FALSE, ],
                                      .mods = "1",
                                      .label = paste( "Reps subset naive",
                                                      .suffix, 
                                                      sep = "" ) )
  
  cat("\n\n------------- NAIVE, REPS ONLY:\n")
  print(.naive.reps.only)
  
  ### Naive, meta-analysis subset
  .naive.meta.only = fit_subset_meta( .dat = .dat[ .dat$isMeta == TRUE, ],
                                      .mods = "1",
                                      .label = paste( "Meta subset naive",
                                                      .suffix, 
                                                      sep = "" ) )
  
  cat("\n\n------------- NAIVE, META ONLY:\n")
  print(.naive.meta.only)
  
  ### Naive, both sources
  .naiveRes = fit_mr( .dat = .dat,
                      .label = paste( "naive",
                                      .suffix, 
                                      sep = "" ),
                      .mods = "isMeta",
                      .write.to.csv = TRUE,
                      .write.table = TRUE,
                      .simple.return = FALSE )
  
  cat("\n\n------------- NAIVE, BOTH SOURCES:\n")
  print(.naiveRes)
  
  ### Moderated, both sources 
  .modRes = fit_mr( .dat = .dat,
                    .label = paste( "mod",
                                    .suffix, 
                                    sep = "" ),
                    .mods = modsS,
                    .write.to.csv = TRUE,
                    .write.table = TRUE,
                    .simple.return = FALSE )
  
  cat("\n\n------------- MODERATED, BOTH SOURCES:\n")
  print(.modRes)
}




################################ FNS FOR BOOTSTRAPPING ################################

# this is Davison & Hinkley's recommendation 
#  see section 3.8 in "Further Topics" chapter
cluster_bt = function(.dat,
                      .clustervar){
  .dat$cluster = .dat[[.clustervar]]
  
  # resample clusters, leaving observations intact
  #https://stats.stackexchange.com/questions/46821/bootstrapping-hierarchical-multilevel-data-resampling-clusters
  # see answer by dash2
  cluster.ids = data.frame(cluster = sample(.dat$cluster, replace = TRUE))
  datb = .dat %>% inner_join(cluster.ids, by = 'cluster')
  return(datb)
}

# # sanity check
# d %>% group_by(study_id) %>% 
#   summarise( mean(yi) )
# # since we're resampling clusters, below should have different clusters but same 
# #  mean for a given cluster as the above:
# b = cluster_bt(.dat = d, "study_id")
# b %>% group_by(study_id) %>% 
#   summarise( mean(yi) )

# get BCa-bootstrapped confidence intervals for a vector of parameters (not just one)
# boot.res: an object from boot()
# type: what kind of bootstrapping to use? (as in boot.ci())
# n.ests: how many parameters were estimated?
get_boot_CIs = function(boot.res,
                        type = "bca",
                        n.ests) {
  bootCIs = lapply( 1:n.ests, function(x) safe_boot_ci(x, boot.res, type) )
  
  return(bootCIs)
}

# for use by get_boot_CIs
# x is an index
safe_boot_ci = function(x, boot.res, type) {
  tryCatch({
    # get CIs for just this index (x)
    CIs = boot.ci(boot.res, type = type, index = x)
    # the middle index "4" on the bootCIs accesses the stats vector
    # the final index chooses the CI lower (4) or upper (5) bound
    c( CIs[[4]][4],
       CIs[[4]][5] )
    
  }, error = function(err){
    # in case CI isn't estimable
    c(NA, NA)
  })
}  # end fn used lapply


# ~ DENSITY PLOT -----------------------------------------

# density plot of age stratified by source
age_densities = function(.dat) {
  # choose axis scaling
  #summary(.dat$mean_agec_mos)
  xmin = -24
  xmax = 30
  tickJump = 6  # space between tick marks
  
  ggplot( data = .dat,
          aes( x = mean_agec_mos,
               fill = studyTypePretty,
               color = studyTypePretty ) ) +
    
    # ensemble estimates shifted to Z=0
    geom_density(alpha = 0.3) +
    
    theme_bw() +
    
    xlab("Mean age (centered; months)") +
    scale_x_continuous( limits = c(xmin, xmax),
                        breaks = seq(xmin, xmax, tickJump)) +
    
    ylab("Density") +
    
    scale_color_manual( values = rev(colors), name = "Source" ) +
    scale_fill_manual( values = rev(colors), name = "Source" ) +
    
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
}


# MODIFIED PUBLICATION BIAS FNS -----------------------------------------

# fns are taken straight from PublicationBias package,
#  but redefine affirmative indicator to use the variable reportedSignif
#  search "EDITED" to find differences from R package
corrected_meta_2 = function( yi,
                             vi,
                             reportedSignif,  # EDITED
                             eta,
                             clustervar = 1:length(yi),
                             model,
                             selection.tails = 1,
                             favor.positive,
                             alpha.select = 0.05,
                             CI.level = 0.95,
                             small = TRUE ) {
  
  # stop if eta doesn't make sense
  if ( eta < 1 ) stop( "Eta must be at least 1.")
  
  # number of point estimates
  k = length(yi)
  
  # calculate alpha for inference on point estimate
  alpha = 1 - CI.level
  
  # warn if clusters but user said fixed
  nclusters = length( unique( clustervar ) )
  if ( nclusters < k & model == "fixed" ) {
    warning( "Clusters exist, but will be ignored due to fixed-effects specification. To accommodate clusters, instead choose model = robust.")
  }
  
  ##### Flip Estimate Signs If Needed #####
  # if favor.positive == TRUE, then we don't need to fit a naive meta-analysis or do anything
  if ( favor.positive == TRUE ) {
    # keep track of whether we flipped for reporting at the end
    flipped = FALSE
    yif = yi
  } else {
    flipped = TRUE
    yif = -yi
  }
  
  # 2-sided p-values for each study even if 1-tailed selection
  pvals = 2 * ( 1 - pnorm( abs(yif) / sqrt(vi) ) )
  
  # affirmative indicator based on selection tails
  # EDITED
  if ( selection.tails == 1 ) A = (reportedSignif == 1) & (yif > 0)
  if ( selection.tails == 2 ) A = (reportedSignif == 1)
  
  k.affirmative = sum(A)
  k.nonaffirmative = k - sum(A)
  
  # eta condition is because svalue automatically calls corrected_meta with eta=1
  #  to fit the naive meta-analysis, and in that case it doesn't matter
  # bm
  if ( eta > 1 ) {
    if ( k.affirmative == 0 | k.nonaffirmative == 0 ) {
      stop( "There are zero affirmative studies or zero nonaffirmative studies. Model estimation cannot proceed.")
    }
  }
  
  
  dat = data.frame( yi, yif, vi, A, clustervar )
  
  
  ##### Fixed-Effects Model #####
  if ( model == "fixed" ) {
    
    # FE mean and sum of weights stratified by affirmative vs. nonaffirmative
    strat = dat %>% group_by(A) %>%
      summarise( nu = sum( 1 / vi ),
                 ybar = sum( yi / vi ) )
    
    # components of bias-corrected estimate by affirmative status
    ybarN = strat$ybar[ strat$A == 0 ]
    ybarS = strat$ybar[ strat$A == 1 ]
    nuN = strat$nu[ strat$A == 0 ]
    nuS = strat$nu[ strat$A == 1 ]
    
    # corrected pooled point estimate
    est = ( eta * ybarN + ybarS ) / ( eta * nuN + nuS )
    
    # inference
    var = ( eta^2 * nuN + nuS ) / ( eta * nuN + nuS )^2
    se = sqrt(var)
    
    # z-based inference
    if ( small == FALSE ) {
      lo = est - qnorm( 1 - (alpha/2) ) * sqrt(var)
      hi = est + qnorm( 1 - (alpha/2) ) * sqrt(var)
      z =  abs( est / sqrt(var) )
      pval.est = 2 * ( 1 - pnorm( z ) )
    }
    
    # t-based inference
    if ( small == TRUE ) {
      df = k - 1
      lo = est - qt( 1 - (alpha/2), df = df ) * sqrt(var)
      hi = est + qt( 1 - (alpha/2), df = df ) * sqrt(var)
      t =  abs( est / sqrt(var) )
      pval.est = 2 * ( 1 - pt( t, df = df ) )
    }
  } # end fixed = TRUE
  
  ##### Robust Independent and Robust Clustered #####
  if ( model == "robust" ) {
    
    # weight for model
    weights = rep( 1, length(pvals) )
    weights[ A == FALSE ] = eta
    
    # initialize a dumb (unclustered and uncorrected) version of tau^2
    # which is only used for constructing weights
    meta.re = rma.uni( yi = yi,
                       vi = vi)
    t2hat.naive = meta.re$tau2
    
    # fit weighted robust model
    meta.robu = robu( yi ~ 1,
                      studynum = clustervar,
                      data = dat,
                      userweights = weights / (vi + t2hat.naive),
                      var.eff.size = vi,
                      small = small )
    
    est = as.numeric(meta.robu$b.r)
    se = meta.robu$reg_table$SE
    lo = meta.robu$reg_table$CI.L
    hi = meta.robu$reg_table$CI.U
    pval.est = meta.robu$reg_table$prob
    eta = eta
  } # end robust = TRUE
  
  return( data.frame( est,
                      se,
                      lo,
                      hi,
                      pval = pval.est,
                      eta = eta,
                      k.affirmative,
                      k.nonaffirmative,
                      signs.recoded = flipped ) )
}



svalue_2 = function( yi,
                     vi,
                     reportedSignif, # EDITED
                     q,
                     clustervar = 1:length(yi),
                     model,
                     alpha.select = 0.05,
                     eta.grid.hi = 200,
                     favor.positive,
                     CI.level = 0.95,
                     small = TRUE,
                     return.worst.meta = FALSE ) {
  
  
  # stop if eta doesn't make sense
  if ( eta.grid.hi < 1 ) stop( "eta.grid.hi must be at least 1.")
  
  # number of point estimates
  k.studies = length(yi)
  
  alpha = 1 - CI.level
  
  # warn if clusters but user said fixed
  nclusters = length( unique( clustervar ) )
  if ( nclusters < k.studies & model == "fixed" ) {
    warning( "You indicated there are clusters, but these will be ignored due to fixed-effects specification. To accommodate clusters, instead choose model = robust.")
  }
  
  # fit uncorrected model
  # no need to pass alpha.select here
  m0 = corrected_meta( yi = yi,
                       vi = vi,
                       eta = 1,
                       model = model,
                       clustervar = clustervar,
                       selection.tails = 1,
                       favor.positive = favor.positive,
                       CI.level = CI.level,
                       small = small )
  
  # stop if q is on wrong side of null
  if ( m0$est > 0 & q > m0$est ) stop( paste( "The uncorrected pooled point estimate is ", round2(m0$est),
                                              ". q must be less than this value (i.e., closer to zero).",
                                              sep = "" ) )
  if ( m0$est < 0 & q < m0$est ) stop( paste( "The uncorrected pooled point estimate is ", round2(m0$est),
                                              ". q must be greater than this value (i.e., closer to zero).",
                                              sep = "" ) )
  
  ##### Flip Estimate Signs If Needed #####
  
  # if favor.positive == TRUE, then we don't need to fit a naive meta-analysis or do anything
  if ( favor.positive == TRUE ) {
    # keep track of whether we flipped for reporting at the end
    flipped = FALSE
  } else {
    flipped = TRUE
    yi = -yi
    q = -q
  }
  
  # 2-sided p-values for each study even if 1-tailed selection
  pvals = 2 * ( 1 - pnorm( abs(yi) / sqrt(vi) ) )
  
  # affirmative indicator under 1-tailed selection
  # EDITED 
  A = (reportedSignif == 1) & (yi > 0)
  
  k.affirmative = sum(A)
  k.nonaffirmative = k.studies - sum(A)
  
  if ( k.affirmative == 0 | k.nonaffirmative == 0 ) {
    stop( "There are zero affirmative studies or zero nonaffirmative studies. Model estimation cannot proceed.")
  }
  
  dat = data.frame( yi, vi, A, clustervar )
  
  
  ##### Fixed-Effects Model #####
  if ( model == "fixed" ) {
    
    if (k.nonaffirmative > 1){
      # first fit worst-case meta
      meta.worst = rma.uni( yi = yi,
                            vi = vi,
                            data = dat[ A == FALSE, ],
                            method = "FE" )
      
      
      est.worst = as.numeric(meta.worst$b)
      lo.worst = meta.worst$ci.lb
    }
    
    if (k.nonaffirmative == 1) {
      est.worst = dat$yi[ A == FALSE ]
      lo.worst = dat$yi[ A == FALSE ] - qnorm(0.975) * sqrt(dat$vi[ A == FALSE ])
    }
    
    # FE mean and sum of weights stratified by affirmative vs. nonaffirmative
    strat = dat %>% group_by(A) %>%
      summarise( nu = sum( 1 / vi ),
                 ybar = sum( yi / vi ) )
    
    # components of bias-corrected estimate by affirmative status
    ybarN = strat$ybar[ strat$A == 0 ]
    ybarA = strat$ybar[ strat$A == 1 ]
    nuN = strat$nu[ strat$A == 0 ]
    nuA = strat$nu[ strat$A == 1 ]
    
    # S-value for point estimate
    sval.est = ( nuA * q - ybarA ) / ( ybarN - nuN * q )
    
    # S-value for CI (to shift it to q)
    # match term names used in Wolfram Alpha
    a = ybarN
    b = ybarA
    c = nuN
    d = nuA
    
    if ( small == FALSE ) k = qnorm( 1 - (alpha/2) )
    if ( small == TRUE ) {
      df = k.studies - 1
      k = qt( 1 - (alpha/2), df = df )
    }
    
    # # version directly from Wolfram
    # termA = a^2 * d * k^2 - (2 * a * c * d * k^2 * q) +
    #           b^2 * c * k^2 -
    #           (2 * b * c * d * k^2 * q) +
    #           c^2 * d * k^2 * q^2 +
    #           c * d^2 * k^2 * q^2 -
    #           c * d * k^4
    
    # manually simplied version
    termA = k^2 * ( a^2 * d -
                      (2 * c * d * q) * (a + b) +
                      b^2 * c +
                      q^2 * (c^2 * d + d^2 * c) -
                      c * d * k^2 )
    
    termB = -a*b + a*d*q + b*c*q - c*d*q^2
    
    termC = a^2 - 2*a*c*q + c^2*q^2 - c*k^2
    
    sval.ci = ( -sqrt(termA) + termB ) / termC
    if ( sval.ci < 0 ) sval.ci = ( sqrt(termA) + termB ) / termC
    
    # # sanity check by inversion
    # # corrected CI limit
    # eta = sval.ci
    # termD = (eta * a + b) / (eta * c + d)
    # termE = k * sqrt( (eta^2 * c + d) / (eta * c + d)^2 )
    # expect_equal( termD - termE,
    #               q )
    # # WORKS!!!
    
  } # end fixed = TRUE
  
  
  ##### Robust Independent and Robust Clustered #####
  if ( model == "robust" ) {
    
    ##### Worst-Case Meta to See if We Should Search at All
    
    if (k.nonaffirmative > 1){
      # first fit worst-case meta to see if we should even attempt grid search
      # initialize a dumb (unclustered and uncorrected) version of tau^2
      # which is only used for constructing weights
      meta.re = rma.uni( yi = yi,
                         vi = vi)
      t2hat.naive = meta.re$tau2
      
      # fit model exactly as in corrected_meta
      meta.worst =  robu( yi ~ 1,
                          studynum = clustervar,
                          data = dat[ A == FALSE, ],
                          userweights = 1 / (vi + t2hat.naive),
                          var.eff.size = vi,
                          small = small )
      
      est.worst = as.numeric(meta.worst$b.r)
      lo.worst = meta.worst$reg_table$CI.L
    }
    
    # robumeta above can't handle meta-analyzing only 1 nonaffirmative study
    if (k.nonaffirmative == 1) {
      est.worst = dat$yi[ A == FALSE ]
      lo.worst = dat$yi[ A == FALSE ] - qnorm(0.975) * sqrt(dat$vi[ A == FALSE ])
    }
    
    ##### Get S-value for estimate
    if ( est.worst > q ) {
      sval.est = "Not possible"
    } else {
      
      # define the function we need to minimize
      # i.e., distance between corrected estimate and the target value of q
      func = function(.eta) {
        # EDITED TO CALL CORRECTED_META_2
        est.corr = corrected_meta_2( yi = yi,
                                     vi = vi,
                                     reportedSignif = reportedSignif,
                                     eta = .eta,
                                     model = model,
                                     clustervar = clustervar,
                                     selection.tails = 1,
                                     favor.positive = TRUE,  # always TRUE because we've already flipped signs if needed
                                     alpha.select = alpha.select,
                                     CI.level = CI.level,
                                     small = small )$est
        return( abs(est.corr - q))
      }
      
      opt = optimize( f = func,
                      interval = c(1, eta.grid.hi),
                      maximum = FALSE )
      sval.est = opt$minimum
      
      # discrepancy between the corrected estimate and the s-value
      diff = opt$objective
      
      # if the optimal value is very close to the upper range of grid search
      #  AND we're still not very close to the target q,
      #  that means the optimal value was above eta.grid.hi
      if ( abs(sval.est - eta.grid.hi) < 0.0001 & diff > 0.0001 ) sval.est = paste(">", eta.grid.hi)
    }
    
    # do something similar for CI
    if ( lo.worst > q ) {
      sval.ci = "Not possible"
      
    } else {
      # define the function we need to minimize
      # i.e., distance between corrected estimate and the target value of q
      func = function(.eta) {
        # EDITED TO CALL CORRECTED_META_2
        lo.corr = corrected_meta( yi = yi,
                                  vi = vi,
                                  reportedSignif = reportedSignif,
                                  eta = .eta,
                                  model = model,
                                  clustervar = clustervar,
                                  selection.tails = 1,
                                  favor.positive = TRUE, # always TRUE because we've already flipped signs if needed
                                  alpha.select = alpha.select,
                                  CI.level = CI.level,
                                  small = small )$lo
        return( abs(lo.corr - q))
      }
      
      opt = optimize( f = func,
                      interval = c(1, eta.grid.hi),
                      maximum = FALSE )
      sval.ci = opt$minimum
      
      # discrepancy between the corrected estimate and the s-value
      diff = opt$objective
      
      # if the optimal value is very close to the upper range of grid search
      #  AND we're still not very close to the target q,
      #  that means the optimal value was above eta.grid.hi
      if ( abs(sval.ci - eta.grid.hi) < 0.0001 & diff > 0.0001 ) sval.ci = paste(">", eta.grid.hi)
    }
    
  }
  
  # s-values less than 1 indicate complete robustness
  # is.numeric is in case we have a "< XXX" string instead of a number
  if ( is.numeric(sval.est) & !is.na(sval.est) & sval.est < 1) sval.est = "Not possible"
  if ( is.numeric(sval.ci) & !is.na(sval.ci) & sval.ci < 1) sval.ci = "Not possible"
  
  # m0 was fit BEFORE flipping signs
  # but q has now been flipped in the latter case in "or" statement below
  if ( (m0$est > 0 & m0$lo < q) | (m0$est < 0 & m0$hi > -q) ) {
    # important: Shiny website assumes that this exact string ("--") for CI can be interpreted as
    #  the naive CI's already containing q
    sval.ci = "--"
    message("sval.ci is not applicable because the naive confidence interval already contains q")
  }
  
  # meta.worst might not exist if, for example, there is only 1 nonaffirmative study
  #  or obviously if the user did not ask for it
  # in those cases, set it to NULL so that return structure can stay the same
  if ( return.worst.meta == FALSE | !exists("meta.worst") ) meta.worst = NULL
  
  
  # if user wanted meta.worst, but we don't have it
  if ( return.worst.meta == TRUE & is.null(meta.worst) ) message("Not returning a worst-case meta-analysis because there were fewer than 2 nonaffirmative studies.")
  
  # # meta.worst might not exist if, for example, there is only 1 nonaffirmative study
  # if ( return.worst.meta == TRUE & exists("meta.worst") ) {
  
  return( list( stats = data.frame( sval.est,
                                    sval.ci = sval.ci,
                                    k.affirmative,
                                    k.nonaffirmative,
                                    signs.recoded = flipped ),
                meta.worst = meta.worst ) )
  
}

