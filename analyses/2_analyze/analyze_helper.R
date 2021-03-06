
################################ MISCELLANEOUS ################################


my_ggsave = function(name,
                     width,
                     height,
                     .results.dir = results.dir,
                     .overleaf.dir = overleaf.dir) {
  
  setwd(.results.dir)
  ggsave( name,
          width = width, 
          height = height)
  
  setwd(.overleaf.dir)
  ggsave( name,
          width = width, 
          height = height)
}

# wrapper for update_result_csv to easily at stat and confidence interval
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
                              .section = section,  # defaults to global var
                              value = NA,
                              print = FALSE,
                              .results.dir = results.dir,
                              .overleaf.dir = overleaf.dir ) {
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
  if ( exists("overleaf.dir") ) {
    setwd(.overleaf.dir)
    write.csv( res, 
               "stats_for_paper.csv",
               row.names = FALSE,
               quote = FALSE )
  }
  
  if ( print == TRUE ) {
    View(res)
  }
}


# stands for "wipe results"
wr = function(){
  setwd(results.dir)
  if( "stats_for_paper.csv" %in% list.files() ) system("rm stats_for_paper.csv")
  setwd(overleaf.dir)
  if( "stats_for_paper.csv" %in% list.files() ) system("rm stats_for_paper.csv")
}

# stands for "view results"
vr = function(){
  setwd(results.dir)
  View( read.csv("stats_for_paper.csv") )
}



################################ FIT META-REGRESSION AND SAVE RESULTS ################################

# for a given set of moderators:
#  - fit meta-regression
#  - get estimate for meta-analysis and for replications
#  - estimate the 3 differences of interest (source coef, difference in Phats)
#  and optionally return them
#  - write a table with all the meta-regression estimates (optionally)
#  - write certain desired stats directly to results csv (optionally)

# .simple.return: should fn return only the stats to be bootstrapped (as a numeric vector)
#  or a more informative dataframe?

# Phats for MA and reps are based on SUBSET models
# .dat is an argument only to allow for bootstrapping
# .mods SHOULD include isMeta
# in returned stats, differences are always meta-analysis - replications
fit_mr = function( .dat,
                   .label = NA,  # name of the analysis
                   .mods, 
                   .write.table = FALSE,
                   .write.to.csv = FALSE,
                   .simple.return = TRUE ) {
  
  # #@TEST ONLY
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
    if ( length(.mods) > 1 ) {
      mR = robu( eval( parse( text = string3 ) ), 
                 data = .dat %>% filter(isMeta == FALSE), 
                 studynum = as.factor(study_id),
                 var.eff.size = vi,
                 modelweights = "HIER",
                 small = TRUE)
      
      calibR = conditional_calib_ests(mR)$calib.shift
    } else {
      calibR = calib_ests( yi = .dat$yi[ .dat$isMeta == FALSE ],
                           sei = .dat$sei[ .dat$isMeta == FALSE ] ) 
    }

    Phat0.rep = mean( calibR > 0 )
    Phat0.2.rep = mean( calibR > 0.2 )
    

    # meta-analysis only
    if ( length(.mods) > 1 ) {
      mM = robu( eval( parse( text = string3 ) ), 
                 data = .dat %>% filter(isMeta == TRUE), 
                 studynum = as.factor(study_id),
                 var.eff.size = vi,
                 modelweights = "HIER",
                 small = TRUE )
      
      calibM = conditional_calib_ests(mM)$calib.shift
    } else {
      calibM = calib_ests( yi = .dat$yi[ .dat$isMeta == TRUE ],
                           sei = .dat$sei[ .dat$isMeta == TRUE ] ) 
    }
   
    Phat0.ma = mean( calibM > 0 )
    Phat0.2.ma = mean( calibM > 0.2 )
    
    Phat0.diff = Phat0.ma - Phat0.rep
    Phat0.2.diff = Phat0.2.ma - Phat0.2.rep
    
    
    ##### 5. Write Results #####
    # row = data.frame( avgDiff = est.ma - est.rep,
    #                   Phat0Diff = runif(n=1, -1,1), # ***obviously fake
    #                   Phat0.2Diff = runif(n=1, -1,1) ) # ***obviously fake
    # return(row)
    
    
    
    if ( .simple.return == FALSE ) {
      .res$avgDiff = meta$b.r[ meta$labels == "isMetaTRUE"]
      .res$avgDiffLo = meta$reg_table$CI.L[ meta$labels == "isMetaTRUE"]
      .res$avgDiffHi = meta$reg_table$CI.U[ meta$labels == "isMetaTRUE"]
      
      #bm
      .res$Phat0.rep = Phat0.rep
      .res$Phat0.ma = Phat0.ma
      .res$Phat0.2.rep = Phat0.2.rep
      .res$Phat0.2.ma = Phat0.2.ma
      .res$Phat0Diff = Phat0.diff
      .res$Phat0.2Diff = Phat0.2.diff
      
    }
    
    
    if ( .simple.return == TRUE ){
      # return as a numeric vector for compatibility with boot()
      return( c( est.ma - est.rep,
                 
                 Phat0.rep,
                 Phat0.ma,
                 
                 Phat0.2.rep,
                 Phat0.2.ma,
                 
                 Phat0.diff,
                 Phat0.2.diff ) ) 
      
    } else return(.res)
  }
  
  
}







fit_subset_meta = function( .dat,
                            .label = NA,  # name of the analysis for the results csv
                            .mods,
                            .simple.return = FALSE) {
  
  # # TEST ONLY
  # .dat = dma
  # #.mods = mod.sets[[2]]  # @can't do this within just reps or just MA
  # .mods = "1"
  # .label = "Reps subset naive"
  # .write.table = FALSE
  # .write.to.csv = FALSE
  # .simple.return = FALSE
  
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



# marginal Phat (no moderators)
get_and_write_phat = function( .dat,
                               .q,
                               label ) {
  
  
  Phat = prop_stronger( q = .q,
                        tail = "above",
                        dat = .dat,
                        yi.name = "yi",
                        vi.name = "vi",
                        cluster.name = "study_id" )
  
  update_result_csv( name = paste( label, "est" ),
                     value = round( 100 * Phat$est ) )
  
  update_result_csv( name = paste( label, "lo" ),
                     value = round( 100 * Phat$lo ) )
  
  update_result_csv( name = paste( label, "hi" ),
                     value = round( 100 * Phat$hi ) )
  
}


# @AD HOC FOR THE MODERATORS ACTUALLY PRESENT IN EACH MODEL
# conditions on all moderators being set to modes in
conditional_calib_ests = function(.model){
  # get data from robu object
  dat = .model$data
  
  bhat = .model$b.r  # vector of coefficient estimates 
  t2 = .model$mod_info$tau.sq  # heterogeneity estimate
  
  # calculate the linear predictor (minus intercept) for each study
  # would need to be modified for other datasets
  otherBhat = bhat[-1]  # non-intercept terms
  otherBhatVars = .model$labels[-1]
  
  # @specific to this set of moderators
  
  # for each bhat, check if it's actually in the model
  # important because some subsets (e.g., boot resamples) could be 
  #  homo on any given moderator
  linpred = bhat[1]
  
  for ( i in 1:length(otherBhatVars) ) {
    
    b = otherBhatVars[i]
    if ( b == "mean_agec" ) {
      linpred = linpred + dat$mean_agec * otherBhat[i]
      
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
  
  # # threshold, shifted to set effect modifiers to 0
  # # would need to be modified for other datasets
  # q.shift = q - ( bhat[2]*z[1] + bhat[3]*z[2] )  
  # # note: below assumes we are considering effects ABOVE the threshold
  # Phat = mean( calib.shift > c(q.shift) )
  
}


# ##### sanity check #1:
# x = conditional_calib_ests(cond.reps.only)
# mod = cond.reps.only  # see analyze.R
# dat = dr
# t2 = mod$mod_info$tau.sq
# 
# myLinpred = mod$b.r[1] + mod$b.r[2] * dat$mean_agec +
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
# myLinpred = mod$b.r[1] + mod$b.r[2] * dat$mean_agec +
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

