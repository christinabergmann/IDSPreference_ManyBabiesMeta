
################################ MISCELLANEOUS ################################

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

# .simple.return: should fn return only the 3 stats to be bootstrapped (as a numeric vector)
#  or a more informative dataframe?
fit_mr = function( .dat,
                   .label = NA,  # name of the analysis
                   .mods, 
                   .write.table = FALSE,
                   .write.to.csv = FALSE,
                   .simple.return = TRUE ) {
  
  # # TEST ONLY
  # .dat = d
  # .mods = mod.sets[[1]]
  # .label = "naive"
  # .write.table = FALSE
  # .write.to.csv = FALSE
  # .simple.return = FALSE
  
  
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
  
  # save this one separately since it will ultimately be returned
  est.rep = est[ meta$labels == "X.Intercept."]
  
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
  if ( .simple.return == FALSE ) .res = data.frame( est.rep = est.rep,
                                                    est.rep.lo = mu.lo[ meta$labels == "X.Intercept." ], 
                                                    est.rep.hi = mu.hi[ meta$labels == "X.Intercept." ] )
  
  
  
  ##### 2. Write Meta-Regression Table #####
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
  }
  
  ##### 3. Get Estimated Mean for Meta-Analysis #####
  # could use linear combo of coefficients and vars, but df are complicated
  #  for robumeta
  # so easiest way is to just refit the model, reversing coding of study type variable
  .mods2 = .mods
  .mods2[ .mods == "isMeta" ] = "isRep"  # now intercept will be for meta-analysis
  linpred.string2 = paste( .mods2, collapse=" + ")
  string2 = paste( "yi ~ ", linpred.string2, collapse = "")
  ( meta2 = robu( eval( parse( text = string2 ) ), 
                  data = d, 
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
  
  ##### 4. Get Phat for Meta, Replications, and Differencee #####
  # ***will fill later if tau^2 > 0
  
  
  ##### 5. Write Results #####
  # row = data.frame( avgDiff = est.ma - est.rep,
  #                   Phat0Diff = runif(n=1, -1,1), # ***obviously fake
  #                   Phat0.2Diff = runif(n=1, -1,1) ) # ***obviously fake
  # return(row)
  
  # sanity check
  expect_equal( est.ma - est.rep, meta$b.r[ meta$labels == "isMetaTRUE"] )
  
  .res$avgDiff = meta$b.r[ meta$labels == "isMetaTRUE"]
  # .res$avgDiffLo = meta$reg_table$CI.L[ meta$labels == "isMetaTRUE"]
  # .res$avgDiffHi = meta$reg_table$CI.U[ meta$labels == "isMetaTRUE"]
  # .res$avgDiff
  
  .res$Phat0Diff = runif(n=1, -1,1) # ***obviously fake
  # .res$Phat0DiffLo = runif(n=1, -1,1) # ***obviously fake
  # .res$Phat0DiffHi = runif(n=1, -1,1) # ***obviously fake
  
  .res$Phat0.2Diff = runif(n=1, -1,1) # ***obviously fake
  
  
  if ( .simple.return == TRUE ){
    # return as a numeric vector for compatibility with boot()
    return( c( est.ma - est.rep,
               runif(n=1, -1,1), # ***obviously fake
               runif(n=1, -1,1) ) ) # ***obviously fake
    
  } else return(.res)
  
}




fit_subset_meta = function( .dat,
                            .label = NA,  # name of the analysis
                            .mods ) {
  
  # # TEST ONLY
  # .dat = dma
  # #.mods = mod.sets[[2]]  # @can't do this within just reps or just MA
  # .mods = "1"
  # .label = "Reps subset naive"
  # .write.table = FALSE
  # .write.to.csv = FALSE
  # .simple.return = FALSE
  
  
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
  

  return(meta)
}



#bm
get_and_write_phat = function( .dat,
                               .q,
                               prefix ) {
  
  
  Phat = prop_stronger( q = .q,
                        tail = "above",
                        dat = .dat,
                        yi.name = "yi",
                        vi.name = "vi",
                        cluster.name = "study_id" )
  
  update_result_csv( name = paste( prefix, "est" ),
                     value = round( 100 * Phat$est ) )
  
  update_result_csv( name = paste( prefix, "lo" ),
                     value = round( 100 * Phat$lo ) )
  
  update_result_csv( name = paste( prefix, "hi" ),
                     value = round( 100 * Phat$hi ) )
  
}





