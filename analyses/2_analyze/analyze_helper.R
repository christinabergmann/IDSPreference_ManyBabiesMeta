
################################ MISCELLANEOUS ################################

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
fit_mr = function( .dat,
                   .label = NA,  # name of the analysis
                   .mods, 
                   .write.table = FALSE,
                   .write.to.csv = FALSE ) {
  
  # # TEST ONLY
  # .dat = d
  # .mods = mod.sets[[1]]
  # .label = "naive"
  # .write.table = TRUE
  # .write.to.csv = TRUE
  
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
  ests = round( est, 2 )
  pvals2 = format_stat(pval)
  # get rid of scientific notation; instead use more digits
  #pvals2[ pval < 0.01 ] = round( pval[ pval < 0.01 ], 3 )
  
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
  
  ##### 4. Get Phat for Meta, Replications, and Differencee #####
  # ***will fill later if tau^2 > 0
  
  
  ##### 5. Write Results #####
  # row = data.frame( avgDiff = est.ma - est.rep,
  #                   Phat0Diff = runif(n=1, -1,1), # ***obviously fake
  #                   Phat0.2Diff = runif(n=1, -1,1) ) # ***obviously fake
  # return(row)
  
  # return as a numeric vector for compatibility with boot()
  return( c( est.ma - est.rep,
                runif(n=1, -1,1), # ***obviously fake
                runif(n=1, -1,1) ) ) # ***obviously fake
  
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
  bootCIs = lapply( 1:n.ests, function(x) boot.ci(boot.res, type = type, index = x) )
  
  # list with one entry per estimate
  # the middle index "4" on the bootCIs accesses the stats vector
  # the final index chooses the CI lower (4) or upper (5) bound
  bootCIs = lapply( 1:n.ests, function(x) c( bootCIs[[x]][[4]][4],
                                             bootCIs[[x]][[4]][5] ) )
}


################################### MODIFIED FROM BOOT PACKAGE ################################### 

# see section called "MM additions"
# I minimally modified the function so that it can proceed even if some of the bootstrap iterates run into errors
# (in this case, Fisher convergence issues) because the boot package version gets confused about dimension mismatches

# source internal boot package functions
# source("bootfuns.R")

my_boot = function (data, statistic, R, sim = "ordinary", stype = c("i", 
                                                                    "f", "w"), strata = rep(1, n), L = NULL, m = 0, weights = NULL, 
                    ran.gen = function(d, p) d, mle = NULL, simple = FALSE, ..., 
                    parallel = c("no", "multicore", "snow"), ncpus = getOption("boot.ncpus", 
                                                                               1L), cl = NULL) 
{
  
  
  call <- match.call()
  stype <- match.arg(stype)
  if (missing(parallel)) 
    parallel <- getOption("boot.parallel", "no")
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") 
      have_snow <- TRUE
    if (!have_mc && !have_snow) 
      ncpus <- 1L
    loadNamespace("parallel")
  }
  if (simple && (sim != "ordinary" || stype != "i" || sum(m))) {
    warning("'simple=TRUE' is only valid for 'sim=\"ordinary\", stype=\"i\", n=0', so ignored")
    simple <- FALSE
  }
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  n <- NROW(data)
  if ((n == 0) || is.null(n)) 
    stop("no data in call to 'boot'")
  temp.str <- strata
  strata <- tapply(seq_len(n), as.numeric(strata))
  t0 <- if (sim != "parametric") {
    if ((sim == "antithetic") && is.null(L)) 
      L <- empinf(data = data, statistic = statistic, stype = stype, 
                  strata = strata, ...)
    if (sim != "ordinary") 
      m <- 0
    else if (any(m < 0)) 
      stop("negative value of 'm' supplied")
    if ((length(m) != 1L) && (length(m) != length(table(strata)))) 
      stop("length of 'm' incompatible with 'strata'")
    if ((sim == "ordinary") || (sim == "balanced")) {
      if (isMatrix(weights) && (nrow(weights) != length(R))) 
        stop("dimensions of 'R' and 'weights' do not match")
    }
    else weights <- NULL
    if (!is.null(weights)) 
      weights <- t(apply(matrix(weights, n, length(R), 
                                byrow = TRUE), 2L, normalize, strata))
    if (!simple) 
      i <- index.array(n, R, sim, strata, m, L, weights)
    original <- if (stype == "f") 
      rep(1, n)
    else if (stype == "w") {
      ns <- tabulate(strata)[strata]
      1/ns
    }
    else seq_len(n)
    t0 <- if (sum(m) > 0L) 
      statistic(data, original, rep(1, sum(m)), ...)
    else statistic(data, original, ...)
    rm(original)
    t0
  }
  else statistic(data, ...)
  pred.i <- NULL
  fn <- if (sim == "parametric") {
    ran.gen
    data
    mle
    function(r) {
      dd <- ran.gen(data, mle)
      statistic(dd, ...)
    }
  }
  else {
    if (!simple && ncol(i) > n) {
      pred.i <- as.matrix(i[, (n + 1L):ncol(i)])
      i <- i[, seq_len(n)]
    }
    if (stype %in% c("f", "w")) {
      f <- freq.array(i)
      rm(i)
      if (stype == "w") 
        f <- f/ns
      if (sum(m) == 0L) 
        function(r) statistic(data, f[r, ], ...)
      else function(r) statistic(data, f[r, ], pred.i[r, 
      ], ...)
    }
    else if (sum(m) > 0L) 
      function(r) statistic(data, i[r, ], pred.i[r, ], 
                            ...)
    else if (simple) 
      function(r) statistic(data, index.array(n, 1, sim, 
                                              strata, m, L, weights), ...)
    else function(r) statistic(data, i[r, ], ...)
  }
  RR <- sum(R)
  res <- if (ncpus > 1L && (have_mc || have_snow)) {
    if (have_mc) {
      parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
    }
    else if (have_snow) {
      list(...)
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", 
                                             ncpus))
        if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parLapply(cl, seq_len(RR), fn)
        parallel::stopCluster(cl)
        res
      }
      else parallel::parLapply(cl, seq_len(RR), fn)
    }
  }
  else lapply(seq_len(RR), fn)
  #t.star <- matrix(, RR, length(t0))  # ~~~ MM commented out
  
  # ~~~~~ MM added
  # number of non-NULL elements of the results vector
  RR = length(unlist(res))
  nulls = sapply( res, is.null)
  res = res[ !nulls ]
  t.star <- matrix(, RR, length(t0))

  # without this, boot.CI gets confused about number of replicates
  R = RR
  # ~~~~~ end of MM additions
  # print(length(res))
  # print(length(RR))
  # if ( length(res) != length(RR) ) browser()
  
  for (r in seq_len(RR)) t.star[r, ] <- res[[r]]
  if (is.null(weights)) 
    weights <- 1/tabulate(strata)[strata]
  boot.return(sim, t0, t.star, temp.str, R, data, statistic, 
              stype, call, seed, L, m, pred.i, weights, ran.gen, mle)
}
