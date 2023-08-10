#' Calculates interaction scores between all reference spectra and spectral sets.
#' The FSE module is one way to produce likely feature shapes supported by statistical
#' association of spectral points across dataset spectra. Singlets contain
#' very little specific information, can match to any feature shape and are
#' therefore excluded.
#'
#' TINA attempts to reduce the number of duplicate features, including those identified
#' in different places on the spectral axis (i.e., misaligned).
#'
#' These features are then matched to full-resolution pure compound reference spectra
#' (e.g. GISSMO, or any other full-res source). The matching reference regions are
#' termed ref-feats, and indeed are hypothesized subsignatures. To test whether these
#' hypothetical subsignatures are present in the actual data, these are then back-fit
#' to the original spectral data using the feature position and quantification.
#'
#' Then, for each backfit (ref-feat to dataset spectrum), a backfit feasibility score
#' is calculated (see filter.matches() for details). This roughly indicates whether
#' or not the fit of the ref-feat is feasible (doesn't grossly exceed the spectral
#' data in any places). Perfectly feasible is 1, not feasible is 0.
#'
#' A good weighted rmse and rval therefore indicates a good shape match between the
#' feature and the ref-feat. A high bff score for a ref-feat and a given dataset
#' spectrum indicates that reference spectrum region has a feasible fit to the data.
#' These two pieces of information provide a basis for a believable association
#' between the reference spectrum and the dataset spectrum.
#'
#' However, there can be numerous versions of each ref-feat, and numerous ref-feats
#' associating a given reference spectrum with a particular dataset spectrum. In
#' fact, the ideal situation is numerous, highly feasible ref-feats associating the
#' two and accounting for 100% of a reference spectrum signature. Thus, for every
#' point along the reference spectrum, we take the highest bff score associated with
#' that point (from any ref-feat that covers a region including that point), with
#' the default being 0. We could sum all these values for a given reference spectrum
#' (again, focusing specifically on a single dataset spectrum). However, small peaks
#' are downweighted by multiplying this composite "point-wise best bff score" vector
#' by the % total reference spectrum intensity of each point. Thus we account for
#' the best-case feasibility of association of each ref spectral point, scaled by
#' that point's overall relevance to the spectral signature. Summing these values
#' across the score vector gives a total association score between a reference
#' spectrum and a dataset spectrum.
#' 
#' Update 27MAY2023: also records best ref feature for each ss x ref pair, to limit
#' evidence presentation to only those which contribute to scores. 
#' Upate 27JUL2023: no longer scoring for span of ref feature, but only for the 
#' non-NA points (exclude NA gaps during scoring. To do this, added feature compression.
#' Ref vectors could also be vastly compressed (individually or otherwise) to save
#' on memory if needed. 
#'
#' @param pars A list of parameters specifying input/output directories and parallelization settings.
#'
#'
#' @return A table of scores for all possible reference spectra and spectral sets.
#' Also returns rfs.used (ref feats actually contributing to each pair's scoref)
#' Combined in a list
#'
#' @export
#'
#'
#' @import parallel
#' @import foreach
#' @import yaml
#' @importFrom magrittr %>%
#' @import pbapply
#'
pair_score_summation <- function(pars, refmat){
  
emptyScore <- function(){
  list(
         bfs.used.tot = NA, # keep these around; they're the best evidence
         bfs.used.res = NA, # keep these around; they're the best evidence
         bfs.used.rmseb = NA, # keep these around; they're the best evidence
         bfs.used.min.score = NA, # keep these around; they're the best evidence
         pair.scores = data.frame(ss.spec = NA,
                                  ref = NA,
                                  score.tot = 0, 
                                  score.res = 0,
                                  score.rmseb = 0,
                                  score.min = 0)
  )
}

## Parallelize the ss.spec - reference pair score summation? ####      
      message('Reading in data for pair score summation...')

        tmpdir <- pars$dirs$temp
        this.run <- paste0(tmpdir)

      # split backfits and combinations by reference spectra used
      # report as pairs with coords in matrix
      
        ss.ref.pairs <- readRDS(paste0(this.run, "/ss.ref.pairs.RDS"))
        # refmat <- readRDS(paste0(this.run, "/temp_data_matching/ref.mat.RDS")) %>% t
      
      # # Normalize refs 
        refmat <- refmat / Rfast::rowsums(refmat, na.rm = T)
        
      # Read in feature data and compress
        feature <- readRDS(paste0(this.run, "/feature.final.RDS"))
        features.compressed <- feature$position %>% compress_stack
        rm(feature)
        
      message('Splitting data for parallelization...')
      # Pre-split the table and refs (reduce overhead) ####
        by.ref <- pbapply::pblapply(unique(ss.ref.pairs$ref), function(r) 
        {
          ref.pairs <- ss.ref.pairs[which(ss.ref.pairs$ref == r),]
          # We don't need the data for the fits, just the descriptive stats
          
          f.nums <- ref.pairs$feat %>% unique
          selection <- features.compressed %>% cstack_selectRows(f.nums)
          
          if (nrow(ref.pairs) > 0){
            
            return(list(ref.num = r,
                        refspec = refmat[r, ],
                        feat.models = selection,
                        feat.nums = f.nums,
                        ref.pairs = ref.pairs))
          }  
          return(NULL)
        })
        
      rm(refmat)
      
      # Good idea to do some load-balancing here... ####
        workloads <- lapply(by.ref, function(x) x$ref.pairs %>% nrow) %>% unlist
        plan <- distribute(workloads, across = pars$par$ncores, messages = F)
        
        chonks <- lapply(1:pars$par$ncores, function(core.number){
          chonk <- by.ref[plan$core.ids == core.number]
          return(chonk)
        })
        
        rm(by.ref)
      
      # Look for specific compound in iterations
      # cmpd <- 'R-Lactate'
      # 
      # refnum <- which(cmpd.names %in% cmpd) %>% .[1]
      # lapply(chonks, function(chonk){
      #   # chonk <- chonks[[1]]
      #   refs <- lapply(chonk, function(r.list){
      #     # r.list <- chonk[[1]]
      #     r.list$ref.num
      #   }) %>% unlist
      #   which(refs %in% refnum)
      # })

      message('Computing scores...')

          
    # Calculate the scores in parallel ####
    
        t1 <- Sys.time()
        
            
            score.list <- mclapply(chonks, mc.cores = pars$par$ncores, 
                                           FUN = function(r.list)                             
            {
              # chonk <- chonks[[3]]                               
              # Go through this list of by.refs
                                                                   
              lapply(chonk, function(r.list){
                tryCatch(
                  {
                        
                        # For each ref:
                        # r.list <- chonk[[209]]
                        
                        # Make sure there are actually spectra in the interactions list ####
                          if(nrow(r.list$ref.pairs) < 1){return(emptyScore())}
                        
                        # Extract the interactions info from the r.list object for this ref ####
                          r.num <- r.list$ref.num
                          refspec <- r.list$refspec #/sum(r.list$refspec, na.rm = T) # already normed
                          ref.pairs <- r.list$ref.pairs
                          feat.models <- cstack_expandRows(r.list$feat.models) %>% is.na %>% "!"()
                          # how to use feat.num to index feat.models?
                            ref.pairs$feat <- ref.pairs$feat %>% factor(levels = r.list$feat.nums) %>% as.integer
                            ref.pairs$rmse.biased <- 1-ref.pairs$rmse.biased
                            # NOTE: this is only temporary and ref.pairs$feat should ONLY be used 
                            # to pull from feat.models after this point in this function!!!
                            # MTJ double-checked this on 28JUL2023
                            
                        # Compute ss interactions (combinations) for this ref ####
                          comb <- expand.grid(unique(ref.pairs$ref), unique(ref.pairs$ss.spec))
                            colnames(comb) <- c("ref", "ss.spec")
        
                        # Set up cumulative score vector objs (will be recycled each sample iteration) ####
                      
                          bff.tot <- list(score.name = 'bff.tot',
                                          scores = rep(0, length(refspec)),
                                          rf.ids = rep(0, length(refspec)))
                          bff.res <- bff.tot
                            bff.res$score.name = 'bff.res'
                          rmseb <- bff.tot
                            rmseb$score.name = 'rmse.biased'
                          min.score <- bff.tot
                            min.score$score.name = 'min.score'
        
                          ref.pairs$min.score <- Rfast::rowMins(cbind(ref.pairs$bff.res,
                                                                      ref.pairs$bff.tot,
                                                                      ref.pairs$rmse.biased),value = TRUE)
                          
                        # Loop through ref-ss.spec combinations and calculate scores ####
                        # - score objects will be re-copied for each lapply iteration, so 
                        #   no need to reset between samples. 
                        # - only the updated score objects are retained for each ss - ref pair
                        # lapply(a, function(x) x$pair.scores) %>% do.call(rbind,.)
                        # lapply(a, function(x) x$bfs.used.tot %>% length) %>% unlist
                          # a <- pblapply(1:nrow(comb), 
                          lapply(1:nrow(comb),
                                 function(i){
                                   
                            ref <- comb$ref[i]
                            ss.spec <- comb$ss.spec[i]   
                            
                            # Compute the scores for this combination
                            tryCatch(
                              {
                                  
                                # Select relevant backfits and matches ####
                                  rp.rows <- which(ref.pairs$ss.spec == ss.spec)
            
                                # Loop through the matches associated with this ref - ss.spec pair ####
                                # Update the bff values in v with any higher bff at that point ####
                                        ## Plotting to show in action: ####
                                          # reg <- refspec %>% trim_sides(out = "inds")
                                          # rs <- refspec
                                          # rs[is.na(rs)] <- 0
                                          # # reg <- 19000:27000
                                          # jpeg(file=paste0("rf.scoring",'_ref_',ref,".jpeg"), width=1200, height=700)
                                          #   simplePlot(xvect = reg, ymat = rs[reg], linecolor = "blue", xdir = 'normal')
                                          # dev.off()
                                          # 
                                          # jpeg(file=paste0("rf.scoring",'_spec_',ss.spec,".jpeg"), width=1200, height=700)
                                          #   simplePlot(xvect = reg, ymat = xmat[ss.spec, reg], linecolor = "blue", xdir = 'normal')
                                          # dev.off()
                                          # 
                                          # j <- 0
                                          # idx <- 0
                                          # 
                                          # # 0th iteration
                                          # jpeg(file=paste0("rf.scoring",j,".jpeg"), width=1200, height=700)
                                          #   scattermore::scattermoreplot(x = reg, y = bff.tot$scores[reg], ylim=c(0,1), cex = 1)
                                          # dev.off()
                                  # Actual loop ####
                                  for (j in rp.rows){
                                    # Description of what this does ####
                                    # j indexes the rows of rp.rows (ss.ref.pairs belonging to this chunk = this ref #)
                                    # that match this dataset spectrum (ss.spec). These could include many matches. A 
                                    # match could also give rise to many rp.rows, because of the different spectra. rp.rows
                                    # could be derived from multiple matches as long as they involve the same ref and 
                                    # have backfits to ss.spec. 
                                        ## Turn on for plotting an example: ####
                                          # idx <- idx + 1
                                          # j <- rp.rows[idx]
                                    # Get the ref range for the matched ref-feat ####
                                      # Where is the ref range? It's in the backfit, now. Need this to calculate pct ref accounted,
                                      # but that is done after we decide which ref points will be included. For now we're just keeping
                                      # track of which ref points have been accounted for. 
                                        feat.model <- feat.models[ref.pairs$feat[j],] %>% .[ref.pairs$feat.start[j] : ref.pairs$feat.end[j]]
                                        ref.pts <- (ref.pairs$ref.start[j] : ref.pairs$ref.end[j]) %>% .[feat.model]
            
                                    # Update the scores objs for each type for this ref feature-ss match ####
                                      # bff total 
                                        
                                        bff.tot <- update_scoreObj(ref.pairs[j, ], bff.tot, ref.pts)
                                        
                                      # bff (resonance) 
                                      
                                        bff.res <- update_scoreObj(ref.pairs[j, ], bff.res, ref.pts)
                                        
                                      # rmse biased
                                      
                                        rmseb <- update_scoreObj(ref.pairs[j, ], rmseb, ref.pts)
                                        
                                      # Or: pick the worst score from each of them
                                      
                                        min.score <- update_scoreObj(ref.pairs[j, ], min.score, ref.pts)
                                      
                                        ## More plotting...can set a delay to update each iteration in slowed time if you want. ####
                                        ## I just print all of them to jpeg and animate the series in ppt.
                                        # Sys.sleep(.2)
                                        # jpeg(file=paste0("rf.scoring",j,".jpeg"), width=1200, height=700)
                                        #   scattermore::scattermoreplot(x = reg, y = bff.tot$scores[reg], ylim=c(0,1), cex = 1)
                                        # dev.off()
                                        
                                  }
                                
                                # Calculate summed score and report as data.frame of ss.spec-reference pairs ####
                                
                                  bff.tot <- sum_score(bff.tot, refspec)
                                  bff.res <- sum_score(bff.res, refspec)
                                  rmseb <- sum_score(rmseb, refspec)
                                  min.score <- sum_score(min.score, refspec)
                                  # if(any(c(bff.tot$scores.tot, bff.res$scores.tot, rmseb$scores.tot, min.score$scores.tot) > 0.9)){browser()}
                                  
                                # Return a list of score information for this sample-ref pair: ####
                                # same format as emptyScore()
                                  list(
                                         bfs.used.tot = bff.tot$rf.ids.tot, # keep these around; they're the best evidence
                                         bfs.used.res = bff.res$rf.ids.tot, # keep these around; they're the best evidence
                                         bfs.used.rmseb = rmseb$rf.ids.tot, # keep these around; they're the best evidence
                                         bfs.used.min.score = min.score$rf.ids.tot, # keep these around; they're the best evidence
                                         pair.scores = data.frame(ss.spec = ss.spec,
                                                                  ref = r.num,
                                                                  score.tot = bff.tot$scores.tot, 
                                                                  score.res = bff.res$scores.tot,
                                                                  score.rmseb = rmseb$scores.tot,
                                                                  score.min = min.score$scores.tot)
                                  )
                            
                              }, 
                              # Return 0 scores if failed in any way
                              error = function(cond){
                                  emptyScore()
                              }
                            )
                            
                        })
                            
                  }, 
                  # Return NA inds and 0 scores if failed in any way (single iteration)
                  error = function(cond){
                    emptyScore()
                  })  
              }) %>% unlist(recursive = F)   
              
            }) %>% unlist(recursive = F, use.names = F)
      
        print(Sys.time() - t1)
     
      # Extract out the scores data.frame rows
        message('\nextracting out the scores data.frame rows...')
        ss.ref.pair.scores <- score.list %>% lapply(function(x) {
          if (is.atomic(x$pair.scores)){
            warning('atomics in pair.scores in pair.score.summation(). excluding.')
            return(NULL)
          }
          x$pair.scores
        }) %>% do.call(rbind,.)
        
      # Each row of ss.ref.pair.scores data frame gets a list of rfs that contributed (to each score, separately).
      # Only included the ss.ref pairs which actually had rfs to add up
     
        rfs.used.tot.rfs <- score.list %>% lapply(function(x) {
          x$bfs.used.tot
        })
        
        rfs.used.res.rfs <- score.list %>% lapply(function(x) {
          x$bfs.used.res
        })
        
        rfs.used.rmseb.rfs <- score.list %>% lapply(function(x) {
          x$bfs.used.rmseb
        })
        
        rfs.used.min.rfs <- score.list %>% lapply(function(x) {
          x$bfs.used.min.score
        })
        
        
        
      # Get the coordinates for each score in the samples x refs matrix
        rfs.used.score.coord <- score.list %>% lapply(function(x) 
          
            data.frame(ss.spec = x$pair.scores$ss.spec,
                       ref = x$pair.scores$ref)
            
        ) %>% do.call(rbind,.)

           
        # Put in a single list
        
          rfs.used <- list(tot = rfs.used.tot.rfs,
                           res = rfs.used.res.rfs,
                           rmseb = rfs.used.rmseb.rfs,
                           min = rfs.used.min.rfs,
                           score.mat.coords = rfs.used.score.coord)
        
      message('\nwriting scores to file...')
      saveRDS(ss.ref.pair.scores, paste0(this.run, "/ss.ref.pair.scores.RDS"))
      saveRDS(rfs.used, paste0(this.run, "/rfs.used.RDS"))
      
      message('Pair score summation complete.')
      # return(ss.ref.pair.scores, bfs.used)
}


#########################################################################################################
#' Update ref-sample scoring vectors (accessory function for pair_score_summation())
#' @param rp ref.pair row 
#' @param score.obj list containing
#'                    - score.name (must be a column name in rp), 
#'                    - score.vect (vector of length(refspec) indexed by ref.pts) containing current scores for each ref point
#'                    - rf.vect (vector of length(refspec) indexed by ref.pts) containing corresponding rf.id for rf which set the current score at each point
#' @param ref.pts indices of score.vect and rf.vect corresponding to this rf; those we are currently updating
#'
#' @return score object with relevant points updated in both vectors  
#' Example: bff.res <- update_scoreObj(ref.pairs[j, ], bff.res, ref.pts)
#'
#' @export
update_scoreObj <- function(rp, score.obj, ref.pts){
  # Get score from refpair row rp
    rf.score <- as.numeric(rp[score.obj$score.name])
  
  # Compare to scores vect (already subsetted for ref.pts)
    replace <- (rf.score > score.obj$scores[ref.pts])
  
  # For points which had a score lower than this rf's score, replace:
    score.obj$scores[ref.pts[replace]] <- rf.score  
    
  # Also record which rf gave rise to those updated scores:
    score.obj$rf.ids[ref.pts[replace]] <- rp$match # record which backfit (match #) it came from

  return(score.obj)
}

#########################################################################################################
#' Sum the score in the scoresObj
#' @param rp ref.pair row 
#' @param score.obj list containing
#'                    - score.name (must be a column name in rp), 
#'                    - score.vect (vector of length(refspec) indexed by ref.pts) containing current scores for each ref point
#'                    - rf.vect (vector of length(refspec) indexed by ref.pts) containing corresponding rf.id for rf which set the current score at each point
#' @param ref.pts indices of score.vect and rf.vect corresponding to this rf; those we are currently updating
#'
#' @return score object with relevant points updated in both vectors  
#' Example: bff.res <- update_scoreObj(ref.pair[j, ], bff.res, ref.pts)
#'
#' @export
sum_score <- function(score.obj, refspec){

  use <- score.obj$scores > 0
  score.obj$scores.tot <- sum(score.obj$scores[use] * refspec[use], na.rm = T)
  score.obj$rf.ids.tot <- unique(score.obj$rf.ids[use])
  
  return(score.obj)
}
   
