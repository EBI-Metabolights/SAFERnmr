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

## Parallelize the ss.spec - reference pair score summation? ####      
      message('Reading in data for pair score summation...')

        tmpdir <- pars$dirs$temp
        this.run <- paste0(tmpdir)

      # split backfits and combinations by reference spectra used
      # report as pairs with coords in matrix
      
        ss.ref.pairs <- readRDS(paste0(this.run, "/ss.ref.pairs.RDS"))
        refmat <- readRDS(paste0(this.run, "/temp_data_matching/ref.mat.RDS")) %>% t
      
      # # Normalize refs 
        refmat <- refmat / Rfast::colsums(refmat, na.rm = T)

      message('Splitting data for parallelization...')
      # Pre-split the table and refs (reduce overhead) ####
        by.ref <- pbapply::pblapply(unique(ss.ref.pairs$ref), function(r) 
        {
          ref.pairs <- ss.ref.pairs[which(ss.ref.pairs$ref == r),]
          # We don't need the data for the fits, just the descriptive stats
          
          if (nrow(ref.pairs) > 0){
            
            return(list(ref.num = r,
                        refspec = refmat[r, ],
                        ref.pairs = ref.pairs))
          }  
          return(NULL)
        })
        
      rm(refmat)
      
      # Good idea to do some load-balancing here...
        workloads <- lapply(chunks, function(x) nrow(ref.pairs)) %>% unlist
        plan <- distribute(workloads, across = pars$par$ncores, messages = F)
        
        chonks <- lapply(1:ncores, function(core.number){
          chonk <- by.ref[plan$core.ids == core.number]
          return(chonk)
        })
        
        rm(by.ref)
      
      
      message('Computing scores...')

          
    # Calculate the scores in parallel ####
    
        t1 <- Sys.time()
        
            
            score.list <- mclapply(chonks, mc.cores = pars$par$ncores, 
                                           FUN = function(r.list)                             
            {
              # chonk <- chonks[[1]]                               
              # Go through this list of by.refs
              chonk.scores <- lapply(chonk, function(r.list){
                  
                # For each ref:
                
                # r.list <- chonk[[1]]
                
                # Make sure there are actually spectra in the interactions list ####
                  if(nrow(r.list$ref.pairs) < 1){return(NULL)}
                 
                # Extract the interactions info from the r.list object for this ref ####
                  r.num <- r.list$ref.num
                  refspec <- r.list$refspec #/sum(r.list$refspec, na.rm = T) # already normed
                  ref.pairs <- r.list$ref.pairs
                
                # Compute ss interactions (combinations) for this ref ####
                  comb <- expand.grid(unique(ref.pairs$ref), unique(ref.pairs$ss.spec))
                    colnames(comb) <- c("ref", "ss.spec")
                    # v.zeros <- rep(0, length(refspec))
                
                # Loop through ref-ss.spec combinations and calculate scores ####
                  lapply(1:nrow(comb), function(i){
                      # print(i)
                    # Select the pair ####
                      ref <- comb$ref[i]
                      ss.spec <- comb$ss.spec[i]
                      
                    # Select relevant backfits and matches ####
                      rp.rows <- which(ref.pairs$ss.spec == ss.spec)
  
                    # Make a tmp score and best.rf.index.used vector of length(refspec) for each score ####
                      
                      cum.bff.tot <- cum.bff.res <- cbt.best <- cbr.best <- rep(0, length(refspec))
  
                    # Loop through the matches associated with this ref - ss.spec pair ####
                    # Update the bff values in v with any higher bff at that point ####
                      for (j in rp.rows){
                        # j indexes the rows of rp.rows (ss.ref.pairs belonging to this chunk = this ref #)
                        # that match this dataset spectrum (ss.spec). These could include many matches. A 
                        # match could also give rise to many rp.rows, because of the different spectra. rp.rows
                        # could be derived from multiple matches as long as they involve the same ref and 
                        # have backfits to ss.spec. 
                        # 
                        # i <- i + 1
                        # j <- rp.rows[i]
                        # Get the ref range for the matched ref-feat
                          # Where is the ref range? It's in the backfit, now.
            
                          ref.range <- ref.pairs$ref.start[j] : ref.pairs$ref.end[j]
                          
                          bff <- ref.pairs[j, ]$bff.tot
                          # print(bff)
                          replace <- (bff > cum.bff.tot[ref.range])
                          cum.bff.tot[ref.range[replace]] <- bff
                          cbt.best[ref.range[replace]] <- ref.pairs$match[j] # record which backfit (match #) it came from 
  
                          bff <- ref.pairs[j, ]$bff.res
                          replace <- (bff > cum.bff.res[ref.range])
                          cum.bff.res[ref.range[replace]] <- bff
                          cbr.best[ref.range[replace]] <- ref.pairs$match[j] # record which backfit (match #) it came from
  
                      }
                      
                    # Calculate score and report as data.frame of ss.spec-reference pairs ####
                    use <- cum.bff.tot > 0
                    score.tot <- sum(cum.bff.tot[use] * refspec[use], na.rm = T)
                    score.res <- sum(cum.bff.res[use] * refspec[use], na.rm = T)
                    
                    list(
                           bfs.used.tot = unique(cbt.best[cbt.best > 0]), # keep these around; they're the best evidence
                           bfs.used.res = unique(cbr.best[cbr.best > 0]), # keep these around; they're the best evidence
                           pair.scores = data.frame(ss.spec = ss.spec,
                                                    ref = r.num,
                                                    score.tot = score.tot, 
                                                    score.res = score.res)
                    )
                    
                })
              }) %>% unlist(recursive = F)   
            }) %>% unlist(recursive = F, use.names = F)
      
        print(Sys.time() - t1)
     
      # Extract out the scores data.frame rows
        message('\nextracting out the scores data.frame rows...')
        ss.ref.pair.scores <- score.list %>% lapply(function(x) {
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
        
      # Coordinate in scores matrix between the two scores; lists of rfs differ
        rfs.used.score.coord <- score.list %>% lapply(function(x) 
          
            data.frame(ss.spec = x$pair.scores$ss.spec,
                       ref = x$pair.scores$ref)
            
          ) %>% do.call(rbind,.)

           
        # Put in a single list
        
          rfs.used <- list(tot = rfs.used.tot.rfs,
                           res = rfs.used.res.rfs,
                           score.mat.coords = rfs.used.score.coord)
        
      message('\nwriting scores to file...')
      saveRDS(ss.ref.pair.scores, paste0(this.run, "/ss.ref.pair.scores.RDS"))
      saveRDS(rfs.used, paste0(this.run, "/rfs.used.RDS"))
      
      message('Pair score summation complete.')
      # return(ss.ref.pair.scores, bfs.used)
}
      