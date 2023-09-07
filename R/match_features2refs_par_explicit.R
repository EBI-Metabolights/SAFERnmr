#' match.features2refs.par.explicit
#'
#' This function matches features to reference spectra using a parallelized loop. 
#' Feature chunking strategy. 
#' Do you match features to refs, or refs to feature? If using batman-style fit, 
#' since refs are more pure, perhaps the latter makes more sense. In other words,
#' fit the cleaner signals to the dirtier ones if using batman. For less biased,
#' use least squares.
#'
#' @param pars A list of parameters containing directories, parallel settings, and matching settings.
#' @return matches.RDS file containing match.info, fits, and peak quality data.
#' @export
#'
#' @import parallel
#' @import foreach
#' @import doParallel
#' @importFrom magrittr %>%
#' 
#' @author MTJ
match_features2refs_par_explicit <- function(pars){
  
  message('\n-----------------------------------------------------------------')
  message('-------------------  Parallel Matching to PRCSs -------------------')
  message('-------------------------------------------------------------------')
  printTime()
  # mem_snapshot(paste0(tmpdir, '/', Sys.time(), '.txt'))

################ Read parameters file ##################
  message('Reading data...\n')
  tmpdir <- pars$dirs$temp

  pad.size <- readRDS(paste0(tmpdir, "/temp_data_matching/pad.size.RDS"))
  
  f.stack.split <- readRDS(paste0(tmpdir, "/temp_data_matching/f.stack.split.RDS"))
    f.stack.split %>% test_nullish('f.stack.split')
    message('Matching ', lapply(f.stack.split, ncol) %>% unlist %>% sum, ' features...\n')
    
  split.scheme <- readRDS(paste0(tmpdir, "/temp_data_matching/split.scheme.RDS"))
    split.scheme %>% test_nullish('split.scheme')
    
  ref.mat <- readRDS(paste0(tmpdir, "/temp_data_matching/ref.mat.RDS")) %>% cstack_expandRows %>% t
    # At some point, we should probably compress these and compute the FTs on the fly...
    # or, better, condense the ref mats and their fts by thresholding and chunking into padded 
    # regions
    # apply(ref.mat, 2, function())
    ref.mat %>% test_nullish('ref.mat')
    
  r.mat <- readRDS(paste0(tmpdir, "/temp_data_matching/rmat.RDS"))
    r.mat %>% test_nullish('r.mat')
    
  # mem_snapshot(paste0(tmpdir, '/', Sys.time(), '.txt'))
    # Par setup ####
      message("Setting up parallel cluster...\n\n")
      ncores <- pars$par$ncores
      my.cluster <- parallel::makeCluster(ncores, type = pars$par$type)
      doParallel::registerDoParallel(cl = my.cluster)
      
      if(foreach::getDoParRegistered()){
        message('\tparallel cluster started on ', foreach::getDoParWorkers(),' cores...\n\n')
      } else {stop('Matching: parallel pool could not be started.')}
      
      
      
# Run the matching loop ####
#       f.subset = split.scheme$f.num,
        # f.stack = f.stack.split,
        # f.mat = f.mat.split,
      t1 <- Sys.time()
      matches.all <- foreach(chunk = 1:length(split.scheme), 
                             split.grp = split.scheme,
                             f.stack = f.stack.split,
                             .combine='c', .multicombine=TRUE,
                             .packages = c('fftw', 'foreach'),
                             .errorhandling="pass") %do%

      {
        
        message('Chunk ', chunk)
        
        # For each feature batch:
            
             # Set up ####
              
              # f.stack <- readRDS(paste0(tmpdir, "/temp_data_matching/f.stack.",chunk,".RDS"))
              # f.mat <- readRDS(paste0(tmpdir, "/temp_data_matching/f.mat.",chunk,".RDS"))
              
              # chunk <- 3
              # split.grp <- split.scheme[[chunk]]
              
              f.subset <- split.grp$f.subset
              # f.stack = f.stack.split[[chunk]]
              # f.mat = f.mat.split[[chunk]]
              
              f.stack <- f.stack %>% matrix(ncol = length(f.subset))
              # f.mat <- f.mat %>% matrix(ncol = length(f.subset))
              # simplePlot(trim_sides(t(f.stack)))
        
              # Loop through all features in the current batch ####
                matches.chunk <-  foreach(f.num = f.subset,           # what actually gets recorded as feat number
                                          f.ind = 1:length(f.subset), # feat number within the chunk
                                          feat = f.stack,             # actual feature profile
                                          .combine='c', .multicombine=TRUE,
                                          .errorhandling="pass") %do%
                {
                  tryCatch(
                        expr = {
                        # mem_snapshot(paste0(tmpdir, '/', Sys.time(), '_feature_',f.num,'.txt'))
                        # Make the fft and conj'd feature on the fly, since it's only ever used once:
                        #  in the next foreach statement. 
                        
                          feat <- feat %>% scale_between # just for good measure
                          padded.feat <- feat %>% c(rep(0, nrow(r.mat) - length(feat)))
                          padded.feat[is.na(padded.feat)] <- 0
      
                          feat.ft <- Conj(fftw::FFT(padded.feat))
                          
                          # f.ind = 2
                          # f.num = f.subset[f.ind]
                          # feat = f.stack[,f.ind, drop = F]
                          # feat.ft = f.mat[,f.ind, drop = F]
                          # simplePlot(trim_sides(feat))
          
                      message("Matching feature: ", f.ind,"...") # not same as f.num
                      message("    - cross-correlating to refs...")
                      
                      # Locate best positions in all available refs  ####
                        allmatches.feat <-  foreach(r.num = 1:ncol(r.mat),
                                                    ref = ref.mat,
                                                    ref.ft = r.mat,
                                                    .combine = 'rbind',
                                                    .errorhandling="pass") %do%
                        {
                          # r.num = 1
                          # ref = ref.mat[,r.num, drop = F]
                          # ref.ft = r.mat[,r.num, drop = F]
                          
                          # simplePlot(feat %>% t %>% trim_sides(out = "inds") %>% feat[.])
                          # simplePlot(ref %>% t %>% trim_sides(out = "inds") %>% ref[.])
                          
                          # Cross-correlate to find locations and scores:
                            matches <- feature_match2ref_slim(f.num, r.num,
                                                              feat, ref,
                                                              pad.size = pad.size,
                                                              feat.ft, ref.ft,
                                                              max.hits = pars$matching$max.hits,
                                                              r.thresh = pars$matching$r.thresh,
                                                              p.thresh = pars$matching$p.thresh)
                            
                            return(matches)
                        }
                      
                        rm(feat.ft)
                      
                      # Escape and return nothing if there are no matches to evaluate: ####
                          if (is.null(allmatches.feat)){return(NA)}
          
                      ######################################################################
          
                      # Evaluate fits for top-scoring positions (regardless of ref) ####
                        message("    - calculating ", nrow(allmatches.feat), " fits...")
                        
                        allmatches.fits <- lapply(1:nrow(allmatches.feat), function(m)
                        {
                            
                          # Get f and r indices for this row
                            f <- allmatches.feat[m, 'feat']
                            r <- allmatches.feat[m, 'ref']
                            feat.pos <- allmatches.feat[m, c('feat.start','feat.end')] %>% as.numeric %>% fillbetween
                            ref.pos <- allmatches.feat[m, c('ref.start','ref.end')] %>% as.numeric %>% fillbetween
          
                          # Get spectral signatures which matched
                            ref <- ref.mat[,r,drop = F] #%>% scale_between
                            # ref.sum <- sum(ref, )
                          # Fit
                            
                            fit <- fit_leastSquares(feat[feat.pos], ref[ref.pos], plots = F, scale.v2 = T)
                              # fit$plot
                            # fit$wasserstein.score <- score_wasserstein(fit$feat.fit, fit$spec.fit)
                            
                            return(fit)
                        })
                        
                        
                        # Extract out minimal fit data ####
                          fit.data <- lapply(allmatches.fits, function(f) f$fit %>% as.numeric) %>% do.call(rbind,.)
                          allmatches.feat[,"fit.intercept"] <- fit.data[,1]
                          allmatches.feat[,"fit.scale"] <- fit.data[,2]
          
                        # Add some different scores from the fits ####
                          message("    - calculating additional scores...")
                          # allmatches.feat[,'wasserstein.score'] <-
                          #   lapply(allmatches.fits, function(x) x$wasserstein.score) %>% unlist
                          allmatches.feat[,'sum.residuals'] <-
                            lapply(allmatches.fits, function(x) x$sum.residuals) %>% unlist
                          allmatches.feat[,'rmse'] <-
                            lapply(allmatches.fits, function(x) x$rmse) %>% unlist
          
                      # Re-weight the rmses using peak quality (for this db) ####
                        message("    - estimating peak quality, weighting RMSEs...\n")
                        peak.quality <- bad_peaks(pos.res.pct.feat = lapply(allmatches.fits,
                                                    function(x) x$overshoot.pct.feat
                                                  ) %>% do.call(rbind,.),
                                                      scale.exp = 2)
                        allmatches.feat[,'rmse.weighted'] <- lapply(1:nrow(allmatches.feat), function(m)
                        {
                            match <- allmatches.fits[[m]]
                            v1 <- match$feat.fit
                            v2 <- match$spec.fit
                            resids <- v1-v2
                            use <- !is.na(resids+peak.quality)
                            rmse <- sum((resids[use] * (1-peak.quality[use])) ^ 2)/sum(use)
                            return(rmse)
                        }) %>% unlist
          
                      # Format results ####
          
                        return(list(matches = allmatches.feat,
                                    peak.quality = peak.quality))
      
                    
                  }, 
                        error = function(cond){
                        return(list(NA))
                  }
                  )
                }
                # The result is a list with the following fields:
                # - matches
                #   - list of data.frames giving the match information
                # - peak.quality
                #   - list of vectors giving the dataset-specific usefulness of each point in feature (down-ranks never-fit points)
                
            # if (is.null(matches.chunk)){message('\tfailed')} else {message('\tsucceeded'); return(matches.chunk)}
          
      }
        print(Sys.time() - t1)

        parallel::stopCluster(my.cluster)
        
        # mem_snapshot(tmpdir, '/', paste0(Sys.time(), '.txt'))
        
############ Format results and save ########################################################################################
    # Compile and save match results ####
        message('Saving match results...\n\n\n')
        saveRDS(matches.all, paste0(tmpdir, "/matches.RDS"))
        unlink(paste0(tmpdir, "/temp_data_matching/rmat.RDS"))
        # matches <- readRDS(paste0(tmpdir, "/matches.RDS"))
        # ####
  message('\n-----------------------------------------------------------------')
  message('-------------------  Parallel Matching Complete -------------------')
  message('-------------------------------------------------------------------')

}
    
    
    
    
        
  

  
  