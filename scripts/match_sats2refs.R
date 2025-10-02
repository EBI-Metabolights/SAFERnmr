  
s <- sat.list[[i]]
plot_protofeature(p = data.frame(driver = s$peak),
          half.window = half.window, ppm = data$ppm,
          xmat = xmat[s$subset,],
          # bgplot = 'stack', line.shape = 'covar', line.color = "corr",
          bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
          showPeaks = FALSE, ref.mask = s$ref.idx, show.mask.bounds = TRUE)

pexp <- expand_protofeature(p = data.frame(driver = s$peak), 
                               xmat[s$subset,], ppm, half.window)

simplePlot(ymat = pexp$cv, xvect = ppm[pexp$specRegion.inds])

featureStack <- pexp$cv %>% length %>% matrix(NA, nrow = 1, ncol = .)
mask <- (pexp$specRegion.inds %in% s$ref.idx)
featureStack[, mask] <- pexp$cv[mask]
# plot(featureStack%>%t)

lib.data <- readRDS(pars$files$lib.data)

ppm.range <- range(pexp$ppmRegion)
ppm.tol <- 1
ref.range <- c(ppm.range[1]-ppm.tol, ppm.range[2]+ppm.tol)

lib.data.processed <- prepRefs_for_dataset(lib.data,
                                           ppm.dataset = ppm,
                                           ref.sig.SD.cutoff = pars$matching$ref.sig.SD.cutoff,
                                           ppm.range = ref.range,
                                           n.cores = pars$par$ncores
)


match_features <- function(
    f.stack
    ){
  
}

# Scale the feature matrix rows ####
  f.stack <- featureStack %>% apply(1, scale_between) # this will also transpose it, so no need to do later

# Put ref spectra in a matrix ####
  ref.mat <- lapply(lib.data.processed, function(x)
    {
      x$mapped$data.compressed %>% expand_stacklist(which.stacks = 'data') %>% .[[1]]
    }
  ) %>% do.call(rbind, .)

  ref.mat[1,] %>% simplePlot
  
  
  # Pad the ref spectra to feature size
  pad.size <- nrow(f.stack) - 1
  r.mat <- ref.mat %>% padmat(use = 0, col.by = pad.size)
  r.mat[is.na(r.mat)] <- 0
  
  # Transpose original matrix so columns are spectra, save and remove it to clear memory ####
    ref.mat <- ref.mat %>% compress_stack(sparse.val = NA)
  
  # List-format the matrices to facilitate parallel
    message('\tSplitting ref matrix to lists...')
    r.mat <- lapply(1:nrow(r.mat), function(r) r.mat[r,])
  
  # Loop through spec matrix, compute fftw::fft()
    message('\tReference matrix fft...')
    r.mat <- mclapply(r.mat, function(ref) fftw::FFT(ref), mc.cores = pars$par$ncores) %>% do.call(cbind,.)
    
    
    

# At this point, we have:
# - features:         f.stack (vertical)
# - references:       ref.mat (vertical; compressed)
# - padded refs:      r.mat
# - fft'd references: r.mat   (vertical)
  
    f.ind = 1
    feat = f.stack[,f.ind, drop = F]
    
    padded.feat <- feat %>% c(rep(0, nrow(r.mat) - length(feat)))
    padded.feat[is.na(padded.feat)] <- 0

    feat.ft <- Conj(fftw::FFT(padded.feat))

    message("    - cross-correlating to refs...")
    
    ref.mat <- ref.mat %>% cstack_expandRows %>% t
    
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
        # 
        # simplePlot(feat %>% t %>% trim_sides(out = "inds") %>% feat[.])
        # simplePlot(ref %>% t %>% trim_sides(out = "inds") %>% ref[.])
        
        # Cross-correlate to find locations and scores:
          matches <- feature_match2ref_slim(f.num, r.num,
                                            feat, ref,
                                            pad.size = pad.size,
                                            feat.ft, ref.ft,
                                            max.hits = pars$matching$max.hits,
                                            r.thresh = .5,#pars$matching$r.thresh,
                                            p.thresh = .01)#pars$matching$p.thresh)
          
                 feat.ft.c <- feat.ft

          
          return(matches)
      }
    

    # Escape and return nothing if there are no matches to evaluate: ####
        # if (is.null(allmatches.feat)){return(NA)}
    
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
            ref <- ref.mat[,r,drop = F]
            
          # Fit
            
            fit <- fit_leastSquares(feat[feat.pos], ref[ref.pos], plots = F, scale.v2 = T)
              # fit$plot
            # fit$wasserstein.score <- score_wasserstein(fit$feat.fit, fit$spec.fit)
            
            return(fit)
        })
     
        # i <- 0
        # i <- i+1
        # plot_fit(allmatches.fits[[i]],type = "auc")
        
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
        matching.time <- Sys.time() - t1
        print(matching.time)

        parallel::stopCluster(my.cluster)
          
          