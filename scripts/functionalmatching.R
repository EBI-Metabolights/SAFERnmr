## Matching Features using Functions

refs <- lib_to_refmat(lib.data.processed)

  refs <- xmat

  f.stack <- stack_sats(sat.list, xmat, ppm, half.window)
  
  feat <- f.stack[f.num,,drop=FALSE]
    simplePlot(feat)
  
# Prep the features and refs ####

  feature.width <- nrow(feat)
  refs.padded.ft <- prep_refs(refs, feature.width)
  refs <- refs %>% t
  
  padded.ref.length <- nrow(refs.padded.ft)
  feat.padded.ft.c <- prep_features(feat, padded.ref.length)
  # match.pack <- prep_features_and_refs(feats, refs)
  
# Do the matching ####
  
  pars$matching$max.hits <- 10
  pars$matching$r.thresh <- .7
  pars$matching$p.thresh <- .01
  
  allmatches.feat <- match_feature(feat, feat.padded.ft.c,
                                   refs, refs.padded.ft)
  
  allmatches.fits <- fit_matches(allmatches.feat, feat, refs)

# Calculate feature specificity score ####
  n.passing <- sum(allmatches.feat$rval >= pars$matching$r.thresh)
  n.refs <- ncol(refs)
  specificity.score <- 1-(n.passing/n.refs)
  
  allmatches.fits$rval %>% sort %>% plot
  
  message('Specificity Score: ', round(specificity.score, 2))

  # Best score (1) should be where feature only binds strongly to true location
  # Worst score (0) is where feature binds everywhere. 

  lib_to_refmat <- function(lib.data.processed){
    refs <- lapply(lib.data.processed, function(x)
      {
        x$mapped$data.compressed %>% expand_stacklist(which.stacks = 'data') %>% .[[1]]
      }
    ) %>% do.call(rbind, .)
  }
  
  stack_sats <- function(sat.list, xmat, ppm, half.window){
    mclapply(sat.list, function(s){
      
      p <- data.frame(driver = s$peak)
      driver <- p$driver
  
      # Driver locates the index, everything else can be built around it
      
      p.abs <- driver - p
      
      fullView <- (driver - half.window):(driver + half.window)
      
      in.bounds <- !(fullView < 1 | fullView > length(ppm))
      
      specreg.inds <- fullView[in.bounds]

      cv <- cr <- rep(NA, length(fullView))
      cv[in.bounds] <- s$covar[in.bounds]
      cr[in.bounds] <- s$corr[in.bounds]
      
      # Make mask
      # If this isn't a protofeature, but just a driver: let "peak" bounds = spec Region
      if (is.null(p.abs$primary.lower)){
        p.abs$primary.lower <- min(specreg.inds)
        p.abs$primary.upper <- max(specreg.inds)
        p.abs$secondary.lower <- p.abs$primary.lower
        p.abs$secondary.upper <- p.abs$primary.upper
      }

      primary.bounds <- c(p.abs$primary.lower, p.abs$primary.upper) %>% sort
      secondary.bounds <- c(p.abs$secondary.lower, p.abs$secondary.upper) %>% sort
      
      
      pk.mask <- rep(0, length(fullView))
      pk.mask[which(fullView %in% primary.bounds) %>% fillbetween] <- 1
      pk.mask[which(fullView %in% secondary.bounds) %>% fillbetween] <- 2
    
      # s <- sat.list[[1]]
      # plot_protofeature(p = data.frame(driver = s$peak),
      #           half.window = half.window, ppm = data$ppm,
      #           xmat = xmat[s$subset,],
      #           # bgplot = 'stack', line.shape = 'covar', line.color = "corr",
      #           bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
      #           showPeaks = FALSE, ref.mask = s$ref.idx, show.mask.bounds = TRUE)
      
      # pexp <- expand_protofeature(p = data.frame(driver = s$peak), 
      #                                xmat[s$subset,], ppm, half.window)
      # pexp$specRegion <- NULL
      # simplePlot(ymat = cv, xvect = ppm[specreg.inds])
      
      feat <- cv %>% length %>% matrix(NA, nrow = 1, ncol = .)
      mask <- (specreg.inds %in% s$ref.idx)
      feat[, mask] <- cv[mask]
      return(feat)
    }, mc.cores = 8) %>% do.call(rbind, .)
  }


  match_features <- function(){
    
  }
  
  prep_features <- function(feat, ref.length){
    # feat = f.stack[,f.num, drop = F]
    # plot(feat)
  
    padded.feat <- feat %>% c(rep(0, ref.length - length(feat)))
    padded.feat[is.na(padded.feat)] <- 0
  
    feat.ft <- Conj(fftw::FFT(padded.feat))
    return(feat.ft)
  }
  
  
  prep_features_and_refs <- function(feats, refs){
    # This provides a packaged var that plugs into the matching function.
    feature.width <- nrow(feats)
    refs.padded.ft <- prep_refs(refs, feature.width)
    
    padded.ref.length <- nrow(refs.padded.ft)
    feats.padded.ft.c <- prep_features(feats, padded.ref.length)
    return(list(feature.width = feature.width,
                feats=feats,
                refs=refs,
                feats.padded.ft.c=feats.padded.ft.c,
                refs.padded.ft=refs.padded.ft))
  }
  
  prep_refs <- function(refs, feature.width){
    # Pad the ref spectra to feature size
    message('\tPadding refs by feature.width - 1...')
    pad.size <- feature.width- 1
    r.mat <- refs %>% padmat(use = 0, col.by = pad.size)
    r.mat[is.na(r.mat)] <- 0
    
    # List-format the matrices to facilitate parallel
      message('\tSplitting ref matrix to lists...')
      r.mat <- lapply(1:nrow(r.mat), function(r) r.mat[r,])
    
    # Loop through spec matrix, compute fftw::fft()
      message('\tReference matrix fft...')
      r.mat <- mclapply(r.mat, function(ref) fftw::FFT(ref), mc.cores = pars$par$ncores) %>% do.call(cbind,.)
      
    return(r.mat)
  }
  
  match_feature <- function(feat, feat.padded.ft.c,
                             refs, refs.padded.ft){
    
    # Locate best positions in all available refs
    
      allmatches.feat <-  foreach(r.num = 1:ncol(refs.padded.ft),
                                  ref = refs,
                                  ref.ft = refs.padded.ft,
                                  .combine = 'rbind',
                                  .errorhandling="pass") %do%
      {
        # r.num = 1
        # ref = refs[,r.num, drop = F]
        # ref.ft = r.mat[,r.num, drop = F]
        # 
        # simplePlot(feat %>% t %>% trim_sides(out = "inds") %>% feat[.])
        # simplePlot(ref %>% t %>% trim_sides(out = "inds") %>% ref[.])
        
        # Cross-correlate to find locations and scores:
          matches <- feature_match2ref_slim(f.num, r.num,
                                            feat, ref,
                                            pad.size = length(feat)-1,
                                            feat.padded.ft.c, ref.ft,
                                            max.hits = 5,#pars$matching$max.hits,
                                            r.thresh = .7,#pars$matching$r.thresh,
                                            p.thresh = .01)#pars$matching$p.thresh)
          
          return(matches)
      }
  
      return(allmatches.feat)
  }
  
  fit_matches <- function(allmatches.feat, feat, ref.mat){
    fits <- lapply(1:nrow(allmatches.feat), function(m)
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
  
        return(fit)
    })
    
    # Extract out minimal fit data ####
      fit.data <- lapply(fits, function(f) f$fit %>% as.numeric) %>% do.call(rbind,.)
      allmatches.feat[,"fit.intercept"] <- fit.data[,1]
      allmatches.feat[,"fit.scale"] <- fit.data[,2]
    
    # Add some different scores from the fits ####
      message("    - calculating additional scores...")
      allmatches.feat[,'sum.residuals'] <-
        lapply(fits, function(x) x$sum.residuals) %>% unlist
      allmatches.feat[,'rmse'] <-
        lapply(fits, function(x) x$rmse) %>% unlist
  
    return(allmatches.feat)
  }


