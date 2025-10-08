## Matching Features using Functions
  Sys.setenv(MallocStackLogging = "0")
  lib.data.processed <- load_lib_data(pars)

  ref.stack <- lib_to_refmat(lib.data.processed) # refs on rows
    rm(lib.data.processed)
  # ref.stack <- xmat # spectra on rows
  
  # Bind feature and ref ranges
  
  feature.stack <- stack_sats(sat.list, xmat, ppm, half.window) # feats on rows

# Prep the features and refs ####
  
  # +/- 1 ppm from each 
  roi <- c(2.25,2.5)
    tol <- 0.1
    roi[1] <- roi[1]-tol
    roi[2] <- roi[2]+tol
  match.pack <- coprep_features_and_refs(feature.stack, ref.stack, ppm, roi, downsampling.factor=4)
  
# Do the matching ####

  matches <- match_features(match.pack)
  
# Functions ####

  downsample_inds <- function(v, dec.factor=4){
    inds <- seq(1,length(v), by=dec.factor)
    return(inds)
  }

  load_lib_data <- function(pars){
  # lib.data <- readRDS(pars$files$lib.data)
      lib.data.700 <- readRDS('/Users/mjudge/Documents/ftp_ebi/gissmo/data.list_700MHz.RDS')
      gissmo.cmpds <- readxl::read_xlsx('/Users/mjudge/Documents/ftp_ebi/gissmo/gissmo2chebi_2024.xlsx')
      
      source('/Users/mjudge/Documents/GitHub/MARIANA_setup_chron/R/add_chebiIDs.R') # on "no-zip" branch
      lib.data.700 <- add_chebiIDs(lib.data = lib.data.700, key = gissmo.cmpds)

  lib.data <- lib.data.700
  
  # ppm.range <- range(pexp$ppmRegion)
  # ppm.range <- c(.5,9.5)
  # ppm.tol <- 1
  # ref.range <- c(ppm.range[1]-ppm.tol, ppm.range[2]+ppm.tol)
  ref.range <- range(ppm)
  
  lib.data.processed <- prepRefs_for_dataset(lib.data,
                                             ppm.dataset = ppm,
                                             ref.sig.SD.cutoff = 0,#pars$matching$ref.sig.SD.cutoff,
                                             ppm.range = ref.range,
                                             n.cores = pars$par$ncores
  )
  return(lib.data.processed)
}

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
    }, mc.cores = pars$par$ncores) %>% do.call(rbind, .)
  }

  match_features <- function(mp){
    # mp <- match.pack
    # pars$matching$max.hits <- 10
    # pars$matching$r.thresh <- .7
    # pars$matching$p.thresh <- .01
    
    # Par setup
    
    my.cluster <- safer_makeCluster(par, nfeats=length(mp$f.numbers))
        
    # Do matching (all refs, per feature):
    
    allmatches.feats <-  foreach(f.num = mp$f.numbers,
                                feat = mp$features,
                                feat.padded.ft.c = mp$features.padded.ft.c,
                                .combine = 'rbind',
                                .errorhandling="pass") %dopar%
    
    {
      # i <- 1
      # f.num<-mp$f.numbers[i]
      # feat = mp$features[,i]
      # feat.padded.ft.c = mp$features.padded.ft.c[,i]
      
      refs = mp$refs
      refs.padded.ft = mp$refs.padded.ft
        
      allmatches.feat <- match_feature(f.num, feat, feat.padded.ft.c,
                                       mp$refs, mp$refs.padded.ft)
      
      if (is.null(nrow(allmatches.feat))){
        allmatches.feat <- NULL
        specificity.score <- Inf
      } else {
        allmatches.feat <- fit_matches(allmatches.feat, feat, mp$refs)
      
  
        # Calculate feature specificity score ####
            n.passing <- sum(allmatches.feat$rval >= pars$matching$r.thresh)
            n.refs <- ncol(mp$refs)
            specificity.score <- n.passing/n.refs
            
            # allmatches.fits$rval %>% sort %>% plot
            
            # message('Specificity Score: ', round(specificity.score, 2))
        
          # Best score (1) should be where feature only binds strongly to true location
          # Worst score (0) is where feature binds everywhere. 
      }
      return(list(matches = allmatches.feat,
                  specificity = specificity.score))
    
    }
    
    parallel::stopCluster(my.cluster)
    message('...parallel cluster closed.')
    return(allmatches.feats)
  }
  
  safer_makeCluster <- function(par, nfeats){
        message("Setting up parallel cluster...\n\n")
        # we want at least one feature per core, but at least one core.
        ncores <- min(pars$par$ncores, nfeats) %>% c(.,1) %>% max
        my.cluster <- parallel::makeCluster(ncores, type = pars$par$type)
        doParallel::registerDoParallel(cl = my.cluster)
        
        if(foreach::getDoParRegistered()){
          message('\tparallel cluster started on ', foreach::getDoParWorkers(),' cores...\n\n')
        } else {stop('Matching: parallel pool could not be started.')}
        return(my.cluster)
  }
  
  prep_features <- function(features, ref.length){

    message('\tPadding features by ref.length (', ref.length, ') - length(feat) (', ncol(features), ')...')
    message('\tPadding features by ref.length - length(feat)...')
    pad.size <- ref.length - ncol(features)
    fsp <-
    mclapply(1:nrow(features), function(f){
      feat <- features[f, ]
      padded.feat <- feat %>% c(rep(0, pad.size),.)
      padded.feat[is.na(padded.feat)] <- 0
      feat.p.ft.c <- Conj(fftw::FFT(padded.feat))
      return(feat.p.ft.c)
    }, mc.cores = pars$par$ncores) %>% do.call(rbind,.) %>% t
    return(fsp)
  }
  
  coprep_features_and_refs <- function(feature.stack, ref.stack, ppm, roi, downsampling.factor=1){
    
    # Assume feature and ref stacks have the same ppm axis (on the columns)
    # roi is in ppm 
    
      # Cut down feature stack to relevant region

      roi <- roi %>% vectInds(., ppm) # ppm
      
      in.range <- lapply(sat.list, function(s){
        # s <- sat.list[[1]]
        !all(is.na(range_intersect(roi, s$finalRegion)))
      
      }) %>% unlist %>% which
    
      # Thin it out some
      selected <- seq(1, length(in.range), length.out=8) %>% in.range[.]
      feature.stack <- feature.stack[selected,]
      
      # Scale the features
      feature.stack <- lapply(1:nrow(feature.stack), function(x){
        feature.stack[x,] %>% scale_between %>% c
      }) %>% do.call(rbind, .)
      # stackplot(feature.stack, vshift = 10)
      
    # Downsample (if doing that)
      message('\t downsampling refs and features by a factor of ',downsampling.factor)
      
      ds.inds.ref <- downsample_inds(ppm, downsampling.factor)
      ref.stack <- ref.stack[,ds.inds.ref]
      ds.ppm <- ppm[ds.inds.ref]
      # stackplot(ref.stack[1:10,], vshift = 10, xvect = ds.ppm)
      
      ds.inds.feat <- seq(1,ncol(feature.stack)) %>% downsample_inds(4)
      feature.stack <- feature.stack[,ds.inds.feat]

      # Adjust ROI for downsampling
        
        roi <- (roi/downsampling.factor) %>% round

    
    # Move on to processing ref.stack
      feature.width <- ncol(feature.stack)
      
      refs.padded.ft <- prep_refs(ref.stack, feature.width, roi)
      refs <- ref.stack %>% t
      
    # Move on to processing feature.stack
      padded.ref.length <- nrow(refs.padded.ft)
      
      features.padded.ft.c <- prep_features(feature.stack, padded.ref.length)
      features <- feature.stack %>% t
      
    return(list(feature.width = feature.width,
                f.numbers=selected,
                features=features,
                refs=refs,
                features.padded.ft.c=features.padded.ft.c,
                refs.padded.ft=refs.padded.ft,
                feature_downsampled_inds=ds.inds.feat,
                ref_downsampled_inds=ds.inds.ref))
    
  }
  
  prep_refs <- function(refs, feature.width, roi){
    # Trim the ref stack to relevant region only:
    ref.stack <- ref.stack[,fillbetween(roi)]
    
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
  
  match_feature <- function(f.num, feat, feat.padded.ft.c,
                             refs, refs.padded.ft){
    
    # Locate best positions in all available refs
    
      allmatches.feat <-  foreach(r.num = 1:ncol(refs.padded.ft),
                                  ref = refs,
                                  ref.ft = refs.padded.ft,
                                  .combine = 'rbind',
                                  .errorhandling="pass") %do%
      {
        r.num = 1
        ref = refs[,r.num, drop = F]
        ref.ft = refs.padded.ft[,r.num, drop = F]
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


