## Matching Features using Functions
  # Sys.setenv(MallocStackLogging = "0")
  lib.data.processed <- load_lib_data(pars)
  xmat <- data$xmat
      ppm <- data$ppm
      tmpdir <- pars$dirs$temp
      
  xmat.lin <- xmat %>% t %>% c
  
  ref.stack.lib <- lib_to_refmat(lib.data.processed) # refs on rows
    rm(lib.data.processed)
    # ref.stack<- ref.stack.lib #
  ref.stack <- xmat # spectra on rows
  # ref.stack <- xmat.lin
  # ppm.lin <- 
  
  # Bind feature and ref ranges
  
  feature.stack <- stack_sats(sat.list, xmat, ppm, half.window) # feats on rows

# Prep the features and refs ####
  
  # +/- 1 ppm from each 
  roi <- c(2.25,2.5)
    tol <- 0.1
    roi[1] <- roi[1]-tol
    roi[2] <- roi[2]+tol
    pars$par$ncores <- 4
  match.pack <- coprep_features_and_refs(feature.stack, ref.stack, ppm, roi, downsampling.factor=8)
  
# Do the matching ####
  
  matches <- match_features(match.pack, fitting = TRUE)
  
  matches <- lapply(1:nrow(matches), function(m){
      # Calculate feature specificity score ####
          matches[m, ]$matches
  }) 
  
  all$rval[1]
  all <- matches %>% do.call(rbind,.)
  any(is.na(all))
  
    # all_clean <- all %>% na.omit()
    
    df_out <- a %>%
      na.omit %>%
      group_by(feat,ref) %>%
      mutate(rval_norm = rval / max(rval)) %>%
      summarise(
        specificity = sum(rval_norm),  # or sum(rval_norm) / n()
        # n_refs_with_hits = n(),         # diagnostic
        .groups = "drop"
      )
    %>% 
      group_by(feat) %>% 
      summarise(
        specificity = mean(specificity)
      )
    
    a <- all %>% na.omit %>% filter(feat==2424)
    
# df_specificity <- all %>%
    #   na.omit() %>%
    # 
    #   # Stage 1: normalize matches within each (feat,ref)
    #   group_by(feat, ref) %>%
    #   mutate(rval_norm = rval / max(rval)) %>%
    # 
    #   # Stage 2: compute ambiguity for each (feat,ref)
    #   summarise(
    #     # number of matches within this reference for this feature
    #     n_matches = n(),
    #     # sorted normalized rvals
    #     top = max(rval_norm),
    #     second = ifelse(n_matches >= 2,
    #                     sort(rval_norm, decreasing = TRUE)[2],
    #                     0),
    #     ambiguity_r = 1 - second,
    #     .groups = "drop"
    #   ) %>%
    # 
    #   # Stage 3: compute feature-level specificity across refs where it exists
    #   group_by(feat) %>%
    #   summarise(
    #     specificity = mean(ambiguity_r, na.rm = TRUE),
    #     n_refs_with_hits = n(),   # optional diagnostics
    #     .groups = "drop"
    #   )

  f.numbers <- match.pack$f.numbers
  
  # keep in mind these are scores applied to 
  plot.sats.grid(sat.list, xmat, ppm, selected = df_specificity$feat, title.strs = df_specificity$specificity)
  
  which.plots <- seq(1,nrow(matches),by=100)
  all.plots <- lapply(, function(x){
    i <- i + 1
    x <- which.plots[i]
    match <- matches[x,]
    
    # Plot them
    
      match <- matches[x,]
      f <- which(match.pack$f.numbers==matches$feat[x])
        
      feat <- match.pack$features[,f]
      
      ref <- match.pack$refs[, match$ref, drop = F]
      fit <- match[c("fit.intercept","fit.scale")]
      plot_match(match, feat, ref, ppm.margin = .25)
      
          # g <- plot_match(match, feat, ref, ppm.margin = .25)
      g <- 
      return(g)
  })
  
  
  
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

  match_features <- function(mp, fit.matches = FALSE){
    # mp <- match.pack
    
    # Par setup
    my.cluster <- safer_makeCluster(par, nfeats=length(mp$f.numbers))
        
    # Do matching (all refs, per feature):
    
    allmatches.feats <-  foreach(f.num = mp$f.numbers,
                                feat = mp$features,
                                feat.padded.ft.c = mp$features.padded.ft.c,
                                .combine = 'rbind',
                                .errorhandling="pass") %dopar%
    
    {
      i <- 16
      f.num<-mp$f.numbers[i]
      feat = mp$features[,i]
      simplePlot(feat)
      feat.padded.ft.c = mp$features.padded.ft.c[,i]
      #
      refs = mp$refs
      refs.padded.ft = mp$refs.padded.ft
      
      allmatches.feat <- match_feature(f.num, feat, feat.padded.ft.c,
                                       mp$refs, mp$refs.padded.ft)
      
      if (fit.matches){
        if (is.null(nrow(allmatches.feat))){
          allmatches.feat <- NULL
          specificity.score <- Inf
        } else {

          allmatches.feat <- fit_matches(allmatches.feat, feat, mp$refs)

        }
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
    features <- lapply(1:nrow(features), function(x) features[x,])
    fsp <-
      mclapply(features, function(feat){
        padded.feat <- feat %>% c(rep(0, pad.size),.)
        padded.feat[is.na(padded.feat)] <- 0
        feat.p.ft.c <- Conj(fftw::FFT(padded.feat))
        return(feat.p.ft.c)
      }, mc.cores = pars$par$ncores) %>% do.call(rbind,.) %>% t
    return(fsp)
  }
  
  coprep_features_and_refs <- function(feature.stack, ref.stack, ppm, roi, downsampling.factor=1){
    # coprep_features_and_refs(feature.stack, ref.stack, ppm, roi, downsampling.factor=8)
    # Assume feature and ref stacks have the same ppm axis (on the columns)
    # roi is in ppm 
    # browser()
      # Cut down ref stack to relevant region
      roi <- roi %>% vectInds(., ppm) # ppm
      reg <- roi %>% fillbetween
      ref.stack <- ref.stack[,reg]
      ppm <- ppm[reg]
      
      # Cut feature stack to relevant features
      in.range <- lapply(sat.list, function(s){
        # s <- sat.list[[1]]
        !all(is.na(range_intersect(roi, s$finalRegion)))
      
      }) %>% unlist %>% which
    
      # Thin it out some
      selected <- in.range#  seq(1, length(in.range), length.out=8) %>% in.range[.]
      feature.stack <- feature.stack[selected,]
      
      # Scale the features
      feature.stack <- lapply(1:nrow(feature.stack), function(x){
        feature.stack[x,] %>% scale_between %>% c
      }) %>% do.call(rbind, .)
      
      # simplePlot(feature.stack[1,])
      # simplePlot(ref.stack[1,], xvect=ppm)
      # stackplot(feature.stack[1:10], vshift = 10)
      
    # Downsample (if doing that)
      message('\t downsampling refs and features by a factor of ',downsampling.factor)
      
      ds.inds.ref <- downsample_inds(ppm, downsampling.factor)
      ref.stack <- ref.stack[,ds.inds.ref]
      ds.ppm <- ppm[ds.inds.ref]
      # simplePlot(ref.stack[1,], xvect = ds.ppm)
      # stackplot(ref.stack[1:10,], vshift = 10, xvect = ds.ppm)
      
      ds.inds.feat <- seq(1,ncol(feature.stack)) %>% downsample_inds(downsampling.factor)
      feature.stack <- feature.stack[,ds.inds.feat]
      # simplePlot(feature.stack[1,])
    
    # Move on to processing ref.stack
      feature.width <- ncol(feature.stack)
      
 ### Experiment with NCC  ###  ###  ###  ###  ###  ###  ###  ###  ### 
 

 ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
  
      refs.padded.ft <- prep_refs(ref.stack, feature.width)
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
                ref_downsampled_inds=ds.inds.ref,
                ppm=ppm))
    
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
  
  match_feature <- function(f.num, feat, feat.padded.ft.c,
                             refs, refs.padded.ft){
    
    # Locate best positions in all available refs
    
      allmatches.feat <-  foreach(r.num = 1:ncol(refs.padded.ft),
                                  ref = refs,
                                  ref.ft = refs.padded.ft,
                                  .combine = 'rbind',
                                  .errorhandling="stop") %do%
      {

        # r.num = 100
        # ref = refs[,r.num, drop = F]
        # ref.ft = refs.padded.ft[,r.num, drop = F]
        # message(r.num)
        # simplePlot(feat)
        #
        # simplePlot(ref %>% c)
        # plotly::plot_ly(data = data.frame(x=mp$ppm[mp$ref_downsampled_inds], y=c(ref)),
        #         x = ~x,
        #         y = ~y,
        #         type = "scatter",
        #         mode = "lines")
        
        # Cross-correlate to find locations and scores:
          # matches <- feature_match2ref_slim(f.num, r.num,
          matches <- feature_match2ref_pcc(f.num, r.num,
                                            feat, ref,
                                            max.hits = 100,#pars$matching$max.hits,
                                            r.thresh = .8)#pars$matching$r.thresh)
          

        # return a 1-row NA-filled placeholder
        if (is.null(matches) || nrow(matches) == 0){
          return(data.frame(
              feat = f.num,
              ref  = r.num,
              lag  = NA,
              rval = NA,
              pval = NA,
              pts.matched = NA,
              pts.feat = NA,
              feat.start = NA,
              feat.end = NA,
              ref.start = NA,
              ref.end = NA
          ))
        }
    
          return(matches)
      } 
      
      return(allmatches.feat)
  }
  
  fit_matches <- function(allmatches.feat, feat, ref.mat){

            # ref.mat <- mp$refs

    fits <- lapply(1:nrow(allmatches.feat), function(m)
    {
        message(m)
      # Get f and r indices for this row
        f <- allmatches.feat[m, 'feat']
        r <- allmatches.feat[m, 'ref']
        feat.pos <- allmatches.feat[m, c('feat.start','feat.end')] %>% as.numeric %>% fillbetween
        ref.pos <- allmatches.feat[m, c('ref.start','ref.end')] %>% as.numeric %>% fillbetween
  
      # Get spectral signatures which matched
        ref <- ref.mat[,r,drop = F] %>% c
        # simplePlot(c(ref))
        
        
      # Fit
        fit <- fit_leastSquares(feat[feat.pos] , ref[ref.pos], plots = FALSE, scale.v2 = TRUE)#;fit$plot
        
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

  plot_match <- function(match, feat, ref, ppm.margin = 1){
                plt_range <- function(roi.inds, ref.ppm, ppm.margin=1){

                  plot.start <- vectInds(ref.ppm[ref_start] + ppm.margin, ref.ppm) 
                  plot.end <- vectInds(ref.ppm[ref_end] - ppm.margin, ref.ppm)
                  
                  
                  reg <- plot.start:plot.end
                  return(reg)
                }
                
                browser()
                
                ref_start <- match$ref_start
                ref_end <- match$ref_end
                roi.inds <- c(ref_start, ref_end)
                
                ref.ppm <- mp$ppm[mp$ref_downsampled_inds]
              
              # Need to account for feat_start and ref_start alignment.
              # Use feat inds to get relative ROI inds, then regenerate actual ROI:
              
                feat.inds.matched <- match$feat_start
                
                offset.start <- (1-match$feat_start) # need relative inds for the edges
                offset.end <- length(feat) - match$feat_end

                roi.inds <- (ref_start + offset.start):(ref_end + offset.end)
                
              # Plot the ref-sized, NA-filled feature
              
                feat.scaled <- match$fit.intercept + match$fit.scale*feat
                
                feat.vect <- rep(NA, length(ref))
                
                feat.vect[roi.inds] <- feat.scaled
              
              # Then limit the range for both together:
                reg <- plt_range(roi.inds, ref.ppm, ppm.margin)
              
                colors.lines <- c("gray", rgb(0, 0, 1, alpha = 0.5))
                
                simplePlot(rbind(c(ref[reg]), c(feat.vect[reg])), xvect=ref.ppm[reg],linecolor = colors.lines)
                
            }

map_ref_xcorr <- function(xc.res, ref) {
  
  N     <- length(xc.res$xcorr)   # = f.len + r.len - 1
  r.len <- xc.res$ref_len
  f.len <- xc.res$feat_len
  
  # FFT linear correlation zero-lag alignment:
  # lag 0 corresponds to index f.len
  ref_start <- f.len
  ref_end   <- f.len + r.len - 1
  
  # Create padded full-length arrays
  vals <- matrix(NA, 3, N)
  rownames(vals) <- c("corr", "ref", "feat")
  
  # NCC already aligned to this indexing
  vals["corr", ] <- xc.res$xcorr
  
  # Place the REF
  vals["ref", ref_start:ref_end] <- as.vector(ref)
  
  # 'use' will be filled later in map_feat_xcorr()
  
  list(
    f.len      = f.len,
    r.len      = r.len,
    N          = N,
    inds = list(
      use        = NULL,
      feat_start = NA,
      feat_end   = NA,
      ref_start  = ref_start,
      ref_end    = ref_end,
      lag        = NA
    ),
    vals = vals
  )
}

map_feat_xcorr <- function(mapped, feat, lag) {
  
  N     <- mapped$N
  f.len <- mapped$f.len
  r.len <- mapped$r.len
  
  ref_start <- mapped$inds$ref_start  # = f.len
  
  vals <- mapped$vals
  
  # --------------------------
  # Place FEAT according to lag
  # lag = -(f.len-1):(r.len-1)
  #
  # feat_start = f.len + lag
  # feat_end   = f.len + lag + f.len - 1
  # -------------------------
  
  feat_start <- ref_start + lag       # = f.len + lag
  feat_end   <- feat_start + f.len - 1
  
  # validity check: must lie in 1..N
  if (feat_start < 1 || feat_end > N) {
    return(NULL)     # should never happen for valid FFT lags
  }
  
  # clear old feat
  vals["feat", ] <- NA
  
  # insert new feat
  vals["feat", feat_start:feat_end] <- feat
  
  # overlap: indices where both exist
  mapped$inds$use <- !(is.na(vals["feat", ]) | is.na(vals["ref", ]))
  
  mapped$inds$feat_start <- feat_start
  mapped$inds$feat_end   <- feat_end
  mapped$inds$lag        <- lag
  mapped$vals            <- vals
  mapped$inds$peak_loc   <- lag + f.len
  
  mapped
}

