## Matching Features using Functions

  sats <- sat.list
  dataset.spectra <- xmat
  downsample.factor <- 8
  tol <- 0.5
  fit.matches <- TRUE

  # For each SAT, add the ref region ####
  
    sats.withranges <- lapply(sats, function(s){
      driver.ppm <- s$driver.initial %>% ppm[.]
      ref.range <- c(-tol, tol) + driver.ppm
        s$ref.range.tol <- vectInds(ref.range, ppm) # this also keeps them in bounds
      return(s)
    })
  
  # Match features and refs ####
  
    sats <- sats.withranges
    dataset.spectra <- xmat %>% t
    fit.matches <- TRUE

    # For each feature:
    spec.data.all <- mclapply(sats, function(s){
      # s <- sats[[500]]
      
      # Downsampling 
        
        f.num <- s$id
        feat.ds <- s_to_feat_ds(s, downsample.factor)
        ref.data <- refs_ds(s, dataset.spectra, ppm, downsample.factor)

        empty.result <- list(specificity = Inf,
                             matches     = NA)
        
        # with downsampling, we lose some tiny feats
        if ( sum( !is.na(feat.ds) )  <4 ){ return(empty.result) }
        
      # Match to ref regions using PCC
        
        matches.feat <- match_feature(f.num, feat.ds, ref.data$refs.ds)

      # Calculate specificity score
      
        specificity.score <- calc_specificity_feat_fast(matches.feat)
        if ( is.nan(specificity.score) ){ return(empty.result) }
      
      # Get fits 
      
        if (fit.matches) {
          
          matches.feat <- fit_matches(matches.feat, feat.ds, ref.data$refs.ds) # downsampled
          # m <- m+1
          # plot_match(matches[m,], feat.ds, refs.ds, ppm.ds, ppm.margin = tol)
        }
        
      # convert ref inds back to actual refmat inds, not ss.inds. (only after fitting, or breaks ref inds)
      matches.feat$ref <- ref.data$ss[matches.feat$ref] 
      
        return(list(
          specificity = specificity.score,
          matches     = matches.feat))      
        
    }, mc.cores = 6)
    saveRDS(spec.data.all, "spec.data.all.RDS")
    # lapply(spec.data, function(s){s$specificity}) %>% unlist %>% as.numeric
  
# Do the matching ####
  
    matches <- lapply(spec.data, function(s){
        s$matches
    }) %>% do.call(rbind, .) %>% na.omit
  
    # matches$rval %>% sort %>% plot(type="l")
  
    matches.sorted <- matches$rval %>% order(decreasing = TRUE) %>% matches[., ]
    
      # m <- 0
      m <- m + 10
      
      
      # plot_match(allmatches.feat[i,], feat, ref, ref.ppm, ppm.margin = 1)
            match <- matches.sorted[m, ]
            f <- match$feat
            s <- sats[[f]]
            
            feat.ds <- s_to_feat_ds(s, downsample.factor)
            r.data <- refs_ds(s, dataset.spectra, ppm, downsample.factor)
            
          # Get spectral signatures which matched
            r <- which(r.data$ss == match$ref)
            ref.ds <- r.data$refs.ds[,r,drop=F] #%>% as.double
            # fit <- fit_leastSquares(feat.ds, ref.ds %>% as.numeric, scale.v2 = FALSE, plots = TRUE); fit$plot
            # match$fit.scale <- 
            
            # match <- fit_matches(matches = match, feat.ds, ref.mat = ref.ds)
            
            message(spec.data.all[[f]]$specificity)
            plot_match(match, feat.ds, ref.ds, r.data$ppm.ds, ppm.margin = 1)
  
    
    # Sort features by specificity
      specificity <- df_out$specificity[match$feat == df_out$feat]


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
  
  match_features <- function(mp, fit.matches = FALSE){
    # mp <- match.pack
    
    # Par setup
    my.cluster <- safer_makeCluster(par, nfeats=length(mp$f.numbers))
        
    # Do matching (all refs, per feature):
    
    allmatches.feats <-  foreach(f.num = mp$f.numbers,
                                feat = mp$features,
                                # feat.padded.ft.c = mp$features.padded.ft.c,
                                .combine = 'rbind',
                                .errorhandling="pass") %dopar%
    
    {
      # i <- 16
      # f.num<-mp$f.numbers[i]
      # feat = mp$features[,i]
      # simplePlot(feat)
      # feat.padded.ft.c = mp$features.padded.ft.c[,i]
      # #
      # refs = mp$refs
      # refs.padded.ft = mp$refs.padded.ft
      
      allmatches.feat <- match_feature(f.num, feat, mp$refs)

      if (fit.matches) {
      
        if (is.null(nrow(allmatches.feat))) {
          allmatches.feat <- NULL
          specificity.score <- Inf
      
        } else {
      
          # always initialize columns so combine() never breaks
          allmatches.feat$fit.intercept <- NA_real_
          allmatches.feat$fit.scale     <- NA_real_
          allmatches.feat$rmse          <- NA_real_
      
          # rows that have a real match and should be fit
          valid_rows <- !is.na(allmatches.feat$rval)
      
          if (any(valid_rows)) {
            fitted <- fit_matches_vectorized(allmatches.feat[valid_rows, ],
                                             feat, mp$refs)
            
            # write fitted values back into original data frame
            allmatches.feat$fit.intercept[valid_rows] <- fitted$fit.intercept
            allmatches.feat$fit.scale[valid_rows]     <- fitted$fit.scale
            allmatches.feat$rmse[valid_rows]          <- fitted$rmse
          }
        }
      }
      
      specificity.score<-NA
      
      # ref.ppm <- mp$ppm[mp$ref_downsampled_inds]
      # i <- 0
      # i <- i + 1
      
      # plot_match(allmatches.feat[i,], feat, ref, ref.ppm, ppm.margin = 1)
      
          #   f <- allmatches.feat[m, 'feat']
          #   r <- allmatches.feat[m, 'ref']
          #   # feat.pos <- allmatches.feat[m, c('feat.start','feat.end')] %>% as.numeric %>% fillbetween
          #   # ref.pos <- allmatches.feat[m, c('ref.start','ref.end')] %>% as.numeric %>% fillbetween
          # # Get spectral signatures which matched
          #   ref <- ref.mat[,r,drop = F] %>% c
          #   # simplePlot(c(ref))
          # # Fit
          #   
          #   fit <- fit_leastSquares(feat[feat.pos] , ref[ref.pos], plots = TRUE, scale.v2 = FALSE);fit$plot
          #   match<- allmatches.feat[m, ]
          #   match$fit.intercept <- fit$fit[1]
          #   match$fit.scale <- fit$fit[2]
          #   plot_match(match, feat, ref, ref.ppm, ppm.margin = 1)
      
      return(allmatches.feat)
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
  
  match_feature <- function(f.num, feat, refs){
    
    # Locate best positions in all available refs
    
      allmatches.feat <-  foreach(r.num = 1:ncol(refs),
                                  ref = refs,
                                  # ref.ft = refs.padded.ft,
                                  .combine = 'rbind',
                                  .errorhandling="stop") %do%
      {

        # r.num = 1
        # ref = refs[,r.num, drop = F]  %>% plot(type="l")
        # message(r.num)
        # feat %>% plot(type="l")
        #
        # simplePlot(ref %>% c)
        # plotly::plot_ly(data = data.frame(x=mp$ppm[mp$ref_downsampled_inds], y=c(ref)),
        #         x = ~x,
        #         y = ~y,
        #         type = "scatter",
        #         mode = "lines")
        
        # Cross-correlate to find locations and scores:
          
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

  fit_matches <- function(matches, feat, ref.mat){
    
    if (is.null(nrow(matches))) {
      return(NULL)
    } else {
  
      # always initialize columns so combine() never breaks
      matches$fit.intercept <- NA_real_
      matches$fit.scale     <- NA_real_
      matches$rmse          <- NA_real_
  
      # rows that have a real match and should be fit
      valid_rows <- !is.na(matches$rval)
      
      if (any(valid_rows)) {
        fitted <- fit_matches_vectorized(matches[valid_rows, ],
                                         feat, ref.mat)
        
        # write fitted values back into original data frame
        matches$fit.intercept[valid_rows] <- fitted$fit.intercept
        matches$fit.scale[valid_rows]     <- fitted$fit.scale
        matches$rmse[valid_rows]          <- fitted$rmse
      }
    }
    return(matches)
  }
  
  fit_matches_vectorized <- function(matches, feat, ref.mat) {
    
    M <- nrow(matches)
    if (M == 0) return(matches)
  
    f.len <- length(feat)
  
    # Preallocate
    a_unscaled <- b_unscaled <- rep(NA_real_, M)
    rmse       <- rep(NA_real_, M)
  
    # -----------------------------
    # Precompute v1 stats ONCE
    # -----------------------------
    # Raw feature
    v1 <- feat
    # scaled feature for RMSE comparison
    fr1 <- range(v1, na.rm = TRUE)
    v1s <- (v1 - fr1[1]) / diff(fr1)
  
    for (j in seq_len(M)) {
  
      r  <- matches$ref[j]
      rs <- matches$ref.start[j]
      re <- matches$ref.end[j]
  
      if (is.na(r) || is.na(rs) || is.na(re) || re < rs)
        next
  
      v2 <- ref.mat[rs:re, r]
  
      # Determine usable overlap region
      use <- !(is.na(v1) | is.na(v2))
  
      if (sum(use) < 3)
        next
  
      x  <- v1[use]
      y  <- v2[use]
  
      # ---------------------------
      # Unscaled regression
      # ---------------------------
      xm <- mean(x)
      ym <- mean(y)
      dx <- x - xm
      dy <- y - ym
  
      var_x <- sum(dx * dx)
      if (!is.finite(var_x) || var_x == 0) next
  
      b <- sum(dx * dy) / var_x
      a <- ym - b * xm
  
      a_unscaled[j] <- a
      b_unscaled[j] <- b
  
      # ---------------------------
      # Scaled RMSE
      # ---------------------------
      # scale v2 on its window
      fr2 <- range(v2, na.rm = TRUE)
      if (!is.finite(diff(fr2)) || diff(fr2) == 0) next
  
      y2s <- (v2 - fr2[1]) / diff(fr2)
      v1s_use <- v1s[use]
  
      # compute scaled regression for RMSE
      ym2s <- mean(y2s[use])
      dy2s <- y2s[use] - ym2s
  
      dx1s_use <- v1s_use - mean(v1s_use)
      var_x1s  <- sum(dx1s_use * dx1s_use)
      if (!is.finite(var_x1s) || var_x1s == 0) next
  
      b_s <- sum(dx1s_use * dy2s) / var_x1s
      a_s <- ym2s - b_s * mean(v1s_use)
  
      fit_vals <- a_s + b_s * v1s_use
      rmse[j] <- sqrt(mean((y2s[use] - fit_vals)^2))
    }
  
    matches$fit.intercept <- a_unscaled
    matches$fit.scale     <- b_unscaled
    matches$rmse          <- rmse
  
    matches
  }

  plot_match <- function(match, feat, ref, ref.ppm, ppm.margin = 1) {
    
    # match<- matches[m,]
    ref.start <- match$ref.start
    ref.end   <- match$ref.end
  
    # 1) Full-length NA vector for plotting feature in ref coordinates
    feat.scaled <- match$fit.intercept + match$fit.scale * feat
    feat.vec <- rep(NA_real_, length(ref))
    feat.vec[ref.start:ref.end] <- feat.scaled
  
    # 2) Determine plotting window in ppm units
    plot_start_ppm <- ref.ppm[ref.start] + ppm.margin
    plot_end_ppm   <- ref.ppm[ref.end]   - ppm.margin
  
    # Determine ppm direction
    ascending <- ref.ppm[1] < ref.ppm[length(ref.ppm)]
    
    # Clamp ppm request to valid range
    lo <- min(ref.ppm)
    hi <- max(ref.ppm)
    
    plot_start_ppm_clamped <- max(lo, min(hi, plot_start_ppm))
    plot_end_ppm_clamped   <- max(lo, min(hi, plot_end_ppm))
    
    # Build region
    if (!ascending) {
      reg <- which(ref.ppm >= plot_start_ppm_clamped &
                   ref.ppm <= plot_end_ppm_clamped)
    } else {
      reg <- which(ref.ppm <= plot_start_ppm_clamped &
                   ref.ppm >= plot_end_ppm_clamped)
    }
    
    # If nothing, fall back to the matched region
    if (length(reg) == 0) {
      reg <- ref.start:ref.end
    }
    
    simplePlot_x(
      rbind(ref[reg], feat.vec[reg]),
      xvect = ref.ppm[reg],
      linecolor = c("gray", rgb(0,0,1,0.5)),
      linewidth = c(.5, 1),
      title.str = paste("Feature", match$feat, " x Ref", match$ref),
      st.str = paste("r =", round(match$rval,3)),
    )
  
  }

  calc_specificity_feat <- function(matches){
    # matches <- matches.feat
    matches %>%
      group_by(feat,ref) %>%
      mutate(rval_norm = rval / max(rval)) %>%
      summarise(
        sum_r_in_ref = sum(rval_norm),  # or sum(rval_norm) / n()
        # n_refs_with_hits = n(),         # diagnostic
        .groups = "drop"
      ) %>%
      # group_by(feat) %>%
      summarise(
        specificity = mean(sum_r_in_ref, na.rm = TRUE),
        n_refs = n(),   # optional diagnostic
        .groups = "drop"
      ) %>% .["specificity"] %>% as.numeric
    
  }

  calc_specificity_feat_fast <- function(matches) {
  
    # Drop NA rows immediately
    matches <- matches[!is.na(matches$rval), ]
    if (!nrow(matches)) return(Inf)
    
    # 1. Compute max rval per ref (vectorized)
    max_rval_per_ref <- tapply(matches$rval, matches$ref, max, na.rm=TRUE)
  
    # 2. Normalize rval using lookup table (no group_by)
    rval_norm <- matches$rval / max_rval_per_ref[matches$ref]
  
    # 3. Sum normalized rval per ref
    sum_r_per_ref <- tapply(rval_norm, matches$ref, sum)
  
    # 4. Specificity = mean across refs
    mean(sum_r_per_ref)
  }
  
  calc_specificity <- function(matches){
    matches %>%
      group_by(feat,ref) %>%
      mutate(rval_norm = rval / max(rval)) %>%
      summarise(
        sum_r_in_ref = sum(rval_norm),  # or sum(rval_norm) / n()
        # n_refs_with_hits = n(),         # diagnostic
        .groups = "drop"
      ) %>%
      # group_by(feat) %>%
      summarise(
        specificity = mean(sum_r_in_ref, na.rm = TRUE),
        n_refs = n(),   # optional diagnostic
        .groups = "drop"
      )
    
  }
    
  s_to_feat_ds <- function(s, downsample.factor){
    # Feature
    # s <- sats[[14]]
    ss <- s$subset
    feat <- rep(NA, length(s$covar))
    feat[ s$pass ] <- s$covar[ s$pass ]
    ds.inds.feat <- downsample_inds(feat, downsample.factor)
    feat.ds <- feat[ds.inds.feat]
    
    return(feat.ds)  
    
  }
    
  refs_ss_ds <- function(s, refs, ppm, downsample.factor){
    
    # Refs
    
    ref.reg <- s$ref.range.tol %>% fillbetween
    ppm.reg <- ppm[ref.reg]
    refs.ss <- refs[ref.reg, s$subset]
    ds.inds.ref <- downsample_inds(refs.ss[,1], downsample.factor)
    refs.ds <- refs.ss[ds.inds.ref, ]# %>% .[,1] %>% plot(type="l")
    ppm.ds <- ppm.reg[ds.inds.ref]# %>% plot(type="l")  
    
    return(
            list(refs.ds = refs.ds,
                 ppm.ds = ppm.ds,
                 ss = s$subset)
    )
    
  }
