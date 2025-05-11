#' log_storm_core: 
#' Locally Optimized Global STORM
#' 
#' Run modified STORM on the provided spectral region and ref shape.
#' Built for accepting corrPocketPairs results. Notes:
#'
#' STORM: Joram Posma's STORM has been adapted and optimized to accept these
#   protofeatures (corrPocketPairs) in the following ways:
# - first, since many of the protofeatures are noise, we provide failure modes
#   and reporting for the following cases:
#   "empty subset",          # empty subset (no spectrum contains signature)
#   "subset degenerated",    # 1-3 spectra in the subset (not enough spectra to
#                              get a reliable correlation)
#   "reference degenerated", # signature degenerates to include < 3 points (not
#                              meaningful to correlate shapes)
#   "did not converge"       # subset continues to change after 24 iterations
#
# - additionally, the correlation r and p-value cutoff q are both used during
#   both the subset selection and reference update steps.
# - we also remove any regions of the reference for which there are fewer than
#   minpeak values after r and p value thresholding. This helps avoid noise.
#
# STORM extracts meaningful features using protofeatures to define the region of
# interest and a rough sketch of the feature shape highly correlated with each
# spectral point. In the future, HCA could be used to cluster potential starting
# feature shapes correlated with each driver, or the nonoptimal subset for each
# point could be re-STORMed to detect any other feature shapes present. It is
# also perfectly reasonable to combine feature shapes from different STORM runs
# for a given dataset, as these comprise a list of somewhat independently tested
# feature shapes, and duplication is not an issue.

#'
#'
#' @param xmat A matrix of spectral data (rohws are spectra, columns are spectral points)
#' @param ppm A vector of the spectral points in ppm (optional, default is all columns of xmat)
#' @param b An integer giving the expansion parameter for the reference peak
#' @param corrthresh A numeric giving the minimum correlation value to be considered for inclusion (for both subset AND reference optimization)
#' @param q A numeric giving the p-value threshold for correlation significance (both subset AND reference optimization)
#' @param minpeak An integer giving the minimum number of points allowed in a run of significant points in the reference
#' @param refSpec A vector of spectral data to use as the initial reference
#' @param ref.idx A vector of the spectral points (columns of xmat) to use as the initial reference
#' @param range.limit A vector of the spectral points (columns of xmat) to use as the initial reference
#'
#' @return A list with components "reconstructed" and "status". "reconstructed" is a matrix
#' containing the reconstructed metabolite concentrations (rohws are samples, columns are metabolites).
#' "status" is a character string indicating whether the method converged successfully or failed.
#'
#' @export log_storm_core
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_path geom_line geom_vline geom_hline ggtitle xlab ylab scale_y_continuous scale_x_continuous
#' @importFrom stringr str_pad
log_storm_core=function(p=NULL, xmat=NULL, ppm=NULL, half.window = 200, corrthresh = .8,
                        q=0.05, minpeak = 10, range.limit=400, plots=FALSE){

############ Setup ##################################################  

    p <- pf
    half.window = 200
    corrthresh = .8
    q=0.05
    minpeak = 10
    range.limit=400
    plots=TRUE
                        
  # Select protofeature
  
    # i <- 11
    # p <- protofeatures[i,]
    # minpeak <- noiseWidth * noise.width.multiple
    hws <- half.window
    
########################################################################################################################
    # Expand protofeature
    
    pexp <- expand_protofeature(p, xmat, ppm, half.window)
    driver <- pexp$driver
    
    wind <- pexp$specRegion.inds
      
      specRegion = pexp$specRegion
      ppmRegion = pexp$ppmRegion
  
      
    # Why is the plot reversing the peaks?
    mask <- pexp$peak.mask>0
    
    ref.idx <- wind[mask]
    refSpec <- pexp$cv[mask]

########################################################################################################################
    
  # Window cannot exceed ppm region
    
    if (is.null(ppm)){ppm <- 1:ncol(xmat)}
    indLimits <- c(1,length(ppm))
  
  # Set ref based on refSpec index and ref.idx
  
    if (length(refSpec) == 1){
      # Assume this is a row index of xmat
        ref.init <- xmat[refSpec, ref.idx]
    } else {
      # Assume this is a ref shape, use as-is (but without NAs)
        ref.init <- refSpec[!is.na(refSpec)]
    }
    
    ref <- ref.init
    
    # For the first driver, use pkMax
    if (is.null(driver)){
      driver.init <- next_driver(ref.profile = ref, current.driver = NULL,
                                 ref.idx = ref.idx, behavior = 'maxPk') %>% .$idx
      
    } else {driver.init <- driver}
      ref.max <- driver.init

############ Initialize for the loop ########################################
  
    fullstack <- 1:nrow(xmat)
    subset.current= fullstack # cannot be same as .previous or loop won't run
    subset.previous <- 0
 
    corr <- NULL
    covar <- NULL

    # The wind (window) vect hold the full spectral region inds, but ref.idx is just for passing vals #####
      
      ref.pass <- rep(FALSE, length(wind))
      ref.pass[wind %in% ref.idx] <- TRUE

    # Set up exit status modes #####
      status <- "succeeded"
      fail.opts <- list("empty subset",          # empty subset
                        "subset degenerated",    # 1-3 spectra in the subset
                        "reference degenerated", # reference < 3 points
                        "did not converge")      # itlimit hit
    
    i=1        
    itlimit = 25
    
    plots.baseName <- 'sat_evolution_'
    
    if (plots){
      plot_protofeature(p, 
                        half.window = half.window, ppm = ppm, xmat=xmat, 
                        bgplot = 'stack', line.shape = 'covar', line.color = "corr", 
                        showPeaks = TRUE, ref.mask = ref.idx, show.mask.bounds = TRUE)
    }

############ Run storm loop ###################################################################

  # Run storm loop until the subset contains exactly the same spectra as the previous one.
  # Or, run while not all previous subset spectra are included in the current subset.
  # - usually, this means that, in the previous iteration, subset.current did not change 
  #   from subset.previous
  #   subset.current is always smaller unless subset.previous is reset to fullstack
    # original: while(length(which(!(subset.previous %in% subset.current)))>0){
    
    while( !all(subset.previous %in% subset.current) & i < itlimit){ 
      
  ## Update the subset ########################################################################
        
    # Update subset.previous to keep track of this loop's starting point #########
      
      subset.previous <- subset.current
      
    # Pull out the data for ref points in the previous subset #############
      
        # xmat[, wind %>% range %>% fillbetween] %>% stackplot
        xmatr=xmat[fullstack, ref.idx]

    # Pull out the subset of spectra which appear to contain the ref ##########################
    
      
          # Calculate similarity between ref shape and shapes in previous subset #####
            
            r=cor(t(xmatr), ref)     # try correlating ref shape to spectra
              # this may help with the stats...
              #   use = "pairwise.complete.obs", method = "pearson"
            
            
            
          # Get the inds of the significantly positively correlated spectra to the ref #####

            a=-abs(r * sqrt((length(r)-2)/(1-r^2)))
            pval=2*pt(a,(length(r)-2))
            
            # A pval AND rval threshold is used. The rval threshold is necessary 
            # to make sure the ref shape is represented faithfully in the subset
            # spectra. If this value increases, the ref shape will be more faithfully
            # preserved through the iterations (although it can grow; points not
            # incorporated into the ref before won't be used in the subset selection).
            
            sspass <- (pval<q & r>corrthresh) %>% which
            
            
              if(length(sspass) < 3) # Failure modes 1 and 2
              {
                  plotrng <- c(min(ref.idx),max(ref.idx))
                  plotreg <- c(min(ref.idx)-length(ref.idx)*1,max(ref.idx)+length(ref.idx)*1)
                  ref.max <- NA
                  corr <- rep(NA, length(plotrng %>% fillbetween))
                  covar <- corr
                  covar[(plotrng %>% fillbetween) %in% ref.idx] <- ref
                  ref.pass <- rep(TRUE, length(corr))
                  status <- fail.opts[[length(sspass)+1]]
                break
              }
            
            
            
          # Subset from the full spectral matrix stack #####
            # subset.current = subset.previous[sspass] # keep the subset of spectra positively correlated with the ref
            subset.current = fullstack[sspass] # keep the subset of spectra positively correlated with the ref
            # xmat[subset.current, ref.idx %>% range %>% fillbetween] %>% simplePlot(xvect = ref.idx %>% range %>% fillbetween)
            # ref %>% simplePlot(xvect = ref.idx)
            # xmat[subset.current, ref.idx %>% range %>% fillbetween] %>% stackplot(xvect = ref.idx %>% range %>% fillbetween)
            
          
  ## Update the ref ###########################################################################          
      
    # Identify the new driver ########
        # if (is.na(ref.max) | is.null(ref.max)){browser()}
        
        plot(x = ref.idx, y = ref)
        abline(v = ref.max)
        ref.max <- next_driver(ref.profile = ref, current.driver = ref.max, 
                               ref.idx = ref.idx, behavior = 'samePk') %>% .$idx 
            
        # ref %>% 
        #   simplePlot(xvect = ref.idx) + 
        #   geom_vline(xintercept = ref.max) + 
        #   geom_vline(xintercept = c(min(wind), max(wind)))
        # xmat[subset.current, ref.idx %>% range %>% fillbetween] %>% 
        #   simplePlot(xvect = ref.idx %>% range %>% fillbetween) + 
        #   geom_vline(xintercept = ref.max) + 
        #   geom_vline(xintercept = c(min(wind), max(wind)))

    # Expand the window for reference identification ##############
          # Center on new max, hws points in either direction
        
        wind <- expandRef_simple(wind, ref.max, hws, ppm)
        
    # STOCSY the new driver within subset.current and the widened window to get new ref ##############
      
        corr<-cor(xmat[subset.current, wind], xmat[subset.current,ref.max])
        covar=cov(xmat[subset.current, wind], xmat[subset.current,ref.max])
        
        # plot_protofeature(p = data.frame(driver = ref.max),
        #                   half.window = hws, ppm = ppm,
        #                   xmat = xmat[subset.current,],
        #                   bgplot = 'stack', line.shape = 'covar', line.color = "corr",
        #                   showPeaks = FALSE, ref.mask = ref.idx)
        
        # plot(corr); abline(h = corrthresh); abline(v = which((wind %>% range %>% fillbetween) == ref.max))
        
    # Clean up the ref with pval, rval, and runlength filtering #######################
      
      # Determine which corrs are significant 
      
        a=-abs(corr * sqrt((length(corr)-2)/(1-corr^2)))
        pval=2*pt(a,(length(corr)-2))

        # The new reference becomes the covariance of the new subset derived above. 
        # This entails thresholding the correlation profile, extracting the 
        # corresponding covariance profile points, and updating the ref inds. The
        # STOCSY correlation to max doesn't have to be super high - just positive. 
        
      # Filter using pval and correlation 
      
        ref.pass <- (pval<q & corr>corrthresh)
        # plot(ref.pass %>% as.integer)
        
      # Remove any runs that are < minpeak. This helps control for expansion
      # by a bunch of noise peaks.
       
        ref.pass <- (ref.pass %>% as.integer %>% runs.labelBy.lengths) > minpeak
        # plot(ref.pass %>% as.integer)
        
      # Check to make sure the ref is valid
      # - contains at least 3 valid points (absolute minimum for a meaningful peak shape)
      # - should there be a contiguous point requirement here?
        
          if ((ref.pass %>% sum(na.rm = TRUE)) < 3){
            plotrng <- c(min(ref.idx),max(ref.idx))
            plotreg <- c(min(ref.idx)-length(ref.idx)*1,max(ref.idx)+length(ref.idx)*1)
            ref.max <- NA
            corr <- rep(NA, length(plotrng %>% fillbetween))
            covar <- corr
            ref.pass <- rep(TRUE, length(corr))
            status <- fail.opts[[3]]
            break
          } # Failure mode 3
        
    # Extract the new ref shape from the thresholded covariance profile #################
       
        ref <- covar[ref.pass] # Update the ref shape using passing ref vals
        ref.idx <- wind[ref.pass] # Also update the ref indices to match new ref
        
      # Plot
      if (plots){
        plot_protofeature(p = data.frame(driver = ref.max),
                          half.window = hws, ppm = ppm,
                          xmat = xmat[subset.current,],
                          bgplot = 'stack', line.shape = 'covar', line.color = "corr",
                          showPeaks = FALSE, ref.mask = ref.idx, show.mask.bounds = TRUE)
        # simplePlot(ref, xvect = ref.idx)
      }
      # Finish the loop by updating the counter
        i <- i+1
    }
  
############ Finish up and return results ########################################################
# Finish up and return results
    # Set variables relevant to output
    # - handle failure mode cases
    # - Use case of SATs is 
    #   - plotting (needs to be indexable on xmat)
    #   - matching (shapes available without use of xmat)
    #     - ppms are important
    #     - gaps are important - assuming that NAs are handled in matching
      last.driver <- ref.max
      ref.max = next_driver(ref.profile = ref, current.driver = ref.max, 
                               ref.idx = ref.idx, behavior = 'samePk') %>% .$idx
      
      wind <- wind %>% range # undo with wind %>% fillbetween
      # Just keep the whole corr and covar
      corr <- corr # [ref.pass %>% which %>% range %>% fillbetween]
      covar <- covar # [ref.pass %>% which %>% range %>% fillbetween]
      ref.pass <- which(ref.pass) # in wind; undo with wind %>% fillbetween %>% .[ref.pass]
      peak <- ref.max # not a driver, but peak
      last.driver
      
      if(i == (itlimit-1)){status <- fail.opts[[4]]}
      
      if (plots){
        # Print into video or grid
        # 
      }
      
  return(list(protofeature = p,
              subset = subset.current,
              finalRegion = wind,
              ref.idx = ref.idx, # ppm inds for ref
              ref.vals = ref,    # ref covariance shape with NAs
              corr = corr,      # passed from last update
              covar = covar,    # passed from last update
              peak = ref.max,   # index in wind
              pass = ref.pass,  # indices in wind
              driver.initial = driver.init,
              status = status,  # see fail.opts
              iterations = i-1) # (completed iterations only) 
         )
 #################################################
}

expandRef_simple <- function(wind, ref.max, hws, ppm, recenter = FALSE){
    # Default is just return the same window - no expansion/recentering
    if (recenter){
      wind <- (ref.max - hws):(ref.max + hws)
    
      # On the ends of the spectra, adjust back in frame
      if (any(wind < 1)){
        wind <- 1:hws
      } else {
        if (any(wind > length(ppm))){
          wind <- (length(ppm)-hws):length(ppm)
        }
      }
    }
    return(wind)
}
