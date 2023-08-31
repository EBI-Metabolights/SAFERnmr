#' My attempt at a pure STORM as described in the manuscript.
#'
#' @param xmat A matrix of spectral data (rows are spectra, columns are spectral points)
#' @param ppm A vector of the spectral points in ppm (optional, default is all columns of xmat)
#' @param b An integer giving the expansion parameter for the reference peak
#' @param corrthresh A numeric giving the minimum correlation value to be considered for inclusion (for both subset AND reference optimization)
#' @param q A numeric giving the p-value threshold for correlation significance (both subset AND reference optimization)
#' @param refSpec A vector of spectral data to use as the initial reference
#' @param ref.idx A vector of the spectral points (columns of xmat) to use as the initial reference
#'
#' @return A list with components "reconstructed" and "status". "reconstructed" is a matrix
#' containing the reconstructed metabolite concentrations (rows are samples, columns are metabolites).
#' "status" is a character string indicating whether the method converged successfully or failed.
#'
#' @export storm_pairplay
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_path geom_line geom_vline geom_hline ggtitle xlab ylab scale_y_continuous scale_x_continuous
#' @importFrom stringr str_pad
storm_mtj=function(xmat=NULL, ppm=NULL, b=30, corrthresh = .8,
                    q=0.05, refSpec=NULL, ref.idx=NULL,
                    driver = NULL){

############ Setup ##################################################  
    
  # Provide checks for input parameters
    # xmat:
    # - must be a matrix with ppm on columns, spectra on rows
    # - must not contain NAs
    # - ideally not scaled
    # 
    # ppm: 
    # - vector of ppm values corresponding to columns of xmat
    # 
    # b: 
    # - integer
    # - number of points to expand the window by on either side of the driver
    # 
    # corrthresh:
    # - correlation threshold for both 
    # 
    # How is the region selected?
    # - ref index is supplied
    # 
  
############ Initialize for the loop ########################################
  
    fullstack <- 1:nrow(xmat)
    subset.current= fullstack # cannot be same as .previous or loop won't run
    subset.previous <- 0
 
    corr <- rep(1,length(ref))
    covar <- ref
    ref.pass <- rep(TRUE,length(ref))
     
    
    # Make a ref.expanded object to hold ref information for extractPeaks #####
      
      ref.expanded <- list(wind = ref.idx %>% range %>% fillbetween)
      ref.pass <- rep(FALSE, length(ref.expanded$wind) )
      ref.pass[ref.expanded$wind %in% ref.idx] <- TRUE
      
      defwidth <- b
        
    # Set up exit status modes #####
      status <- "succeeded"
      fail.opts <- list("empty subset",          # empty subset
                        "subset degenerated",    # 1-3 spectra in the subset
                        "reference degenerated", # reference < 3 points
                        "did not converge")      # itlimit hit
    
    i=1        
    itlimit = 25
    
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
      
        # xmat[, ref.idx %>% range %>% fillbetween] %>% simplePlot
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
            subset.current = fullstack[sspass] # keep the subset of spectra positively correlated with the ref

  ## Update the ref ###########################################################################          
      
    # Identify the new driver ########
        ref.max <- next_driver(ref.profile = ref, current.driver = ref.max, 
                               ref.idx = ref.idx, behavior = 'simpleMax') %>% .$idx 
            
      # Alt: Use Static maxwidth (expansion amount, number of points = peak width in initial corrpocketpair) #######           
          maxwidth <- b
          bounds.sig.new <- ref.expanded$wind %>% range    # consider region = current region span +/- maxwidth points

    # Expand the window for reference identification ##############
        
        ref.expanded$newWind <- expand_window(window = bounds.sig.new,  # (bounds or vect are okay)
                                              within = seq_along(ppm),  # by = b+1) # other option; comment out lower line
                                              by = maxwidth*2,            # allow expansion by at least 1, 
                                              keep.nas = F)             # don't keep NAs

        ref.expanded$wind <- ref.expanded$newWind # "newWind" = window, but re-centered (as below in original storm)
        
        
    # STOCSY the new driver within subset.current and the widened window to get new ref ##############
      
        corr<-cor(xmat[subset.current, ref.expanded$wind], xmat[subset.current,ref.max])
        covar=cov(xmat[subset.current, ref.expanded$wind], xmat[subset.current,ref.max])

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
       
        ref=covar[ref.pass]
        
      # Also update the ref indices to match new ref
        
        ref.idx= ref.expanded$wind # the below all stems from this, which are ppm inds
        ref.idx=ref.idx[ref.pass]
        ref.expanded$wind <- ref.idx %>% range %>% fillbetween
        
    # (Plotting) #################
        plotrng <- range(ref.idx)
        plotreg <- c(min(ref.idx)-length(ref.idx)*1,max(ref.idx)+length(ref.idx)*1)
      
      # Finish the loop by updating the counter
        i <- i+1
    }
  
############ Finish up and return results ########################################################
# Finish up and return results
     # Set variables relevant to output
    # - handle failure mode cases
    # - ensure index ranges match up
    
      ref.max = next_driver(ref.profile = ref, current.driver = ref.max, 
                               ref.idx = ref.idx, behavior = 'samePk') %>% .$idx
      
      finalreg <- plotrng %>% fillbetween
      corr <- corr[ref.pass %>% which %>% range %>% fillbetween]
      covar <- covar[ref.pass %>% which %>% range %>% fillbetween]

      if(i == (itlimit-1)){status <- fail.opts[[4]]}
      
  return(list(subset = subset.current,
              finalRegion = finalreg,
              plotRegion = plotreg %>% fillbetween,
              ref.idx = ref.idx,
              ref.vals = ref,
              corr = corr,
              covar = covar,
              peak = ref.max,   # index in ref.expanded$wind
              pass = ref.pass,  # indices in ref.expanded$wind
              driver.initial = driver.init,
              status = status,  # see fail.opts
              iterations = i-1) # (completed iterations only) 
         )
 #################################################
}
