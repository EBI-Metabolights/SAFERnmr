# Statistal Annotation Tag Extraction

# Follows:
# - setup
# - load_data
# - protofeatures

# Parameter setup ####
    # Override for now:
    only.region.between <- pars$corrpockets$only.region.between
    # only.region.between <- pars$corrpockets$only.region.between
      if (is.null(only.region.between))                       # which ppms to run fse between
        {only.region.between <- range(ppm)}                   #   (default is all)
    correlation.r.cutoff <- pars$storm$correlation.r.cutoff   # rvalue cutoff for both subset selection (ref shape) and ref update (STOCSY)
    q <- pars$storm$q                                         # q param from storm (pval cutoff after mhtc)
    b <- pars$storm$b                                         # number of peak widths to expand ref by on each side

  # Plotting
    number.of.plots <- pars$storm$number.of.plots             # pdf of all extracted features will be plotted. Choose
                                                              # only 150 of these (evenly spaced) or suffer the
                                                              # consequences...
                                                    
    plot.location <- paste0(tmpdir,"/plots/")                     # where to put the plot (just dump into run folder)
    dir.create(plot.location, showWarnings = FALSE)

 
# Run code ####

     
        bounds <- vectInds(only.region.between, ppm)
        testregion <- bounds[2]:bounds[1]
        
        colwithPair <- pocketPairs$peakMap %>% is.na %>% "!"(.) %>% 
          t %>% rowSums(na.rm = TRUE) %>% ">"(.,0) %>% which
        regions_subset <- (colwithPair %in% testregion) %>% which
        storm_rnd1 = list()
        
        ncores <- pars$par$ncores
        chunks <- lapply()
        
        message("Running storm on ",length(regions_subset), " provided protofeatures between ",ppm[bounds[1]]," and ",ppm[bounds[2]]," ppm.")
        
        storm_rnd1 <- 
              mclapply(regions_subset,
              # pblapply(regions_subset,
                  function (x) {
                    # TryCatch will 
                    tryCatch(
                      expr = {
                              # Set up the region
                                # x <- regions_subset[1793]
                                
                                driver <- colwithPair[x]
                                peakPos <- pocketPairs$peakMap[,driver] %>% is.na %>% "!"(.) %>% which
                                pair.region <- pocketPairs$regions[peakPos,driver]
          
                              # Do storm_pairplay
                    
                                # Set params
                                  pw <- span(peakPos)/2 %>% ceiling
                                  wind <- pair.region
                                  shift <- range(ppm[pair.region])
                                  
                                # Use original covariance signal within corr bounds as shape seed
                                # (could also use best spectrum index)
                                  shape <- pocketPairs$cov[peakPos,driver]
                                  bestSpec = cor( xmat[ ,wind] %>% t, shape ) %>% which.max
                                  
                                # Do the storm
                                  
                                  res <- storm_pairplay(xmat, ppm,
                                                              b = (pw * b) %>% ceiling, corrthresh = correlation.r.cutoff, q = q,
                                                              minpeak = noisewidth, refSpec = shape, ref.idx = pair.region,
                                                              driver = driver)
                                  res$cpp.driver <- driver
                                  
                                  return(res)
                        
                      },
                      error = function(cond){
                        return('setup error')
                        }
                      )
                  }, mc.cores = pars$par$ncores
              )
        # Note: errors in the loop are captured and passed out as strings.
        # NULL elements are not possible, although parts of an element could be.
        # Those are checked below.  
