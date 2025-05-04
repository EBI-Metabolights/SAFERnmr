# Statistal Annotation Tag Extraction

# Follows:
# - setup
# - load_data
# - protofeatures

# 

#   protofeatures

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

        storm_rnd1 = list()
        
        n.cores <- pars$par$ncores
        
  # Set up multicore
    
    # split up the ppm vector into chunks
    # but first, randomize it for load balancing 
    
    pfs.in.region <- (protofeatures$table$driver <= bounds[1] & protofeatures$table$driver <= bounds[2] ) %>% protofeatures$table[.,]
        
    pfs.rand <- sample(seq_along(rownames(pfs.in.region)))
      unrand <- order(pfs.rand, decreasing = FALSE)
    chunk.size <- ceiling(length(pfs.rand) / n.cores)
    pf.chunk.assignments <- split(pfs.rand, ceiling(seq_along(pfs.rand) / chunk.size))
    
    pfs.split <- pfs.in.region %>% split(seq_along(pfs.rand))
      
    pf.chunks <- lapply(pf.chunk.assignments, function(x) pfs.split[x])

    results <- mclapply(pf.chunks, function(pfs){
      # pfs <- pf.chunks[[1]]
      # 
      lapply(pfs, function(pf){
        # pf <- pfs[[1]]
        
        log_storm_core(xmat=xmat, ppm=ppm, half.window = 200, corrthresh = .8,
                        q=0.05, minpeak = 10, range.limit=400, plots=FALSE)
        
      })
      
    }, mc.cores = n.cores)
      
        message("Running storm on ",length(regions_subset), " provided protofeatures between ",ppm[bounds[1]]," and ",ppm[bounds[2]]," ppm.")
        
        storm_rnd1 <- 
              mclapply(regions_subset,
              # pblapply(regions_subset,
                  function (x) {
                    # TryCatch will 
                    tryCatch(
                      expr = {
                              # Set up the region

                              # Do storm_pairplay
                    
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
