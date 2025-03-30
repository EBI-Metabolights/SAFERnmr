# fse2
params_loc <- '/Users/mjudge/Downloads/safer/results/1739774749/params.yaml'
pars <- yaml::yaml.load_file(params_loc, eval.expr = TRUE)

fse <- function(pars){
  
  tmpdir <- pars$dirs$temp

################ Get MTBLS1 data from RDS ##################

    X_raw <- readRDS(pars$files$spectral.matrix)
    # X_raw <- readRDS(paste0(tmpdir, "/spectral.matrix.RDS"))
      xmat <- X_raw[-1,]          # spectral matrix (each row is a spectrum; 
                                  # doesn't require alignment. normalization ok
                                  # but not necessary. scaling, no.)
      ppm <- X_raw[1,]            # ppm vector (corresponding to cols of xmat)
      digital.res <- ppm %>% diff %>% mean %>% abs # ppm per element
      
      message('digital resolution of xmat is ', digital.res, ' ppm.')
      if (!is.null(pars$opts$npoints)){
        if (pars$opts$npoints < length(ppm)){
          
          # Re-interpolate dataset spectra to a lower number of points to save compute
            rs <- resample_spectra(xmat, ppm, npoints = pars$opts$npoints, cores = pars$par$ncores)
            ppm <- rs$ppm
            xmat <- rs$spectra
            if (nrow(rs$spectra) < pars$tina$min.subset){stop('The minimum number of spectra is: ', pars$tina$min.subset,'. The submitted dataset only has ', nrow(rs$spectra),'. SAFER terminated.')}
            
            rm(rs)
          
          digital.res <- ppm %>% diff %>% mean %>% abs # ppm per element
          message('Using opts:npoints @ ',pars$opts$npoints, ' points. 
                  \nNew xmat digital resolution: ', digital.res)

        }
      }
      # if not set, do nothing (warning is printed)
      
      
################ Set up parameters ##################
    
  # Corr Pocket Pairs 

    half.window <- (pars$corrpockets$half.window / digital.res) %>% ceiling  
        if (half.window > 1000){stop('Window size is too large. Please keep to < 1000 points (~ ', round(1000 * digital.res, 4),' ppm for this dataset).')}
                                # Window for the sliding correlation calculation. 
                                # This x 2 should capture any 2 adjacent resonances
                                # in a multiplet. 
                                # Provided in ppm, converted here to (column) elements

    noise.percentile <- pars$corrpockets$noise.percentile     
                                # noise characterization... For every spectral
                                # point, we calculate a correlation peak. If you
                                # average all of the resulting peak shapes, 99%
                                # of them will be at least n points wide. All of 
                                # the peaks have a width of at least 3 (one point
                                # on either side of the driver, due to the way
                                # a peak is defined). What fraction contain 5
                                # points? 30 points? This parameter determines the 
                                # fraction of all spectral corr peaks that noise
                                # should fit within, and sets the noise width 
                                # accordingly. Higher is more permissive.
                                
     cpp.rcutoff <- pars$corrpockets$rcutoff
                                # correlation cutoff for picking highest secondary
                                # peak in corrpocket pair extraction. Generally
                                # ~ 0.75 should do fine.
       
  # STORM

    only.region.between <- pars$corrpockets$only.region.between
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
    
################ Use corrPocketPairs to extract likely j-pairs ##################
 
correlation_pocket_pairs <-  function(x, ppm, ws, reg = NULL, plotHeatmap = FALSE, wdlimit = 0.99,
                             rcutoff = 0.5){
  
  x <- xmat
  ppm <- ppm
  ws <- 100
  reg <- NULL # indices of ppm
  plotHeatmap <- FALSE
  wdlimit <- 0.95
  plotHeatmap <- FALSE
  rcutoff <- 0.5
  pars$corrpockets$only.region.between <- c(-0.5,10.5)   
      
      
      # ####
################ Run STORM on these corrpairs ##################  

    # Run code    
        bounds <- vectInds(only.region.between, ppm)
        testregion <- bounds[2]:bounds[1]
        
        colwithPair <- pocketPairs$peakMap %>% is.na %>% "!"(.) %>% 
          t %>% rowSums(na.rm = TRUE) %>% ">"(.,0) %>% which
        regions_subset <- (colwithPair %in% testregion) %>% which
        storm_rnd1 = list()
        
        
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
          
################ Report run stats  ######
        
          fmodes <- lapply(storm_rnd1, 
                           function(x) {
                             if (is.character(x)){return(x)} # this will get any errors from setup or storm
                             if (x$status == 'succeeded'){
                               if (any(is_nullish(x))){
                                 # Even if STORM succeeded, it may contain NULLs in some return value elements. 
                                 return(  paste0(  is_nullish(x) %>% which %>% names, " contains NULL")) 
                               }
                             }
                             return(x$status)
                           })
          
          succeeded <- lapply(fmodes, function(x) x == 'succeeded') %>% unlist
          failed <- !succeeded
          
          message(str_c("Failed iterations (count): ", sum(failed), " (",
                        (sum(failed)/length(regions_subset) * 100) %>% round, " %)"))
          
          message(str_c("Succeeded iterations (count): ", sum(succeeded), " (",
                        (sum(succeeded)/length(regions_subset) * 100) %>% round, " %)"))
          
          # Print out breakdown of statuses
            fm <- fmodes %>% plyr::ldply(rbind)
            colnames(fm) <- 'status'
            fm <- fm %>% group_by(status) %>% count %>% as.data.frame
            
            print(fm)
      
            
# ##### Save #####
  
    fse.result <- list(storm_features = storm_rnd1[succeeded],
                       xmat = xmat,
                       ppm = ppm,
                       noisewidth = noisewidth)
            
    fse.result %>% test_nullish    
    message("Saving results...")

    saveRDS(fse.result, paste0(tmpdir, "/fse.result.RDS"))
    # fse.result <- readRDS(paste0(tmpdir, "/fse.result.RDS"))
    
################ Plotting Results #######################

  # THIS DOESN'T WORK WHEN XDATA ARE TOO BIG
  
    # plot.filename <- paste0('features_',
    #                         pars$study$id,
    #                         "_np_",noise.percentile,
    #                         "_r_",correlation.r.cutoff,
    #                         "_b_",b,".pdf")     # what to name the plot file

        # # Plot all the storm results in grid (postage stamp) format
        # 
          
            # everyNth <- every_nth(select = number.of.plots, 
            #                       from = sum(succeeded))
            # 
            # plot_stormRefRegions_grid(xmat, ppm,
            #                           storm_rnd1[succeeded %>% which %>% .[everyNth]], # if not doing a small region
            #                           plotLoc = plot.location,
            #                           filename = plot.filename,
            #                           calcStocsy = FALSE, n_xticks = 4)
       # Plot all in region
          # use <- seq(12810, 12972)
            # plot_stormRefRegions_grid(xmat,ppm,
            #                           storm_rnd1[succeeded %>% which %>% .[use]], # if not doing a small region
            #                           plotLoc = plot.location,
            #                           filename = str_c(plot.filename,".citrate.pdf"),
            #                           calcStocsy = FALSE,n_xticks = 4)
    
    message("\nData written to ", tmpdir, "/fse.result.RDS")
    message("\nFeature Shape Extraction completed.\n\n\n")
    message('-------------------------------------------------------')
    message('-------------------       FSE       -------------------')
    message('-------------------------------------------------------')
    
    fm$status <- paste0(fm$status,'.SATs')
    fm <- rbind(data.frame(status = 'protofeatures', 
                           freq = length(regions_subset)),
                fm)
    fm <- setNames(data.frame(t(fm[,-1])), fm[,1] %>% str_replace(' ', '.'))
    fm$n.samples <- nrow(xmat)
    fm$n.points <- ncol(xmat)
    
    return(fm)
}
