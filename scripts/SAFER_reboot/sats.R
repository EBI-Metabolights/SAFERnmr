# Statistal Annotation Tag Extraction

# Follows:
# - setup
# - load_data
# - protofeatures

# 

#   protofeatures

# Parameter setup ####

    # Passing from previous functions ####
    
    xmat <- data$xmat
    ppm <- data$ppm
    tmpdir <- pars$dirs$temp
    
    # Override for now:
    
      only.region.between <- pars$corrpockets$only.region.between
      # only.region.between <- pars$corrpockets$only.region.between
        if (is.null(only.region.between))                       # which ppms to run fse between
          {only.region.between <- range(ppm)}                   #   (default is all)
      correlation.r.cutoff <- pars$storm$correlation.r.cutoff   # rvalue cutoff for both subset selection (ref shape) and ref update (STOCSY)
      q <- pars$storm$q                                         # q param from storm (pval cutoff after mhtc)
      b <- pars$storm$b                                         # number of peak widths to expand ref by on each side

    # Plotting ####
      number.of.plots <- pars$storm$number.of.plots             # pdf of all extracted features will be plotted. Choose
                                                                # only 150 of these (evenly spaced) or suffer the
                                                                # consequences...
                                                      
      plot.location <- paste0(tmpdir,"/plots/")                     # where to put the plot (just dump into run folder)
      dir.create(plot.location, showWarnings = FALSE)
  
   
# Run code ####

     
        bounds <- vectInds(only.region.between, ppm)

        storm_rnd1 = list()
        
        n.cores <- pars$par$ncores
        
        half.window = half.window <- (pars$corrpockets$half.window / data$digital.res) %>% ceiling
        
  # Set up multicore
    
    # split up the ppm vector into chunks
    # but first, randomize it for load balancing 
    
    pfs.in.region <- ((protofeatures$table$driver <= bounds[1]) & (protofeatures$table$driver >= bounds[2]) ) %>% protofeatures$table[.,]
        
    pfs.rand <- sample(seq_along(rownames(pfs.in.region)))
      unrand <- order(pfs.rand, decreasing = FALSE)
    chunk.size <- ceiling(length(pfs.rand) / n.cores)
    pf.chunk.assignments <- split(pfs.rand, ceiling(seq_along(pfs.rand) / chunk.size))
    
    pfs.split <- pfs.in.region %>% split(seq_along(pfs.rand))
      
    pf.chunks <- lapply(pf.chunk.assignments, function(x) pfs.split[x])

    results <- mclapply(pf.chunks, function(pfs){
      # pfs <- pf.chunks[[1]]
      # 
      statuses <- lapply(pfs, function(pf){
        # pf <- pfs[[143]]
        # message(pf$driver)
        s <- 
        tryCatch(
          expr = {
          
          log_storm_core(p = pf, xmat=xmat, ppm=ppm, half.window = half.window, corrthresh = .8,
                        q=0.05, minpeak = protofeatures$noiseWidth * protofeatures$noise.width.multiple, 
                        min.subset = 6,
                        plots=FALSE)
            # p = pf
            # half.window = 200
            # corrthresh = .8
            # q=0.05
            # minpeak = protofeatures$noiseWidth * protofeatures$noise.width.multiple
            # min.subset = 6
            # plots=TRUE
            
        },warning = function(w){
          # message('iteration ', which(lapply(pfs, function(x) x$driver) %>% unlist == pf$driver))
          return(list(protofeature = pf,
                      warning = w,
                      status = 'warning',
                      details = str_c('pf iteration: ', which(lapply(pfs, function(x) x$driver) %>% unlist == pf$driver))
                      )
                 )
        }, error = function(e){
          # message('iteration ', which(lapply(pfs, function(x) x$driver) %>% unlist == pf$driver))
          return(list(protofeature = pf,
                      warning = e,
                      status = 'error',
                      details = str_c('pf iteration: ', which(lapply(pfs, function(x) x$driver) %>% unlist == pf$driver))
                      )
                 )
        })
        
        s2 <- 
        tryCatch(
          expr = {

            i <- i + 10
            i
            
            p <- pfs[[i]]
            plot_protofeature(p,
                      half.window = half.window, ppm = data$ppm,
                      xmat = data$xmat,
                      bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
                      showPeaks = TRUE, show.mask.bounds = FALSE)
                        
            pexp <- expand_protofeature(p, xmat, ppm, half.window)
            driver <- pexp$driver
            
            wind <- pexp$specRegion.inds
            x <- wind
            specRegion = pexp$specRegion
            
            mean.spec <- colMeans(specRegion)
            
            fittedSpecs <- lapply(1:nrow(specRegion), function(m){
              # m <- m + 1
              fit <- fit_leastSquares(specRegion[m, ] %>% c, mean.spec, plots = FALSE)
              # fit$plot
            })
            
            
            simplePlot(y)

          log_storm_core(p = pf, xmat=xmat, ppm=ppm, half.window = half.window, corrthresh = .8,
                        q=0.05, minpeak = protofeatures$noiseWidth * protofeatures$noise.width.multiple, 
                        min.subset = 6,
                        plots=FALSE)
            # p = pf
            # half.window = 200
            # corrthresh = .8
            # q=0.05
            # minpeak = protofeatures$noiseWidth * protofeatures$noise.width.multiple
            # min.subset = 6
            # plots=TRUE
            
        },warning = function(w){
          # message('iteration ', which(lapply(pfs, function(x) x$driver) %>% unlist == pf$driver))
          return(list(protofeature = pf,
                      warning = w,
                      status = 'warning',
                      details = str_c('pf iteration: ', which(lapply(pfs, function(x) x$driver) %>% unlist == pf$driver))
                      )
                 )
        }, error = function(e){
          # message('iteration ', which(lapply(pfs, function(x) x$driver) %>% unlist == pf$driver))
          return(list(protofeature = pf,
                      warning = e,
                      status = 'error',
                      details = str_c('pf iteration: ', which(lapply(pfs, function(x) x$driver) %>% unlist == pf$driver))
                      )
                 )
        })
        
        return(s)

      })
      
      
    }, mc.cores = n.cores) %>% unlist(recursive = FALSE)
      
    ## Recombine into one list
    
    statuses <- results %>% lapply(function(x) x$status)
    
    
##  ####
        # message("Running storm on ",length(regions_subset), " provided protofeatures between ",ppm[bounds[1]]," and ",ppm[bounds[2]]," ppm.")
        # 
        # storm_rnd1 <- 
        #       mclapply(regions_subset,
        #       # pblapply(regions_subset,
        #           function (x) {
        #             # TryCatch will 
        #             tryCatch(
        #               expr = {
        #                       # Set up the region
        # 
        #                       # Do storm_pairplay
        #             
        #                           res$cpp.driver <- driver
        #                           
        #                           return(res)
        #                 
        #               },
        #               error = function(cond){
        #                 return('setup error')
        #                 }
        #               )
        #           }, mc.cores = pars$par$ncores
        #       )
        # # Note: errors in the loop are captured and passed out as strings.
        # # NULL elements are not possible, although parts of an element could be.
        # # Those are checked below.  

## Report ####
fmodes <- lapply(statuses, 
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
              (sum(failed)/length(pfs.split) * 100) %>% round, " %)"))

message(str_c("Succeeded iterations (count): ", sum(succeeded), " (",
              (sum(succeeded)/length(pfs.split) * 100) %>% round, " %)"))

# Print out breakdown of statuses
  fm <- fmodes %>% plyr::ldply(rbind)
  fm$.id <- NULL
  colnames(fm) <- 'status'
  fm <- fm %>% group_by(status) %>% count %>% as.data.frame
  
  print(fm)
  
  
## Plotting ####

  # return(s$status)
  # s$subset
  # s$finalRegion
  # s$ref.idx
  # s$ref.vals
  # s$covar
  
  simplePlot(xmat[,only.region.between %>% vectInds(ppm) %>% fillbetween])
  sat.list %>% lapply(function(x) x$peak) %>% unlist %>% sort %>% plot(y = 1:length(sat.list), x=.)
  
  sat.list <- results[succeeded]
  
  i <- i + 1
  s <- sat.list[[i]]
  
  plot_protofeature(p = data.frame(driver = s$peak),
            half.window = half.window, ppm = data$ppm,
            xmat = xmat[s$subset,],
            bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
            showPeaks = FALSE, ref.mask = s$ref.idx, show.mask.bounds = TRUE)

  
  i <- i + 10
  i
  plot_protofeature(pfs.in.region[i, ],
            half.window = half.window, ppm = data$ppm,
            xmat = data$xmat,
            bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
            showPeaks = TRUE, show.mask.bounds = FALSE)
  