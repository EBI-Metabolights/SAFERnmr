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
      # only.region.between <- c(4.5,5)
      # only.region.between <- c(2,3)
      only.region.between <- c(0,10)
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
        
        half.window <- (pars$corrpockets$half.window / data$digital.res) %>% ceiling
        
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
            # p = pf
            half.window = 200
            corrthresh = .9
            q=0.05
            minpeak = protofeatures$noiseWidth * protofeatures$noise.width.multiple
            min.subset = 6
            plots=FALSE
    
    results.all.cores <- mclapply(pf.chunks, function(pfs){
      # pfs <- pf.chunks[[1]]
      # 
      results.core <- lapply(pfs, function(pf){
        # pf <- pfs[[143]]
        # pf <- pfs[[1]]
        # message(pf$driver)
        s <- 
        tryCatch(
          expr = {

          log_storm_core(p = pf, data = data, half.window = half.window, corrthresh = corrthresh,
                        q=q, minpeak = protofeatures$noiseWidth * protofeatures$noise.width.multiple, 
                        min.subset = min.subset,
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
        # ####
        # s2 <-
        # tryCatch(
        #   expr = {
        # 
        #     # i <- i + 10
        #     i<- 70
        # 
        #     p <- pfs[[i]]
        #     plot_protofeature(p,
        #               half.window = half.window, ppm = data$ppm,
        #               xmat = data$xmat,
        #               bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
        #               showPeaks = TRUE, show.mask.bounds = FALSE)
        # 
        #     pexp <- expand_protofeature(p, xmat, ppm, half.window)
        #     driver <- pexp$driver
        # 
        #     wind <- pexp$specRegion.inds
        #     x <- wind
        #     specRegion = pexp$specRegion
        # 
        #     simplePlot(xmat[,wind], ppm[wind])
        #     simplePlot(specRegion, ppm[wind])
        # 
        # 
        #   s <- log_storm_core(p = p, data=data, half.window = half.window, corrthresh = .8,
        #                 q=0.05, minpeak = protofeatures$noiseWidth * protofeatures$noise.width.multiple,
        #                 min.subset = 6,
        #                 plots=FALSE, local.fits = fits)
        # 
        #     # p = pf
        #     # half.window = 200
        #     # corrthresh = .8
        #     # q=0.05
        #     # minpeak = protofeatures$noiseWidth * protofeatures$noise.width.multiple
        #     # min.subset = 6
        #     # plots=TRUE
        #     plot_protofeature(p = data.frame(driver = s$peak),
        #               half.window = half.window, ppm = data$ppm,
        #               xmat = xmat[s$subset,],
        #               bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
        #               showPeaks = FALSE, ref.mask = s$ref.idx, show.mask.bounds = TRUE)
        # 
        # },warning = function(w){
        #   # message('iteration ', which(lapply(pfs, function(x) x$driver) %>% unlist == pf$driver))
        #   return(list(protofeature = pf,
        #               warning = w,
        #               status = 'warning',
        #               details = str_c('pf iteration: ', which(lapply(pfs, function(x) x$driver) %>% unlist == pf$driver))
        #               )
        #          )
        # }, error = function(e){
        #   # message('iteration ', which(lapply(pfs, function(x) x$driver) %>% unlist == pf$driver))
        #   return(list(protofeature = pf,
        #               warning = e,
        #               status = 'error',
        #               details = str_c('pf iteration: ', which(lapply(pfs, function(x) x$driver) %>% unlist == pf$driver))
        #               )
        #          )
        # })
        
        return(s)

      })
      
      return(results.core)
      
    }, mc.cores = n.cores) %>% unlist(recursive = FALSE)
      
    ## Recombine into one list
    # results.all.cores <- results.core
    statuses <- results.all.cores %>% lapply(function(x) x$status)
    
    
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
  
  sat.list <- results.all.cores[succeeded]
  
  sat.list <- lapply(1:length(sat.list), function(x){
    s <- sat.list[[x]]
    s$id <- x
    return(s)})
  
  # simplePlot(xmat[,only.region.between %>% vectInds(ppm) %>% fillbetween])
  # stackplot(xmat[,only.region.between %>% vectInds(ppm) %>% fillbetween])
  sat.list %>% lapply(function(x) x$peak) %>% unlist %>% sort %>% plot(y = 1:length(sat.list), x=.)
  sortOrder <- sat.list %>% lapply(function(x) x$peak) %>% unlist %>% order()
  sat.list <- sat.list[sortOrder]
  
  i <- 0
  i <- i + 100
  s <- sat.list[[i]]
  
  plot_protofeature(p = data.frame(driver = s$peak),
            half.window = half.window, ppm = data$ppm,
            xmat = xmat[s$subset,],
            # bgplot = 'stack', line.shape = 'covar', line.color = "corr",
            bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
            showPeaks = FALSE, ref.mask = s$ref.idx, show.mask.bounds = TRUE)

  # i <- i + 100
  # i
  # plot_protofeature(pfs.in.region[i, ],
  #           half.window = half.window, ppm = data$ppm,
  #           xmat = data$xmat,
  #           bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
  #           showPeaks = TRUE, show.mask.bounds = FALSE)
  # 
  # plot_sat(p, half.window, ppm, xmat, bgplot='overlayed', line.shape='covar', line.color='corr', showPeaks=TRUE, ref.mask = NULL, show.mask.bounds=FALSE)
  
  # Imagine that you scroll across, seeing the feature shapes that correspond to the spectral point you're on.
  # When you find the shape, you click or press enter to trigger matching, etc. for it. 
  
plot.sats.grid <- function(sat.list, xmat, ppm, selected, title.strs=NA, include.spectra = F){
  # selected <- seq(from=1, to=length(sat.list), by = 10)
  sat.list <- sat.list[selected]
  sat.list <- lapply(1:length(sat.list), function(s){
    
    this.title <- title.strs[s]
    
    if (is.na(this.title)){
      this.title <- ''
    }
    
    sat.list[[s]]$title <- this.title
    sat.list[[s]]
  })
  
  # plots <- pbapply::pblapply(selected, function(sat.index){
  plots <- mclapply(sat.list, function(s){
    # sat.index <- selected[1]
    # s <- sat.list[[sat.index]]
    
    pexp <- expand_protofeature(p = data.frame(driver = s$peak), 
                                xmat[s$subset,], data$ppm, half.window)
      
    cv <- pexp$cv
    ppms <- pexp$ppmRegion
    range(ppms)
    ref.mask <- s$ref.idx
    ref.mask.region <- pexp$specRegion.inds %in% ref.mask
    
    g1 <- simplePlot(cv, xvect=ppms,linecolor = 'gray')
    
    cv[!ref.mask.region] <- NA
    ppms[!ref.mask.region] <- NA
    
    cv.fit <- fit_leastSquares(cv, colMeans(pexp$specRegion), plots = TRUE, scale.v2 = FALSE)
      # cv.fit$plot

    if (include.spectra){
      colors.lines <- c(rep("gray", nrow(pexp$specRegion)), 'blue')
      g1 <- simplePlot(rbind(pexp$specRegion,
                       cv.fit$feat.fit),
                 pexp$ppmRegion,
                 linecolor = colors.lines)
    } else {
      colors.lines <- c('gray', 'blue')
      g1 <- simplePlot(rbind(colMeans(pexp$specRegion), 
                       cv.fit$feat.fit), 
                 pexp$ppmRegion, 
                 linecolor = colors.lines)
    }
    
    
    
    # g2 <- simplePlot(cv, xvect=ppms,linecolor = 'red')
    
    # g <- g1+g2
    
    g1 + ggtitle(s$title)
    
    
    
    
    
    
    
    
    # g1 <- simplePlot(covar.filtered, xvect = ppm.vals.filtered, n_xticks = 4)
    # Add the peak bounds
      # g1 <- g1 + 
      #   geom_vline(xintercept = ppm[pexp$primary.bounds], linetype = 2, col = "black") +
      #   geom_vline(xintercept = ppm[pexp$secondary.bounds], linetype = 2, col = "black")
      # g1 <- 
      
    # plot_protofeature(p = data.frame(driver = s$peak),
    #           half.window = half.window, ppm = data$ppm,
    #           xmat = xmat[s$subset,],
    #           # bgplot = 'stack', line.shape = 'covar', line.color = "corr",
    #           bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
    #           showPeaks = FALSE, ref.mask = s$ref.idx, show.mask.bounds = TRUE)
    
  }, mc.cores = min(pars$par$ncores, length(selected))) 
  
  plots %>% grid_pdf(plotLoc=tmpdir, filename="/sats.pdf")
  
  # UX Idea:
    # Look at grid plot
    
    # When feature is selected, display its plot_protofeature(overlay)
    # if switch is flipped, display its plot_protofeature(stackplot)                                                                                     
  grid_pdf <- function(plots=NULL, plotLoc="./", filename="grid_plot.pdf"){
  # How big to make the page? 2 inches for each plot, and grid will be square.
    dim <- 3*round(sqrt(length(plots)))
    pdf(file = str_c(plotLoc,filename),   # The directory you want to save the file in
        width = dim, # The width of the plot in inches
        height = dim)
    
    gridExtra::grid.arrange(grobs = plots)
    
    dev.off()  
  }
}
  
  
  