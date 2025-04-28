# Derive protofeatures for a spectral matrix

# Follows:
# - setup
# - load data
# For surviving protofeatures:
# data.frame(primary.lower,
#            primary.upper,
#            secondary.lower,
#            secondary.upper,   # relative driver position
#            res.center,        # relative driver position
#            index = x$index,
#            row.names = NULL)
# Viewing protofeatures:
  # p <- protofeatures.split[[i]] 
  ## or
  i <- i + 10
  p <- protofeatures[i, ]
  plot_protofeature(p, 301, ppm, xmat, bgplot='overlayed')


protofeatures <- function(pars, xmat, noise.width.multiple = 2, top.n.peaks = 5, n.cores = 6){
  
  # Get correlation pockets (protofeatures)
  
  ################ Set up parameters ##################
    
    plot.location <- pars$dirs$temp
  
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
         
################ Use corrPocketPairs to extract likely j-pairs ##################
 
    # Run corrpocketPairs on everything
      
      pocketPairs <- correlation_pocket_pairs(xmat, ppm, ws = half.window, plotHeatmap = FALSE,
                                        wdlimit = noise.percentile, # **** **** **** #
                                        noise.width.multiple = 2,
                                        top.n.peaks = 5,
                                        rcutoff = cpp.rcutoff, n.cores = 6)
      
      pocketPairs %>% debug_write("pocketPairs.RDS", pars)
      # pocketPairs <- readRDS(paste0(pars$dirs$temp, "/debug_extra.outputs", "/pocketPairs.RDS"))

      
    # Clean up the result
      # Remove nulls
      pocketPairs$peakBounds <- pocketPairs$peakBounds[!is.null(pocketPairs$peakBounds)]
      
      
    # Unlist into individual pairs
    
      protofeatures <- mclapply(pocketPairs$peakBounds, function(x) {
        
          null.ones <- is.null(x$secondary)
          
          lapply(x$secondary[!null.ones], function(y) {
            
            if (!is.null(y)){
              data.frame(primary.lower = x$primary['lower'],
                         primary.upper = x$primary['upper'],
                         secondary.lower = y[['lower']],
                         secondary.upper = y[['upper']],
                         res.center = x$res.center, 
                         index = x$index,
                         row.names = NULL)
            }
            else {
              NULL
            }
          }) %>% do.call(rbind,.)

      }, mc.cores = 6) %>% do.call(rbind,.)
    
    # Convert to relative inds
      drivers <- protofeatures$index
      protofeatures <- protofeatures$res.center - protofeatures
      protofeatures$driver <- drivers
      protofeatures$index <- NULL
      protofeatures$res.center <- NULL
      
    # Expand protofeature
    
          # i <- 0
          # 
          # i <- i + 1
          # p <- protofeatures[i,]
          # i
          # plot_protofeature(p, 200, ppm, xmat, bgplot='stack', line.shape = 'covar', line.color = 'corr')
          # 
          # debug_write(protofeatures, 'protofeatures.RDS', pars)
          # 
    # Report number of pairs
      
      numPairs <- nrow(protofeatures)
      window.index <- (-half.window):(half.window)
      
      pdf(file = paste0(plot.location, "corrpeak_distribution.pdf"),   # The directory you want to save the file in
          width = 4, # The width of the plot in inches
          height = 4) # The height of the plot in inches
      
        pocketPairs$noiseDist %>% plot(x = window.index, ylab="Fraction of peaks including index", xlab="Window index")
        noisewidth <- sum(pocketPairs$noiseDist >= noise.percentile)
        abline(h = noise.percentile, col="red")
        title(ylab = "", main = "Average Extracted Diagonal Peak Shape (pre-filtering)")
      
      dev.off()
      # here's a thought: if you filter all peaks based on the noise feature shape,
      # the relative prominence of true signal using those boundaries is going to be
      # minimal because signal is locally pretty flat, while noise will mostly be
      # captured within those bounds. Peaks could be classified based on the % of their
      # actual signal (using their actual bounds) captured by the n% cutoff bounds.

      # ####
    
  return(list(protofeatures = protofeatures,
              noiseWidth = noiseWidth,
              noise.width.multiple = noise.width.multiple))
}