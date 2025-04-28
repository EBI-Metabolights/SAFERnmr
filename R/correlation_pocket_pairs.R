#' correlation_pocket_pairs
#'
#' @param x a matrix with numerical values.
#' @param ppm a numeric vector of the same length as \code{ncol(x)} specifying the ppm for each column.
#' @param ws an integer specifying half the size of the sliding window used to calculate correlations.
#' @param reg an optional vector specifying the column indices to consider in \code{x}.
#' @param plotHeatmap a logical indicating whether to plot a heatmap of the correlation matrix.
#' @param wdlimit a numeric (0,1) specifying the minimum fraction of sliding windows that need to contain a pocket to consider a point noise. Like asking "what 1-the minimum fraction of points you think might be noise?". If at least 5% of the points are noise, choose 0.95 
#' @param rcutoff a numeric (0,1) specifying the correlation coefficient cutoff to use when extracting significant peaks.
#'
#' @return a list with the following elements:
#' \item{regions}{a matrix specifying the indices of the peaks in each column of the original matrix, filtered to exclude peaks that are too small or only have unidirectional interactions.}
#' \item{corr}{a matrix of correlation coefficients between columns of \code{x}.}
#' \item{cov}{a matrix of covariance values between columns of \code{x}.}
#' \item{peakMap}{a matrix specifying the indices of the peaks in each column of the correlation matrix.}
#' \item{noiseDist}{the proportion of pockets containing each window index.}
#'
#' @import magrittr
#'
#' @export
correlation_pocket_pairs <-  function(x, ppm, ws, reg = NULL, plotHeatmap = FALSE, wdlimit = 0.99,
                                      noise.width.multiple = 2, top.n.peaks = 5, rcutoff = 0.5, n.cores = 10){
  
  # assuming that the matrix is the full matrix, and ppm inds are the columns
  # x <- xmat
  # ppm <- ppm
  # ws <- 100
  # reg <- NULL # indices of ppm
  # plotHeatmap <- FALSE
  # wdlimit <- 0.95
  # plotHeatmap <- FALSE
  # rcutoff <- 0.5
  # pars$corrpockets$only.region.between <- c(-0.5,10.5)
  # noise.width.multiple = 2
  # top.n.peaks = 5
  # include upper and lower bounds for peak width?
  
  if (is.null(reg)) {reg <- seq_along(ppm)}

##############################################################################################################     
  # Local correlation/covariance calculation
  
  message("Computing local correlations across columns of x...")
  res <- slidingCorr(x = x[,reg], 
                     ws = ws,
                     extractPockets = TRUE, 
                     plotting = FALSE, vshift = 10,
                     ppm = ppm[reg], n.cores = 5) # ppm only used for plotting
 
##############################################################################################################     
  # Peak Extraction
  
  message("Extracting peaks from local correlations...")
  cc <- res$corr_compact
  cv <- res$cov_compact
  colnames(cc) <- reg
  colnames(res$cov_compact) <- reg
  cc[cc<=0] <- 0
  
  # Use center peak for each column ####
    centers <- res$isPocket
    windowDist <- ( centers %>% rowSums(na.rm = TRUE) ) / ncol(cc)
    noiseWidth <- sum(windowDist > wdlimit)
    
  # We can exclude columns altogether which don't pass this threshold
    notNoise <- which(colSums(centers, na.rm = TRUE) >= noise.width.multiple*noiseWidth)
    # these are the inds of cc columns that pass
    
    # Pull out the n highest non-center peaks that pass the noiseWidth threshold ####
    
    # Only want the correlations that have >= noiseWidth correlation
    cc.split <- lapply(notNoise, function(i) list(corrs = cc[,i],
                                                  covar = cv[,i],
                                                  col = i))
    res.center <- res$center

    cc.peaks <- parallel::mclapply(cc.split, function(col.info){
    # cc.peaks <- lapply(cc.split, function(col.info){
        
        # j <- lapply(cc.split, function(ci) {ci$col == 18501}) %>% unlist %>% which
        # col.info <- cc.split[[1765]]
        cc.col <- col.info$corrs
        # cv.col <- col.info$covar
        # NOTE: Everything in here is in window indices
        
        
        # It might be better to get covariance peaks whose maxima also have correlation > r.thresh
        # Or the covariance peak which contains the highest corr
        # Really, we don't know if we can trust the covariance to provide good peakshape
        # A sufficiently wide correlation peak is an indication that there is something underneath in at least a subset of the data
        # but the shape of the correlation profile does not preserve the relative intensities of the peaks.
        # What if the correlation peak shapes were used, but were scaled by their respective covariance profiles (batman fit or mean-match)
        
        # Pull vect and peaks
          peaks <- extractPeaks_corr(cc.col, plots = FALSE)
          primary.peak <- (peaks$peaks %in% res.center)
          secondary.peaks <- which(!primary.peak)
          useful.points <- cc.col > rcutoff
          
          bigEnough <- lapply(secondary.peaks, function(x) peaks$bounds[[x]] %>% 
                                unlist %>% fillbetween %>% useful.points[.] %>% sum) %>% unlist >= noiseWidth*noise.width.multiple
          
          # If no peaks worth extracting, then skip this column
            if (!any(bigEnough)){return(NULL)}
            # browser()
          # Which secondary peaks are wide enough?
            pk.idxs <- secondary.peaks[bigEnough]
            pk.locs.cc.col <- peaks$peaks[pk.idxs]
            
          # Pick the tallest n of those
            pk.maxima <- cc.col[pk.locs.cc.col]
            pk.ranks <- order(pk.maxima, decreasing = TRUE)
            bestPeaks.idx <- pk.ranks[1:min(top.n.peaks, length(pk.ranks))] %>% secondary.peaks[.]
          # bounds <- peaks$bounds[bestPeaks.idx] %>% do.call(cbind,.)
          
          # Package up:
          result <- list(primary = peaks$bounds[primary.peak]%>%unlist,
                      secondary = bestPeaks.idx %>% secondary.peaks[.] %>% pk.idxs[.] %>% peaks$bounds[.],
                      index = col.info$col,
                      res.center = res.center)
        
          # If any secondary peaks were null, remove them. Not sure why this happens yet.
            result$secondary <- result$secondary[!is.null(result$secondary)] # need to follow up on these cases!
            
          # Development/Debugging:
            # i <- 0
            # 
            # i <- i + 1
            # 
            # p <- result
            # driver <- p$index
            # # peak.inds <- c(p$primary.lower:p$primary.upper, p$secondary.lower:p$secondary.upper)
            # primary.peak.inds <- c(p$primary %>% fillbetween)
            # secondary.peak.inds <- c(p$secondary[[i]] %>% unlist %>% fillbetween)
            # 
            # shape <- cv.col
            # plot(x = 1:length(shape), y = shape, type = 'l')
            #   lines(x = primary.peak.inds, shape[primary.peak.inds], col='blue', lwd=2)
            #   lines(x = secondary.peak.inds, shape[secondary.peak.inds], col='blue', lwd=2)
            #   abline(v=p$secondary[[i]]%>%unlist)
            # 
            # plot_protofeature(driver, ws, ppm, xmat, res, bgplot='overlayed')
            
            ####

            # shape <- cv.col
            # plot(x = 1:length(shape), y = shape, type = 'l')
            #   lines(x = primary.peak.inds, shape[primary.peak.inds], col='blue', lwd=2)
            #   lines(x = secondary.peak.inds, shape[secondary.peak.inds], col='blue', lwd=2)
            #   abline(v=p$secondary[[i]]%>%unlist)

            ####
            
            # shape <- rep(NA, length(cc.col))
            #   
            # # Scale primary peak
            # 
            #   primary.peak.scaled <- fit_batman(feat = cc.col[primary.peak.inds],
            #              spec = cv.col[primary.peak.inds],
            #              exclude.lowest = 0.5, plots = FALSE) #TRUE
            #   shape[primary.peak.inds] <- primary.peak.scaled$feat.fit
            #   
            # # Scale secondary peak
            # 
            #   secondary.peak.scaled <- fit_batman(feat = cc.col[secondary.peak.inds],
            #              spec = cv.col[secondary.peak.inds],
            #              exclude.lowest = 0.5, plots = FALSE) #TRUE
            #   shape[secondary.peak.inds] <- secondary.peak.scaled$feat.fit
            #   
              
            ## ---------
            # Track down the cases where secondary appears
            # if (any(is.null(result$secondary))){
            #   browser()
            # }
            # 
          return(result) # index = index in "reg", which indexes cc.corr
                # This way of storing the peaks allows us to use:
                # cc.peak <- cc.peaks[[1]]
                # bounds <- cc.peak$index - (res.center-cc.peak$primary)
                # but remember - these are mainly to index the corr and cov mats.
# })
    }, mc.cores = n.cores)
    # }, mc.cores = pars$par$ncores)
    
  
  # non.null <- cc.peaks %>% lapply(function(x) !is.null(x)) %>% unlist %>% which
  # cc.peaks <- cc.peaks[non.null]
  # 
  # cc.peaks[[1]]
    
##############################################################################################################     
  # Filtering
  
  # message("Filtering results (no peaks < size of noise; only bidirectional relationships)...")
  
  # # At this point, we need to remove peaks that don't have a partner...
  # # This enforces complete connectivity in a correlation cluster. This 
  # # may not make sense for 1' and 2' correlation peaks, however, because of the
  # # 'love triangle' that could exist within a multiplet.

  # Get the inds of the non-NA elements, convert to ppm inds
    # indsmat (from slidingCorr()) is just the column of reg.
        
      # regions = lapply(1:ncol(res$indsmat), function(j) (res$indsmat[,j]- 1 + min(reg)) %>% range(na.rm = TRUE) )
      # regions <- res$indsmat - 1 + min(reg)
      # regions[!pks] <- NA
      
      # *** You need: the peak bounds (assume first peak is primary) and the inds for the window.
      
       message("corrPocketPairs() finished.")
        
    return(list(corr = cc,
                cov = res$cov_compact,
                peakBounds = cc.peaks,
                noiseDist = windowDist,
                noiseWidth = noiseWidth,
                noise.width.multiple = noise.width.multiple)) # % of pockets containing each windowInd
}
