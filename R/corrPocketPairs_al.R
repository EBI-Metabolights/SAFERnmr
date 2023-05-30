#' corrPocketPairs
#'
#' @param x a matrix with numerical values.
#' @param ppm a numeric vector of the same length as \code{ncol(x)} specifying the ppm for each column.
#' @param ws an integer specifying half the size of the sliding window used to calculate correlations.
#' @param reg an optional vector specifying the column indices to consider in \code{x}.
#' @param plotHeatmap a logical indicating whether to plot a heatmap of the correlation matrix.
#' @param wdlimit a numeric (0,1) specifying the minimum fraction of sliding windows that need to contain a pocket to consider a point noise.
#' @param rcutoff a numeric (0,1) specifying the correlation coefficient cutoff to use when extracting significant peaks.
#'
#' @return a list with the following elements:
#' \item{regions}{a matrix specifying the indices of the peaks in each column of the original matrix, filtered to exclude peaks that are too small or only have unidirectional interactions.}
#' \item{corr}{a matrix of correlation coefficients between columns of \code{x}.}
#' \item{cov}{a matrix of covariance values between columns of \code{x}.}
#' \item{peakMap}{a matrix specifying the indices of the peaks in each column of the correlation matrix.}
#' \item{noiseDist}{the proportion of pockets containing each window index.}
#'
#' @importFrom magrittr %>%
#'
#' @export
corrPocketPairs_al <-  function(x, ppm, ws, reg = NULL, plotHeatmap = FALSE, wdlimit = 0.99,
                             rcutoff = 0.75){
  
  if (is.null(reg)) {reg <- seq_along(ppm)}

##############################################################################################################     
  # Local correlation/covariance calculation
  
  message("Computing local correlations across columns of x...")
  res <- slidingCorr(x = x[,reg], 
                     ws = ws,
                     extractPockets = TRUE, 
                     plotting = FALSE, vshift = 10,
                     ppm = ppm[reg])
 
##############################################################################################################     
  # Peak Extraction
  
  message("Extracting peaks from local correlations...")
  cc <- res$corr_compact
  pts <- matrix(data = FALSE, nrow = nrow(cc), ncol = ncol(cc))
  pks <- pts
  iscenter <- pts
  pkID <- matrix(NA, nrow = nrow(cc), ncol = ncol(cc))
  
  
  # Use center peak for each column ####
    centers <- res$isPocket
    windowDist <- ( centers %>% rowSums(na.rm = TRUE) ) / ncol(cc)
    noiseWidth <- sum(windowDist > wdlimit)
    
  # We can exclude columns altogether which don't pass this threshold
    notNoise <- which(colSums(centers, na.rm = TRUE) >= noiseWidth)
  
  
  
  # Now get the highest non-center peak that passes the noiseWidth threshold
    for (i in notNoise){
        pkID[centers[,i],i] <- i
    }
    
    
    iscenter <- !is.na(pkID)
    
    # Next id will be one more than the max in pkID
      pid <- max(notNoise) + 1
      
    # Pull out the highest non-center peak that passes the noiseWidth threshold
      for (i in notNoise){
        # Pull vect and peaks
          peaks <- extractPeaks_corr(cc[,i], plots = FALSE)
          notcenter <- which(!(peaks$peaks %in% res$center))
          # bigEnough <- lapply(notcenter, function(x) peaks$bounds[[x]] %>% 
          #                       unlist %>% 
          #                       diff) %>% unlist >= noiseWidth
          
          bigEnough <- lapply(notcenter, function(x) peaks$bounds[[x]] %>% 
                                unlist %>% fillbetween %>% 
                                 cc[.,i] %>% ">" (.,rcutoff) %>% sum) %>% unlist >= noiseWidth
          
          # If no peaks worth extracting, then skip this column
            if (!any(bigEnough)){next}
          
          bestPeak <- which.max( cc[peaks$peaks[notcenter[bigEnough]],i] ) %>% notcenter[.]
          inds <- peaks$bounds[bestPeak] %>% unlist %>% fillbetween
          # pks[inds,i] <- TRUE
          pkID[inds,i] <- pid
          # iscenter[inds, i] <- TRUE
      }
      
      pks <- !is.na(pkID)
      
      
##############################################################################################################     
  # Filtering
  
  message("Filtering results (no peaks < size of noise; only bidirectional relationships)...")
  
  # At this point, we need to remove peaks that don't have a partner...
        # In other words, remove column if there is no center + noncenter peak
        nonCenterPks <- pks & (!iscenter)
        
        rmcols <- !apply(nonCenterPks, 2, any)
        pkID[,rmcols] <- NA
        pks[,rmcols] <- FALSE

  # Get the inds of the non-NA elements, convert to ppm inds
    # indsmat (from slidingCorr()) is just the column of reg.
      regions <- res$indsmat - 1 + min(reg)
      regions[!pks] <- NA
    
  # For each column, record the interactions as a pairs of spectral points between the center and its best peak
        # Interactions between spectral points means for each noncenter peak,
        # - the ppm index of each one (xmat column number) is derived from the regions inds
        # - the column number in cc gives the xmat column number of the center peak it belongs to
        # because the inds of cc correspond to the relative ind in regions, no need to translate
        # indices. 

        # Get column number in ppm (driver index) 
            
            colInPPM <- regions[nonCenterPks]
            
        # Use ind2sub to get the column number in cc (interaction with that )
            ccind <- nonCenterPks %>% which
            subs <- ind2subR(ccind, m = nrow(cc)) # get the columns
            colInCC <- subs$cols
        
        # Bind into pairs (associating ppms with ppms) and sort each one, convert into df
          
          prs <- sortPairs(rbind(colInCC, colInPPM)) %>% matrix(., ncol = 2) %>% data.frame
          prs <- cbind(prs, ccind) # keep track of the pos in cc explicitly.
                                   # if this gets a count of 1 later, we'll simply
                                   # remove that ind from the peak mask later. It 
                                   # may be that we want to remove the col altogether. 
                                   
            colnames(prs) <- c("pt1","pt2", "ccind")
            
        # Find the pairs which are found at least twice (necessary condition for bi-directionality of the correlation)
          uniquePairs <-prs %>%
            dplyr::group_by(pt1,pt2) %>%
            dplyr::mutate(Count = n()) %>%
            dplyr::ungroup() %>%
            dplyr::distinct() %>%
            dplyr::filter(., Count < 2)

        # Remove peaks with only majority unidirectional associations
        
          # cp <- matrix(0, nrow(cc), ncol(cc))
          # cp[nonCenterPks] <- 2
          # cp[uniquePairs$ccind] <- 5
          
          # heatmap(cp, Colv = NA, Rowv = NA, scale="none")
          
          # For each column of the filter, determine the ratio of unique points to bidirectional points
            f <- matrix(NA, nrow(cc), ncol(cc))
              f[nonCenterPks] <- 1
              f[uniquePairs$ccind] <- 0
            bidirs <- f %>% t %>% rowSums(na.rm = TRUE)
            allassoc <- nonCenterPks %>% t %>% rowSums(na.rm = TRUE)
            
            ratios <- bidirs/allassoc
            
         # Filter out columns whose responder peaks are mostly unidirectional points.
         # This effectively symmetrizes the matrix without transforming it/modifying data.  
         
            fpks <- pks
            fpks[,ratios < 0.5] <- FALSE 
        
       
            # fpks <- pks
            # fpks[uniquePairs$ccind] <- FALSE
            
            pkID[!fpks] <- NA
      
      # Plot result
            if (plotHeatmap){
              message("Plotting heatmap...")
              cp <- cc
              cp[is.na(pkID)] <- NA
              heatmap(cp, Colv = NA, Rowv = NA, scale="none")
            }
          
       message("corrPocketPairs() finished.")
        
    return(list(regions = regions,     # these are unfiltered
                corr = cc,
                cov = res$cov_compact,
                peakMap = pkID,        # %>% is.na can be used as a filter for all other outputs. 
                noiseDist = windowDist)) # % of pockets containing each windowInd
}