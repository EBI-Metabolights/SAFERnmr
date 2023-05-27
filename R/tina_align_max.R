#' Aligns the features in the input matrix based on their maximum peaks
#'
#' @param feature a list object containing a matrix of feature values and a matrix of their corresponding positions
#' @param scaling a logical value indicating whether to scale the feature values to between 0 and 1 based on the highest peak's intensity
#'
#' @return a list object containing
#'
#' @export
align.max <- function(feature, scaling = TRUE){
  
  # Get the lags (-index of max peak not including bounds)
    # max.inds <- feature$stack %>% apply(1, which.max)
    
    peaks <- lapply(1:nrow(feature$stack), function(x) {
      
      # Extract the peaks from the feature:
        thisrow <- feature$stack[x,]
        pks <- extractPeaks_corr(thisrow)
        
      # Take max true peak index (not detected on boundaries):
        passpks <- pks$peaks[pks$truePeak]
        maxpk <- thisrow[passpks] %>% c(-Inf,.) %>% which.max - 1 # this is the index in passpks
          # c(-Inf,.) => allows empties in which.max with default = 1. 
          # -1 => if no true peaks, returns 0 (i.e. no lag; desirable behavior)
        
      # Return peak ind (or 0 if no peak found)
        return(list(maxpk = max(c(0,passpks[maxpk])),
                    pks = pks,
                    pkind.max = max(c(0,maxpk))))
          # likewise, default is 0
          
    })
    
    # # Filter out features without true max peaks (left and right bounds, true local max, see extractPeaks_corr)
      max.inds <- lapply(peaks, function(x) x$maxpk) %>% unlist

    # Subset everything with a true max peak
      # keep <- !(max.inds == 0) 
      peaks <- peaks#[keep]
      featureStack <- feature$stack#[keep,]
      positionStack <- feature$position#[keep,]
      driver.relative <- feature$driver.relative#[keep]
      max.inds <- max.inds#[keep]
      lags <- -max.inds

    # Scale center peak between 1 and 0 to 
      
      if (scaling){
          scaledFeatures <- lapply(1:nrow(featureStack), function(x) {
               
            # Extract the peaks from the feature:
                v <- featureStack[x, ,drop = FALSE]
                
              # Focus on the highest (soon-to-be center) peak
                pks <- peaks[[x]]$pks
                maxpk <- peaks[[x]]$pkind.max
                bounds <- pks$truePeak %>% which %>% .[maxpk] %>% pks$bounds[[.]]
              
              # Use the intensity of the bound with the lower valley 
                smallerBound <- bounds %>% unlist %>% v[.] %>% min
                
              # Use the peak max intensity
                pkmax <- peaks[[x]]$pks$peaks[maxpk] %>% v[.]
                
              # Scale so peak is between 0 and 1
                v <- (v - smallerBound)/(pkmax-smallerBound) # 
                      
            # Return scaled feature (or 0 if no peak found)
                return(v)
        })
          featureStack <- do.call(rbind, scaledFeatures)
      }
    
  # Align the features to the max peak  
    
    # Get the inds in the widened matrix that will hold the max-aligned features
      s.inds <- outer(1:ncol(featureStack), lags, "+") %>% t
      s.inds <- s.inds - min(s.inds) + 1
    
    # Build a new feature and position matrix with aligned features
      fs.ma <- ps.ma <- matrix(NA, nrow = nrow(featureStack), ncol = s.inds %>% span)

      for (i in 1:nrow(featureStack)){
        fs.ma[i, s.inds[i,]] <- featureStack[i,]
        ps.ma[i, s.inds[i,]] <- positionStack[i,]
      }
      
      # feature$stack %>% heatmap(Rowv = NA, Colv = NA, scale = "row")
      # fs.ma %>% heatmap(Rowv = NA, Colv = NA, scale = "row")
     
      driver.relative <- lapply(1:length(driver.relative), function(i){
        # i <- 1
        driver <- positionStack[i, driver.relative[i]]
        d.r <- which(ps.ma[i, ] %in% driver)
        if (length(d.r) == 1){return(d.r)} else {return(NA)}
      }) %>% unlist
      
    # Update feature object for new inds
      feature$stack <- fs.ma
      feature$position <- ps.ma
      feature$subset$ss.all <- feature$subset$ss.all#[keep,]
      feature$subset$sizes <- feature$subset$sizes#[keep]
      feature$region$sizes <- feature$region$sizes#[keep]
      # feature$driver.relative <- lapply(1:sum(keep),
      #                                    function(x) which.max( feature$stack[x, ] )
      #                                    ) %>% unlist
      feature$driver.relative <- driver.relative
      

    # These are no longer relevant; remove (only for adjacent feature comparison)
      feature$subset$diffs <- NULL
      feature$subset$overlaps <- NULL
      feature$subset$fraction <- NULL
      feature$region$overlaps <- NULL
      feature$region$fraction <- NULL
  
  return(feature)
}