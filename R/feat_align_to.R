#' Align features in "align" to a specified feature
#' Uses fft convolution. Concatenates features with padding and does a single fft calc.
#' Note: don't do this for huge feature numbers, or will cost in memory.
#'
#' @param align matrix of feature profiles (NAs allowed, but will be zero-filled)
#' @param to row index (of matrix) for feature to align to
#' @param max.hits number of top-ranked conv peaks to report
#' @param max.lag mag lag to consider
#'
#' @return lag table 
#'
#' @importFrom magrittr %>%
#' @importFrom fftw FFT
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr group_by
#' @importFrom dplyr slice
#' 
#' @export
feat_align_to <- function(align, to, max.hits = 1, max.lag = Inf){
  
  
  # align and to should be same size (NA padded if needed)
  
  # stackplot(align, vshift = 3)
  # align <- align[100,,drop = F]
  
  
  
  if (is.vector(align)){feature.stack <- matrix(align, nrow = 1)} else {feature.stack <- align}
  if (is.vector(to)){feat <- matrix(to, nrow = 1)} else {feat <- to}
  feat.nas <- is.na(feat)
  
  # Just use R's ccf function (slow): ####
    # align[is.na(align)] <- 0
    # feat[is.na(feat)] <- 0
    # r <- ccf(c(feat),
    #          c(align), lag.max = length(feat)-1)
    # r$acf %>% c %>% simplePlot
    # browser()
    # simplePlot(feat)
    # simplePlot(align)
  
    
  ############ Initial calculations  ####
     
  # No NAs allowed - zeros will work for now
  
    featureStack <- feature.stack
    featureStack[is.na(featureStack)] <- 0
    
    feat[feat.nas] <- 0
    
  # Need to pad the feature matrix so they comparisons can't overlap
  
    pad.size <- ncol(featureStack)
    padding <- rep(0,pad.size)

    featureStack <- apply(featureStack, 1, function(feat) c(feat, padding)) # matrix is transposed!
    
    feat <- c(feat, padding)
    
    padded.feat.width <- nrow(featureStack)

  # Build concatenated feature vector 
  #   The features are padded in the front, so add a pad at the start 
    allfeats <- c(padding, c(featureStack))
  
  # Pre-compute the Conj(FFT()) for the vector
    
    af.ft.c <- Conj(fftw::FFT(allfeats))
    # af.ft.c <- fftw::FFT(rev(allfeats))
    
    # nonzero <- allfeats != 0
        # rm(allfeats)
        

    fill.len <- length(af.ft.c) - padded.feat.width

  ############ Do the vectorized cross-correlations ####
        
            # Cross-correlate to find locations and scores: ####

              # Do the FFT-based convolution ####
                
                feat.ft <- feat %>%                 # take this feature
                  c(., rep(0, fill.len) ) %>%     # zero-pad the end of the feature to match the big vect length
                  fftw::FFT(.)          # compute its FT
    
                # r <- ccf(feat %>% c(., zeros(1, fill.len) ), 
                #          allfeats, lag.max = length(feat)-1)

                r <- (feat.ft*af.ft.c) %>%          # multiply the big vectors
                  fftw::FFT(.,inverse = TRUE) %>%   # inverse FT the result
                  Re %>% c %>%                      # take the real part and unlist
                  rev 
                
                # r <- r / ((length(allfeats)) * (length(feat) - fill.len))
                
                  # # Plot: 
                  #   scattermore::scattermoreplot(x = 1:length(r) - pad.size , y = r)
                  

                  rm(feat.ft)
                

              # Get convolution maxima (candidate lags) ####
                
                lags <- localMaxima(r)
                  # plot(r, x = 1:length(r), xlab = 'index in big vect')
                  

              # Pick the highest lag for each feature range ####
                 
                # Label each lag with the feature range it came from
                
                  shiftedLags <- lags - pad.size                        # Account for padding in front of each feature. 
                                                                        # Remember that each feature has padding in front of it, 
                                                                        # of length(biggest feature). Thus, if a conv peak is in the 
                                                                        # second half of a feature on our huge vector, it actually 
                                                                        # points to a match in the next feature. We can account 
                                                                        # for this by shifting all lags by one pad size. 
                    # plot(r, x = (1:length(r))-pad.size, xlab = 'shift by pad.size')
                  
                  f2.with.hit <- (lags / padded.feat.width) %>%         # bin each lag into a feature range
                                    ceiling                             # round up
                  
                # Get the lag within the matched feature
                
                  lag.f2 <- shiftedLags -                               # start with the lags in the big vector (minus front pad)
                              padded.feat.width*(f2.with.hit - 1)       # subtract the starting point for the matched feature
                    
                  
                  
                    rm(shiftedLags)
                    
                    
            
                
                # Now that each conv peak is annotated, compile all the hits into a data.frame
                  hits <- data.frame(f1 = 1,
                                     f2 = f2.with.hit, 
                                     lag.in.f2 = lag.f2,
                                     pos.big = lags,
                                     val = r[lags])
                  
                  # hits <- hits[hits$f1 != hits$f2,] # don't keep cases where feature matches itself
                  hits <- hits[hits$f2 <= ncol(featureStack), ]
                  hits <- hits[hits$lag.in.f2 <= max.lag, ]
                
                    rm(r)
                    rm(lag.f2)
                    rm(lags)
                    rm(f2.with.hit)
                  
                  
              # This allows us to select the top conv peak within each feature in the big vector. 
              #   The resulting lag.in.f2 for each is the lag that f1 needs relative to f2. 
              
                # hits <- hits %>% group_by(f2) %>%
                #   slice(which.max(val)) %>%           # take only the top conv peak in each f2 region
                #   as.data.frame
              
              # New way allows selection of top [max.hits] conv peaks per f2. 
                hits <- hits %>% 
                  arrange(desc(val)) %>%                       # sort all rows by val first
                  group_by(f2) %>%                             # group the hits by f2 index (i.e., what are the hits for each feature 2?)
                  slice(                                       # for each f2 slice...
                          1:min( 
                                  c(  length(val),             # if there are more hits for this f2 than [max.hits], limit to
                                      max.hits                 # the top [max.hits] rows (already sorted by val)
                                    ) 
                                )                         
                        ) %>% as.data.frame

              # browser()
              # # Plot the hit to show it actually works. (for debugging) ####
                # n <- 0
                # n <- n + 1
                # hit <- hits[n, ]
                # 
                #   fs <- feature.stack[c(hit$f1,hit$f2), ]
                #   pair.lag <- c(0,hit$lag.in.f2)
                # 
                #   indsmat <- outer((1:ncol(fs)), pair.lag, "-") %>% t
                #   indsmat <- indsmat - min(indsmat) + 1
                # 
                #   linds <- sub2indR(rows = c(1,2), cols = indsmat, m = 2)
                # 
                #   tmpmat <- matrix(NA, 2, max(indsmat))
                #   tmpmat[linds] <- fs
                # 
                #   tmpmat %>% trim_sides %>% simplePlot
            # Export hit table ####
             
        
    return(hits)
        
}