#' TINA Setup
#'
#' Given a list of peak regions and a matrix of spectra, creates a matrix to house
#' all the features and returns it in a list format along with additional 
#' information about the subsets and regions of interest.
#'
#' @param a A list of peak regions.
#' @param xmat A matrix of spectra.
#'
#' @return A list containing the feature stack, position stack, subsets, and 
#' region information.
#'
#' @examples
#' # Create a mock input list of peak regions and a matrix of spectra
#' peaks <- list(list(finalRegion = c(1, 2, 3), ref.vals = c(1, 2, 3), 
#'                    ref.idx = c(1, 2, 3), subset = c(1, 2)),
#'               list(finalRegion = c(4, 5, 6), ref.vals = c(4, 5, 6), 
#'                    ref.idx = c(4, 5, 6), subset = c(3, 4)))
#' spectra <- matrix(1:12, ncol = 4)
#' 
#' # Run the function
#' tina_setup(a = peaks, xmat = spectra)
#'
#' @export
tina_setup <- function(a, xmat){
  
      message("TINA Setup...")
    
      # Build a matrix to house all the features
        # Extract final region inds with regard to peak
          relinds <- lapply(1:length(a),
                            function(x) a[[x]]$finalRegion %in% a[[x]]$ref.idx %>% which)

          featureStack <- relinds %>% unlist %>% max %>% matrix(NA, nrow = length(a), ncol = .)
          positionStack <- featureStack
          subsets <- matrix(FALSE, nrow(featureStack), nrow(xmat))
          
        # Fill in the values
          for (i in 1:length(a)){
            featureStack[i,relinds[[i]]] <- a[[i]]$ref.vals
            positionStack[i,relinds[[i]]] <- a[[i]]$ref.idx
            subsets[i, a[[i]]$subset] <- TRUE
          }

        # Look at the subset overlap
                  # rowSums(subsets) %>% plot(xlab = "Feature (sorted)",
                  #                                    ylab = "Subset size")
                  # title(main = "Subset Size Distribution pre-afc")
                  # abline(h = 5)
          subset.diffs <- subsets %>% diff %>% abs %>% rowSums
          subset.overlaps <- (subsets[-1,] & subsets[-nrow(subsets),]) %>% rowSums
          subset.sizes <- subsets %>% rowSums
                  # (subset.diffs / subset.overlaps) %>% sort %>% plot(xlab = "Feature Comparison",
                  #                                                    ylab = "Subsets diff/overlap")
                  # subset.ratios <- (subset.diffs / subset.overlaps) 
                  # subset.ratios[subset.ratios < .5] %>% sort %>% plot(xlab = "Feature Comparison",
                  #                                                     ylab = "Subsets diff/overlap")
                  
          

          # If the comparison is missing > 0.5 of the spectra in either subset
            subset.frac.1 <- subset.overlaps / subset.sizes[-nrow(subsets)]
            subset.frac.2 <- subset.overlaps / subset.sizes[-1]
            subset.frac <- subset.frac.1
            nds <- subset.frac.1 > subset.frac.2
            subset.frac[nds] <- subset.frac.2[nds] 
            
            
            # subset.frac %>% sort %>% plot(xlab = "Feature Comparison",
            #                               ylab = "Feature Subset Overlap (lesser)")
            #   abline(h = 0.5, col = "red")
              
        # Look at the position region overlap (speed up later)
            reg.frac <- rep(0, nrow(positionStack)-1)
            reg.overlap <- reg.frac
            for (i in 1:length(reg.frac)){
              r1 <- positionStack[i,] %>% range(na.rm = TRUE) %>% fillbetween
              r2 <- positionStack[i+1,] %>% range(na.rm = TRUE) %>% fillbetween
              reg.overlap[i] <- sum(r1 %in% r2)
              reg.frac.1 <- reg.overlap[i] / length(r1)
              reg.frac.2 <- reg.overlap[i] / length(r2)
              reg.frac[i] <- min(c(reg.frac.1, reg.frac.2))
            }
            reg.sizes <- lapply(relinds, length) %>% unlist
            # reg.frac %>% sort %>% plot(xlab = "Feature Comparison",
            #                            ylab = "Feature Region Overlap (lesser)")
            #   abline(h = 0.5, col = "red")
                    
    return(list(stack = featureStack,
                position = positionStack,
                subset = list(ss.all = subsets,
                              diffs = subset.diffs,
                              sizes = subset.sizes,
                              overlaps = subset.overlaps,
                              fraction = subset.frac),
                region = list(overlaps = reg.overlap,
                              fraction = reg.frac,
                              sizes = reg.sizes)
                )
           )
}