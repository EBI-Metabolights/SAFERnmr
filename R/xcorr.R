#' Calculate cross-correlation between reference and shifted spectra
#'
#' Wrote my own cross-correlation function for more consistent behavior in the
#' context of comparing ref signatures to spec data.
#' This takes a (wider) spec vector and slides it across the ref signature
#' lags. The lags are used to extract a window of spec for each shift (lag).
#' These windows are collected as rows of a matrix, which is then used to compute
#' all correlations to the ref. MTJ 2023
#'
#' @param spec Numeric matrix. A matrix of shifted spectra.
#' @param ref.vals Numeric vector. The values of the reference spectrum.
#' @param ref.idx Integer vector. The indices of the reference spectrum.
#' @param lags Integer vector. The range of lags to consider.
#' @param min.overlap Integer. The minimum number of overlapping data points between reference
#'   and shifted spectra to calculate correlation.
#'
#' @return A list containing the following elements:
#'   \item{rvals}{Numeric vector. The correlation scores for each shifted spectrum.}
#'   \item{pvals}{Numeric vector. The p-values for each correlation score.}
#'   \item{overlaps}{Numeric vector. The number of overlapping data points between reference
#'     and each shifted spectrum.}
#'   \item{not.used}{Logical vector. Indicates which shifted spectra did not meet the minimum
#'     overlap requirement and were not used in calculating correlation scores or p-values.}
#'   \item{lags}{Integer vector. The lags used for each shifted spectrum.}
#'   \item{bestFit}{List. Contains information on the best-matching shifted spectrum. The list
#'     has the following elements:
#'     \item{which.lag}{Logical vector. Indicates which lags were the best match.}
#'     \item{r}{Numeric. The correlation score for the best match.}
#'     \item{lag}{Integer. The lag corresponding to the best match.}
#'     \item{pval}{Numeric. The p-value for the correlation score of the best match.}
#'     \item{overlaps}{Integer. The number of overlapping data points between the reference and
#'       the best-matching shifted spectrum.}
#'     \item{specInds.matched}{Integer vector. The indices of the data points in the shifted
#'       spectrum that correspond to the overlapping region with the reference spectrum. This can
#'       be used to plot the overlapping regions.}}
#'
#' @export
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' ref <- rnorm(10)
#' ref.idx <- 1:10
#' spec <- matrix(rnorm(1000), nrow = 10)
#' lags <- -5:5
#'
#' # Calculate cross-correlation
#' xcorr(spec, ref, ref.idx, lags)
xcorr <- function(spec, ref.vals, ref.idx, lags, min.overlap = 3){
  
  # Build a matrix where each row is the spec data extraction resulting from a 
  # different lag with respect to ref.
    
    indmat <- outer(ref.idx %>% range %>% fillbetween, lags, FUN = "+")
    
    spec.allshifts <- (spec[indmat] %>% pracma::Reshape(., nrow(indmat), ncol(indmat)))
   
  # Format ref so missing vals (indicated by gaps in ref.idx) are NAs:
    ref <- rep(NA, ref.idx %>% span)
    ref[ref.idx-min(ref.idx)+1] <- ref.vals
    
    
  # (optional) Don't calc corr is not enough overlap
      # Determine length of overlap
    
        # calculate # of mutually non-na elements between ref and each row of spec.allshifts
        # More explicit:
        # overlaps <- apply(spec.allshifts, 2, function(shiftedSpec) ((ref %>% is.na) | 
        #                                                             (shiftedSpec %>% is.na)) %>% "!"(.) %>% sum)
        # 
        # Same speed but simpler (addition with an NA defaults to NA):
          
          # Non-NAs overlap:
            
            overlaps <- apply(spec.allshifts, 2, function(x) (ref + x) %>% 
                                is.na %>% "!"(.) %>% sum)
            
            # sd.na.or.zero <- apply(spec.allshifts, 1 , function(x) sd(x , na.rm = T) == 0)
            # sd.na.or.zero[is.na(sd.na.or.zero)] <- T
            # 
            
            dontCalc <- overlaps < min.overlap 
              
          
            if (sum(dontCalc) == length(overlaps)){return(earlyOut(dontCalc,overlaps))}
            
            if (sum(!dontCalc) == 0){return(earlyOut(dontCalc,overlaps))}
            
            lags <- lags[!dontCalc]
            overlaps <- overlaps[!dontCalc]
            spec.allshifts[,dontCalc] <- NA
            
          # NA-agnostic region overlap:
            #overlaps <- rep(length(ref),ncol(spec.allshifts))
      
      
    # heatmap(spec.allshifts %>% t, Rowv = NA, Colv = NA, scale = "none")
    # spec %>% simplePlot
    
    # Potentially much faster - but how does it work?
      # convscores <- outer(spec,ref,"*") %>% rowSums(na.rm = TRUE)
      # bestMatch <- (convscores/overlaps) %>% which.max
      # Calculate corr for the best match instead of all of them
      
    # heatmap(conv, Rowv = NA, Colv = NA, scale = "none")

  # Calculate the correlation between the ref and all rows of spec.allshifts
    # if (is_empty(spec.allshifts[,-dontCalc])){browser()}
    
    # spec.allshifts[is.na(spec.allshifts)] <- NaN
    
    # if(isempty(spec.allshifts[ , !dontCalc])){browser()}
    scores <- suppressWarnings(
      cor(spec.allshifts[ , !dontCalc],
                  ref,
                  use = "pairwise.complete.obs",
                  method = "pearson"))
    
    scores[is.infinite(scores)] <- NA
            
      # Calculate pvalue using # elements available to compare (overlaps, above)
      # and the correlation from the comparison:
          a=-abs(scores * sqrt( (overlaps-2) /(1-scores^2)))
          pvals=2*pt(a,(overlaps-2))
          # pvals %>% log %>% plot
          bestFit.idx <- which.min(pvals)
          bestFit.score <- scores[bestFit.idx]
          bestFit.pval <- pvals[bestFit.idx]
      
      # Convert back to lag 
        # bestFit.idx is the index of the row of indmat that was the best fit
        # that, in turn, corresponds to an index of lag
          
          bestFit.lag <- lags[bestFit.idx]
          indmat <- indmat[,!dontCalc]
          matchingInds <- indmat[,bestFit.idx]
          
          # **** Export rvals, pvals, etc. as NA if not passing dontCalc - this
          # allows them to be readily compiled at the same lengths for when in batches
          full.rvals <- full.pvals <- full.overlaps <- full.lags <- rep(NA, ncol(spec.allshifts))
          
          full.rvals[!dontCalc] <- scores
          full.pvals[!dontCalc] <- pvals
          full.overlaps[!dontCalc] <- overlaps
          full.lags[!dontCalc] <- lags
          
  return(list(rvals = full.rvals, 
              pvals = full.pvals,
              overlaps = full.overlaps,
              not.used = dontCalc,
              lags = full.lags,
              bestFit = list(which.lag = which(full.lags %in% bestFit.lag),
                             r = bestFit.score,
                             lag = bestFit.lag,
                             pval = bestFit.pval,
                             overlaps = overlaps[bestFit.idx],
                             specInds.matched = matchingInds) # matchinginds used to extract spec for plotting
              )
         )
}

#' Return early output for `xcorr` if no comparisons should be computed
#' 
#' This function returns a list of NA values for `rvals`, `pvals`, `lags`, and `bestFit`, 
#' as well as the `overlaps` values and `dontCalc` boolean vector that caused the 
#' function to exit early.
#'
#' @param dontCalc A boolean vector indicating which comparisons are not calculated
#' @param overlaps A numeric vector of the overlaps between the reference and each 
#'   shift of the spectrum
#' 
#' @return A list containing NA values for `rvals`, `pvals`, `lags`, and `bestFit`, 
#'   as well as the `overlaps` values and `dontCalc` boolean vector.
#' 
#'
#' @export
earlyOut <- function(dontCalc,overlaps){
  dummy <- rep(NA, length(overlaps))
  return(list(rvals = dummy, 
                pvals = dummy,
                overlaps = overlaps,
                not.used = dontCalc,
                lags = dummy,
                bestFit = list(which.lag = NA,
                               r = NA,
                               lag = 0,
                               pval = NA,
                               overlaps = NA,
                               specInds.matched = NA)
              )
         )
         
}