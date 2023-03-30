#' Detect baseline effect in feature profile (also indicates monotonic features, which are easy to get with STORM)
#' Goes through a couple of checks:
#' - is the peak prominence (max vs. lower bound) of at least one peak > prom.ratio x the entire feature intensity range
#' - do the valleys and peaks in the feature profile correlate strongly (r > cutoff.corr)? If so, likely a baseline effect.
#' 
#' @param feat A numeric vector containing feature data.
#' @param cutoff.corr The maximum allowed correlation between predicted and
#'   observed peak intensities.
#' @param prom.ratio The minimum required ratio of peak prominence to feature
#'   intensity range.
#' @return A list containing two logical values. The first indicates whether
#'   the feature passed the peak prominence filter, and the second indicates
#'   whether the feature passed the peak fit filter.
#'
#' @export
detect.baseline.effect <- function(feat, cutoff.corr = 0.99, prom.ratio = 0.3){
      
        # feat <- feature$stack[1, , drop = FALSE] %>% trim.sides
        # Extract peaks from feature
          pks <- extractPeaks_corr(feat)
          # pks <- extractPeaks_corr(feat %>% trim.sides, plots = TRUE)
          
            # True peaks have 2 bounds < peak max
              truepks <- pks$truePeak %>% which
              # # Also ensure that neither of the adjacent peak maxima are greater
              #     lmax <- feat[pks$peaks] %>% localMaxima
              #     true.lmax <- truepks %in% lmax
              
            # Chuck it if singlet or no peaks
              if (sum(truepks)<=1){return(list(pass.prom = FALSE,
                                               pass.fit = FALSE))}
              
        # Extract local prominences from feature
          proms <- lapply(1:length(truepks), function(x){
            # Take the diff between the peak max and the lower of the two bounds.
                                                  pkmax <- feat[pks$peaks[x]]
                                                  lbound <- min(feat[pks$bounds[[x]]%>%unlist])
                                                  return(pkmax - lbound)
                                            }) %>% unlist
        # Get feature intensity range
          irange <- range(feat, na.rm = TRUE) %>% diff
      
        # Get the relative peak sizes and compare to cutoff
          ratios <-  proms / irange
          pass.prom <- any(ratios > prom.ratio)

        # See if bounds track peaks (requires > 1 peak)
          bnds <- pks$bounds %>% unlist %>% as.numeric
          bnd.vals <- feat[bnds]
          peaks <- pks$peaks
          
          # Fit quadradic model
            
            bndmodel <- lm(bnd.vals ~ poly(bnds,2))
            
          # Interpolate peaks points based on model
            pkvals.pred <- predict(bndmodel, data.frame(bnds=pks$peaks)) %>% as.numeric
            # confint(bndmodel, level=0.95)
            # plot(bnds, bnd.vals, col='deepskyblue4')
            #   lines(bnds,feat[bnds],col='firebrick1',lwd=3)
            #   lines(x = pks$peaks, y = feat[pks$peaks], col = "blue")
            
          # Strong correlation (implies?) strong baseline effect - remove feature
             pass.fit <- cor(pkvals.pred, feat[pks$peaks]) < cutoff.corr
          
        # If there are any ratios above the threshold, pass filter
          return(list(pass.prom = pass.prom,
                      pass.fit = pass.fit))
 
}