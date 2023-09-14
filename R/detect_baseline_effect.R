#' Detect baseline effect in feature profile (also indicates monotonic features, which are easy to get with STORM)
#' Goes through a couple of checks:
#' -> linearly interpolate across NA gaps
#' - is the feature a singlet?
#' - is the peak prominence (max vs. lower bound) of at least one peak > prom.ratio x the entire feature intensity range
#' - do the valleys and peaks in the feature profile correlate strongly (r > cutoff.corr)? If so, likely a baseline effect.
#'
#' @param feat A numeric vector containing feature data.
#' @param cutoff.corr The maximum allowed correlation between predicted and
#'   observed peak intensities.
#' @param prom.ratio The minimum required ratio of peak prominence to feature
#'   intensity range.
#'   
#' @importFrom magrittr %>%
#'   
#' @return A list containing two logical values. The first indicates whether
#'   the feature passed the peak prominence filter, and the second indicates
#'   whether the feature passed the peak fit filter.
#'
#' @export
detect_baseline_effect <- function(feat, cutoff.corr = 0.99, prom.ratio = 0.3){
  
      # Get (simple) peaks, check if singlet
        # Connect noncontiguous feature bits with linear interp
      
        # feat <- feature$stack[1, , drop = FALSE] %>% trim_sides
        # Extract peaks from feature
          
          feat <- pracma::interp1(x = 1:length(feat), y = c(feat))
          pks <- extractPeaks_corr(feat)

            # True peaks have 2 bounds < peak max
              truepks <- pks$truePeak %>% which
              
              res <- erode2(feat)
              # plot_spec(spec = df$r2, ppm = df$iteration)

              not.singlet = length(truepks)>1
              
              if (is_nullish(not.singlet)){not.singlet <- FALSE}
             
        # Extract local prominences from feature ####
          
          if (length(truepks) > 0){ # can't calculate prominences of nonexistent peaks
            
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
              
          } else {pass.prom <- FALSE}

          if (is_nullish(pass.prom)){pass.prom <- FALSE}
            
        # See if bounds track peaks (requires > 1 peak) ####
            if (sum(truepks)>1){ # need > 1 peaks to do this
              
              bnds <- pks$bounds %>% unlist %>% as.numeric
              bnd.vals <- feat[bnds]
              peaks <- pks$peaks
              
                # Fit quadradic model
                  
                  bndmodel <- lm(bnd.vals ~ poly(bnds,2))
                  
                # Interpolate peaks points based on model
                  # Include edge points, peak bounds, and peak maxima
                    allpoints.with.edges <- unique(sort(c(pks$peaks, bnds, 1, length(feat))))
                    
                  # Predict the reduced profile if the peaks were removed
                    pkvals.pred <- predict(bndmodel, data.frame(bnds=allpoints.with.edges)) %>% as.numeric
                  # confint(bndmodel, level=0.95)
                  # plot(bnds, bnd.vals, col='deepskyblue4')
                  #   lines(bnds,feat[bnds],col='firebrick1',lwd=3)
                  #   lines(x = pks$peaks, y = feat[pks$peaks], col = "blue")
                  
                # Strong correlation (implies?) strong baseline effect - remove feature
                   pass.fit <- cor(pkvals.pred, feat[allpoints.with.edges]) < cutoff.corr
                 
              } else {pass.fit <- F}
              
            if (is_nullish(pass.fit)){pass.fit <- F}
            
        # Return results ####
          return(list(pass.prom = pass.prom,
                      pass.fit = pass.fit,
                      not.singlet = not.singlet))
 
}