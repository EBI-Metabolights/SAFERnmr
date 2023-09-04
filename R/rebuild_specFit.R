# rebuild_specFit ############################################################################
#' 
#' ...
#'
#' @param sf spec feature information (includes backfit information)
#' @param feature features object for study
#' @param xmat spectral matrix for study
#' @param refmat reference spectral matrix for study (spectra on rows)
#' 
#' @return just the fit data
#' 
#' @importFrom magrittr %>%
#' 
#' @export
rebuild_specFit <- function(sf, feature, xmat, refmat){
  
  # sf <- sfs[x, ]
      
  # Get feature (the profile being fit) ####
  
      # Get the feature model
        
        # Starting from feature and match info:
        
          f.matched.pts <- sf$feat.start:sf$feat.end
          feat.gaps <- feature$position[sf$feat, f.matched.pts] %>% is.na
    
      # Get the ref segment (rf)
      
        # Starting from refmat and match info and feature model:
        
          rf <- refmat[ sf$ref, sf$ref.start:sf$ref.end ]
        
        # NA-fill the feature gaps
            
          rf[feat.gaps] <- NA

  # Get spec ####
  
    # Get the spectrum data
      
      spec.reg <- c(sf$spec.start,sf$spec.end) %>% expand_window(within = c(1,ncol(xmat)))
        window.loc <- sf$spec.start - min(spec.reg) # store this for now
  
      spec.segment <- xmat[sf$ss.spec, spec.reg]
        # simplePlot(spec.segment)
      window.pos <- 1:length(rf) + window.loc
      
      
  # Apply the fit 
  
            # lil.spec <- spec.segment[window.pos] 
            lil.ref <- rf
            
            #   plot_fit(fit)
            
          # Build expanded ref vect
          
            ref.segment <- rep(NA, length(spec.segment))
            ref.segment[window.pos] <- lil.ref
      
          ref.segment <- ref.segment * sf$fit.scale + sf$fit.intercept
          
    return(list(feat = ref.segment,
                spec = spec.segment))
}
