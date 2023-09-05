# opt_specFit_slim ############################################################################
#' 
#' Given a specFit row (match info plus ss.spec information), expand the matched
#' spec region to span-1 on each side of the specfit. For each internal lag, do 
#' batman-style fitting and assess based on climbing to the first local maximum 
#' spec signal accounted for (starting from current fit position). Return the 
#' updated specFit row. 
#' 
#' Note: this is meant to be run during backfitting. Using after that will not 
#' update some fields calculated there, e.g.:
#' 
#'    "bffs.res" "bffs.tot" "rmse" "rmse.biased" "pct.ref"
#'
#' @param sf spec feature information (includes backfit information)
#' @param feature features object for study
#' @param spec spectral matrix for study
#' @param refmat reference spectral matrix for study (spectra on rows)
#' 
#' @return updated specFit row (sf), with optimized fit coefficients and spec start/end 
#' 
#' @importFrom magrittr %>%
#' 
#' @export
opt_specFit_slim <- function(sf, feat.model, refspec, spec){
  
  # sf <- sfs[x, ]
  #   sf <- match.info[x, ] %>% [propagated to ss.specs]
  # feat.model <- feature$position[sf$feat, ]
  # refspec <- refmat[ sf$ref, ]
  # spec <- xmat[ sf$ss.spec, ]
     
  # Get feature (the profile being fit) ####
  
      # Get the feature model
        
        # Starting from feature and match info:
        
          f.matched.pts <- sf$feat.start:sf$feat.end
          feat.gaps <- feat.model[ f.matched.pts ] %>% is.na
    
      # Get the ref segment (rf)
      
        # Starting from refmat and match info and feature model:
        
          rf <- refspec[ sf$ref.start:sf$ref.end ]
        
        # NA-fill the feature gaps
            
          rf[feat.gaps] <- NA

  # Get spec ####
  
    # Get the spectrum data
      
      spec.reg <- c(sf$spec.start,sf$spec.end) %>% expand_window(within = c(1,length(spec)))
        window.loc <- sf$spec.start - min(spec.reg) # store this for now
  
      spec.segment <- spec[spec.reg]
        # simplePlot(spec.segment)
      window.pos <- 1:length(rf) + window.loc
      
      
  # Optimize location by batmat fitting ####
      
      res <- batman_local_opt(feat = rf, spec.segment, 
                              window.pos,
                              exclude.lowest = 0.5, 
                              opt.on = 'fraction.spec.accounted')
      # plot_fit(res$fit)
      
  # Update match information using lag ####
  
    sf$fraction.spec.accounted <- res$fit$fraction.spec.accounted
    sf$fit.rval <- res$fit$rval
    sf$fit.intercept <- res$fit$intercept
    sf$fit.scale <- res$fit$ratio
    sf$spec.start <- sf$spec.start - res$lag
    sf$spec.end <- sf$spec.end - res$lag
    
    return(list(sf = sf,
                feat = res$fit.ref,
                spec = res$fit.spec))
}
