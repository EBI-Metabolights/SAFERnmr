# opt_specFit ############################################################################
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
#' @param xmat spectral matrix for study
#' @param refmat reference spectral matrix for study (spectra on rows)
#' 
#' @return updated specFit row (sf), with optimized fit coefficients and spec start/end 
#' 
#' @importFrom magrittr %>%
#' 
#' @export
opt_specFit <- function(sf, feature, xmat, refmat){
  
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
      
      
  # Optimize location by batmat fitting ####
      
      # res <- batman_fit_in(feat = rf, spec.segment, window.pos)
      res <- batman_local_opt(feat = rf, spec.segment, 
                              window.pos,
                              exclude.lowest = 0.5, 
                              opt.on = 'fraction.spec.accounted')
      # plot_fit(res$fit)
      
  # Update match information using lag ####
  
    sf$score <- res$fit$fraction.spec.accounted * res$fit$rval
      if (sf$score > 1){warning('somehow in opt_specFit, score exceeded 1!')}
    sf$fit.intercept <- res$fit$intercept
    sf$fit.scale <- res$fit$ratio
    sf$spec.start <- sf$spec.start - res$lag
    sf$spec.end <- sf$spec.end - res$lag
    
    return(list(sf = sf,
                feat = res$fit.ref,
                spec = res$fit.spec))
}
