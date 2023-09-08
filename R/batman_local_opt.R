# batman_local_opt ############################################################################
#' 
#' Batman-fit a shorter profile (feat) to a longer segment (spec.segment) at 
#' each internal lag. Pick the optimal lag based on the fraction.spec.accounted
#' by the fit (where the two signals overlap, excluding NAs). This works because
#' the Batman fitting does not allow feature signal to exceed spec.segment intensity
#' (excluding the exclude.lowest fraction of points in feature). The rationale is
#' that a feature which matches well will account for more of the relevant spec 
#' signal once fit. Since batman fitting is just calculation of ratios, this goes
#' pretty quick. 
#' 
#' Used to optimize after convolution-based matching during matching, or during 
#' backfitting. 
#'
#' @param sf spec feature information (includes backfit information)
#' @param feature features object for study
#' @param xmat spectral matrix for study
#' @param refmat reference spectral matrix for study (spectra on rows)
#' 
#' @return list: fit optimal fit information
#'               lag optimal lag value
#'               fit.ref (same size as spec.segment)
#'               fit.spec (same size as spec.segment)
#' 
#' @importFrom magrittr %>%
#' 
#' @export
batman_local_opt <- function(feat, spec.segment, window.pos, exclude.lowest = 0.5,
                             opt.on = 'fraction.spec.accounted' # or 'ratio'
                             ){
  
          # Loop through lags and assess the best fit based on % spectral signal
          # accounted for.
          # - NOTE: constrained by simply climbing from lag = 0 to local fit 
          #         quality maximum

            # which lags can we look at?
            # - must be valid inds
              
            max.lag <- tryCatch(
              {
                
                min(
                  c(pk_bounds(spec.segment) %>% unlist %>% matrix(nrow = 2) %>% diff %>% max,
                  pk_bounds(feat) %>% unlist %>% matrix(nrow = 2) %>% diff %>% max,
                  length(feat)-1),
                  na.rm = TRUE
                )
              },
              error = function(cond)
              {
                floor((length(feat)-1)/2) # want to put a limit on this!
              }
            )
              
            lags <- (-max.lag:max.lag)
            
            # Run the loop 
              
              fit.vals <- lapply(lags, function(lag){
                
                tryCatch(
                  {
                    lil.spec <- spec.segment[window.pos - lag]
                    lil.ref <- feat
                    # if (length(lil.spec) != length(lil.ref)){browser()}
                    fit <- fit_batman(feat = lil.ref, spec = lil.spec, exclude.lowest = 0.5)
                      data.frame(fsa = fit$fraction.spec.accounted,
                                 ratio = fit$ratio) # if optimization is constrained to peak, use ratio to opt

                  },
                  # warning = function(cond){
                  #     data.frame(fsa = 0,
                  #                ratio = 0)
                  # },
                  error = function(cond){
                      data.frame(fsa = 0,
                                 ratio = 0)
                  })
                  
              }) %>% rbindlist
            
            # Optimize to local max. (otherwise you risk being outside of 
            # the signal which produced the correlation to begin with!)
              
              if (opt.on == 'fraction.spec.accounted'){
                opt.vals <- fit.vals$fsa
              } else {
                if (opt.on == 'ratio'){
                  opt.vals <- fit.vals$ratio
                } else {
                  opt.vals <- fit.vals$fsa
                }
              }
              
              opt.lag <- climb(start.pt = max.lag + 1, v = opt.vals) %>% lags[.]
              # opt.lag <- pct.acct %>% which.max %>% lags[.]
          
          # Calculate final fit using optimal lag
            lil.spec <- spec.segment[window.pos - opt.lag] 
            lil.ref <- feat
            
            fit <- fit_batman(feat = lil.ref, spec = lil.spec, exclude.lowest = 0.5)

            fit$rval <- cor(fit$feat.fit, fit$spec.fit, use = "pairwise.complete.obs")
            
            #   plot_fit(fit, type = 'auc')
            
          # Build expanded ref vect
          
            ref.segment <- rep(NA, length(spec.segment))
            ref.segment[window.pos - opt.lag] <- fit$feat.fit
              # plot_fit(list(feat.fit = ref.segment,
              #               spec.fit = spec.segment), type = 'auc')
              # simplePlot(ref.segment)
              # simplePlot(spec.segment)
            
        return(list(fit = fit,
                    lag = opt.lag,
                    fit.ref = ref.segment,
                    fit.spec = spec.segment))
  
}

##################################################################################

#' Don't optimize, just fit
#'
#' @param sf spec feature information (includes backfit information)
#' @param feature features object for study
#' @param xmat spectral matrix for study
#' @param refmat reference spectral matrix for study (spectra on rows)
#' 
#' @return list: fit optimal fit information
#'               lag optimal lag value
#'               fit.ref (same size as spec.segment)
#'               fit.spec (same size as spec.segment)
#' 
#' @importFrom magrittr %>%
#' 
#' @export
batman_fit_in <- function(feat, spec.segment, window.pos, exclude.lowest = 0.5
                      ){
  
          # Loop through lags and assess the best fit based on % spectral signal
          # accounted for.
          # - NOTE: constrained by simply climbing from lag = 0 to local fit 
          #         quality maximum

              
            
              opt.lag <- 0
              # opt.lag <- pct.acct %>% which.max %>% lags[.]
          
          # Calculate final fit using optimal lag
            lil.spec <- spec.segment[window.pos - opt.lag] 
            lil.ref <- feat
            
            fit <- fit_batman(lil.ref, lil.spec, exclude.lowest = 0.5)
            fit$rval <- cor(fit$feat.fit, fit$spec.fit, use = "pairwise.complete.obs")
            
            #   plot_fit(fit, type = 'auc')
            
          # Build expanded ref vect
          
            ref.segment <- rep(NA, length(spec.segment))
            ref.segment[window.pos - opt.lag] <- fit$feat.fit
              # plot_fit(list(feat.fit = ref.segment,
              #               spec.fit = spec.segment), type = 'auc')
              # simplePlot(ref.segment)
              # simplePlot(spec.segment)
            
        return(list(fit = fit,
                    lag = opt.lag,
                    fit.ref = ref.segment,
                    fit.spec = spec.segment))
  
}
