#'#' Spec-feature extraction
#'
#'  Using feature shape from filtered fse results, attempt to align the feature to each 
#'  spectrum in the dataset (within the spectral region from which the feature was 
#'  initially extracted). Accept the adjustment if correlation to the initial feature 
#'  improves and passes the provided r value threshold. Then recalculate the feature
#'  profile using STOCSY of the original driver. Calculate the least squares fit of
#'  the feature to each spectrum which passed the correlation threshold. 
#'  
#'  Think of this as: try alignment of spectral region to extracted feature shape. 
#'  If improved, accept, on a spectrum-by-spectrum basis. Then follow up with a 
#'  final round of STORM. 
#'  
#'  Used after filtering features, as it can take a while. 
#' 
#' - no need to test range contradictions
#' - just do fft-based alignment to current profile for ALL spectra
#' - allow lags of +/- (feature width-1) 
#' - least squares fits, rvals, rmses calculated. Use expand.fit() to get vals from fit. 
#' - updated subset (ss) inds, best lag for each ss spec, and updated feature object returned.
#' - in feat, profile and ss will be updated. MUST calculate feature positions for
#'   each ss spectrum after this: 
#'   
#'   feat <- list(profile = profile,
#'                position = f.pos,
#'                ss = ss,
#'                driver.relative = driver)
#'                
#'   al.fts <- apply_lags_feat(feat = feat,
#'                             xmat = xmat, 
#'                             lags = lags$lag.in.f2)
#'   
#' 
#' @param feature feature object (all feature)
#' @param f.ind feature index to do sfe on
#' @param xmat spectral matrix
#' @param ppm ppm vector
#' @param r.thresh correlation threshold
#' 
#' @return feat = aligned$feat,
#' lags for the aligned subset spectra
#' fits (fit feature profile to each spectrum)
#' rmses for fits
#' rvals for each subset spectrum (in case you want to redo filtering)
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @author MTJ

sfe <- function(feature, f.ind, xmat, ppm, r.thresh = 0.8){
  # Example inputs ####
    # feature <- feature.ma
    # f.ind <- 500
    # xmat <- fse.result$xmat
    # ppm <- fse.result$ppm
  
  # Clean up the feature (trim NAs on sides) ####
    profile <- feature$stack[f.ind, , drop = F]
    
    f.pos <- feature$position[f.ind, ,drop = F]
      cols <- f.pos %>% trim_sides(out = "inds")
      driver <- (cols %in% feature$driver.relative[f.ind]) %>% which
      f.pos <- f.pos[cols]
      profile <- profile[cols]
    ss <- feature$subset$ss.all[f.ind,] %>% which
  
  # Co-optimize alignment and feature profile ####
    # Includes correlation thresholding for feature shape after aligning
      
      feat <- list(profile = profile,
                   position = f.pos,
                   ss = ss,
                   driver.relative = driver)
      
      
      aligned <- tryCatch(
        expr = {
          align_spec2feat(feat = feat, xmat = xmat, r.thresh = r.thresh)
          # stackplot(aligned$valsmat)

        }, 
        error = function(cond){
          NULL
        }
      )
      
      test_nullish(aligned)
      
  # Fit the feature to each passing spectrum ####
    
    fits <- lapply(1:nrow(aligned$valsmat), function(m) {
      
      # * fit now has tryCatch included, returns NA-filled fit obj if failed
        fit <-  fit_leastSquares(v1 = aligned$feat$profile,
                                 v2 = aligned$valsmat[m,], 
                                 ppm = ppm[aligned$feat$position],
                                 plots = F,
                                 scale.v2 = F)

        return(list(fit = fit$fit,
                    rmse = fit$rmse))
    })
    
    # Note: this will fail if plots are null. Just returning vals.
    fits %>% test_nullish
    # Cannot get a null value out of this
    
    failed.specs <- lapply(fits, function(fit) is.na(fit$rmse)) %>% unlist
    succeeded <- !failed.specs
      if (any(succeeded)){
        fits <- fits[succeeded]
        aligned$lags <- aligned$lags[succeeded]
        aligned$rvals <- aligned$rvals[succeeded]
      }
    
    rmses <- lapply(fits, function(fit) fit$rmse) %>% unlist
    fits <- lapply(fits, function(fit) fit$fit)
    
  
  # Return updated feature info ####

        return(list(feat = aligned$feat,
                    lags = aligned$lags,
                    fits = fits,
                    rmses = rmses,
                    rvals = aligned$rvals)
             )
}