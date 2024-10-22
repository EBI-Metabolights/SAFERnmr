#' resample_spectra
#'
#' Linearly interpolates (or decimates) a ppm vector and spectral matrix to a 
#' (typically lower) number of spectral points.
#' Preserves spectral range and intensities. Should be replaced with binning or 
#' moving average at some point.
#'
#' @param xmat spectral matrix (spectra on rows)
#' @param ppm ppm vector to match xmat columns
#' @param npoints number of points to include
#' @param cores number of cores for parallel processing (default 1)
#'
#' @return spectral matrix and ppm value linearly interpolated to new number of points
#'
#' @importFrom stats approx
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply
#'
#' @export
resample_spectra <- function(xmat, ppm, npoints, cores = 1){
   
  if (ppm[2] < ppm[1]){
    ppm.resampled <- seq(from = max(ppm), to = min(ppm), length.out = npoints)
  } else {
    ppm.resampled <- seq(from = min(ppm), to = max(ppm), length.out = npoints)
  }
  
  xl <- split(xmat, row(xmat))
  
  X.resampled <- mclapply(xl, function(spec)
    {
      
      # Linear interpolation to dataset ppm axis ####
        # spec <- xl[[1]]
        # roi <- c(2.4,6.8) %>% sort
        # inds <-  vectInds(roi, ppm) %>% fillbetween
        #  simplePlot(spec[inds], xvect = ppm[inds])
        
        spec.resampled <- approx(ppm, spec, 
                                 ppm.resampled, # get values for these ppms (can be selective, or if goes oob, then NA)
                                 method="linear", rule=1, f=0) %>% .$y
        
        # inds <-  vectInds(roi, ppm.resampled) %>% fillbetween
        # simplePlot(spec.resampled[inds], xvect = ppm.resampled[inds])
        # simplePlot(spec[inds], xvect = ppm.resampled[inds])
        
      # Return a copy of the list element with mapped data added ####
        return(spec.resampled)
        
    }, mc.cores = cores
  ) %>% do.call(rbind,.)

  return(list(ppm = ppm.resampled,
              spectra = X.resampled))
}
