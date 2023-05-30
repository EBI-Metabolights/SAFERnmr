#' Full Width at Half Maximum (FWHM) calculation
#'
#' Given an isolated peak, calculates its Full Width at Half Maximum (FWHM),
#' i.e., the width of the peak at half of its maximum height.
#' 
#' @param isolatedpk a numeric vector representing an isolated peak
#' @return a numeric value representing the FWHM of the peak
#'
#' @importFrom magrittr %>%
#'
#'
#' @export
fwhm <- function(isolatedpk){
  hm <- max(isolatedpk)/2
  hm.pts <- localMinima(abs(isolatedpk - hm)) %>% sort %>% .[1:2]
  return(diff(hm.pts))
}