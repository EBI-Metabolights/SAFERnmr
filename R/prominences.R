#' Calculate the prominence of each peak
#'
#' The prominence of a peak is the height of the peak relative to the lowest contour line that encloses
#' the peak and no higher peak. This function calculates the prominence of each peak in a vector.
#'
#' @param peaks A list of peak information as output from the \code{peakdet} function.
#' @param v2 A vector of the y-values of the peaks.
#'
#' @return A vector of the prominences of each peak.
#'
#' @importFrom magrittr %>%
#' 
#' @export
prominences <- function(peaks, v2){
  maxima <- v2[peaks$peaks]
  proms <- lapply(1:length(maxima), function(p) {
    v.hts <- peaks$bounds[[p]] %>% as.numeric %>% v2[.]
    return(min(maxima[p] - v.hts, na.rm = T))
  }) %>% unlist
  return(proms)
}