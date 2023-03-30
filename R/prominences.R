#' Calculate the prominence of each peak
#'
#' The 'prominence' of a peak is the height of the peak relative to the lower bound
#' in this case.
#'
#' @param peaks A list of peak information as output from the \code{extractPeaks_corr} function.
#' @param v2 A vector of the y-values of the peaks.
#'
#' @return A vector of the prominences of each peak.
#' @export
#'
#' 
#' 
#' @importFrom magrittr %>%
#' 
prominences <- function(peaks, v2){
  maxima <- v2[peaks$peaks]
  proms <- lapply(1:length(maxima), function(p) {
    v.hts <- peaks$bounds[[p]] %>% as.numeric %>% v2[.]
    return(min(maxima[p] - v.hts, na.rm = T))
  }) %>% unlist
  return(proms)
}