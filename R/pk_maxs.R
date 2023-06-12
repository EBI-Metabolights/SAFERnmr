#' Get indices of maxima of true peaks in a vector. Wrapper for extractPeaks_corr()
#'
#' This function takes a numeric vector \code{v} and a logical mask \code{mask}, and returns the indices of the local maxima in \code{v} that are not masked by \code{mask}.
#'
#' @param v A numeric vector
#' @param mask A logical vector of the same length as \code{v}, indicating which values of \code{v} should be masked (i.e., excluded from consideration as local maxima)
#' @return A numeric vector containing the indices of the local maxima in \code{v} that are not masked by \code{mask}
#' @importFrom magrittr %>%
#' @examples
#' v <- c(1, 2, 3, 2, 1)
#' mask <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
#' pk_maxs(v, mask)
#' # Returns: 3
#' @export
pk_maxs <- function(v, mask){
  pks <- extractPeaks_corr(v, mask = mask)
  return(pks$peaks[pks$truePeak] %>% unlist)
}