#' Get bounds for local maxima in a vector
#'
#' This function takes a numeric vector \code{v} and returns the bounds of the regions containing the local maxima in \code{v}. Updated to include only true peaks (those with a lower value on either side; no endpoints).
#'
#' @param v A numeric vector
#' @return A list of numeric vectors, where each element corresponds to the bounds of a region containing a local maximum in \code{v}
#' @examples
#' v <- c(1, 2, 3, 2, 1, 2, 2, 1, 2)
#' pk.bounds(v)
#' # Returns: list(c(1, 3), c(4, 5))
#' @export
pk_bounds <- function(v){
  pks <- extractPeaks_corr(v)
  return(pks$bounds[pks$truePeak])
}