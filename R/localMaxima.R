#' Find the indices of the local maxima of a numeric vector
#'
#' Given a numeric vector, this function returns the indices of the local maxima.
#' A local maximum is defined as a data point that is greater than its adjacent
#' points. The function adds a minimum value to each end of the vector to ensure that the
#' endpoints can also be identified as local maxima, if appropriate.
#'
#' @param v A numeric vector
#'
#' @return A numeric vector containing the indices of the local maxima of \code{v}
#' @importFrom magrittr %>%
#'
#' @examples
#' localMaxima(c(1, 2, 3, 2, 1))
#'
#' @export
localMaxima <- function(v){
  # MTJ 2022
  # Credit https://stackoverflow.com/a/6836583/8104821  
  return(((c(min(v),v,min(v)) %>% diff %>% sign %>% diff) == -2) %>% which)
}