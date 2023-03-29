#' Find local minima in a vector
#'
#' Given a numeric vector, returns the indices of local minima in the vector.
#'
#' @param v a numeric vector
#' @return a numeric vector containing the indices of local minima in \code{v}.
#' @export
localMinima <- function(v){
  return( v %>% diff %>% ">"(0) %>% c(FALSE,.,TRUE) %>% diff %>% ">"(0) %>% which )
}