#' Find local minima in a vector
#'
#' Given a numeric vector, returns the indices of local minima in the vector. Endpoints
#' can be included. Specifically, a point is considered a local minimum if there
#' are no lesser points adjacent.
#'
#' @param v a numeric vector
#' @return a numeric vector containing the indices of local minima in \code{v}.
#'
#' @importFrom magrittr %>%
#'
#' @export
localMinima <- function(v) {
  return(v %>% diff() %>% ">"(0) %>% c(FALSE, ., TRUE) %>% diff() %>% ">"(0) %>% which())
}
