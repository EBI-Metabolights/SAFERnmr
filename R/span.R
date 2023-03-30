#' Calculates the span of a numeric vector
#'
#' Calculates the number of elements needed to contain fillbetween(range(v)), where v is an
#' integer array of indices.
#'
#' @param v A numeric vector
#' @return The span of the vector
#''
#' @export
span <- function(v){
  return(v %>% range(., na.rm = TRUE) %>% diff +1)
}