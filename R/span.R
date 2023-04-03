#' Calculates the span of a numeric vector
#'
#' Calculates the span of a numeric vector, which is defined as the range of the values
#' plus one.
#'
#' @param v A numeric vector
#' @return The span of the vector
#''
#' @export
span <- function(v) {
  return(v %>% range(., na.rm = TRUE) %>% diff() + 1)
}
