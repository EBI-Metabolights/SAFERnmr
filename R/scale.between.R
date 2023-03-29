#' Scale values in a vector to be between a specified range.
#'
#' This function takes a numeric vector and scales the values to be between a specified lower and upper bound.
#' 
#' @param data A numeric vector to be scaled.
#' @param lower A numeric value representing the lower bound of the output range. Default is 0.
#' @param upper A numeric value representing the upper bound of the output range. Default is 1.
#'
#' @return A numeric vector with values scaled between the lower and upper bounds.
#'
#' @examples
#' x <- c(1, 3, 5, 7, 9)
#' scale.between(x, 0, 10)
#' # Returns: 0 3 5 7 10
#'
#' @export
scale.between <- function(data, lower = 0, upper = 1){
  return((data - min(data,na.rm = TRUE))/(max(data,na.rm = TRUE)-min(data,na.rm = TRUE)) # make positive and between 0 and 1
          * (upper-lower)                         # stretch to fit abs(new range)
          + lower)                                # slide down to lower bound
}
