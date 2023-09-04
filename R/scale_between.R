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
#' scale_between(x, 0, 10)
#'
#' @export
scale_between <- function(data, lower = 0, upper = 1){
  if (any(is_nullish(data)) | all(is.na(data))){return(NA)}
  return((data - min(data,na.rm = TRUE))/(max(data,na.rm = TRUE)-min(data,na.rm = TRUE)) # make positive and between 0 and 1
          * (upper-lower)                         # stretch to fit abs(new range)
          + lower)                                # slide down to lower bound
}
# upper = 1
# lower = 0
# data = vector
# bottom <- min(data, na.rm = TRUE)
# top <- max(data, na.rm = TRUE)
# data.range <- top-bottom
# new.range <- upper - lower
# data <- (data + lower) / (new.range/data.range) + bottom  
# # Scaling, broken down
#   plot(data)
#   
#   data <- data - bottom
#     plot(data)
#   
#   # Scale data to new range size
#   data <- data * (new.range/data.range) # combine range sizes, then multiply to vector
#     plot(data)
#     
#   # Put bottom of data @ bottom of new range
#   data <- data - lower
#     plot(data)

# # Undo scaling for data to show it works
#   data <- (data + lower) / (new.range/data.range) + bottom
#   all(data == vector)