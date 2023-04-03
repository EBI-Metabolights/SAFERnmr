#' sortPairs.R
#'
#' Given a two-row matrix where each column contains a pair of numbers,
#' this function sorts the pairs such that each column's first element is less
#' than its second element. The function assumes that pairs are already
#' sorted, and quickly checks if any pairs need to be switched.
#'
#' @param pairs A two-row matrix where each column contains a pair of numbers
#' @return A two-row matrix with pairs sorted in ascending order
#' @examples
#' pairs <- matrix(c(1, 3, 2, 5, 4, 3), nrow = 2)
#' sortPairs(pairs)
#' # Output:
#' #      [,1] [,2] [,3]
#' # [1,]    1    2    3
#' # [2,]    3    5    4
#'
#' @export
sortPairs <- function(pairs) {
  # Fast sort
  # assume no switch
  spairs <- pairs

  # Find columns that need switching: is second element < first element?
  doswitch <- pairs[2, ] < pairs[1, ]

  # For columns that need switching, reverse the order
  spairs[2, doswitch] <- pairs[1, doswitch]
  spairs[1, doswitch] <- pairs[2, doswitch]

  return(spairs)
}
