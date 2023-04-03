#' Assigns a unique label to each run of 1s in a binary vector.
#'
#' This function takes a binary vector and assigns a unique integer label to each
#' run of 1s in the vector. The labels start at 1 and increment by 1 for each new
#' run. Runs of 0s are assigned a label of 0.
#'
#' @param v A binary vector.
#' @return A numeric vector with the same length as `v` containing the labels
#'   assigned to each run of 1s.
#' @examples
#' v <- c(0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1)
#' run.labels(v)
#'
#' @export
run.labels <- function(v) {
  breaks <- (c(0, v) %>% diff()) > 0
  labs <- breaks %>% cumsum()
  labs[!v] <- 0
  return(labs)
}


#' Computes the lengths of runs of 1s in a binary vector.
#'
#' This function takes a binary vector and computes the lengths of all runs of
#' consecutive 1s in the vector. The output is a numeric vector, with length equal
#' to the number of 1s in the input vector, where each point in each run is labeled
#' by its run length.
#'
#' @param v A binary vector.
#' @return A numeric vector containing the lengths of all runs of consecutive
#'   1s in the input vector.
#' @examples
#' v <- c(0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1)
#' runs.labelBy.lengths(v)
#'
#' @export
runs.labelBy.lengths <- function(v) {
  if (is.logical(v)) {
    v2 <- rep(1, length(v))
    v2[!v] <- 0
    v <- v2
  }
  runs <- rle(v)
  runs$values <- runs$values * runs$lengths
  return(runs %>% inverse.rle())
}
