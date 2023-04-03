#' Pad matrix with rows and columns of specified value
#'
#' This function takes a matrix and pads it with rows and columns of a specified
#' value. The number of rows and columns to add can also be specified. The new
#' rows and columns will have the same value.
#'
#' @param x A matrix to pad
#' @param use The value to fill new rows and columns with (default is NA)
#' @param row.by The number of rows to add to each side of the matrix
#' @param col.by The number of columns to add to each side of the matrix
#'
#' @return A padded matrix
#' @export
#'
#' @examples
#' m <- matrix(1:4, ncol = 2)
#' padmat(m, use = 0, row.by = 1, col.by = 1)
#'
#' m2 <- matrix(1:6, ncol = 2)
#' padmat(m2, row.by = 2, col.by = 3)
#'
#' m3 <- matrix(1:9, ncol = 3)
#' padmat(m3, use = -1, row.by = 1, col.by = 2)
#'
#' @seealso \code{\link{sub2indR}}
#'
#' @keywords matrix manipulations
padmat <- function(x, use = NA, row.by = 0, col.by = 1) {
  if (!is.matrix(x)) {
    x <- as.matrix(x) %>% t()
  }
  coords <- expand.grid(row = 1:nrow(x), col = 1:ncol(x))
  rows.in.newmat <- nrow(x) + row.by * 2
  xp <- matrix(data = use, rows.in.newmat, ncol(x) + col.by * 2)
  linds <- sub2indR(
    rows = coords$row + row.by,
    cols = coords$col + col.by,
    rows.in.newmat
  )
  xp[linds] <- x
  return(xp)
}
