fit_leastSquares_fast <- function(v1, v2) {
  
  # Use-mask
  use <- !(is.na(v1) | is.na(v2))
  if (sum(use) < 3)
    return(list(a = NA_real_, b = NA_real_, rmse = NA_real_))
  
  x <- v1[use]
  y <- v2[use]
  
  # Closed-form linear regression: y = a + b*x
  xm <- mean(x)
  ym <- mean(y)
  dx <- x - xm
  dy <- y - ym
  
  b <- sum(dx * dy) / sum(dx * dx)
  a <- ym - b * xm
  
  # RMSE (fast, minimal)
  fit_vals <- a + b * x
  rmse <- sqrt(mean((y - fit_vals)^2))
  
  list(a = a, b = b, rmse = rmse)
}
