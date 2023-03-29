#' Calculate Pearson correlations and covariance using STOCSY for a subset of variables.
#' 
#' 
#' @param X A matrix of spectral data with variables in the columns.
#' @param driverInd A vector of indices for the variables that will be used as drivers.
#' 
#' @return A list with two elements: \code{cc}, a vector of Pearson correlations 
#'         between each variable in \code{X} and the variables in \code{driverInd}; and
#'         \code{cv}, a vector of covariance values between each variable in \code{X} 
#'         and the variables in \code{driverInd}.
#' 
#' @importFrom coop pcor covar
#'
#' @export
stocsy_slim <- function(X, driverInd){
  res <- list()

  l <- X[, driverInd]
  res$cc <- apply(X, 2, function(x, ll = l) {coop::pcor(x, ll)})
  res$cv <- apply(X, 2, function(x, ll = l) {coop::covar(x, ll)})
  
  return(res)
}