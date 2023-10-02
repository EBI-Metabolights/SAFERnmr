#' Expand a correlation window around a given peak index to the nearest correlation minima outside of the window
#' bounds.
#'
#' @param xmat A matrix of NMR spectra.
#' @param ppm A vector of ppm values corresponding to the columns of xmat.
#' @param peakInd The index of the peak around which to expand the window.
#' @param wind A vector defining the window (column inds) around the peak.
#'
#' @importFrom magrittr %>%
#'
#' @return A list with the following elements:
#'   \item{peak}{The index of the peak used for window expansion.}
#'   \item{intcorr}{The integral of the correlation within the expanded window.}
#'   \item{data}{A data frame containing the correlation values, covariance values, and corresponding ppm values within the expanded window.}
#'   \item{sr}{The STOCSY object for the expanded window.}
#'   \item{wind}{The original window vector.}
#'   \item{newWind}{The expanded window vector.}
#'   \item{corrRbound}{The index of the right boundary of the expanded window in the \code{wind} vector.}
#'   \item{corrLbound}{The index of the left boundary of the expanded window in the \code{wind} vector.}
#'

#' @importFrom magrittr %>%
#' @importFrom ggplot2 geom_line geom_vline geom_point scale_x_reverse
#' @export
corr_expand_window <- function(xmat, ppm, peakInd, wind){
  # setClass("stocsyObject", slots=list(r="numeric", cov="numeric", driver= "numeric"))
  
  stocsy <- function(mat, sel.reg, driver,
                     plotting = FALSE){
    
    driver.ind <- vectInds(driver, sel.reg)
    
    # sr <- new("stocsyObject", 
    #           r = cor(mat, mat[, driver.ind]), 
    #           cov = cov(mat, mat[, driver.ind]), 
    #           driver = driver,
    #           driver.relative = driver.ind)
    
    sr <- list(
      r = cor(mat, mat[, driver.ind]),
      cov = cov(mat, mat[, driver.ind]),
      driver = driver,
      driver.relative = driver.ind
    )
    
    return(sr)
  }
  
  # peakInd <- peaks.init[i]
  # Select a region that's close by and not out of bounds
  # no bigger than max(ppm), no smaller than min(ppm), within n * 2 indices
  
    wind <- max((min(wind)),1):min((max(wind)),length(ppm)) # stay in bounds of ppm vector
    offset <- peakInd-min(wind) # (this might change if wind is out of bounds)
    
  # Local STOCSY within window

    sr <- stocsy(xmat[,wind], wind, peakInd,
                 plotting = FALSE)

  # Get local mins in correlation vector around peak
    localMins <- localMinima(sr$r)
  
    corrRbound <- (localMins > offset+1) %>%  # Take local mins to the right of driver peak    --- inds in wind
      which() %>% localMins[.] %>% min() %>%  # get the ind of the leftmost one --- inds in wind
      c(.,length(wind)) %>% min()   # take lesser between that and right bound of wind. ----  
    # #%>% wind[.] #%>% ppm[.] # transform to ppm inds, and then to ppm vals --- inds in ppm -> ppm values 
    
    corrLbound <- (localMins < offset+1) %>%  # Take the local mins < peak
      which() %>% localMins[.] %>% max() %>%  # get the ind of the rightmost one to left of peak
      c(.,1)            %>% max()   # make sure it's not less than the left wind bound 
  # if you wanted ppm inds or vals:  %>% wind[.] #%>% ppm[.]
    # browser()
  # Check extracted window by plotting. Color in these STOCSYs is abs(cc)
                  # df <- data.frame(r = sr$r,inds = wind)
                  # g <- ggplot(data = df, aes(x = inds, y = r)) +
                  #   geom_line() +
                  #   geom_vline(xintercept = peakInd, linetype = 2, col = "grey") +
                  #   geom_vline(xintercept = wind[c(corrLbound, corrRbound)], linetype = 1, col = "grey") +
                  #   geom_point(data = data.frame(xs = wind[localMins], ys = sr$r[localMins]),
                  #          aes(x = xs, y = ys)) +
                  #   scale_x_reverse()
                  # plot(g)
  
  # Store data within list
    vals <- list()
    vals$peak <- peakInd
    vals$intcorr <- sum(sr$r[corrLbound:corrRbound])
    vals$data <- rbind(corr = sr$r[corrLbound:corrRbound],
                       cov = sr$cov[corrLbound:corrRbound],
                       inds = wind[corrLbound:corrRbound])
    vals$sr <- sr
    vals$wind <- wind
    vals$newWind <- wind[corrLbound:corrRbound]
    vals$corrRbound <- corrRbound
    vals$corrLbound <- corrLbound
  return(vals)
}