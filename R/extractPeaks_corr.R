#' Extract local maxima and minima from a correlation vector
#'
#' This function takes a correlation vector and returns the indices of the local maxima and minima, and the bounds for the corresponding windows around each maximum.
#'
#' @param corr A correlation vector.
#' @param mask A boolean mask for excluding points.
#' @param plots Whether to plot the resulting peaks and bounds.
#'
#' @return A list containing the peak indices, the bounds of the windows around each peak, and the indices of any masked regions, and whether or not each peak is a true peak (i.e., require a lower point on each side)
#'
#' @importFrom graphics plot
#' @importFrom ggplot2 aes geom_line geom_vline geom_point scale_x_reverse
#' @importFrom magrittr %>%
#'
#' @examples
#' corr <- rnorm(100)
#' peaks <- extractPeaks_corr(corr, plots = TRUE)
#'
#' @export
extractPeaks_corr <- function(corr, mask = NULL, plots = FALSE) {
  # peakInds <- peaks.init[i]
  # for each peak ind, expand to closest corr min
  # for the entire list of bounds, report the range.
  # all inds within window, then converted to ppm inds

  if (is.matrix(corr)) {
    corr <- as.numeric(corr)
  }
  if (is.null(mask)) {
    mask <- !is.na(corr)
  }

  # Get inds of local maxima using a mask for statistical cutoff (note: these don't include bounds)
  excluded <- (!mask) %>% which()
  corr[excluded] <- min(corr)

  # Ensure the region bounds can be maxs or mins

  localMaxs <- localMaxima(corr)
  localMins <- localMinima(corr)

  # Include mask bounds in local minima (excluding mask == False) , and ensure no duplicates
  diffMask <- mask %>%
    as.integer() %>%
    diff()
  maskBounds <- diffMask %>%
    "!="(., 0) %>%
    which()
  maskBounds <- maskBounds + as.integer(diffMask[maskBounds] > 0)

  localMins <- c(localMins, maskBounds) %>%
    sort() %>%
    unique()

  # Get the local minima about each point
  # For each max point, extend bounds in the left and right directions
  # Stop until you hit a window bound or local min, whichever is first

  minlist <- lapply(
    localMaxs,
    function(x) corr_expand(peak = x, localMins = localMins, vRange = c(1, length(corr)))
  )

  # Check extracted window by plotting. Color in these STOCSYs is abs(cc)
  if (plots) {
    df <- data.frame(r = corr, inds = seq_along(corr))
    g <- ggplot(data = df, aes(x = inds, y = r)) +
      geom_line(na.rm = TRUE) +
      geom_vline(xintercept = unlist(minlist), linetype = 1, col = "grey") +
      geom_point(
        data = data.frame(xs = localMaxs, ys = corr[localMaxs]),
        aes(x = xs, y = ys), na.rm = TRUE
      ) +
      # Indicate excluded points
      geom_point(
        data = data.frame(xs = excluded, ys = corr[excluded]),
        aes(x = xs, y = ys),
        shape = 4, na.rm = TRUE
      ) +
      scale_x_reverse()
    plot(g)
  }

  return(list(
    peaks = localMaxs,
    bounds = minlist,
    maskBounds = maskBounds,
    truePeak = !(localMaxs %in% c(1, length(corr), maskBounds))
  ))
}
