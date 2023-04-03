#' Plot the STORM results
#'
#' Plot the storm results (s) from storm_play, etc.
#' Plot light gray ref region for optimal subset.
#' use in lapply().
#' Pass in s:
#'   nothing outside
#'
#' @param xmat numeric matrix of spectral data
#' @param ppm numeric vector of length equal to ncol(xmat)
#' @param s STORM results from storm_play or similar function
#' @param bgplot a character vector indicating how the plot is rendered:
#'        "overlayed" or "stack". Default is "overlayed".
#' @param vshift numeric scalar indicating the vertical shift between spectra when bgplot is "stack". Default is 1.
#' @param hshift numeric scalar indicating the horizontal shift between spectra when bgplot is "stack". Default is 0.05.
#' @param calcStocsy a logical indicating whether to calculate STOCSY plot. Default is TRUE.
#' @param n_xticks an integer indicating the number of ticks on the x-axis.
#'        Default is 4.
#'
#' @return a ggplot object with the plot of the spectral data, the reference region, and the STOCSY plot (if applicable).
#'
#' @importFrom ggplot2 geom_path geom_vline scale_colour_gradientn
#' @importFrom ggnewscale new_scale_color
#' @importFrom gridExtra grid.arrange
#' @importFrom matrixStats rowSds
#' @importFrom colorRamps matlab.like2
#'
#'
#' @export
#'
#' @examples
#' plot_storm_refRegions(xmat = xmat, ppm = ppm, s = stormResults[[1]])
#'
#' @seealso \code{\link{plot_stormRefRegions_grid}}, \code{\link{plot_stormRefRegions_summarize}}
plot_storm_refRegions <- function(xmat, ppm, s, bgplot = "overlayed", vshift = 1, hshift = 0.05, calcStocsy = TRUE, n_xticks = 4) {
  # Plot the storm results (s) from storm_play, etc.
  # Plot light gray ref region for optimal subset.
  # use in lapply().
  # Pass in s:
  #   nothing outside
  #
  # MTJ 2022

  specRegion <- xmat[
    s$subset,
    s$finalRegion
  ]
  ppmRegion <- ppm[s$finalRegion]

  if (bgplot == "overlayed") {
    g <- simplePlot(specRegion, ppmRegion, n_xticks = n_xticks)
  }

  if (bgplot == "stack") {
    g <- stackplot(specRegion, ppmRegion, vshift = vshift, hshift = hshift)
  }


  ##############
  # Plot the reference in black, bold, floating in middle of spectra
  # Scale the profile to the spectra (average ppm-wise mean matching)
  refInReg <- s$finalRegion %in% s$ref.idx %>% which()
  ref <- s$ref.vals * mean(matrixStats::rowSds(specRegion[, refInReg])) / matrixStats::rowSds(t(s$ref.vals)) # normalize
  ref <- ref + (median(t(specRegion[, refInReg])) - median(ref)) # median-center

  # Make a copy of the ref over the plot (final) region with NAs in the gap
  stretchedRef <- rep(NA, length(s$finalRegion))
  stretchedRef[refInReg] <- ref

  spr <- t(specRegion[, refInReg])

  # New scaling (as of 26AUG22)
  spr <- spr - min(spr)
  stretchedRef <- stretchedRef - min(stretchedRef, na.rm = TRUE)
  stretchedRef <- stretchedRef * max(spr) / max(stretchedRef, na.rm = TRUE) + min(specRegion)
  stretchedRef <- stretchedRef * 0.5


  # Make the plot
  g <- g + geom_path(
    data = data.frame(
      refvals = stretchedRef,
      refppms = ppm[s$finalRegion]
    ),
    mapping = aes(x = refppms, y = refvals),
    colour = "black", linewidth = .5,
    lineend = "round",
    linejoin = "round",
    linemitre = "5"
  )


  ##############

  if (calcStocsy) {
    ce <- corr_expand_window(
      xmat = specRegion, ppm = seq_along(s$finalRegion),
      peakInd = s$finalRegion %in% s$peak %>% which(),
      wind = seq_along(s$finalRegion)
    )
    corr <- ce$sr@r
    stocsyLine <- ce$sr@cov
    bounds <- c(
      ce$corrLbound,
      ce$corrRbound
    )
  } else {
    # Use corr/cov profile in s

    # Add ref sections drivers and bounds for each corr peak
    mask <- (refInReg %>% range() %>% fillbetween()) %in% refInReg
    refPeaks <- extractPeaks_corr(s$corr[refInReg %>%
      range() %>%
      fillbetween()],
    mask = mask,
    plots = FALSE
    )

    # Find the mask bounds (ends inclusive)
    bounds <- c(
      mask %>% which() %>% min(),
      refPeaks$maskBounds,
      mask %>% which() %>% max()
    ) %>% unique()

    # Finalregion is sometimes smaller than corr. All corrs are reported from
    # the last window in storm, but finalregion is the filled inds for
    # range(ref.idx). Make sure we're only using the corrs from within finalregion

    scorr <- s$corr[mask %>%
      which() %>%
      range() %>%
      fillbetween()]
    scov <- s$cov[mask %>%
      which() %>%
      range() %>%
      fillbetween()]

    # Any gaps should be filled with NAs, just like the ref.
    corr <- rep(NA, length(scorr))
    stocsyLine <- corr
    corr[refInReg] <- s$corr[refInReg]
    stocsyLine[refInReg] <- s$covar[refInReg]
  }
  spr <- t(specRegion[, refInReg])

  # New scaling (as of 26AUG22)
  spr <- spr - min(spr)
  stocsyLine <- stocsyLine - min(stocsyLine, na.rm = TRUE)
  stocsyLine <- stocsyLine * max(spr) / max(stocsyLine, na.rm = TRUE) + min(specRegion)


  # Plot the stocsy line over the current plot


  g <- g + ggnewscale::new_scale_color() +
    geom_line(
      data = data.frame(
        refvals = stocsyLine,
        refppms = ppm[s$finalRegion],
        corr = corr
      ),
      mapping = aes(x = refppms, y = refvals, colour = corr),
      linewidth = 1.25,
      lineend = "round",
      linejoin = "round",
      linemitre = "5"
    ) +
    scale_colour_gradientn(
      colours = colorRamps::matlab.like2(10),
      limits = c(-1, 1)
    )

  # Plot corr boundaries as vert lines
  g <- g +
    geom_vline(xintercept = ppm[s$peak], linetype = 2, col = "grey")

  return(g)
}
