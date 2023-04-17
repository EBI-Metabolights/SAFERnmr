#' Plot a correlation window with local minima
#'
#'
#' @param sr An object of class Spectrum containing correlation information
#' @param lbound An integer specifying the left boundary of the correlation window to be plotted
#' @param rbound An integer specifying the right boundary of the correlation window to be plotted
#' @param showplot A logical indicating whether the plot should be displayed on the screen (default = TRUE)
#' @param restrictTo A numeric vector of length two specifying the limits of the x-axis in the plot
#'
#' @return A ggplot2 object of the correlation window with local minima highlighted
#' @export
#'
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom scales breaks_pretty
#' @importFrom utils str_c
#' @importFrom colorRamps matlab.like2
plot_stocsy_corrWindow <- function(sr, lbound, rbound, showplot = TRUE,
                                   restrictTo = NULL) {
  # Modified from metabom8 stocsy function to show local minima
  # Intended use is to visualize the results from corr_expand()
  # lbound and rbound are indices
  # MTJ 2022

  if (is.null(restrictTo)) {
    restrictTo <- c(1, length(sr@ppm))
  }
  plotregion <- restrictTo[1]:restrictTo[2]

  df <- data.frame(cc = sr@r[plotregion], cv = sr@cov[plotregion], ppm = sr@ppm[plotregion])

  csc_lab <- paste("r(d=", sr@driver, ", X)", sep = "")
  cap <- paste0("Sample size: n=", nrow(sr@X))

  g1 <- ggplot(df, aes_string(
    x = "ppm",
    y = "cv",
    colour = "abs(cc)"
  )) +
    geom_line() +
    scale_x_reverse(breaks = breaks_pretty()) +
    scale_colour_gradientn(
      colours = colorRamps::matlab.like2(10),
      limits = c(0, 1),
      name = csc_lab
    ) +
    labs(
      x = expression(delta ~ {}^1 * H ~ (ppm)),
      y = gsub("^r", "cov", csc_lab),
      title = str_c("Peak @ ", sr@driver),
      caption = cap
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(colour = "black"),
      panel.grid.minor.x = element_blank()
    )

  g1 <- g1 + geom_vline(xintercept = sr@driver, linetype = 2, col = "grey")
  g1 <- g1 + geom_vline(xintercept = c(lbound, rbound), linetype = 1, col = "grey")

  if (showplot) {
    plot(g1)
  }

  return(g1)
}
