#' Plot stacked spectra with ggplot2 and ggridges
#'
#' This function plots stacked spectra using ggplot2 and ggridges, which
#' creates a density plot of the spectra that can help to visualize overlaps
#' and peaks.
#'
#' @param ymat A numeric matrix of spectral intensities.
#' @param xvect A numeric vector of ppm values.
#' @param vshift A numeric value for shifting spectra vertically (spectra are scaled down to this value).
#' @param hshift A numeric value for shifting spectra horizontally (as % of actual range).
#' @param shift_by A numeric vector for shifting spectra by row number.
#' @param xdir A character string for x-axis direction, either "reverse" or "forward".
#' @param show A logical value for whether or not to display the plot.
#'
#' @import ggplot2
#' @import ggridges
#' @import dplyr
#'
#' @return A ggplot2 object.
#'
#' @examples
#' data(metab)
#' stackplot(metab[, 1:10], xvect = metab[, 11])
#'
#' @export
stackplot <- function(ymat = NULL, xvect = NULL,
                      vshift = 1, hshift = 0,
                      shift_by = NULL,
                      xdir = "reverse",
                      show = TRUE) {
  # Requires
  require(ggplot2)
  require(ggridges)
  require(dplyr)
  # Handle single vector case:
  if (is.vector(ymat)) {
    ymat <- t(matrix(ymat))
  }

  # Shift by row number
  if (is.null(shift_by)) {
    shift_by <- 1:nrow(ymat)
  }

  # Default create xvect
  if (is.null(xvect)) {
    xvect <- 1:ncol(ymat)
  }

  # Modify the intensities using vshift
  #   Turns out, the yshift per se is set. So what this is actually doing is adjusting
  #   the effective spacing between spectra by just making the intensities lower
  #   for each spectrum. No vshift needed.

  # Calculate actual hshift value to be fraction of on ppm range
  actualH <- xvect %>%
    range() %>%
    diff() %>%
    abs() %>%
    "*"(hshift)


  # Melt into df for plotting in ggplot
  df <- as.data.frame(t(ymat))
  colnames(df) <- 1:ncol(df)
  df$ppm <- xvect

  d <- reshape2::melt(df, id.vars = "ppm")
  colnames(d) <- c("ppm", "specNumber", "Spectral Intensity")
  d$color <- rep(alpha("gray", 0.6))
  d$shiftBy <- NA
  d <- d %>% filter(!is.na(`Spectral Intensity`))

  # Scale the spectra down to range of shift_by
  d$`Spectral Intensity` <- d$`Spectral Intensity` %>% scale.between(.,
    lower = 0,
    upper = (max(shift_by) - min(shift_by)) / vshift
  )

  # Calculate horizontal shifts for ppm and label the shift-by column
  d$shiftBy <- shift_by[d$specNumber[1:nrow(d)]]
  d$ppm <- d$ppm - shift_by[d$specNumber[1:nrow(d)]] * actualH


  ##################################################################################
  library(ggridges)
  g <- ggplot(
    d,
    aes(
      x = `ppm`,
      y = `shiftBy`,
      height = `Spectral Intensity`,
      group = `shiftBy`
    )
  ) +
    ggridges::geom_ridgeline(na.rm = TRUE)

  if (xdir == "reverse") {
    g <- g + scale_x_reverse(breaks = breaks_pretty())
  } else {
    g <- g + scale_x_continuous(breaks = breaks_pretty())
  }

  g <- g +
    # scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_ridges(
      grid = FALSE,
      center_axis_labels = TRUE
    ) +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none"
    )

  if (show) {
    plot(g)
  }
  return(g)
}
