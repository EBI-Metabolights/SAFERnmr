#' Plot correlation window for a given pocket
#'
#' @param vals a list containing values for the correlation window including sr, wind, corrRbound, and corrLbound
#'
#' @return a ggplot object displaying the correlation window
#'
#' @importFrom ggplot2 theme element_blank
#' @import plot_stocsy_corrWindow
#'
#'
#' @export
plot_correlation_pocket <- function(vals) {
  sr <- vals$sr
  wind <- vals$wind

  # Plotting Ripped and modified from stocsy.R (Torbin's code)

  span <- vals$corrRbound - vals$corrLbound + 1

  plotRegion <- c(vals$corrLbound - span, vals$corrRbound + span)

  g <- plot_stocsy_corrWindow(sr, wind[vals$corrRbound], wind[vals$corrLbound], showplot = FALSE, restrictTo = plotRegion)
  g <- g + theme(
    title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) # ,panel.border =  element_blank())

  return(g)
}
