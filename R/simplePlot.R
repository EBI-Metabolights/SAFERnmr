#' Create a simple plot using ggplot2
#'
#' This function creates a simple plot using ggplot2, based on a matrix or dataframe input.
#' It can handle both single vector and multiple vectors input, and can generate x-axis ticks
#' with a specified number of ticks or using a default breaks_pretty() function.
#'
#' @param ymat A matrix or dataframe, with each column being a spectrum to plot.
#' @param xvect A numeric vector for the x-axis values.
#' @param n_xticks An integer specifying the number of x-axis ticks to generate.
#' @param xdir A character string indicating the direction of the x-axis, either "reverse" or "forward".
#' @param linecolor A character string indicating the color of the plotted lines.
#' @param opacity A numeric value specifying the opacity of the plotted lines.
#' @param linewidth A numeric value specifying the width of the plotted lines.
#'
#' @return A ggplot object.
#'
#' @import ggplot2 scales
#'
#' @examples
#' data(mtcars)
#' plotmat <- scale(mtcars[, -1], center = TRUE, scale = TRUE)
#' simplePlot(plotmat, linecolor = "black", opacity = 0.6, linewidth = 0.5)
#'
#' simplePlot(plotmat[, 1], linecolor = "red", opacity = 0.6, linewidth = 0.5)
#'
#' simplePlot(plotmat, xvect = mtcars$wt, n_xticks = 5, linecolor = "black", opacity = 0.6, linewidth = 0.5)
#'
#' @export
simplePlot <- function(ymat = NULL, xvect = NULL, n_xticks = NULL, xdir = "reverse",
                       linecolor = "gray", opacity = 0.6, linewidth = 0.5){
  
  
  # Handle single vector case:
    
    if (is.vector(ymat))
      {ymat <- ymat %>% c %>% as.matrix %>% t}
  
  # Default create xvect
    if (is.null(xvect)){xvect <- 1:ncol(ymat)}
    # xvect <- range(xvect, na.rm = T) %>% fillbetween
  
  # Melt into df for plotting in ggplot
    df <- as.data.frame(t(ymat))
    colnames(df) <- 1:ncol(df)
    df$ppm <- xvect
    
    d <- reshape2::melt(df, id.vars="ppm")
    colnames(d) <- c("ppm", "specNumber", "Spectral Intensity")
    # d <- d %>% filter(!is.na(`Spectral Intensity`))
    d$color <- rep(alpha(linecolor, opacity))
  
  ##############    
  g <- ggplot2::ggplot(d) + 
    ggplot2::geom_line(linetype = "solid",
                       aes(ppm, `Spectral Intensity`, col = specNumber), 
                       na.rm=TRUE, linewidth = linewidth) + 
    ggplot2::scale_color_manual(values=d$color)
  
  if (is.null(n_xticks)){
    if (xdir == "reverse"){
      g <- g + ggplot2::scale_x_reverse(breaks = scales::breaks_pretty())
    } else {
      g <- g + ggplot2::scale_x_continuous(breaks = scales::breaks_pretty())
    }
  }else{
    if (xdir == "reverse"){
      g <- g + ggplot2::scale_x_reverse(breaks = scales::breaks_extended(n = n_xticks))
    } else {
      g <- g + ggplot2::scale_x_continuous(breaks = scales::breaks_extended(n = n_xticks))
    }
  }
    
  g <- g +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = element_text(colour = "black",size = 12), 
                  legend.position = "none",
                  axis.text.y = ggplot2::element_blank(),
                  # axis.title.x = element_text(size = 16,vjust = -0.5),
                  axis.title.x = ggplot2::element_blank(),
                  axis.title.y = ggplot2::element_blank(),
                  axis.ticks = ggplot2::element_blank(),
                  # axis.title = ggplot2::element_text(size = 12,vjust = 0.5),
                  panel.border =  ggplot2::element_blank(),
                  panel.grid.minor = ggplot2::element_blank(),
                  # panel.grid.major = ggplot2::element_line(color = "gray",
                  #                                         size = 0.1,
                  #                                         linetype = 1),
                  panel.grid.major = ggplot2::element_blank())
  return(g)
}