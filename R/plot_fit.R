#' Plot the spectral and feature fits
#'
#' @param fit An object of class "fit"
#' @param type A character string specifying the type of plot. Possible values are "simple", "auc", and "color.line".
#' @param ppm A numeric vector containing the ppm values.
#' @param color A numeric vector containing the color values.
#'
#' @return A ggplot object
#'
#' @examples
#' plot_fit(fit, type = "simple", ppm = ppm, color = color)
#'
#' @importFrom ggnewscale new_scale_color
#' @importFrom colorRamps matlab.like2
#'
#'
#' @export
plot_fit <- function(fit, type = "simple", 
                     ppm = NULL,
                     color = NULL){
  
  # Note: need to pass a continuous ppm vector/xaxis. If there are NAs, this function
  # will make strange plots (later, interpolate a linear scale between the provided points).
  
  if (is.null(ppm)){ppm <- which(!is.na(fit$spec.fit))}
  # if (any(is.na(ppm))){complete_indsVect(ppm)}
    
  if (type == "simple"){
    
      g <- simplePlot(fit$spec.fit,
                      xvect = ppm,
                      linecolor = "gray",
                      opacity = .9, 
                      linewidth = 1) # ,-res$residuals

      g <- g + new_scale_color() +
              ggplot2::geom_line(data = data.frame(vals = fit$feat.fit,
                                                   ppm = ppm,
                                                   corr = rep(-1, length(ppm))),   # 
                        mapping = aes(x = ppm, y = vals, colour = corr),
                        linewidth = .5,
                        na.rm = TRUE) +
              scale_colour_gradientn(colours = matlab.like2(10),
                                     limits = c(-1, 1))
  }
  
  if (type == "auc"){
      ff <- fit$feat.fit
        ff[fit$fit.neg] <- 0
      g <- simplePlot(fit$spec.fit,
                      xvect = ppm,
                      linecolor = "black",
                      opacity = 1, 
                      linewidth = .75) # ,-res$residuals
      
      feature.chunks <- run.labels(!is.na(ff))
      
      for (i in unique(feature.chunks[feature.chunks>0])){
        # i <- feature.chunks[2]
        this.chunk <- feature.chunks == i
        g <- g + new_scale_color() +
              ggplot2::geom_area(data = data.frame(vals = ff[this.chunk],
                                                   ppms = ppm[this.chunk]),
                        na.rm = TRUE,
                        mapping = aes(x = ppms, y = vals),
                        fill="pink", alpha=0.6,
                        color = "black",    # line color
                        lwd = 0.25,    # line width
                        linetype = 1)
      }

  }
  if (type == "color.line"){
      
      g <- simplePlot(fit$spec.fit,
                      xvect = ppm,
                      linecolor = "gray",
                      opacity = .9, 
                      linewidth = 1.5) # ,-res$residuals
      
      # g <- g + geom_path(data = data.frame(vals = fit$spec.fit,
      #                                      ppms = ppm),
      #                    na.rm = TRUE,
      #                    mapping = aes(x = ppms, y = vals),
      #                    colour = "black", linewidth = .25)

      g <- g + new_scale_color() +
              ggplot2::geom_line(data = data.frame(vals = fit$feat.fit,
                                          ppm = ppm,
                                          corr = 1-color), 
                        mapping = aes(x = ppm, y = vals, colour = corr),
                        linewidth = 1,
                        na.rm = TRUE) +
              scale_colour_gradientn(colours = matlab.like2(10),
                                     limits = c(-1, 1))
    # plot(g)
  }
      
  return(g)
}