library(wavelets)

xmat=data$xmat
ppm=data$ppm
pexp <- expand_protofeature(p, xmat, ppm, half.window)
driver <- pexp$driver

  plot_protofeature(p,
            half.window = half.window, ppm = data$ppm,
            xmat,
            bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
            showPeaks = TRUE)
  p$driver
  p.width <- 
  abs(p$primary.upper - p$primary.lower) + 
    abs(p$secondary.upper - p$secondary.lower)
  
wind <- pexp$specRegion.inds

specRegion = pexp$specRegion
ppmRegion = pexp$ppmRegion

y <- specRegion[1,]
x <- ppmRegion
simplePlot(y, ppmRegion)

library(plotly)
library(baseline)

plot_ALS_stack_plotly <- function(x, y, ppm, p,
                                  lambda_values = 10^seq(3, 8, by = 1),
                                  p_asym = 0.01,
                                  vspace = 0.001,
                                  show = c("baseline", "detrended", "both")) {
  
  show <- match.arg(show)
  
  # Clean y
  y_clean <- as.numeric(y)
  if (any(!is.finite(y_clean))) {
    warning("Replacing NA/NaN/Inf with 0")
    y_clean[!is.finite(y_clean)] <- 0
  }
  
  # ROI bounds
  window_rng <- p$driver + range(c(p$primary.lower, p$secondary.lower, p$primary.upper, p$secondary.upper))
  window_min_idx <- window_rng[1]
  window_max_idx <- window_rng[2]
  window_min_ppm <- ppm[window_min_idx]
  window_max_ppm <- ppm[window_max_idx]
  
  plt <- plot_ly()
  vshift <- max(y_clean) * vspace
  y_axis_max <- max(y_clean) + vshift * length(lambda_values)
  y_axis_min <- min(y_clean)
  
  for (i in seq_along(lambda_values)) {
    lambda_i <- lambda_values[i]
    y_offset <- vshift * (i - 1)
    
    # Baseline using 'baseline' package (ensure matrix input)
    bl <- baseline(matrix(y_clean, nrow = 1), method = "als",
                   lambda = lambda_i, p = p_asym, maxit = 100)
    
    baseline_est <- getBaseline(bl)[1, ]
    y_corrected <- getCorrected(bl)[1, ]
    
    # Plot original
    plt <- plt %>%
      add_trace(x = x, y = y_clean + y_offset, type = 'scatter', mode = 'lines',
                line = list(color = "gray", width = 1), showlegend = FALSE)
    
    if (show %in% c("baseline", "both")) {
      plt <- plt %>%
        add_trace(x = x, y = baseline_est + y_offset, type = 'scatter', mode = 'lines',
                  line = list(color = "red", width = 2, dash = "dot"), showlegend = FALSE)
    }
    
    if (show %in% c("detrended", "both")) {
      plt <- plt %>%
        add_trace(x = x, y = y_corrected + y_offset, type = 'scatter', mode = 'lines',
                  line = list(color = "black", width = 2), showlegend = FALSE)
    }
  }
  
  # ROI vertical lines
  plt <- plt %>%
    add_segments(x = window_min_ppm, xend = window_min_ppm,
                 y = y_axis_min, yend = y_axis_max,
                 line = list(color = "red", width = 2, dash = "dash")) %>%
    add_segments(x = window_max_ppm, xend = window_max_ppm,
                 y = y_axis_min, yend = y_axis_max,
                 line = list(color = "red", width = 2, dash = "dash"))
  
  plt <- plt %>%
    layout(title = paste("ALS Stackplot (show:", show, ")"),
           xaxis = list(title = "ppm"),
           yaxis = list(title = "Intensity (offset)", range = c(y_axis_min, y_axis_max)),
           showlegend = FALSE)
  
  return(plt)
}



# To explore baseline fits only:
plot_ALS_stack_plotly(x, y, ppm, p, show = "baseline")

# To see original + detrended:
plot_ALS_stack_plotly(x, y, ppm, p, show = "detrended")

# To see all three (original, baseline, detrended):
plot_ALS_stack_plotly(x, y, ppm, p, show = "both")
