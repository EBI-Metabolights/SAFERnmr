library(ggplot2)
library(patchwork)

# input must be individual peaks (given in protofeatures table)

  # p <- protofeatures.split[[i]]
  # plot_protofeature(p, ws, ppm, xmat, bgplot='overlayed')

plot_protofeature <- function(p, ws, ppm, xmat, bgplot='overlayed', line.shape='covar', line.color='corr'){
  
  driver <- p$index
  # peak.inds <- c(p$primary.lower:p$primary.upper, p$secondary.lower:p$secondary.upper)
  
  shape <- switch(line.shape,
                  covar = pocketPairs$cov[, driver],
                  corr = pocketPairs$corr[, driver])

  color.vect <- switch(line.color,
                       covar = pocketPairs$cov[, driver],
                       corr = pocketPairs$corr[, driver])
  # if not dealing with correlations, make sure colors are mapped to range instead of [-1, 1]
  if (any(color.vect < -1 | color.vect > 1)){
    cvals.range <- range(color.vect)
  } else {
    cvals.range <- c(-1,1)
  }


  fullView <- (driver - ws+1):(driver + ws-1)
  
  specreg.inds <- keep_inds_in_bounds(check = fullView, 
                      against = seq_along(ppm))
  
  # specRegion = matrix(NA, nrow=nrow(xmat), ncol=length(specreg.inds))
  
  outside_ppms <- fullView[!(fullView %in% specreg.inds)]
  
  if (length(outside_ppms)>0){
    
  }
  
  specRegion = xmat[,
                    specreg.inds]
  
  ppmRegion = ppm[specreg.inds]
  
  n.colors <- 10
  cmap <- matlab.like2(n.colors)
  darkRed <- cmap[n.colors]

  
  # 1. Original stackplot
  
    primary.bounds <- c(p$primary.lower,p$primary.upper)
    secondary.bounds <- c(p$secondary.lower,p$secondary.upper)
    
    g1 <- switch(bgplot,
                 overlayed = simplePlot(specRegion, ppmRegion, n_xticks = 5),
                 stack = stackplot(specRegion, ppmRegion, vshift = 10, hshift = 0))
    
  # Add the peak bounds
  
  g1 <- g1 + 
    geom_vline(xintercept = ppmRegion[secondary.bounds], linetype = 2, col = "black") +
    geom_vline(xintercept = ppmRegion[primary.bounds], linetype = 2, col = "black") +
    geom_vline(xintercept = ppm[driver], linetype = 2, col = darkRed)
    
  
  # 2. Correlation-colored plot
  df <- data.frame(
    ppms = ppmRegion,
    shape = shape,
    color.vect = color.vect
  )
  
  g2 <- ggplot(df, aes(x = ppms, y = shape, colour = color.vect)) +
    geom_line(linewidth = 1.25) +
    scale_colour_gradientn(colours = cmap, limits = cvals.range) +
    geom_vline(xintercept = ppmRegion[secondary.bounds], linetype = 2, col = "black") +
    geom_vline(xintercept = ppmRegion[primary.bounds], linetype = 2, col = "black") +
    geom_vline(xintercept = ppm[driver], linetype = 2, col = darkRed) + 
    scale_x_reverse() + 
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
  
  # 3. Stack them vertically
  combined_plot <- g1 / g2 + plot_layout(ncol = 1, heights = c(5, 1))  # Adjust heights if needed
  
  # 4. Display
  print(combined_plot)
}
