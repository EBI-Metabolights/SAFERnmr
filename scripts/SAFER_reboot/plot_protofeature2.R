library(ggplot2)
library(patchwork)

# input must be individual peaks (given in protofeatures table)
# recalculates the correlation and covar between driver and xmat for +/- half.window
# half.window could be larger than half.window used in protofeature, but should not be smaller or else peak bounds may not fit.
# If larger, only the peaks in protofeature are reported. 

  # p <- protofeatures.split[[i]] 
  ## or
  # p <- protofeatures[i, ]
  # plot_protofeature(p, half.window, ppm, xmat, bgplot='overlayed')

plot_protofeature <- function(p, half.window, ppm, xmat, bgplot='overlayed', line.shape='covar', line.color='corr', showPeaks=TRUE, ref.mask = NULL, show.mask.bounds=FALSE){
  

  # Decide if the 
    pexp <- expand_protofeature(p, xmat, ppm, half.window)
    
  # Set up line shape and colors
  
    shape <- switch(line.shape,
                    covar = pexp$cv,
                    corr = pexp$cr)
    
    n.colors <- 10
    cmap <- matlab.like2(n.colors)
    darkRed <- cmap[n.colors]
  
    color.vect <- switch(line.color,
                         covar = pexp$cv,
                         corr = pexp$cr)
  
      # if not dealing with correlations, make sure colors are mapped to range instead of [-1, 1]
      if (any(color.vect < -1 | color.vect > 1)){
        cvals.range <- range(color.vect)
      } else {
        cvals.range <- c(-1,1)
      }
  
  # Fix x-axis limits for both plots to ensure alignment
      xlim_fixed <- range(pexp$ppmRegion, na.rm = TRUE) %>% rev
    
  # 1. Original stackplot
  
    g1 <- switch(bgplot,
                 overlay = simplePlot(pexp$specRegion, pexp$ppmRegion, n_xticks = 5),
                 stack = stackplot(pexp$specRegion, pexp$ppmRegion, vshift = 10, hshift = 0))
    
  if (showPeaks){
  # Add the peak bounds
    g1 <- g1 + 
      geom_vline(xintercept = ppm[pexp$primary.bounds], linetype = 2, col = "black") +
      geom_vline(xintercept = ppm[pexp$secondary.bounds], linetype = 2, col = "black")
  }
    
  # add driver
  
    g1 <- g1 + geom_vline(xintercept = ppm[pexp$driver], linetype = 2, col = darkRed)
    
    g1 <- g1 + 
          theme(
            axis.title.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank()
          )
    
  # 2. Correlation-colored plot
  df <- data.frame(
    ppms = pexp$ppmRegion,
    shape = shape,
    color.vect = color.vect
  )
  
  # Set up to have gray where line != ref.mask points
  if (!is.null(ref.mask)) {
    # ref.mask <- pexp$specRegion.inds[pexp$peak.mask>0]
    ref.mask.region <- pexp$specRegion.inds %in% ref.mask
    df$final_color <- ifelse(ref.mask.region, df$color.vect, NA)
  } else {
    df$final_color <- df$color.vect
  }
  
df_lines <- df %>%
  mutate(
    xend = lead(ppms),
    yend = lead(shape),
    color_start = final_color,
    color_end = lead(final_color)
  ) %>%
  filter(!is.na(xend), !is.na(yend), !is.na(color_start), !is.na(color_end)) %>%
  mutate(color_avg = (color_start + color_end) / 2)

g2 <- ggplot(df_lines) +
  geom_segment(aes(x = ppms, xend = xend, y = shape, yend = yend, color = color_avg),
               linewidth = 1.25, lineend = "round") +
  scale_color_gradientn(colours = cmap, limits = cvals.range, na.value = "gray") +
  scale_x_reverse(limits = xlim_fixed, expand = c(0, 0), oob = scales::oob_keep) +
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black", size = 12),
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
  
  # Add driver
    
    g2 <- g2 + geom_vline(xintercept = ppm[pexp$driver], linetype = 2, col = darkRed)
    
  if (!is.null(ref.mask) & show.mask.bounds){
    mask <- pexp$specRegion.inds %in% ref.mask
    diffMask <- mask %>% diff
        maskBounds <- diffMask %>% "!="(.,0) %>% which 
          maskBounds <- pexp$specRegion.inds[maskBounds + as.integer(diffMask[maskBounds] > 0)]
    
    g2 <- g2 + geom_vline(xintercept = ppm[maskBounds], linetype = 1, col = "gray")
  }
  
  # 3. Stack them vertically
  # combined_plot <- g1 / g2 + plot_layout(ncol = 1, heights = c(5, 1))  # Adjust heights if needed
  combined_plot <- g1 / g2 + 
                    plot_layout(ncol = 1, heights = c(5, 1)) & 
                    theme(plot.margin = margin(0, 0, 0, 0), panel.spacing = unit(0, "pt"))

  
  # 4. Display
  return(combined_plot)
}
