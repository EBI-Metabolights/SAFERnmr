# Derive protofeatures for a spectral matrix

# Follows:
# - setup
# - load data
# For surviving protofeatures:
# data.frame(primary.lower,
#            primary.upper,
#            secondary.lower,
#            secondary.upper,   # relative driver position
#            res.center,        # relative driver position
#            index = x$index,
#            row.names = NULL)
# Viewing protofeatures:
  # p <- protofeatures.split[[i]] 
  ## or
  # i <- i + 10
  # p <- protofeatures[i, ]
  # plot_protofeature(p, 301, ppm, xmat, bgplot='overlayed')

# protofeatures <- compute_protofeatures(pars, xmat, noise.width.multiple = 2, top.n.peaks = 5, n.cores = 6)

compute_protofeatures <- function(pars, data, noise.width.multiple = 2, top.n.peaks = 5, n.cores = 6){
  
  # Get correlation pockets (protofeatures)
  
  ################ Set up parameters ##################
    
    plot.location <- pars$dirs$temp
    
    # Corr Pocket Pairs 
  
      xmat <- data$xmat
      ppm <- data$ppm
      
      half.window <- (pars$corrpockets$half.window / data$digital.res) %>% ceiling  
          if (half.window > 1000){stop('Window size is too large. Please keep to < 1000 points (~ ', round(1000 * data$digital.res, 4),' ppm for this dataset).')}
                                  # Window for the sliding correlation calculation. 
                                  # This x 2 should capture any 2 adjacent resonances
                                  # in a multiplet. 
                                  # Provided in ppm, converted here to (column) elements
  
      noise.percentile <- pars$corrpockets$noise.percentile     
                                  # noise characterization... For every spectral
                                  # point, we calculate a correlation peak. If you
                                  # average all of the resulting peak shapes, 99%
                                  # of them will be at least n points wide. All of 
                                  # the peaks have a width of at least 3 (one point
                                  # on either side of the driver, due to the way
                                  # a peak is defined). What fraction contain 5
                                  # points? 30 points? This parameter determines the 
                                  # fraction of all spectral corr peaks that noise
                                  # should fit within, and sets the noise width 
                                  # accordingly. Higher is more permissive.
                                  
       cpp.rcutoff <- pars$corrpockets$rcutoff
                                  # correlation cutoff for picking highest secondary
                                  # peak in corrpocket pair extraction. Generally
                                  # ~ 0.75 should do fine.
         
################ Use corrPocketPairs to extract likely j-pairs ##################
 
    # Run corrpocketPairs on everything
      
      pocketPairs <- correlation_pocket_pairs(xmat, ppm, ws = half.window, plotHeatmap = FALSE,
                                        wdlimit = noise.percentile, # **** **** **** #
                                        noise.width.multiple = 2,
                                        top.n.peaks = 5,
                                        rcutoff = cpp.rcutoff, n.cores = 6)
      
      pocketPairs %>% debug_write("pocketPairs.RDS", pars)
      # pocketPairs <- readRDS(paste0(pars$dirs$temp, "/debug_extra.outputs", "/pocketPairs.RDS"))

      
    # Clean up the result
      # Remove nulls
      pocketPairs$peakBounds <- pocketPairs$peakBounds[!is.null(pocketPairs$peakBounds)]
      
      
    # Unlist into individual pairs
      message('Compiling protofeature table...')
      protofeatures <- mclapply(pocketPairs$peakBounds, function(x) {
        
          null.ones <- is.null(x$secondary)
          
          lapply(x$secondary[!null.ones], function(y) {
            
            if (!is.null(y)){
              data.frame(primary.lower = x$primary['lower'],
                         primary.upper = x$primary['upper'],
                         secondary.lower = y[['lower']],
                         secondary.upper = y[['upper']],
                         res.center = x$res.center, 
                         index = x$index,
                         row.names = NULL)
            }
            else {
              NULL
            }
          }) %>% do.call(rbind,.)

      }, mc.cores = 6) %>% do.call(rbind,.)
    
    # Why are the primary peaks not lower?
      
    # Convert to relative inds
      drivers <- protofeatures$index
      protofeatures <- protofeatures$res.center - protofeatures
      protofeatures$driver <- drivers
      protofeatures$index <- NULL
      protofeatures$res.center <- NULL
      noisewidth <- sum(pocketPairs$noiseDist >= noise.percentile)
      protofeatures$index <- seq_along(1:nrow(protofeatures))
      
    # Expand protofeature
    
          # i <- 0
          # 
          # i <- i + 1
          # p <- protofeatures[i,]
          # i
          # plot_protofeature(p, 200, ppm, xmat, bgplot='stack', line.shape = 'covar', line.color = 'corr')
          # 
          # debug_write(protofeatures, 'protofeatures.RDS', pars)
          # 
    # Report number of pairs
      
      numPairs <- nrow(protofeatures)
      window.index <- (-half.window):(half.window)
      
      pdf(file = paste0(plot.location, "corrpeak_distribution.pdf"),   # The directory you want to save the file in
          width = 4, # The width of the plot in inches
          height = 4) # The height of the plot in inches
      
        pocketPairs$noiseDist %>% plot(x = window.index, ylab="Fraction of peaks including index", xlab="Window index")
        abline(h = noise.percentile, col="red")
        title(ylab = "", main = "Average Extracted Diagonal Peak Shape (pre-filtering)")
      
      dev.off()
      # here's a thought: if you filter all peaks based on the noise feature shape,
      # the relative prominence of true signal using those boundaries is going to be
      # minimal because signal is locally pretty flat, while noise will mostly be
      # captured within those bounds. Peaks could be classified based on the % of their
      # actual signal (using their actual bounds) captured by the n% cutoff bounds.

      # ####
      message('compute_protofeatures finished.')
  return(list(table = protofeatures,
              noiseWidth = noisewidth,
              noise.width.multiple = noise.width.multiple))
}

# Protofeatures are just a sketch of a correlation signature fragment
# The actual shapes, xmat segment,nad 

expand_protofeature <- function(p, xmat, ppm, half.window){
    
    driver <- p$driver

  # Driver locates the index, everything else can be built around it
    
    p.abs <- driver - p
    
    fullView <- (driver - half.window):(driver + half.window)
    
    in.bounds <- !(fullView < 1 | fullView > length(ppm))
    
    specreg.inds <- fullView[in.bounds]
    
    specRegion = xmat[,
                      specreg.inds]
    
    ppmRegion = ppm[specreg.inds]
    
      # simplePlot(specRegion, xvect = ppmRegion)
      
  # Recalculate cov and corr
  
    cv <- cr <- rep(NA, length(fullView))
    cv[in.bounds] <- cov(xmat[,driver], specRegion)
    cr[in.bounds] <- cor(xmat[,driver], specRegion)
    
  # Apply peaks
  
    # If this isn't a protofeature, but just a driver: let "peak" bounds = spec Region
    if (is.null(p.abs$primary.lower)){
      p.abs$primary.lower <- min(specreg.inds)
      p.abs$primary.upper <- max(specreg.inds)
      p.abs$secondary.lower <- p.abs$primary.lower
      p.abs$secondary.upper <- p.abs$primary.upper
    }
    # browser()
    primary.bounds <- c(p.abs$primary.lower, p.abs$primary.upper) %>% sort
    secondary.bounds <- c(p.abs$secondary.lower, p.abs$secondary.upper) %>% sort
    
    
    pk.mask <- rep(0, length(fullView))
    pk.mask[which(fullView %in% primary.bounds) %>% fillbetween] <- 1
    pk.mask[which(fullView %in% secondary.bounds) %>% fillbetween] <- 2
      # plot(specreg.inds, cr)
      # abline(v = c(driver), col='red')
      # abline(v = c(primary.bounds, secondary.bounds))
      # par(new=TRUE)
      # plot(specreg.inds, pk.mask)
    # peak.inds <- specreg.inds[pexp$peak.mask > 0]
    
    
    return(list(driver = driver,
                cv = cv, # covariance between driver and specRegion.inds - mask with [peak.mask>0] for ref.idx in log_storm
                cr = cr, # correlation between driver and specRegion.inds
                specRegion.inds = specreg.inds, # ppm indices for region of interest - mask with [peak.mask>0] for ref.idx in log_storm
                specRegion = specRegion, # xmat[, specRegion.inds] 
                ppmRegion = ppmRegion, # ppm vals in specRegion.inds
                peak.mask = pk.mask, # labels peaks 1, 2 >0 in specRegion.inds
                primary.bounds = primary.bounds, # ppm indices, bounds of peak
                secondary.bounds = secondary.bounds)) # ppm indices, bounds of peak
}

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
  

  # Expand the protofeature
    pexp <- expand_protofeature(p, xmat, ppm, half.window)
    
  # Set up line shape and colors - decide which line to plot
  
    shape <- switch(line.shape,
                    covar = pexp$cv,
                    corr = pexp$cr)
    
    n.colors <- 10
    cmap <- matlab.like2(n.colors)
    darkRed <- cmap[n.colors]
  
  # Decide how the line should be colored
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
  scale_x_reverse(breaks=breaks_pretty(), limits = xlim_fixed, expand = c(0, 0), oob = scales::oob_keep) +
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
    
  # Set tick labels 
    g2 <- g2 + theme(
                      axis.text.x = element_text(size = 12),
                      axis.ticks.x = element_line(),
                      axis.title.x = element_blank(),
                      plot.margin = margin(0, 0, 0, 0)
                    )
    
  if (!is.null(ref.mask) & show.mask.bounds){
    mask <- pexp$specRegion.inds %in% ref.mask
    diffMask <- mask %>% diff
        maskBounds <- diffMask %>% "!="(.,0) %>% which 
          maskBounds <- pexp$specRegion.inds[maskBounds + as.integer(diffMask[maskBounds] > 0)]
    
    g2 <- g2 + geom_vline(xintercept = ppm[maskBounds], linetype = 1, col = "gray")
  }
  
  # Align x axes manually
    # library(grid)
    
    # Align widths
    # max_widths <- grid::unit.pmax(g1$widths, g2$widths)
    # g1grob$widths <- max_widths
    # g2grob$widths <- max_widths
    # 
    # # New page and draw manually
    # grid.newpage()
    # grid.draw(rbind(g1grob, g2grob, size = "last"))
    # 
  # # 3. Stack them vertically
  combined_plot <- g1 / g2 + plot_layout(ncol = 1, heights = c(5, 1))  # Adjust heights if needed
  #   combined_plot <- g1 / g2 +
  #     plot_layout(heights = c(5, 1)) &
  #     theme(
  #       plot.margin = margin(0, 0, 0, 0),
  #       panel.spacing = unit(0, "pt")
  #     )

  
  # 4. Display
  return(combined_plot)
}

# profile <- switch(line.shape,
#                 covar = pexp$cv,
#                 corr = pexp$cr)

profile <- pexp$cv

# color.vect <- switch(line.color,
#                      covar = pexp$cv,
#                      corr = pexp$cr)

color.vect <- pexp$cr

ppms <- pexp$ppmRegion


# 1. Original plot

# basePlot <- switch(bgplot,
#                    overlay = simplePlot(pexp$specRegion, pexp$ppmRegion, n_xticks = 5),
#                    stack = stackplot(pexp$specRegion, pexp$ppmRegion, vshift = 10, hshift = 0))

basePlot <- simplePlot(pexp$specRegion, pexp$ppmRegion, n_xticks = 5)

# What I'm building is a suite of plotting functions which 
# 1) produce a base plot (e.g. stackplot or overlay)
# 2) overlay SAT-related features onto the plot (e.g. correlation lines, peak bounds, etc.)
  
# accessory function: add peak bounds
  # ppm[pexp$primary.bounds]
  # 
  # if (peakBounds){
  # # Add the peak bounds
  #   basePlot <- basePlot + 
  #     geom_vline(xintercept = , linetype = 2, col = "black") +
  #     geom_vline(xintercept = ppm[pexp$secondary.bounds], linetype = 2, col = "black")
  # }
  # g2...
  # if (!is.null(ref.mask) & show.mask.bounds){
  #   mask <- pexp$specRegion.inds %in% ref.mask
  #   diffMask <- mask %>% diff
  #       maskBounds <- diffMask %>% "!="(.,0) %>% which 
  #         maskBounds <- pexp$specRegion.inds[maskBounds + as.integer(diffMask[maskBounds] > 0)]
  #   
  #   g2 <- g2 + geom_vline(xintercept = ppm[maskBounds], linetype = 1, col = "gray")
  
# accessory function: add driver
  # basePlot <- basePlot + geom_vline(xintercept = ppm[pexp$driver], linetype = 2, col = darkRed)
  # 
  # basePlot <- basePlot +
  #       theme(
  #         axis.title.x = element_blank(),
  #         axis.text.x  = element_blank(),
  #         axis.ticks.x = element_blank()
  #       )
  # g2...
  # # Add driver
  #   
  #   g2 <- g2 + geom_vline(xintercept = ppm[pexp$driver], linetype = 2, col = darkRed)

             
plot_addProfile <- function(profile, color.vect, ppms, basePlot){
  
  # Set up line shape and colors - decide which line to plot
  
    n.colors <- 10
    cmap <- matlab.like2(n.colors)
    darkRed <- cmap[n.colors]
  
  # Decide how the line should be colored
  
      # if not dealing with correlations, make sure colors are mapped to range instead of [-1, 1]
      if (any(color.vect < -1 | color.vect > 1)){
        cvals.range <- range(color.vect)
      } else {
        cvals.range <- c(-1,1)
      }
  
  # Fix x-axis limits for both plots to ensure alignment
      xlim_fixed <- range(ppms, na.rm = TRUE) %>% rev
    
    
  # 2. Correlation-colored plot
  df <- data.frame(
    ppms = ppms,
    shape = profile,
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
    scale_x_reverse(breaks=breaks_pretty(), limits = xlim_fixed, expand = c(0, 0), oob = scales::oob_keep) +
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
  
    
  # Set tick labels 
    g2 <- g2 + theme(
                      axis.text.x = element_text(size = 12),
                      axis.ticks.x = element_line(),
                      axis.title.x = element_blank(),
                      plot.margin = margin(0, 0, 0, 0)
                    )
    
    g1 + g2
  
  # # # 3. Stack them vertically
  # combined_plot <- g1 / g2 + plot_layout(ncol = 1, heights = c(5, 1))  # Adjust heights if needed
  # 
  # # 4. Display
  # return(combined_plot)
}

plot_addProfile(profile, color.vect, ppms, basePlot)

