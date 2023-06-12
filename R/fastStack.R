#' Fast alternative to stackplot(), geom_area-based plots, or the like.
#' Uses scattermore or ggplot2, but pre-computes points covered by lower spectra
#' (i.e. lower row #; foreground) and fills with NA to give the appearance of 
#' being layered. 
#'
#' Row 1 is the bottom spectrum in the plot. Note: x is currently  minmax scaled [0,1]. 
#'
#' @param x spectral matrix
#' @param ppm ppm vector
#' @param raster use scattermore? default FALSE gives geom_line based plot. Both return ggplot obj.
#' @param vshift how much to shift each consecutive row upwards (% of standard deviation for x)
#' @param pixels vector of (rows, cols) - used by scattermore (raster = T)
#' @param pointsize ~ linewidth - used by scattermore (raster = T). 0 is fastest (1 pixel)
#' @param interpolate - used by scattermore (raster = T). Nicer result if (default) TRUE.
#' 
#' @return ggplot object; stackplot of spectra
#' @importFrom magrittr %>%
#' @importFrom scales breaks_pretty
#' @import ggplot2
#' @importFrom scattermore geom_scattermost
#' 
#' @export 
fastStack <- function(x, ppm, 
                      raster = F, vshift = 1, pixels = c(512, 512), pointsize = 0, interpolate = T){
    vshift <- vshift * sd(x, na.rm = T)
    # x <- scale.between(x) # catch coefficients? to apply to features?
    xs <- apply(x, 2, function(r) r + (1:nrow(x)) * vshift)
   
    xs <- rm.covered.points(xs)
    
    # Raster or no?
      
      if (raster){
        # xs <- scale.between(xs, 1, nrow(mat))
  
        df.lines <- data.frame(ppm = ppm, int = floor(xs) %>% t %>% c) %>% na.omit
  
        g <- ggplot() +
          ggplot2::scale_x_reverse(breaks = scales::breaks_pretty()) + 
          scattermore::geom_scattermost(xy = df.lines,
                                        interpolate = interpolate,
                                        pointsize = pointsize,
                                        pixels = pixels) +
          theme_clean_nmr() + 
          ggplot2::scale_x_reverse(breaks = scales::breaks_pretty())
        
      } else {
        
        g <- simplePlot(xs, xvect = ppm)
        
      }
    
    g
  
}


# ####

# This is a much, much faster version of project_features.stackplot() that offers scattermore-driven
# rastered plotting, as well as vectorized vector-graphic ggplots. Still working on plotting AUC for 
# features, but this is pretty good for now. 
# 
fastStack.withFeatures <- function(xmat, ppm,
                                   raster = T, 
                                   bfs, plt.pars){
  
   
   
# Put features in an x-sized matrix:
  # For each spectrum, put the features in:
    
    # Expand the ranges (better vis, but keep in bounds of ppm axis): ####
      # exp.by <- ceiling(plt.pars$exp.by / ((ppm %>% range %>% diff) / length(ppm)))
      # exp.ranges <- apply(bfs$fit.positions, 1, function(x){x %>% range(na.rm = TRUE) %>% sort})
      # if(length(exp.ranges) == 0){return(NULL)}
      #   lbound <- 1
      #   ubound <- length(ppm)
      #   exp.ranges[1, ] <- apply(exp.ranges[1,,drop=F] - exp.by, 2, function(x) max(c(x,lbound)))
      #   exp.ranges[2, ] <- apply(exp.ranges[2,,drop=F] + exp.by, 2, function(x) min(c(x,ubound)))
      
    # Get the xmat cols of interest ####
      # cols.x <- exp.ranges %>% range(na.rm = T) %>% fillbetween
      # message(paste(plt.pars$xlim, sep = " "))
      cols.x <- plt.pars$xlim %>% vectInds(ppm) %>% fillbetween
      ss.rows <- bfs$fit.xrow
      x.rows <- unique(ss.rows)
    
    # If simply getting a row for each fit. ####
      f.stack <- lapply(1:length(ss.rows), function(from.ss) {
        # from.ss <- ss.rows[1]

          # Line the rfs up with cols.x
            pos <-  bfs$fit.positions[from.ss, ]
            val <-  bfs$fit.feats[from.ss, ]

            # Make the matrix
            ss.vals <- rep(NA, length(cols.x))

            # Get the positions within the matrix
              # Since the introduction of xlim, these do not match 
              # however:
              #   we can convert the pos vector to pos within cols.x
              #   but: need to check to make sure these inds aren't 
              #   negative or > length(cols.x)
              # * inds, val, and pos are size(feature)
              # * ss.vals can be a different size.
              # * make indices of size(pos) that shoot values into the correct
              #   places in ss.vals. 
              
              inds <- (pos - min(cols.x) + 1)
                not.in.region <- inds < 1 | inds > length(cols.x)
                use <- !(not.in.region | is.na(val))
                
              ss.vals[inds[use]] <- val[use] # don't use NA elements of val

            return(ss.vals)

      }) %>% do.call(rbind,.)

    # Get x data: ####
        
        x <- xmat[x.rows, cols.x]
    
    # Apply the vshifts ####
      
        vshift <- plt.pars$vshift * sd(x, na.rm = T)
        vshifts <- (1:nrow(x)) * vshift # make vector of shifts, simple row # multiple of vshift.
        xs <- apply(x, 2, function(r) r + vshifts)
        # simplePlot(xs)
        
        # Apply the appropriate shift to each bf (according to the row it matches in xs)
          # they're all currently at 0. Not necessarily ordered. But all belonging 
          # to row 3 in xmat, for instance, will jump up that amount. Likewise with each
          # of the other rows in xs the bf.fits could belong to (and share vshift with).
          
        f.stack <- lapply(1:nrow(f.stack), function(r) {
          
          f.stack[r, ] + vshifts[x.rows == bfs$fit.xrow[r]]#-0.5*vshift
          
        }) %>% do.call(rbind,.)
        # simplePlot(f.stack)

    # Apply h shifts
    # ppm.mat = outer(ppm[cols.x], hshifts, "+")
        
    # Work out which points should be removed to give the impression of being covered up? ####
        #   - if value is < any other values lower down in that column, set to NA
        #   - i.e. remove points if they are exceeded by a lower row
        #   - for features, don't compute (if they overlap, we may want to show that)
        
        f.rows.in.x <- lapply(ss.rows, function(ssr) which(x.rows %in% ssr)) %>% unlist
        
        xs <- rm.covered.points(xs)
          # simplePlot(xs)
        # f.stack <- rm.covered.points(xs, f.stack, f.rows.in.x)
          # simplePlot(f.stack)
          # use this for features?
  
    # # Interleave xs and f.stack ####
    #   
    #   # Pad with row of zeros so the first row has something to subtract
    #     xs <- rbind(rep(0, ncol(xs)), xs)
    #   
    #   # For each xs, put xs, then all feats from that f.rows.in.x
    #   # - for the ribbons, compute the ymin as the row beneath it
    #     fullstack <- lapply(2:nrow(xs), function(row.xs) {
    #       
    #       
    #       same.spec <- f.rows.in.x == row.xs
    #       
    #       # If any features for this row of xs, return with labels. If not, return xs.row with label. 
    #       if (any(same.spec)){
    #         
    #         list(data = rbind(f.stack[same.spec, ], xs[row.xs, ]), # want features on top
    #              label = c(rep("feature", sum(same.spec)), "spec"),
    #              mins = repmat(xs[row.xs-1, ], sum(same.spec)+1, 1))
    #         
    #       } else {
    #         
    #         list(data = xs[row.xs, ],
    #              label = "spec",
    #              mins = xs[row.xs-1, ])
    #       }
    #       
    #     }) %>% unlist(recursive = F) %>% split(., names(.))        
    #   
    #   xs <- xs[-1, ]
    #   
    #   stack <- fullstack$data %>% do.call(rbind,. )
    #   mins <- fullstack$mins %>% do.call(rbind,. )
    #     # simplePlot(mins)
    #     # simplePlot(stack)
    #   labels <- fullstack$label %>% unlist(use.names = F)
    #   rownames(stack) <- make.names(labels, unique = TRUE)
    #   rownames(mins) <- rownames(stack)
      
  # Plotting ####
      if (raster){ 
        
          
        df.lines <- data.frame(ppm = ppm[cols.x],
                               int = xs %>% t %>% c,
                               color = alpha('black', alpha = 1)) %>% na.omit
        
        
        # g <- ggplot() +
        #   scattermore::geom_scattermost(xy = df.lines,
        #                                 interpolate = plt.pars$interpolate,
        #                                 pointsize = plt.pars$pointsize,
        #                                 pixels = plt.pars$pixels) +
        #   theme_clean_nmr()

         # Calculate feature fills
          
          # outlines
          df.feats <- data.frame(ppm = ppm[cols.x], 
                                 int = f.stack %>% t %>% c,
                                 color = alpha('blue', alpha = 5)) %>% na.omit
          # df.feats <- data.frame(ppm = df.feats$ppm, 
          #                            int = outer(df.feats$int, vshift*(seq(-10,10,by = 2))/10, '-') %>% c,
          #                            df.feats$color)
          
        # df <- rbind(df.lines, df.feats)
        # g <- g +
        #   scattermore::geom_scattermost(xy = df.feats,
        #                         xlim = c(max(df.feats[, 1]), min(df.feats[, 1])),
        #                         interpolate = plt.pars$interpolate,
        #                         pointsize = plt.pars$pointsize,
        #                         pixels = plt.pars$pixels,
        #                         color = 'blue') 
          # ggplot2::scale_x_reverse(breaks = scales::breaks_pretty())
        
          # yrange <- range(c(df.lines$int, df.feats$int))
          points.has.feature <- df.lines$ppm %in% df.feats$ppm
          yrange <- range(c(df.lines$int[points.has.feature], df.feats$int))
          
          xrange <- range(c(df.lines$ppm, df.feats$ppm)) %>% rev
          
          scattermore::scattermoreplot(
                                        x = df.lines$ppm,
                                        y = df.lines$int,
                                        xlab = 'ppm',
                                        ylab = '',
                                        size = plt.pars$pixels,
                                        cex = .0,
                                        ylim = yrange,
                                        xlim = xrange,
                                        col = 'black',
                                        yaxt="n"
                                      ) 
          par(new=TRUE)
          scattermore::scattermoreplot(
                                        x = df.feats$ppm,
                                        y = df.feats$int,
                                        xlab = 'ppm',
                                        ylab = '',
                                        size = plt.pars$pixels,
                                        cex = .1,
                                        ylim = yrange,
                                        xlim = xrange,
                                        col = alpha('blue', alpha = .2),
                                        yaxt="n"
                                      ) 
          
          # scattermore::scattermoreplot(
          #                           x = df$ppm,
          #                           y = df$int,
          #                           xlab = 'ppm',
          #                           ylab = '',
          #                           xlim = c(max(df$ppm), min(df$ppm)),
          #                           col = df$color
          #                         ) 
          
    
      } else {
        
        # Simpleplot approach:
        # - plot all as lines
        # - features in blue
          
          # Write color for each line, assuming all specs at first: ####
            spec <- list(line.color = 'black',
                         line.width = .5,
                         fill.color = 'white',
                         fill.opacity = 1
                         )
            feat <- list(line.color = 'blue',
                         line.width = 1,
                         fill.color = 'pink',
                         fill.opacity = .25
                         )

            line.colors <- rep(spec$line.color, length(labels))
            # fill.colors <- rep(alpha(spec$fill.color, spec$fill.opacity), length(labels))
            fill.colors <- rep(spec$fill.color, length(labels))
            fill.alphas <- rep(spec$fill.opacity, length(labels))
            line.widths <- rep(spec$line.width, length(labels))

            # then modify the feature colors:
              line.colors[labels == 'feature'] <- feat$line.color
              fill.colors[labels == 'feature'] <- feat$fill.color
              fill.alphas[labels == 'feature'] <- feat$fill.opacity
              line.widths[labels == 'feature'] <- feat$line.width
          

          # Make df, then melt for plotting ####
            df <- as.data.frame(t(stack))
            df$ppm <- ppm[cols.x]
            
            # df.mins <- as.data.frame(t(mins))
            # df.mins$ppm <- ppm[cols.x]
            
            d <- reshape2::melt(df, id.vars="ppm")
            # d.mins <- reshape2::melt(df.mins, id.vars = "ppm")
            colnames(d) <- c("ppm", "label", "Spectral Intensity")
            # d$type <- 'spec'
              feats <- grepl('feature', d$label)
              specs <- !feats
              # d$type[feats] <- 'feature'
              # d$ymin <- d.mins$value
              
          # Make plot in one shot: ####
            g <- ggplot2::ggplot(d) + 
                ggplot2::geom_line(mapping = aes(ppm, `Spectral Intensity`, color = label), 
                                   linetype = "solid",
                                   na.rm=TRUE, 
                                   linewidth = .5, 
                                   show.legend = F) + theme_clean_nmr() +
                ggplot2::scale_color_manual(values=line.colors) + 
                ggplot2::scale_linewidth_manual(values=line.widths) + 
                ggplot2::scale_alpha_manual(values=fill.alphas)
            g
            
            
          # # Ribbons/areas ####
          #     
          #  # ribbon min should be:  max(this.row-1, min(this.row, na.rm = T), na.rm = T)
          #   
           # g <- ggplot2::ggplot(d[specs, ]) +
           #        ggplot2::geom_area(mapping = aes(ppm, `Spectral Intensity`, color = label),
           #                           linetype = "solid",
           #                           na.rm=TRUE,
           #                           linewidth = .5,
           #                           show.legend = F) + theme_clean_nmr() +
           #        ggplot2::scale_color_manual(values = rep(line.spec, length( unique( d$label[specs]) ) ) )
        
          
          # add in features
           # ggplot2::ggplot(d[feats, ], aes(x = ppm, y = `Spectral Intensity`, color = label)) + 
                # geom_line(na.rm = T, show.legend = F)
              # g + 
              #   ggplot2::geom_ribbon(data = d[feats, ],
              #                        aes(x = ppm, ymin = ymin, ymax = `Spectral Intensity`), 
              #                        color = alpha('pink', alpha = 0.5),
              #                        na.rm = TRUE,
              #                        show.legend = F) +
              #    theme_clean_nmr()
            
            
      }

        # return(g)

        
}
  
      


          
        # # Overlay the fit feature AUC ####
        # 
        #   g <- g + geom_ribbon(data = data.frame(ymax = f.rev[s, ],
        #                                          ymin = min(x.rev[s, ], na.rm = T),
        #                                          ppms = f.ppm[s, ]),
        #            # na.rm = TRUE,
        #            mapping = aes(x = ppms, ymin = ymin, ymax = ymax),
        #            colour = "gray",
        #            linewidth = 0.5,
        #            fill="pink",alpha=0.4)

