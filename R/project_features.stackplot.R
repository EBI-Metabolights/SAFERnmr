#' Generate a stacked area plot of the spectra for each feature in \code{xmat}, overlaid with the corresponding best-fit features
#'
#' @param xmat A matrix of spectral data with rows corresponding to spectra and columns to ppm values
#' @param ppm A vector of ppm values corresponding to the columns of \code{xmat}
#' @param label A character string indicating the feature label to be used in the plot title
#' @param bestfits A list containing the best fit data for each feature
#' @param exp.by The amount by which to expand the range of each feature's ppm range in the plot
#' @param vshift The amount by which to vertically shift each spectrum in the plot
#' @param hshift The amount by which to horizontally shift each spectrum in the plot
#' @param sort.rows A logical indicating whether or not to sort the rows of the spectra in each plot by their summed intensity
#'
#' @return A list containing the plot objects and parameters used to create the plot
#' 
#' @import ggplot2
#' @importFrom pracma repmat
#' @importFrom scales breaks_pretty
#' @importFrom reshape2 melt
#' 
#' @export
project_features.stackplot <- function(xmat, ppm,
                                       label,
                                       bestfits,
                                       exp.by = 0.05,
                                       vshift = 2,
                                       hshift = 0.002,
                                       sort.rows = F,
                                       plt.rng = 1:ncol(xmat)){
 
  require(ggplot2)
  require(pracma)
  require(scales)
  require(reshape2)
  
  # Only plot the fits that passed:
    bestfits$fit.feats <- bestfits$fit.feats[bestfits$pass.fit,, drop = F]
    bestfits$fit.positions <- bestfits$fit.positions[bestfits$pass.fit, , drop = F]
    bestfits$fit.xrow <- bestfits$fit.xrow[bestfits$pass.fit]

  # Limit the data to the relevant ppm range
    # Identify the ppm ranges from the spec.feature column data ####
        # ppmranges
       
        
  
        # exp.by <- 0.05 # ppm units on either side
        # vshift <- 2 # fraction of spectrum median
        # hshift <- 0.002 # ppm
        
        # Expand the ranges (better vis)
          
          exp.ranges <- apply(bestfits$fit.positions, 1, function(x){x %>% range(na.rm = TRUE) %>% ppm[.] %>% sort})
          if(length(exp.ranges) == 0){return(NULL)}
            lbound <- min(ppm)
            ubound <- max(ppm)
            exp.ranges[1, ] <- apply(exp.ranges[1,,drop=F] - exp.by, 2, function(x) max(c(x,lbound)))
            exp.ranges[2, ] <- apply(exp.ranges[2,,drop=F] + exp.by, 2, function(x) min(c(x,ubound)))
        
    # If any non-intersecting ppm ranges, then split into multiple plots
        # ntrsxns <- range.intersect.all(exp.ranges) 
        #   diag(ntrsxns) <- 0
        # v <- rowSums(ntrsxns) > 0
        # rl <- run.labels(v)
        #   rl <- rl - min(rl) + 1 # handles solo spec-feature case
        
        # plot.list <- unique(rl)
        
        plot.list <- range.groups(exp.ranges)
        nplots <- length(plot.list)
        
        # Compile info for each plot
          plots <- lapply(plot.list, function(p){
            list(rows = sort(p), # sort because max.cliques gives ~random order
                 range = exp.ranges[, p] %>% c %>% range)
          })
    
    # For each in plot.list (each ppm range)
      plots.all.ranges <- 
            lapply(1:nplots, function(p){
              # Cut the xmat to that region, but retain the ppms
              #   * ppms are the x-axis, not inds
                # p <- 1
                plotn <- plots[[p]]
                
                # Get the xmat rows for this plot 
                #   (i.e. those spec-features with this broader ppm range)
                #   - bestfits$fit.xrow is the xmat rows used for fitting. 
                #   - plotn$rows are the relative rows indices for the fit.feats matrix
                
                rowinds <- bestfits$fit.xrow[plotn$rows]
                colinds <- plotn$range %>% vectInds(., ppm) %>% fillbetween
                
                 # Make little x, shift vals ####
                  
                  if (length(rowinds) > 1){
                    x.rev <- xmat[rev(rowinds), colinds,drop = F]
                  } else {x.rev <- xmat[rowinds, colinds, drop = F]}

                  # Sort by row intensity, if desired
                    if (sort.rows){
                      neworder <- x.rev %>% rowSums(na.rm = TRUE) %>% order(decreasing = T)
                      # Correct the order of the rows for x, rowinds, plotn$rows
                      plotn$rows <- plotn$rows[neworder]
                      rowinds <- rowinds[neworder]
                      x.rev <- x.rev[neworder, ]
                    }
                  # Shift the values according to vshift ####
                    v.adjustments <- rev(1:nrow(x.rev)) * vshift * mean(x.rev, na.rm = T)
                    x.rev <- x.rev + v.adjustments
                      # simplePlot(x.rev, xvect = ppm[colinds])
                      
                  # Also shift the fit features ####
                    if (length(rowinds) > 1){
                      f.rev <- bestfits$fit.feats[rev(plotn$rows), ,drop = F] + v.adjustments
                    } else {f.rev <- bestfits$fit.feats + v.adjustments}

                      na.cols <- is.na(colSums(f.rev))
                      na.vals <- f.rev[, na.cols]
                      
                 # Shift the ppm vals according to hshift ####
                    h.adjustments <- rev(1:nrow(x.rev)) * hshift
                    x.ppm <- pracma::repmat(ppm[colinds], n = nrow(x.rev), m = 1) + h.adjustments
                   
                  # Also reverse the ppmranges for the fit features, and adjust them with hshift
                    if (length(rowinds) > 1){
                      f.ppm <- bestfits$fit.positions[rev(plotn$rows), ,drop = F] %>% 
                        apply(1, complete.indsVect) %>% t %>%
                        apply(., 1, function(x) ppm[x]) %>% t + h.adjustments
                    } else {
                      f.ppm <- bestfits$fit.positions[plotn$rows, ,drop = F] %>% 
                        apply(1, complete.indsVect) %>% t %>%
                        apply(., 1, function(x) ppm[x]) %>% t + h.adjustments
                    }
                      # na.ppms <- f.ppm[, na.cols]

                  # Melt xrev into df for plotting in ggplot ####
                    df <- as.data.frame(t(x.rev))
                    colnames(df) <- 1:ncol(df)
                    df$ppm <- colinds
                    
                    d <- reshape2::melt(df, id.vars="ppm")
                    colnames(d) <- c("ppm", "specNumber", "Spectral Intensity")
                    d$ppm <- c(t(x.ppm))
                    # d$color <- rep(alpha("gray", 0.6))
                    
                  # Make plot ####
                    g <- ggplot2::ggplot(d[d$specNumber==1,])
        
                  # Reverse axis ####
                    g <- g + ggplot2::scale_x_reverse(breaks = scales::breaks_pretty())
        
                  # Make simple ####
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
                    
                # For each spectrum, add plot objs: ####
                  for (s in 1:nrow(x.rev)){
                    # s <- 0
                      
                      # s <- s+1
                              g <- g + geom_area(data = data.frame(vals = x.rev[s,],
                                                                   ppms = x.ppm[s,]),
                                       na.rm = TRUE,
                                       mapping = aes(x = ppms, y = vals),
                                       colour = alpha("black", 0.6),
                                       linewidth = 0.5,
                                       fill="white",alpha=1)
                              
                            # Overlay the fit feature AUC ####

                              g <- g + geom_ribbon(data = data.frame(ymax = f.rev[s, ],
                                                                     ymin = min(x.rev[s, ], na.rm = T),
                                                                     ppms = f.ppm[s, ]),
                                       # na.rm = TRUE,
                                       mapping = aes(x = ppms, ymin = ymin, ymax = ymax),
                                       colour = "gray",
                                       linewidth = 0.5,
                                       fill="pink",alpha=0.4)
                      # plot(g)
                  }
                if (!is.null(label)){f.str <- str_c("Feature ", label, ", ")} else {f.str <- ''}
                g <- g + labs(title = str_c(f.str,
                                            round(min(f.ppm,na.rm = T),2),"-",
                                            round(max(f.ppm,na.rm = T),2)," ppm")) + 
                theme(plot.title = element_text(hjust = 0.5))
               return(list(g.obj = g,
                           row.inds = rowinds))
      
            })
    
    return(list(plots = plots.all.ranges,
                pars = list(exp.by = exp.by,
                            vshift = vshift,
                            hshift = hshift)
                )
           )
    # plot(plots.all.ranges[[1]]$g.obj)
}