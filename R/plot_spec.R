#' Plot spectrum in plotly
#' 
#' plot reference spectrum (or any spectrum) using plotly, with event handler for region selection.
#' 
#' @param spec spectrum (vector or 1d matrix)
#' @param ppm ppm vector for spectrum
#' @param aucs area under the curves (not used currently)
#' @param title plot title (i.e. compound name)
#' @param source.name source event name (for shiny interactive input passing)
#' 
#' @return A plotly figure.
#' 
#' @importFrom plotly plot_ly layout
#' @importFrom RColorBrewer brewer.pal
#' @importFrom magrittr %>%
#'
#' 
#' @export
plot_spec <- function(spec, ppm, aucs = NULL, title = '', source.name = 'refspec'){
# plot umap scores (layout) and label clusters (if given)
# (cluster labels for each umap.obj$layout row)
# MTJ2022

      df <- data.frame(ppm = ppm,
                       intensity = spec)
      df$intensity[is.na(df$intensity)] <- 0

  
     # Since plotly only does select boxes for plots with markers, add two trivial
     # and invisible markers to our plot:
        pts <- seq(1, length(df$ppm), length.out = 2)
        
        fig <- plotly::plot_ly(source = source.name) %>%
          add_trace(x = df$ppm, y = df$intensity, type = 'scatter', mode = 'lines') %>%           # spectrum (what we want to see)
          add_trace(x = df$ppm[pts], y = df$intensity[pts], type = 'scatter', mode = 'markers',   # markers (fully transparent)
                        marker = list(opacity = 0),
                    showlegend = F)                                                               # turn off legend
      

    # Add ref feature aucs (if present):
      if (!is.null(aucs)){
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
      }
                # g <- g + labs(title = str_c("Feature ", label, ", ",
                #                             round(min(f.ppm,na.rm = T),2),"-",
                #                             round(max(f.ppm,na.rm = T),2)," ppm")) + 
                # theme(plot.title = element_text(hjust = 0.5))
               
        # fig <- plotly::plot_ly(data = df,
        #                        x = x, y = y,
        #                        type="scatter", mode = "lines") %>%
        #        plotly::layout(yaxis = list(zeroline = FALSE),
        #                       xaxis = list(zeroline = FALSE) #,showlegend = FALSE
        #                       )
         
         # fig <- ggplotly(g)
         fig <- fig %>% layout(  xaxis = list(autorange = "reversed"),
                                 title = title)
         
         fig

}


