# drawHeatmap ##################################################################
#' 
#' Makes a plotly heatmap for a matrix in the style of browse_evidence().
#'
#' @param mat a numerical matrix of values, can have rownames and colnames.
#' @param dropRowNames whether or not to replace the rownames with 1:nrow(mat)
#' @param clipRowNames number of characters to cut off rownames to (if too long)
#' @param source.name tag for tracking in shiny apps
#' @return interactive plotly heatmap object
#' 
#' @importFrom magrittr %>%
#' @importFrom stringr str_replace
#' @import plotly
#' 
#' @export
drawHeatmap <- function(mat, dropRowNames = F, clipRowNames = NA, source.name = 'heatmap'){
  
  if(dropRowNames){
    rownames(mat) <- NULL
  } else {
  
    if(!is.na(clipRowNames)){ # assume integer
        rownames(mat) <- rownames(mat) %>% stringr::str_trunc(clipRowNames, ellipsis = '')
    }
  
    # If rownames are not unique, modify them ####
      # - what we want to do is give the first n characters
        rn <- data.frame(name = rownames(mat)) %>% 
                dplyr::group_by(name) %>%
                # add row number within each group
                dplyr::mutate(instance = dplyr::row_number()) %>%
                dplyr::ungroup()
        
        rn$unique.name <- lapply(1:nrow(rn), function(r) {
          
            paste(
                    rn$name[r],
                    rep(' ', rn$instance[r]-1),
                    collapse = ''
            )
     
          }) %>% unlist
        rownames(mat) <- rn$unique.name

}

if(nrow(mat) == 1){
      
      fig <- 
            plot_ly(source = source.name, # x = colnames(mat),
                    z = mat, #zmin = min(mat, na.rm = T), zmax = max(mat, na.rm = T),
                    zauto = TRUE,
                    y = rownames(mat),
                    x = colnames(mat),
                    colors = grDevices::colorRamp(c("#ffeda0","#feb24c", "#f03b20")),
                    type = 'heatmap',
                    hovertemplate = paste('<b>Compound</b>: %{y}',
                                          '<br><b>Sample</b>: %{x}',
                                          '<br><b>Score</b>: %{z:.2f}',
                                          '<extra></extra>')
                    ) %>%
            layout(
                    yaxis = list(
                                  # title = rownames(mat),
                                  zeroline = FALSE,
                                  showline = FALSE,
                                  showticklabels = TRUE,
                                  ticks = FALSE,
                                  showgrid = FALSE
                                ),
                    xaxis = list(title = list(text='Sample', font = list(size = 20), standoff = 25)),
                    hoverlabel = list(font=list(size=15)))
                   
} 
else {
  fig <- 
      plot_ly(source = source.name, # x = colnames(mat),
              type = 'heatmap',
              y = rownames(mat),
              x = colnames(mat),
              z = mat, 
              colors = grDevices::colorRamp(c("#ffeda0","#feb24c", "#f03b20")),
              zauto = TRUE, #zmin = 0, zmax = 1, # color scale
              hovertemplate = paste('<b>Compound</b>: %{y}',
                                    '<br><b>Sample</b>: %{x}',
                                    '<br><b>Score</b>: %{z:.2f}',
                                    '<extra></extra>')
              ) %>%
      layout(yaxis = list(title = list(text='Compound', font = list(size = 20), standoff = 25)),
              xaxis = list(title = list(text='Sample', font = list(size = 20), standoff = 25)),
              hoverlabel = list(font=list(size=15))
             )
  
}

fig
}

# drawScatterScores ############################################################
#' 
#' 
#' Makes a plotly scatterplot in the style of browse_evidence().
#'
#' @param mat a numerical matrix of values, can have rownames and colnames.
#' @param dropRowNames whether or not to replace the rownames with 1:nrow(mat)
#' @param source.name tag for tracking in shiny apps
#' @return interactive plotly heatmap object
#' 
#' @importFrom magrittr %>%
#' @importFrom stringr str_replace
#' @import plotly
#' 
#' @export
drawScatterScores <- function(mat, dropRowNames = F, source.name = 'scatter'){
  if(dropRowNames){
    rownames(mat) <- NULL
  }
  
  # Convert to df
    
    df <- ind2subR(1:length(mat), nrow(mat))
    df$score <- c(mat)
    df$Compound <- rownames(mat) %>% .[df$cols]
    df$Sample <- df$cols
    df <- data.frame(df)
  
  

  if(nrow(mat) == 1){
            fig <-
              plot_ly(source = source.name, # x = colnames(mat),
                      type = 'scatter',
                      mode = 'markers',
                      df, x = ~cols, y = ~score,
                      #size = ~score,
                      color = ~score,
                      colors = colorRamp(c("#ffeda0","#feb24c", "#f03b20")),
                      marker = list(  
                                      # colorbar = list(title = "Evidence Score"),
                                      cauto = FALSE,
                                      cmin = 0,
                                      cmax = 1
                                    ),
                      hovertemplate = paste('<br><b>Sample</b>: %{x}',
                                            '<br><b>Score</b>: %{y:.2f}',
                                            '<extra></extra>')
                      ) %>%
              layout(
                      yaxis = list(title = list(text='Score', font = list(size = 20), standoff = 25)),
                      xaxis = list(title = list(text='Sample', font = list(size = 20), standoff = 25)),
                      hoverlabel = list(font=list(size=15))) %>% hide_colorbar()
              # colorbar(limits = c(0, 1),
              #                  title = "Evidence Score")

  } 
    else {
    fig <- plot_ly(source = source.name, # x = colnames(mat),
                   y = rownames(mat) %>% stringr::str_trunc(30),
                      type = 'scatter',
                      mode = 'markers',
                      df, x = ~rows, y = ~cols,
                      size = ~score,
                      color = ~score,
                      colors = colorRamp(c("#ffeda0","#feb24c", "#f03b20")),
                      marker = list(    
                                      # colorbar = list(title = "Evidence Score"),
                                      cauto = FALSE,
                                      cmin = 0,
                                      cmax = 1
                                    ),
                      
                      hovertemplate = paste('<br><b>Sample</b>: %{x}',
                                            '<br><b>Score</b>: %{y:.2f}',
                                            '<extra></extra>')
                      ) %>%
              layout(
                      yaxis = list(title = list(text='Compound', font = list(size = 20), standoff = 25)),
                      xaxis = list(title = list(text='Sample', font = list(size = 20), standoff = 25)),
                      hoverlabel = list(font=list(size=15))
              ) %>% hide_colorbar()
      
              # colorbar(limits = c(0, 1),
              #          title = "Evidence Score")
    
  } 
  fig
}


# plot_spec ############################################################
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

      df <- data.frame(ppm = ppm %>% c,
                       intensity = spec %>% c)
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


