#' Plot UMAP scores
#' 
#' Plots UMAP scores (layout) and labels clusters (if given) using a dataframe input (for shiny app).
#' 
#' @param df A data.frame with columns "x", "y", and "cluster".
#' @param clusters An optional vector of cluster labels for each point in the layout.
#' 
#' @return A plotly figure.
#' 
#' @importFrom plotly plot_ly layout
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export
plot_umap.scores.df <- function(df, clusters = NULL){
# plot umap scores (layout) and label clusters (if given)
# (cluster labels for each umap.obj$layout row)
# MTJ2022

  if (is.null(clusters)){
    
    # Set up standard marker
       std.marker <- list(size = 7.5,
                     color = 'rgba(255, 182, 193, .9)',
                     line = list(color = 'rgba(152, 0, 0, .8)',
                                 width = 2))

        fig <- plotly::plot_ly(data = data.frame(x = umap.obj$layout[,1],
                                                 y = umap.obj$layout[,2]),
                                           x = x, y = y,
                                           type="scatter", mode = "markers",
                                           marker = std.marker) %>%
                           plotly::layout(yaxis = list(zeroline = FALSE),
                                          xaxis = list(zeroline = FALSE) #,showlegend = FALSE
                                          )
         fig
         return(fig)
     
  }else{
    
    # Assign colors based on cluster

      cmap <- rep(RColorBrewer::brewer.pal(12, "Paired"), 
                  length.out = length(unique(df$cluster)))
      
      mycolors <- lapply(1:length(unique(df$cluster)), 
                         function(x) rep(cmap[x], 
                                         sum(df$cluster %in% x)
                                         )
                         ) %>% unlist
   
    # Plot
      fig <- plotly::plot_ly(type = 'scatter',
                             data = df, x = df$x, y = df$y, 
                             color = as.character(df$cluster),
                             colors = mycolors,
                             text = paste("Cluster ", df$cluster),
                             hoverinfo = 'text',
                             mode = 'markers',
                             marker = list(size = 7.5, 
                                           opacity = 0.5,
                                           line = list(width = 2,
                                                       opacity = 1))
                             ) %>%
        
           plotly::layout(yaxis = list(zeroline = FALSE),
                          xaxis = list(zeroline = FALSE),
                          showlegend = FALSE)
  
        fig
        return(fig)
  }
  
}

# Could also calculate ellipses for each HCA cutoff point, and draw them on top of one another like a topo map...
