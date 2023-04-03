#' Plot UMAP Scores
#'
#' This function plots the UMAP scores (layout) and optionally labels the clusters
#' (if given).
#'
#' @param umap.obj A UMAP object generated using the umap package.
#' @param clusters A vector of cluster labels for each row of umap.obj$layout.
#'
#' @return A plotly plot of the UMAP scores, with optional cluster labels.
#'
#'
#' @export
#' @importFrom plotly plot_ly layout
#' @importFrom RColorBrewer brewer.pal
plot_umap.scores <- function(umap.obj, clusters = NULL) {
  # plot umap scores (layout) and label clusters (if given)
  # (cluster labels for each umap.obj$layout row)
  # MTJ2022

  if (is.null(clusters)) {

    # Set up standard marker
    std.marker <- list(
      size = 7.5,
      color = "rgba(255, 182, 193, .9)",
      line = list(
        color = "rgba(152, 0, 0, .8)",
        width = 2
      )
    )

    fig <- plotly::plot_ly(
      data = data.frame(
        x = umap.obj$layout[, 1],
        y = umap.obj$layout[, 2]
      ),
      x = x, y = y,
      type = "scatter", mode = "markers",
      marker = std.marker
    ) %>%
      plotly::layout(
        title = str_c(
          "Features in UMAP - n_neighbors = ", umap.obj$config$n_neighbors,
          " and min_dist = ", umap.obj$config$min_dist
        ),
        yaxis = list(zeroline = FALSE),
        xaxis = list(zeroline = FALSE),
        showlegend = FALSE
      )
    fig
    return(fig)
  } else {

    # Assign colors based on cluster

    df <- data.frame(
      x = umap.obj$layout[, 1],
      y = umap.obj$layout[, 2],
      cluster = clusters
    )
    df$key <- rownames(umap.obj$layout)

    cmap <- rep(RColorBrewer::brewer.pal(12, "Paired"),
      length.out = length(unique(df$cluster))
    )

    mycolors <- lapply(
      1:length(unique(df$cluster)),
      function(x) {
        rep(
          cmap[x],
          sum(df$cluster %in% x)
        )
      }
    ) %>% unlist()

    # Plot
    fig <- plotly::plot_ly(
      type = "scatter",
      data = df, x = df$x, y = df$y,
      color = as.character(df$cluster),
      colors = mycolors,
      text = paste("Cluster ", df$cluster),
      hoverinfo = "text",
      mode = "markers",
      marker = list(
        size = 7.5,
        opacity = 0.5,
        line = list(
          width = 2,
          opacity = 1
        )
      )
    ) %>%
      plotly::layout(
        title = str_c(
          "Features in UMAP - n_neighbors = ", umap.obj$config$n_neighbors,
          " and min_dist = ", umap.obj$config$min_dist
        ),
        yaxis = list(zeroline = FALSE),
        xaxis = list(zeroline = FALSE),
        showlegend = FALSE
      )

    fig
    return(fig)
  }
}
