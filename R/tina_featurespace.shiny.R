#' Run Feature Space Shiny app
#'
#' Runs a Shiny app that displays a UMAP layout with associated clusters, and allows the user to select
#' points and view their associated extracted features.
#'
#' @param study A character string specifying the accession number of the study.
#' @param umap.obj A UMAP object, output of umap() or a related function.
#' @param cluster.labs A vector of cluster assignments for each point in the UMAP layout.
#' @param featureStack A object of extracted features for each point in the UMAP layout.
#'
#' @importFrom shiny fluidPage plotOutput sidebarLayout sidebarPanel mainPanel h1 h2 h4 titlePanel verbatimTextOutput
#' @importFrom plotly event_register plotlyOutput
#' @importFrom plotly plot_ly layout event_data
#' @import RColorBrewer
#' @export
#'
#' @author MTJ2022
#' @export
runFeatureSpace.shiny <- function(study,
                                  umap.obj,
                                  cluster.labs,
                                  featureStack) {

  # Prep the data ----
  # Convert to df outside of the plot_umap.scores function ####
  df <- data.frame(
    x = umap.obj$layout[, 1],
    y = umap.obj$layout[, 2],
    cluster = cluster.labs
  )
  df$key <- rownames(umap.obj$layout)

  # Build the app ----
  # Define UI ----
  ui <- fluidPage(
    titlePanel(h1(stringr::str_c(study, ": Feature Space"))),
    sidebarLayout(
      position = "right",
      sidebarPanel(
        h2("Extracted Features for selected points"),
        plotOutput("featurePlot")
        # verbatimTextOutput("se")
      ),
      mainPanel(
        h2("UMAP Layout with AP clusters projected"),
        h4(stringr::str_c(
          "Features in UMAP - n_neighbors = ", umap.obj$config$n_neighbors,
          " and min_dist = ", umap.obj$config$min_dist
        )),
        h4("use selection box to view features for any points"),
        plotlyOutput("umap.clusters")
      )
    )
  )

  # Define server logic ----
  server <- function(input, output) {
    # Feature Space Plot, record selection in event_register
    output$umap.clusters <- renderPlotly({

      # Plot all the features with their cluster assignments
      plot_umap.scores.df(df, df$cluster) %>%
        layout(dragmode = "zoom") %>%
        event_register(event = "plotly_brushed") # plotly_selected

      # Alternatively, we could simply plot the exemplars, then use clust inds to plot features.
      # ... this would be much quicker, but less data.
      # It would be nice to have the option
      # ... to select individual data points or entire clusters at will.
    })


    # Actually plot using the indices ####
    output$featurePlot <- renderPlot({
      # Identify selected points (can't just get indices at present)
      sel.range <- event_data("plotly_brushed")
      in.x <- df$x >= sel.range$x[1] & df$x <= sel.range$x[2]
      in.y <- df$y >= sel.range$y[1] & df$y <= sel.range$y[2]
      selectedPoints <- which(in.x & in.y)
      # Plot the features for those points (if )
      if (length(selectedPoints > 0)) {
        fs <- featureStack[selectedPoints, , drop = FALSE] %>% trim.sides()
        simplePlot(fs)
      } else {
        NULL
      }
    })
  }

  # Run the app ----
  shinyApp(ui = ui, server = server)
}
