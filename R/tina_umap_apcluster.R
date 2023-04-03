#' Combines and clusters features in a matrix
#'
#'
#' @param featureStack A matrix of features to be clustered.
#' @param doUMAP Logical indicating whether to perform UMAP projection on the features.
#' @param umap.n_neighbors the number of approximate nearest neighbors used to construct
#' the initial high-dimensional graph. It effectively controls how UMAP balances local
#' versus global structure - low values will push UMAP to focus more on local structure by
#' constraining the number of neighboring points considered when analyzing the data in high dimensions,
#' while high values will push UMAP towards representing the big-picture structure while losing fine detail.
#' @param umap.min_dist the minimum distance between points in low-dimensional space. This parameter controls
#'  how tightly UMAP clumps points together, with low values leading to more tightly packed embeddings.
#' Larger values of min_dist will make UMAP pack points together more loosely, focusing instead on the
#' preservation of the broad topological structure.
#' @param umap.n_epochs The number of iterations to run.
#' @param ap.rcutoff rcutoff for filtering out correlations before clustering
#' @param ap.q q parameter is the shared value for the initial preferences s(i,k) for a priori clustering. See apcluster doc.
#' @param ap.max.iter The maximum number of iterations for the ap cluster.
#' @param ap.lam The "self-responsibility" parameter for the affinity propagation clustering algorithm.
#' @param max.plots The maximum number of scatter plots to generate for the clusters.
#' @param plot.loc The directory in which to save the scatter plot PDF.
#' @param plot.name The name of the scatter plot PDF.
#'
#' @return A list containing UMAP and clustering results, cluster assignments, pairwise
#' distances and correlations.
#' @export
tina_combineFeatures <- function(featureStack,
                                 doUMAP = TRUE,
                                 umap.n_neighbors = 20,
                                 umap.min_dist = 0.1,
                                 umap.n_epochs = 500,
                                 ap.rcutoff = 0.99,
                                 ap.q = .95,
                                 ap.max.iter = 1000,
                                 ap.lam = 0.9,
                                 max.plots = 250,
                                 plot.loc = "./",
                                 plot.name = "study_feature_clusters.pdf") {

  ############ Setup ####
  # Scale the features ####
  featureStack <- apply(featureStack, 1, scale.between) %>% t()

  # Calculate pairwise correlations ####
  message("Computing feature correlations (assuming max-aligned)...")
  allcors <- cor(featureStack %>% t(), use = "pairwise.complete.obs", method = "pearson")

  ############ UMAP the features first? ####
  if (doUMAP) {

    # Setup ####
    # Get data, zero-fill
    mat <- featureStack
    mat[is.na(mat)] <- 0

    # Set up config struct
    custom.config <- umap::umap.defaults
    custom.config$n_epochs <- umap.n_epochs
    custom.config$n_neighbors <- umap.n_neighbors
    custom.config$min_dist <- umap.min_dist
    custom.config$random_state <- 123
    cc <- custom.config

    # Run one umap ####

    message(
      "Running umap with min_dist = ", cc$min_dist,
      " and n_neighbors = ", cc$n_neighbors, "..."
    )
    tictoc::tic()
    feats.umap <- umap::umap(
      d = mat,
      config = cc,
      method = "naive"
    )
    tictoc::toc()

    pts <- feats.umap$layout
  } else {
    pts <- featureStack
  }

  ############ Affinity Propagation Clustering based on UMAP layout ####

  # Calculate pairwise euclidean distances between points in the umap projection ####
  d <- apcluster::negDistMat(pts, r = 2)

  # Filter similarity matrix based on pairwise correlations ####
  message("Removing similarities for feature pairs where r < ", ap.rcutoff, " (Pearson).")
  d[allcors < ap.rcutoff] <- -Inf # makes them ineligible (best)

  # Run Affinity Propagation Clustering ####
  #   Ideally, we'd provide a rough estimate of the number of features for the dataset,
  #   and then use the q estimation function to try a couple q's in that neighborhood.

  message("Running affinity propagation clustering on umap'd feature coordinates (a2a)...")
  tictoc::tic()
  apres <- apcluster::apcluster(d,
    maxit = ap.max.iter,
    lam = ap.lam,
    q = ap.q, details = TRUE
  )

  # Print report ####
  message("Affinity propagation clustering complete. Produced ", length(apres@clusters), " clusters.")
  tictoc::toc()


  ############ Plot Clusters ####
  # Produce plots ####
  message("Generating plots. Progress:")
  clusters <- apres@clusters
  everyNth <- c(T, rep(F, max(floor(length(clusters) / max.plots) - 1, 0)))
  plots <- pblapply(clusters[everyNth], function(x) {
    fs <- featureStack[x, , drop = FALSE] %>% trim.sides()
    return(simplePlot(fs))
  })
  # Print plots to file ####
  message("Printing plots to file ", paste0(plot.loc, plot.name), "...")
  dim <- 3 * round(sqrt(length(plots)))
  pdf(
    file = paste0(plot.loc, plot.name), # The directory you want to save the file in
    width = dim, # The width of the plot in inches
    height = dim
  )
  gridExtra::grid.arrange(grobs = plots)
  dev.off()
  message("Complete.")

  ##############################
  # Ideally, we'd now like to associate each of these features back to the featureStack rows,
  # but this can easily be done with the clustering results as row indices.


  return(list(
    umap.results = feats.umap,
    apcluster.results = apres,
    clusters = clusters,
    umap.distances = d,
    corrs = allcors
  ))
}
