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
                                 plot.name = "study_feature_clusters.pdf"){
  
  ############ Setup ####        
    # Scale the features ####
      featureStack <- apply(featureStack, 1, scale.between) %>% t
      
    # Calculate pairwise correlations ####
      message("Computing feature correlations (assuming max-aligned)...")
      allcors <- cor(featureStack %>% t, use = "pairwise.complete.obs", method = "pearson")
      
  ############ UMAP the features first? ####
  if (doUMAP){
    
      # Setup ####
          # Get data, zero-fill
            mat <- featureStack
            # simplePlot(featureStack[1:50,] %>% trim.sides)
            mat[is.na(mat)] <- 0

          # Set up config struct
            custom.config <- umap::umap.defaults
              custom.config$n_epochs <- umap.n_epochs
              custom.config$n_neighbors <- umap.n_neighbors
              custom.config$min_dist <- umap.min_dist
              custom.config$random_state <- 123
              cc <- custom.config

      # Run one umap ####
          
          message("Running umap with min_dist = ", cc$min_dist,
                  " and n_neighbors = ", cc$n_neighbors ,"...")
          tictoc::tic()
            feats.umap <- umap::umap(d = mat,
                           config = cc,
                           method = 'naive')
          tictoc::toc()

      # Plot the umap ####
            # plot(feats.umap$layout[,1],feats.umap$layout[,2])
            
            # title(paste0("Features in UMAP - n_neighbors = ", feats.umap$config$n_neighbors,
            #             " and min_dist = ", feats.umap$config$min_dist))

              # plot_umap.scores(feats.umap)
              
              #,
              #                marker = list(size = 10,
              #                              color = 'rgba(255, 182, 193, .9)',
              #                              line = list(color = 'rgba(152, 0, 0, .8)',
              #                                          width = 2)))
              # fig <- fig %>% layout(title = paste0("Features in UMAP - n_neighbors = ", feats.umap$config$n_neighbors,
              #                               " and min_dist = ", feats.umap$config$min_dist),
              #          yaxis = list(zeroline = FALSE),
              #          xaxis = list(zeroline = FALSE))
              # 
              # fig
      
      # Try to optimize umap params ####
        #   df <- opt_umap_ser(mat, config = cc, 
        #                      try.nns = seq(20,100, by = 20),
        #                      try.mindists = seq(0.05, 1, length.out = 5),
        #                      n.dim = 2)
        #   dims <- df$ax %>% unique %>% length
        #   
        #   plots <- lapply(seq(1, nrow(df), by = dims), function(x) 
        #     {
        #       rows <- seq(x, x + dims-1)
        #       umap.layout <- data.frame(x = df[rows[1], 4:ncol(df)] %>% as.numeric,
        #                                 y = df[rows[2], 4:ncol(df)] %>% as.numeric)
        # 
        #     # Plot
        #       g <- ggplot(umap.layout,aes(x,y)) +
        #             geom_point(shape=1,size=.1) +
        #             ggtitle(paste0("nn = ", df[x, "nn"], ", ",
        #                     paste0("md = ",    df[x, "md"]))) +
        #             theme_classic()
        #       return(g)
        #     }
        #   )
        #   
        # # Rank results by spread
        #   dmats <- lapply(seq(1, nrow(df), by = dims), function(x) 
        #     {
        #       rows <- seq(x, x + dims-1)
        #       umap.layout <- data.frame(x = df[rows[1], 4:ncol(df)] %>% as.numeric,
        #                                 y = df[rows[2], 4:ncol(df)] %>% as.numeric)
        #       
        #     # Calc spread
        #       d <- dist(umap.layout, method = 'euclidean') %>% as.matrix
        #       d.s <- sort(d)
        #       d.s <- d[seq(1, length(d.s), length.out = 1E6)] / max(d.s)
        #       # scattermore::scattermoreplot(seq_along(d.s), sort(d.s))
        #     }
        #   )
        #   
          # df <- opt_umap_par(mat, config = cc, ncores = 3,
          #                    try.nns = seq(5,100, length.out = 5),
          #                    try.mindists = seq(0.05, 1, by = .4))
          
          # try_nns <- seq(10,30, by = 5)
          # # try_mindists <- seq(0.05, 1, by = .2)
          # 
          # feats.umaps <-
          #               lapply(try_nns, function(x) {
          #                 cc <- custom.config
          #                 cc$n_neighbors <- x
          #                 cc$min_dist <- .1
          #                 message("Running umap with min_dist = ", cc$min_dist, " and n_neighbors = ", cc$n_neighbors ,"...")
          #                 tictoc::tic()
          #                 res <- umap::umap(d = mat,
          #                                config = cc,
          #                                method = 'naive')
          #                 tictoc::toc()
          #                 return(res)
          #               })
          # rm(res)
          # plots <- lapply(feats.umaps,
          #                 function(res) {
          #                   # res <- feats.umaps[[1]]
          #                   df <- data.frame(x = res$layout[,1], y = res$layout[,2])
          #                   g <- ggplot(df,aes(x,y)) +
          #                         geom_point(shape=1,size=.1) +
          #                         ggtitle(paste0("n_neighbors = ", res$config$n_neighbors)) +
          #                         # ggtitle(paste0("min_dist = ", res$config$min_dist)) +
          #                         theme_classic()
          #                   return(g)
          #                 })
          # 
            # dim <- 3*round(sqrt(length(plots)))
            # pdf(file = paste0("/Users/mjudge/Documents/tmp_res/","umap_opt.pdf"),   # The directory you want to save the file in
            #     width = dim, # The width of the plot in inches
            #     height = dim)
            # gridExtra::grid.arrange(grobs = plots, ncol = length(df$md %>% unique))
            # dev.off()
            
    # Select the best umap projection
        
            
          
            
    pts <- feats.umap$layout
  } else {pts <- featureStack}
      
  ############ Affinity Propagation Clustering based on UMAP layout ####
        
        # Calculate pairwise euclidean distances between points in the umap projection ####
          # pts <- feats.umap$layout
          # pts <- featureStack
          # pts <- res$scores[,1:2]
          #   plot(pts[,1],pts[,2])
          d <- apcluster::negDistMat(pts, r=2)
          
        
        # Filter similarity matrix based on pairwise correlations ####
          # d[allcors<0.9] <- min(d) # set to minimum (better)
          message("Removing similarities for feature pairs where r < ", ap.rcutoff, " (Pearson).")
          d[allcors< ap.rcutoff] <- -Inf # makes them ineligible (best)
          
        # Run Affinity Propagation Clustering ####
        #   Ideally, we'd provide a rough estimate of the number of features for the dataset,
        #   and then use the q estimation function to try a couple q's in that neighborhood.
        
          message("Running affinity propagation clustering on umap'd feature coordinates (a2a)...")
          tictoc::tic()
            apres <- apcluster::apcluster(d,
                                          maxit = ap.max.iter,
                                          lam = ap.lam,
                                          q = ap.q, details = TRUE)
            
        # Do HCA on top of apcluster result?
          # apcluster::aggExCluster()
            # d[allcors< ap.rcutoff] <- 0
            

        # Do OPTICS + DBSCAN? 
          # see optics_play.R, par grid search for 2 params
                     
        # Louvain clustering?
            
        
        # Print report ####
          # apcluster::plot(apres)
          message("Affinity propagation clustering complete. Produced ", length(apres@clusters), " clusters.")
          tictoc::toc()

        
  ############ Plot Clusters ####
    # Produce plots ####
        message("Generating plots. Progress:")
        clusters <- apres@clusters
        everyNth <- c(T,rep(F, max(floor(length(clusters) / max.plots ) - 1,0)))
        plots <- pblapply(clusters[everyNth], function(x) {
                        fs <- featureStack[x, ,drop = FALSE] %>% trim.sides
                        # fss <- apply(fs, 1, function(thisrow) scale.to.minmax(thisrow, fs)) %>% t
                        # ref <- list(signature = list(idx.wind = 1:ncol(fs),
                        #             vals = compSpec[1,idx]),
                        #             wind = list(inds = compSpec[1,] %>% seq_along))
                        #             ... + plot_addRef(specRegion = ,ref = ,ppm = 1:ncol(fs))
                        return(simplePlot(fs))
                        })
    # Print plots to file ####
        message("Printing plots to file ", paste0(plot.loc, plot.name), "...")          
        dim <- 3*round(sqrt(length(plots)))
        pdf(file = paste0(plot.loc,plot.name),   # The directory you want to save the file in
            width = dim, # The width of the plot in inches
            height = dim)
        gridExtra::grid.arrange(grobs = plots)
        dev.off()
        message("Complete.")      
        
  ##############################
  # Ideally, we'd now like to associate each of these features back to the featureStack rows, 
  # but this can easily be done with the clustering results as row indices. 
        
        
  return(list(umap.results = feats.umap,
              apcluster.results = apres, 
              clusters = clusters,
              umap.distances = d,
              corrs = allcors))
} 