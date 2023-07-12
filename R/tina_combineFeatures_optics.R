#' TINA clustering using OPTICS
#' Tina Is Not Alignment ;)
#' 
#' Why? 
#' - We can't carefully compare all features in all cases, but we can use some rougher
#'   comparisons to guide the process and reduce computational burden. 
#' How?
#' - distance metric is just euclidean distance squared and negated. Computed using
#'   parallelDist::parDist(), which multithreads according to the ncores set in 
#'   our params.yaml file. It's pretty damn efficient, and offers other distance
#'   metrics that we can pass through at some point (although euclidean makes sense).
#' - Clustering using density-based clustering should allow us to find local minima
#'   in the feature comparison landscape, where clusters represent pockets of shape
#'   space. These are the dominant shapes in our spectra. We want to avoid 
#'   overclustering, but still reduce the number of features to a set which has 
#'   multiple instances of each shape (e.g. reliably extracted) and smaller (so
#'   we do fewer match comparisons).
#' Why OPTICS?
#' - It's density-based (and handles different densities), uses few parameters, 
#'   it produces a hierarchical arrangement of features, and it's wicked fast.
#' How to optimize eps?
#' - Set the max to a high value, then compute the hierarchy. Scan a range of eps
#'   values and see how many clusters get picked. I usually see a curve with a 
#'   max. Since we want to avoid overclustering, choose that point as the cutoff. 
#'   NOTE: this only works when minPts > 1.
#'   
#'  Assumes you aligned the features in the input matrix somehow (e.g. based on their maximum peaks)
#' 
#' MTJ 2022
#'   
#' 
#'
#' @param featureStack feature profile matrix
#' @param max.eps (OPTICS) maximum epsilon value (may need adjustment, possibly related to mean dist or #' of points because they're scaled profiles?)
#' @param minPts (OPTICS) minimum number of points needed to define core distance
#' @param eps.stepsize (OPTICS) step size for assessing optimal eps. value
#' @param max.plots max number of plot to generate. Evenly spaced subset will be chosen to meet this requirement
#' @param plot.loc where the plot gets written to file
#' @param plot.name plot filename (example.pdf)
#' @param nfeats a logical value indicating whether to scale the feature values to between 0 and 1 based on the highest peak's intensity
#' @param dist.threads a logical value indicating whether to scale the feature values to between 0 and 1 based on the highest peak's intensity
#'
#' @return a list object containing
#' @importFrom magrittr %>%
#' @import dbscan
#' @import parallelDist
#' @import pbapply
#' @import gridExtra
#'
#' @export
#### Function ####        
tina_combineFeatures_optics <- function(featureStack,
                                         max.eps = 50,
                                         minPts = 2,
                                         eps.stepsize = .01,
                                         max.plots = 600,
                                         plot.loc = "./",
                                         plot.name = "feature_clusters.pdf",
                                         nfeats = 10000,
                                         dist.threads = 1){
  
  
      # feature <- readRDS( paste0(pars$dirs$temp, "/feature.RDS"))
      # feature.dist.pw(feature$stack, 1:100, pars)

  
  ############ Setup ####        
    # Scale the features ####
      featureStack <- apply(featureStack, 1, scale_between) %>% t
      
    # Calculate pairwise correlations ####
      # message("Computing feature correlations (assuming max-aligned)...")
      # allcors <- cor(featureStack %>% t, use = "pairwise.complete.obs", method = "pearson")
      
  ############ Clustering  ####
        
          message("Computing (aligned) feature similarities...")
          t1<-Sys.time()
            
            mat <- featureStack[1:min(nfeats,nrow(featureStack)),]
            mat[is.na(mat)] <- 0 # No NAs allowed
            
            d <- parallelDist::parDist(mat, method = 'euclidean', threads = dist.threads)
              closeAllConnections()
              d <- d^2
              d <- 0-d
              # scattermore::scattermoreplot(x = 1:length(d), y = d %>% as.vector %>% sort)
            # saveRDS(d, './data.d.RDS')
            # d <- readRDS('./data.d.RDS')
          print(Sys.time()-t1)

         # # Do OPTICS ####
         
          t1<-Sys.time()
          message('Doing OPTICS clustering on features using euclidean distances...')
            opres <- dbscan::optics(d, eps = max.eps, minPts = minPts, search = "kdtree") # kdtree only available for euclidean dist
          print(Sys.time()-t1)
          
          message('Evaluating cluster number at a range of eps_cl values...')
            eps.vals <- seq(0, max.eps/minPts, eps.stepsize)
            num.clusters <- pblapply(eps.vals, function(x) {
              dbres <- dbscan::extractDBSCAN(opres, eps_cl = x)
              run.labels(dbres$cluster) %>% max
              # run.labels(zeroes) + 
            }) %>% unlist
            
            eps.val <- num.clusters %>% which.max %>% eps.vals[.]
          
          # Plot
          pdf(file = paste0(plot.loc,'/optics.opt.epsvalue.pdf'),   # The directory you want to save the file in
            width = 3, # The width of the plot in inches
            height = 3)
          
            plot(eps.vals, num.clusters)
            abline(v = eps.val, col = 'red')
            abline(h = max(num.clusters), col = 'red')
          dev.off()
          
            
          message('Using maximum cluster number: ', max(num.clusters), ' (eps_cl = ', eps.val, ')')
            dbres <- dbscan::extractDBSCAN(opres, eps_cl = eps.val)
            clusters <- run.labels(dbres$cluster)
            zeroes <- clusters == 0
            
          message('OPTICS Clustering Results: ')
          message(sum(zeroes), ' features (', round(sum(zeroes)/length(zeroes), 2)*100,'%) were excluded by clustering. The other ', sum(!zeroes), ' features (', round(sum(!zeroes)/length(zeroes), 2)*100,'%) fell into ', max(clusters),' shape clusters...')
          
          cluster.sizes <- lapply(unique(clusters[clusters>0]), function(x) which(clusters == x) %>% length) %>% unlist
          
          pdf(file = paste0(plot.loc,'/optics.clusters.sizeDistribution.pdf'),   # The directory you want to save the file in
            width = 3, # The width of the plot in inches
            height = 3)

            # hist(cluster.sizes, breaks = 100, xlab = 'Cluster Size', ylab = 'Frequency')
            plot(sort(cluster.sizes), xlab = 'Cluster', ylab = 'Cluster Size')

          dev.off()
          
          # # Do HDBSCAN (pretty good) ####
          # t1 <- Sys.time()
          # hdbres <- dbscan::hdbscan(x = mat, minPts = 2, gen_hdbscan_tree = T, gen_simplified_tree = T, verbose = T)
          # print(Sys.time()-t1)
          # # hdbres1 <- hdbres
          # # hdbres <- hdbres1
          # # plot(hdbres)
          # message('Optimizing cut point...')
          # cutpts <- seq(0.0, 10, by = .01)
          # nclusts <- pblapply(cutpts, function(cutpt)
          #   {
          #     clusters <- stats::cutree(hdbres$hc, h = cutpt)
          #     length(unique(clusters))
          #   }) %>% unlist
          # plot(x = cutpts, y = nclusts)
          # ep <- smerc::elbow_point(x = cutpts, y = nclusts)
          # 
          # clusters <- stats::cutree(hdbres$hc, h = ep$x) # take minimum clustering
          #   abline(h = max(clusters))
        
  ############ Plot Clusters ####
    # Produce plots ####
        message("Generating plots. Progress:")
        c.labs <- unique(clusters)
        cluster.list <- lapply(c.labs, function(x) which(clusters == x))
        
        clusters <- cluster.list[c.labs > 0] # remove noise points
        everyNth <- c(T,rep(F, max(floor(length(clusters) / max.plots ) - 1,0)))
        plots <- pblapply(clusters[everyNth], function(x) 
          {
              # Check to see if plot will be empty (should never happen)
                if (all(is.na(featureStack[x, ,drop = FALSE]))){warning('TINA: clustering: no non-NA columns in cluster ', x,'!')}
          
                fs <- featureStack[x, ,drop = FALSE] %>% trim_sides
                return(simplePlot(fs))
          })
        
    # Print plots to file ####
        message("Printing plots to file ", paste0(plot.loc, '/',plot.name), "...")
        dim <- 3*round(sqrt(length(plots)))
        pdf(file = paste0(plot.loc,'/',plot.name),   # The directory you want to save the file in
            width = dim, # The width of the plot in inches
            height = dim)
        gridExtra::grid.arrange(grobs = plots)
        dev.off()
        message("Complete.")      
        
  ##############################
  # Ideally, we'd now like to associate each of these features back to the featureStack rows, 
  # but this can easily be done with the clustering results as row indices. 
        
        
  return(list(optics.result = opres,
              eps.cutoff = eps.val, # gotten using maximum of distribution
              clusters = cluster.list,
              labels = c.labs,
              cluster.sizes = cluster.sizes))
} 
