#' The CLUSTERING part of tina ####
#' Ultimately, it doesn't matter much what clustering method is used, as this is
#' primarily a means of combining highly similar features to reduce the computational
#' burden of pairwise comparisons to reference spectra.
#' Taking this part out into its own function to allow reversion upon failure. 
#'
#' @param pars A list of parameters for the TINA pipeline.
#' @param feature A feature obj
#' @param min.features minimum features required to do OPTICS
#' @param do.clustering actually do clustering. If FALSE, save dummy clusters objs.
#'
#' @return Nothing, saves cluster files and plots to pars$dirs$temp
#'
#' @import pbapply
#' @importFrom gridExtra grid.arrange
#' @importFrom magrittr %>%
#' 
#' @export
cluster_features <- function(pars, feature, min.features = 1000, do.clustering = F){

  dummyClusters <- function(clust.numbers){
        clusters <- list(method = 'none',
                         results = NA,
                         cluster.labs = clust.numbers,
                         groups = as.list(clust.numbers))
        return(clusters)
  }

    this.run <- pars$dirs$temp
    feature$sfe <- NULL
      
    # Check to see if clustering param was provided; if not, turn it off
    
    if (nrow(feature$stack) > min.features & do.clustering){
     # OPTICS-based ####
        
      clusters <- tryCatch(expr = { 
        
          message('\n\t---- OPTICS-based clustering ----')
          printTime()
          t1 <- Sys.time()
          results <- tina_combineFeatures_optics(feature$stack,
                                                 max.eps = 50,
                                                 minPts = 2,
                                                 eps.stepsize = .01,
                                                 max.plots = 600,
                                                 plot.loc = this.run,
                                                 plot.name = "feature_clusters.pdf",
                                                 nfeats = pars$tina$nfeats,
                                                 dist.threads = pars$par$ncores)
    
          # Label the "noise" points as individual clusters
            noiseclust <- which(results$labels == 0)
            stray.labels <- seq_along(results$clusters[[noiseclust]]) + max(results$labels)
            stray.feats <- results$clusters[[noiseclust]] %>% as.list
            
          # What if a cluster is > max size allowed?
    
          # Update the results object
          # (remove noiseclust, add to the end the broken noise clust)
            results$clusters <- c(results$clusters[-noiseclust], stray.feats)
            print(Sys.time() - t1)
            
            list(method = 'optics',
                 results = results,
                 cluster.labs = clusts2labs(results$clusters),
                 groups = results$clusters)
    
            
        }, error = function(cond){
          message('\n\tIn TINA: OPTICS clustering failed. Reverting to no clustering.')
          dummyClusters(1:nrow(feature$stack))
        })
     
      
    } else {
     # SKIP CLUSTERING ####
        message('\n\tClustering on < 1000 features is not recommended. Skipping to avoid artifacts.')
        clusters <- dummyClusters(1:nrow(feature$stack))
    }
    
    # Produced object: clusters. Check for nullish elements:
        
        # clusters %>% test_nullish
        # clusters$results 
        
###############################################################################################################   
    
    # check each cluster's quality
    # - calculate all pairwise alignments + ls fits for cluster features
    #   - returned as pw rmse matrix, lags
    # - identify key feature
    # - align all to key
    # - rmse-based threshold for removing features from cluster
    # - single features (or no clustering) : pass through with dummy vals
    
      message('\n\tcheckClusters: optimizing intra-cluster alignment and selecting representative features...')
    
          # Force garbage collection to slim down workspace
            gc()
          printTime()
          t1 <- Sys.time()
            # *** Note: this is parallelized for ncores - 2
            # *** Note: currently not using rmse cutoff.
            clust.checks <- checkClusters(clusters = clusters, feature = feature, # just needs stack
                                        par.cores = pars$par$ncores, 
                                        par.type = pars$par$type)
          print(Sys.time() - t1)
  
          
          clust.info <- clust.checks$clust.info
          clusters <- clust.checks$clusters
          clusters.prev <- clust.checks$clusters.prev
         
          keys <- lapply(clust.info, function(ci){
            ci$key.feat
          }) %>% unlist
          
        clust.info %>% test_nullish
        
        message('\n\tcheckClusters complete. Saving...\n')     
        
      clusters %>% debug_write("clusters.RDS", pars)
      # clusters <- readRDS(paste0(pars$dirs$temp, "/debug_extra.outputs", "/clusters.RDS"))

      clust.info %>% debug_write("clust.info.RDS", pars)
      # clust.info <- readRDS(paste0(pars$dirs$temp, "/debug_extra.outputs", "/clust.info.RDS"))
      
  ############ Plot the cleaned clusters ####
  tryCatch(expr = {
      if (do.clustering & pars$tina$plots$cleaned.clusters){
      # Produce plots ####
          message("Generating plots. Progress:")
          cluster.list <- clusters$groups
          
          everyNth <- every_nth(select = pars$storm$number.of.plots, 
                                from = cluster.list)
          clust.subset <- seq_along(cluster.list) %>% .[everyNth]
          
          plots <- pblapply(clust.subset, function(x)
            {
            # print(x)
                  feat.inds <- clust.info[[x]]$labels
                  fs <- feature$stack %>% 
                    lag_features(., clust.info[[x]]$lag.table, 
                                 to = clust.info[[x]]$key.feat) %>% 
                    trim_sides 
                  
                  return(simplePlot(fs))
            })
          
      # Print plots to file ####
          message("Printing checked cluster plots to file ", "...")
          dim <- 3*round(sqrt(length(plots)))
          pdf(file = paste0(this.run,'/','checked.clusters.pdf'),   # The directory you want to save the file in
              width = dim, # The width of the plot in inches
              height = dim)
          gridExtra::grid.arrange(grobs = plots)
          dev.off()
          message("Complete.")            
      }
    }, error = function(cond){message('Plotting cleaned clusters failed after checkClusters() in TINA. May want to look into this.')}
  )
###############################################################################################################          
      
# Make output objects that make sense: ####

  # TINA Results: what have we done? ####
  
      # Filtered features
      # - filt
      
      # Produced feature object
      # - feature stack
      # - stack.pre.sfe (for record; not used)
      # - position stack
      # - ss
      #   - matrix 
      #   - diffs (not used)
      #   - sizes (not used)
      #   - overlaps (not used)
      #   - fraction (not used)
      # - region
      #   - sizes (size of the feature profile, not super useful)
      #   - fraction (not used)
      #   - overlaps (not used)
      # - driver.relative (driver position; column of [].stack)

      # Produced feature object (max-aligned)
      # - feature stack 
      # - position stack 
      # - ss
      #   - matrix 
      # - region
      #   - this could be updated to hold the list of regions from
      # - driver.relative (driver position; column of [].stack)

      
  # final clusters object ####
  # - labels (single vector labeling each cluster with its ID, matches feature stack rows)
  # - keys (list of feature stack row/clusters$labels indices of cluster representatives)
  # - info (results of check.clusters)
  #   - label (cluster number)
  #   - features.aligned ()
  #   - key.feat
  #   - lag.table
  #   - profile
  #   - rmse.mean.clust (mean rmse of all features aligned to key feature)
  #   - rmse.pw.fits (pairwise fits between features)
  #   - rmse.to.best (each feature's rmse to the key feature)

      cluster.final <- list(labels = clusters$cluster.labs,
                            clust.results = clusters$results,
                            keys = keys,
                            info = clust.info,
                            groups = clusters$groups,
                            method = clusters$method)
      cluster.final[-which(names(cluster.final) == "clust.results")] %>% # OPTICS results will contain Inf. at times
        test_nullish 
      
      saveRDS(cluster.final, paste0(pars$dirs$temp, "/cluster.final.RDS"))
}