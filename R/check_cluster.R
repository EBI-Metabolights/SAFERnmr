#' Calculate feature cluster information
#' - align features
#' - choose representative (key) feature
#' - calculate profile for cluster
#' 
#' @param clust.feats cluster feature indices (rows in feature.stack)
#' @param feature.stack matrix holding feature profiles (feature$stack)
#'
#' @importFrom magrittr %>%
#'
#' @return clust.info list; labels of features, rmses, profile, key feature, lag.table, other info
#'
#' @export
#' @importFrom magrittr %>%
    check.cluster <- function(clust.feats, feature.stack){
        
        if (length(clust.feats) == 1){
          return(singleFeat(feature.stack[clust.feats, ,drop = F], clust.feats))
        }
      
    # Choose representative shape
      # Calculate all pairwise lsfs and rmses
          
          feature.stack <- feature.stack[clust.feats, ,drop = F] %>% trim_sides
            # simplePlot(feature.stack)
          
      # Align all the features to each other ####
      
        lag.table <- feat_align(feature.stack, max.hits = 1, counter = F)
       
      # Remove duplicate matches ####
        
        lag.table <- cbind(lag.table$f1,lag.table$f2) %>% t %>% sortPairs %>% t %>% duplicated %>% "!"(.) %>% lag.table[.,]
      
      # Build relative alignment matrix ####
        # al.mat <- zeros(nrow(feature.stack))
        #   amlis <- sub2indR(rows = lag.table$f1, cols = lag.table$f2, m = nrow(feature.stack))
        #   rownames(al.mat) <- 1:nrow(feature.stack)
        #   colnames(al.mat) <- rownames(al.mat)
        # al.mat[amlis] <- lag.table$lag.in.f2
          
      # Get the RMSE ####
          lag.table <- lapply(1:nrow(lag.table), function(m) {
            pair <- lag.table[m, ]
            al.fts <- lag_features(feature.stack, pair) %>% trim_sides # al.fts[2,] %>% simplePlot
            fit <- fit_leastSquares(al.fts[1,], al.fts[2,], plots = F, scale.v2 = T)
              # fit$plot %>% plot
            pair$rmse <- fit$rmse
            return(pair)
          }) %>% do.call(rbind,.)
                
      # Pick key feature by minimizing average rmse  ####
        # Build best pairwise rmse matrix (optimal alignment, then rmse for each pair of features in the cluster):
          linds <- sub2indR(rows = lag.table$f1, cols = lag.table$f2, m = nrow(feature.stack))
          rmses <- matrix(0, nrow(feature.stack), nrow(feature.stack))
          rmses[linds] <- lag.table$rmse
            rmses[lower.tri(rmses)] <- rmses %>% upper.tri %>% rmses[.]
            centrality.score <- colMeans(rmses) 
            key.feat <- which.min(centrality.score)
              
        # Align the other features to the key 
        #   - if you trim_sides, you lose any indexing within the pos.stack
        #   - however, that can be recovered using lag.table. Plus, the aligned feats 
        #     are already out of index (the relative shift of the key feature will need
        #     to be derived using trim_sides(out = "inds"), but this can be done if needed)
          al.fts <- lag_features(feature.stack, lag.table, to = key.feat) %>% trim_sides # %>% simplePlot

            # Calc fits ####
              fits.pw <- lapply(1:nrow(al.fts), function(m) {
                fit <- fit_leastSquares(al.fts[m, ],
                                        al.fts[key.feat, ],
                                        plots = F,
                                        scale.v2 = T)
                  # fit$plot %>% plot
                return(fit)
              })
            
              al.fts <- lapply(fits.pw, function(fit) fit$feat.fit) %>% do.call(rbind,.) %>% trim_sides
              rmses.to.best <- lapply(fits.pw, function(fit) fit$rmse) %>% unlist
              rmse.mean.clust <- mean(rmses.to.best)

          # Get final profile using rmse-weighted mean of each column ####
              # Define weights
                
                w <- 1/rmses.to.best
                w[is.infinite(w)] <- max(w[!is.infinite(w)])
                
                prof <- al.fts %>% apply(., 2, function(x) stats::weighted.mean(x, w))
              
          # Reset lag.table and key.feat inds (currently relative within cluster) to actual feature inds 
            
            lag.table$f1 <- lag.table$f1 %>% clust.feats[.]
            lag.table$f2 <- lag.table$f2 %>% clust.feats[.]
            key.feat <- key.feat %>% clust.feats[.]
            
          cluster.info <- list(labels = clust.feats,
                               rmse.mean.clust = rmse.mean.clust,
                               profile = prof,
                               features.aligned = al.fts,
                               key.feat = key.feat,
                               lag.table = lag.table,
                               rmses.pw.fits = rmses,
                               rmses.to.best = rmses.to.best)
          return(cluster.info)
    }

    singleFeat <- function(feat, ind){
      feat <- trim_sides(feat)
      return(
              list(labels = ind,
                   rmse.mean.clust = 0,
                   profile = feat,
                   features.aligned = feat,
                   key.feat = ind,
                   lag.table = data.frame(f1 = ind, f2 = ind, lag.in.f2 = 0, pos.big = NA, val = NA, rmse = 0),
                   rmses.pw.fits = 0,
                   rmses.to.best = 0)
      )
    }