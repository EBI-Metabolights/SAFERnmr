#' Calculate feature cluster information - wrapper for check.cluster that parallelizes it
#' See check.cluster documentation. 
#' No features are lost during this check, but clusters which fail will be split into
#' 
#' 
#' @param clusters cluster data from tina_combineFeatures_optics (and adjustments in tina)
#' @param feature feature object containing data on all the features
#'
#' @importFrom magrittr %>%
#' @import parallel
#' @import foreach
#' @import doParallel
#'
#' @return clust.info.split: (list of cluster information for each cluster, concatenated)
#'         clusters.prev: clustering assignments before check - includes failed
#'         clusters: updated clusters object 
#'          - adjusted.during.check: TRUE is added if any clusters are split. This indicates
#'            that the clusters in clusters$results (from OPTICS) will NOT be the same as
#'            adjusted clusters. Use the clusters$cluster.labs and clusters$groups.
#'         
#'
#' @export
#' @importFrom magrittr %>%
checkClusters <- function(clusters, feature, 
                          par.cores =  1,
                          par.type = 'PSOCK'){
  
  
  # Try max-align and clustering
    # feature <- feature.ma
    # clusters.by.grp <- clusters$results$clusters[clusters$results$labels != 0]
    clusters.by.grp <- clusters$groups
    
    # Serially ####
    # message('Cluster ', i)
    # clust.info <- lapply(1:length(clusters.by.grp), function(i){
    #     message('Cluster ', i)
    #     c.info <- check.cluster(clust.feats = clusters.by.grp[[i]], 
    #                   feature.stack = feature$stack, 
    #                   rmse.cutoff = 0.3)
    #     return(c.info)
    # })
    
    # In parallel ####
      # Split the clusters and features into list
        clust <- lapply(clusters.by.grp, function(clabs){
          list(labels = clabs,
               f.stack = feature$stack[clabs, ,drop = F])
        })
    
        # Sort the clusters by size to optimize par processing order:
          sort.by.size <- order(  lapply(clust, function(x) x$f.stack %>% nrow) %>% unlist  )
          clust <- clust[sort.by.size]
        # *** this should be replaced by distribute()
      
      
    ############ Set up for parallel compute ####

      # Par dependencies ####
        library(doParallel)
        library(parallel)
        library(foreach)

      # Par setup ####

        ncores <- par.cores
        message('Starting parallel pool on ', ncores, ' cores...')
        my.cluster <- parallel::makeCluster(ncores, type = par.type)
        doParallel::registerDoParallel(cl = my.cluster)
        if (foreach::getDoParRegistered()){
          message('\tpool started on ', foreach::getDoParWorkers(), ' cores')
        }
        
      # Par run ####
        # i <- 1
        
        clust.info <-   foreach(cluster = clust,
                                .combine='c',
                                .errorhandling="pass") %dopar% 
                          {
                                  # print(i)
                                  # cluster <- clust[[193]]  
                              c.info <- tryCatch(
                                    expr = {
                                      
                                              c.info <- check_cluster(clust.feats = 1:length(cluster$labels), 
                                                                      feature.stack = cluster$f.stack)
                                              
                                              c.info$labels <- c.info$labels %>% cluster$labels[.]
                                              c.info$key.feat <- c.info$key.feat %>% cluster$labels[.]
                                              c.info$lag.table$f1 <- c.info$lag.table$f1 %>% cluster$labels[.]
                                              c.info$lag.table$f2 <- c.info$lag.table$f2 %>% cluster$labels[.]
                                              
                                              return(c.info)
                                      
                                    },
                                    error = function(cond){
                                      
                                              c.info <- list(
                                                 labels = cluster$labels,
                                                 rmse.mean.clust = NA,
                                                 profile = NA,
                                                 features.aligned = NA,
                                                 key.feat = NA,
                                                 lag.table = NA,
                                                 rmses.pw.fits = NA,
                                                 rmses.to.best = NA)
                                              
                                              return(c.info)
                                    }
                                )
                                  
                                    return(c.info)

                          }
      # Clean up ####
        parallel::stopCluster(my.cluster)
        message('Pool closed.')

      clust.info.split <- split(clust.info, names(clust.info))
      
      clust.info.split <- lapply(1:length(clust.info.split$key.feat), function(i){
        list(labels = clust.info.split$labels[i]$labels %>% unlist(recursive = F),
             features.aligned = clust.info.split$features.aligned[i]$features.aligned,
             key.feat = clust.info.split$key.feat[i]$key.feat,
             lag.table = clust.info.split$lag.table[i]$lag.table,
             profile = clust.info.split$profile[i]$profile,
             rmse.mean.clust = clust.info.split$rmse.mean.clust[i]$rmse.mean.clust,
             rmses.pw.fits = clust.info.split$rmses.pw.fits[i]$rmses.pw.fits,
             rmses.to.best = clust.info.split$rmses.to.best[i]$rmses.to.best)
      })

      # Unsort the clusters:
        clust.info.split <- clust.info.split[order(sort.by.size)]
      
      # If any clusters failed the check, unlist the members to their own clusters.
      
        failed.check <- lapply(clust.info.split, function(x) x %>% is_nullish %>% any) %>% unlist
        
        if (any(failed.check)){
          
          message('\n\t', sum(failed.check), ' clusters failed check and will be split. Their pre-check labels are also saved...')
          new.clusters <- lapply(clust.info.split[failed.check], function(x) {
            
                                  # Loop through the cluster members and create a single-feat cluster
                                    lapply(x$labels, function(f){
                                      
                                      singleFeat(feature$stack[f, ], f)
                                      
                                    })
                                    
                                }) %>% unlist(recursive = FALSE)
        
        
        
          # Delete bad clusters and cat on their split members:
            # First, save a record of previous clusters:
              clusters.prev <- clusters
              clusters.prev$results <- NULL # no need to repeat
              clusters.prev$failed <- failed.check
              
            clust.info.split <- c(clust.info.split[!failed.check], new.clusters)
            
          # Adjust the clusters object to reflect the split clusters, and return
            
            # Get lists of feature labels for each cluster
              groups <- lapply(clust.info.split, function(x) x$labels)
            
            # Get the feature-wise cluster label assignments from groups
              assignments <- lapply(1:length(groups), function(x) 
                
                data.frame(feat.num = groups[[x]], 
                           lab = x)
                
                )%>% do.call(rbind,.) %>% arrange(feat.num)
              
            # Modify clusters obj and add indication that it was adjusted
              clusters$cluster.labs <- assignments$lab
              clusters$groups <- groups
              clusters$adjusted.during.check <- TRUE
        
          return(
                  list(clust.info = clust.info.split,
                       clusters = clusters,
                       clusters.prev = clusters.prev)
                 )
              
        } else {
          # If there were no failed clusters
          return(
                  list(clust.info = clust.info.split,
                       clusters = clusters,
                       clusters.prev = NULL)
                 )
        }
}

  

# Commands to run script in terminal 
# cd /Users/mjudge/Documents/galaxy.instance/annoTateR
# R
# setwd("/Users/mjudge/Documents/GitHub/MARIANA-NMR")
#   source("./deps.R")
# setwd("/Users/mjudge/Documents/galaxy.instance/annoTateR")
# pars <- yaml.load_file("./data/params.yaml", eval.expr=TRUE) # load param file from copy
# feature.stack <- readRDS(paste0(pars$dirs$temp,'/feature.stack.RDS'))
# matches <- feat_align(featureStack = feature.stack, max.hits = 1)



