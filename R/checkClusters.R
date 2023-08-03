#' Calculate feature cluster information - wrapper for check.cluster that parallelizes it
#' See check.cluster documentation. 
#' 
#' @param clusters cluster data from tina_combineFeatures_optics (and adjustments in tina)
#' @param feature feature object containing data on all the features
#'
#' @importFrom magrittr %>%
#' @import parallel
#' @import foreach
#' @import doParallel
#'
#' @return clust.info.split (list of cluster information for each cluster, concatenated)
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
          
      
      
    ############ Set up for parallel compute ####

      message("Setting up parallel cluster...\n\n")
      
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
                                  c.info <- check_cluster(clust.feats = 1:length(cluster$labels), 
                                                          feature.stack = cluster$f.stack)
                                  
                                  c.info$labels <- c.info$labels %>% cluster$labels[.]
                                  c.info$key.feat <- c.info$key.feat %>% cluster$labels[.]
                                  c.info$lag.table$f1 <- c.info$lag.table$f1 %>% cluster$labels[.]
                                  c.info$lag.table$f2 <- c.info$lag.table$f2 %>% cluster$labels[.]
                                  # i <- i + 1
                                  # if(any(is_nullish(c.info))){
                                  #   c.info$
                                  # } else {
                                    return(c.info)
                                  # }
                                  
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
      
            
      return(clust.info.split)
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



