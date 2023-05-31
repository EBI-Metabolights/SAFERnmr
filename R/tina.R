#' Tina Function. Tina Is Not Alignment.
#'
#' This is a bit of a semantic argument, but spectral alignment typically refers to
#' modifying spectra or reposition peaks such that they align with each other from
#' sample to sample. Because we've done FSE, features have much more information
#' which allows them to be detected in multiple regions combined by shape similarity.
#' Instead of modifying spectra, TINA reformulates alignment as a clustering
#' problem in feature (shape) space.
#'
#' There are also several filtering steps employed before this.
#' We cannot eliminate all poor shapes, but there are a couple of useful heuristics
#' which generally reduce unnecessary computation downstream. First, there are many
#' feature shapes which are either quite poor in quality, or do not contain
#' sufficient information to be useful for annotation. We keep features with:
#'   - in defined ppm range (generally [-1, 11])
#'   - long enough runs? (contains runs > noisewidth*3 adjacent points)
#'   - large enough subset? (>= 5 spectra with feature, good correlation reliability)
#'   - has at least 1 true peak > 0.3 of the range of intensities? (not just the
#'     side of a broad peak; not monotonic)
#' Since the same features will often be extracted multiple times (either in the
#' same spectral region, or other regions; i.e. same peak, but misaligned), it is
#' advantageous to reduce this redundancy by clustering feature shapes. We accomplish
#' this using a combination of UMAP projection and Affinity Propagation clustering.
#' Before comparing feature shapes, we align them to their maximum intensity
#' resonance. This is quick, and usually performs well enough.
#' UMAP uses euclidean distance as a default. In practice, UMAP is able to sort
#' feature shapes into tight clusters when low n_neighbors (e.g. 5) and min_dist
#' (e.g. 0.05) are used. This does not capture global relationships as well, but
#' we only use it to identify tight clusters.
#' All pairwise correlations (PCCs) are calculated for the feature shapes. A mask
#' for correlation thresholding is applied to the distance matrix (generated using
#' apcluster::negDistMat(pts, r=2), squared negative euclidean distance) to ensure
#' that clustered features have a high correlation as well. apcluster uses a q
#' parameter to optimize the initial preferences. Higher q -> stricter clusters.
#' Raising the lambda (dampening) parameter helps avoid oscillations which prevent
#' convergence, although raising this too high can make updates too slow to
#' converge within the number of iterations.
#' See
#' https://cran.r-project.org/web/packages/apcluster/vignettes/apcluster.pdf for
#' a full description of affinity propagation parameters
#' Ultimately, it doesn't matter much what clustering method is used, as this is
#' primarily a means of combining highly similar features to reduce the computational
#' burden of pairwise comparisons to reference spectra.

#'
#' @param pars A list of parameters for the TINA pipeline.
#'
#' @return A list containing the results and cluster labels from TINA
#'
#' @export
#' @import ggridges
#' @importFrom umap umap
#' @importFrom apcluster apcluster
#' @importFrom plotly plot_ly
#' @importFrom tictoc tic toc
#' @importFrom gridExtra grid.arrange
#' @importFrom magrittr %>%
#'
tina <- function(pars){
  message('-------------------------------------------------------')
  message('-------------------      TINA       -------------------')
  message('-------------------------------------------------------')
  message('\n\n\n')
  
  
################ Read parameters file ##################
  
  tmpdir <- pars$dirs$temp
  this.run <- paste0(tmpdir)

# Params ####

    fse.result <- readRDS(paste0(this.run,"/fse.result.RDS"))
      xmat <- fse.result$xmat
      ppm <- fse.result$ppm
      
    bounds <- pars$tina$bounds        # only consider signatures within this region (ppm)
    # rcutoff <- 0.99           # can also be set with "elbow" method
    # pcutoff <- 0.001          # cutoff for afc correlations (basically does nothing)
    # overlap.fraction <- 0.5   # don't compare features if they don't overlap by at least half their points
    #                           # additionally, when calculating the combined feature shape, only use points
    #                           # for which this % of features have data
    min.subset <- pars$tina$min.subset          # don't keep features if their subsets are too small (basically does nothing)
    prom.ratio <- pars$tina$prom.ratio

########### Setup ###################################################################################


    # Run tina_setup to set up successful features

        feature <- tina_setup(fse.result$storm_features, fse.result$xmat)

    # If throttling comparisons, keep just a random subset of features
        if (pars$debug$enabled == TRUE) {

          max.features <- pars$debug$throttle_features

        } else {

          max.features <- nrow(feature$stack)

        }

    # Run feature filter function to get the filter:
    #   - in defined ppm range?
    #   - long enough runs?
    #   - large enough subset?
    #   - has at least 1 true peak > 0.3 of the range of intensities?


    # Build feature filter

          filts  <- filterFeatures(feature, ppm = ppm,
                                    ppm.range = bounds, min.runlength = fse.result$noisewidth*2,
                                    min.subset = min.subset, prom.ratio = prom.ratio, give = "filter",
                                    max.features = max.features)

          saveRDS(filts, paste0(tmpdir, "/filts.RDS"))
          # filts <- readRDS(paste0(tmpdir, "/filts.RDS"))

    # Re-run tina_setup on the filtered features

        filt <- filts$all$not.null &
                     filts$all$inbounds.ppm &
                     filts$all$ss.pass &
                     filts$all$not.singlet &
                     filts$all$not.monotonic &
                     filts$all$not.null.driver &
                     filts$all$rand.subset

      # Rerun setup on filtered features
        feature <- tina_setup(fse.result$storm_features[filt], xmat)
       
        # saveRDS(feature, paste0(tmpdir, "/feature.RDS"))
        # feature <- readRDS(paste0(tmpdir, "/feature.RDS"))
        
################ Plotting Filtered features #######################
 
  # Plot the good features that passed: ####
    if (pars$tina$plots$filtered.features){
        message('printing good features to pdf...')

        everyNth <- every.nth(select = pars$storm$number.of.plots,
                              from = sum(filt))

        plot_stormRefRegions_grid(xmat = NULL, ppm,
                                  fse.result$storm_features %>% .[filt] %>% .[everyNth], # if not doing a small region
                                  plotLoc = paste0(tmpdir,'/'),
                                  filename = 'filtered.features.pdf',
                                  calcStocsy = FALSE, n_xticks = 4)
    }
        
  # Plot some of each filter type, if any, to get an idea of what we removed: ####

    if (pars$tina$plots$filtered.out){
      res <- lapply(1:length(filts$all), function(filt.num){
        # Plot the bad ones
        # filt.num <- 6
        bad.ones <- !filts$all[[filt.num]]
        filt.name <- filts$all %>% names %>% .[filt.num]

        if (any(bad.ones)){
          message('\n\n', sum(bad.ones), ' features failed ', filt.name, ' filter. Plotting profiles to pdf...')

              everyNth <- every.nth(select = pars$storm$number.of.plots,
                                    from = sum(bad.ones))
              message('Plotting every ', sum(everyNth), ' / ', length(everyNth), ' features...')
              plot_stormRefRegions_grid(xmat = NULL, ppm,
                                        fse.result$storm_features %>% .[bad.ones] %>% .[everyNth], # if not doing a small region
                                        plotLoc = paste0(tmpdir,'/'),
                                        filename = paste0('failed_filter_', filt.name, '.pdf'),
                                        calcStocsy = FALSE, n_xticks = 4)
        }
      })
      rm(res)
    }
      
      

########### SFE  ################################################################################## 

        
    # Do spec-feature extraction for all features
        message('\n\nsfe running...')
        t1 <- Sys.time()
            features.specd <- parallel::mclapply(1:nrow(feature$stack),
                                                 FUN = function(i){
            # features.specd <- lapply(1:nrow(feature$stack),
            #                          FUN = function(i){
                # print(i)
                # i <- 1
                sfe(feature, i,
                         xmat,
                         ppm,
                         r.thresh = 0.8)
                # try({})
                # if (!is.empty(res)){return(res)}else{return(i)}

              }, mc.cores = pars$par$ncores)

            message('Parallel sfe done on ', length(features.specd), ' features.')
        print(Sys.time() -t1)
        
        saveRDS(features.specd, paste0(tmpdir, "/features.specd.RDS"))
        # features.specd <- readRDS(paste0(tmpdir, "/features.specd.RDS"))

########### Plot sfe results  ################################################################################## 
        
    # Plot: single spec: ####
      # res <- features.specd[[1]]
      #   
      #   feat.fits <- lapply(1:length(res$feat$ss), function(i){
      # 
      #     ss.spec <- res$feat$ss[i]
      #     fit.feat <- apply.fit.specfeat(profile = res$feat$profile,
      #                                    fit.coefficients = res$fits[[i]]) %>% c
      # 
      #     full.pos <- res$feat$position %>% range(na.rm = T) %>% fillbetween
      #     spec.data <- fse.result$xmat[ss.spec, full.pos - res$lags[[i]]]
      # 
      #     fit <- list(spec.fit = spec.data,
      #                 feat.fit = fit.feat)
      # 
      #     # g <- plot.fit(fit = fit, type = 'auc', ppm = fse.result$ppm[full.pos])
      # 
      #     return(list(fit.feats = fit.feat,
      #                 fit.xrow = ss.spec,
      #                 fit.positions = full.pos,
      #                 plot = NULL
      #                 ))
      #   }) %>% unlist(recursive = F)
      # #   
      #   feat.fits <- split(feat.fits, names(feat.fits))
      #     feat.fits$fit.feats <- feat.fits$fit.feats %>% lapply(unlist) %>% do.call(rbind,.)
      #     feat.fits$fit.positions <- feat.fits$fit.positions %>% lapply(unlist) %>% do.call(rbind,.)
      #     feat.fits$fit.xrow <- feat.fits$fit.xrow %>% unlist
      #     feat.fits$plot <- feat.fits$fit.positions %>% unlist(recursive = F)
      #     feat.fits$pass.fit <- T
      # 
          
    # Plot: all specs: ####
      #   
        # Compute the feature stackplots

          # p <- project_features.stackplot(fse.result$xmat, fse.result$ppm,
          #                                  label = NULL,
          #                                  bestfits = feat.fits,
          #                                  exp.by = 0.05, # ppm
          #                                  vshift = 1, # pass as slide bar?
          #                                  hshift = 0, # pass as slide bar?
          #                                  sort.rows = F)
          # p$plots[[1]] %>% plot
          # 
        # Without the areas highlighted
          # i <- 1
          # i <- i + 500
          # res <- features.specd[[i]]
          # x.reg <- apply.lags.feat(feat = res$feat,
          #                          xmat = fse.result$xmat,
          #                          lags = res$lags)
          # simplePlot(x.reg$vals)
          # x.reg$cols.range
          
        
        
        
########### Adjust feature object with sfe results  ############################################################
    # Apply/add new feature info to old features #####
      # profile stack must be updated
      # - needs to allow for resizing?
      # - sfe will return trimmed profiles and undo any alignment of features
      
        # feature$stack.pre.sfe <- feature$stack
        # 
        # Get the width of the new feature mat ####
          n.cols <- lapply(features.specd, function(res){
            res$feat$profile %>% length
          }) %>% unlist %>% max

          # lags <- lapply(features.specd, function(x) x$lags) %>% unlist
          # hist(lags[abs(lags) < 50], breaks = 100)

          # if (n.cols != ncol(feature$stack.pre.sfe))
          # {
          #   message('feature stack now contains ',
          #           n.cols, ' columns, instead of ',
          #           ncol(feature$stack.pre.sfe))
          # }

        # Make new feature mat ####
          feature$stack <- lapply(features.specd, function(res){
            needed <- n.cols - length(res$feat$profile)
            return(  c(res$feat$profile, rep(NA, needed))  )
          }) %>% do.call(rbind, .)


          # feature$position stays the same (base feature)
          # feature$subset gets updated
          # feature$region is a question. I think range(ranges) makes sense.

          # res$feat
          # res$lags
          # res$fits
          # res$rmses
          # res$rvals
    
    # Align to feature maximum ####
 
        message('Aligning features...')
        feature.ma <- align.max(feature, scaling = FALSE)

        saveRDS(feature.ma, paste0(tmpdir, "/feature.ma.RDS"))
        # feature.ma <- readRDS(paste0(tmpdir, "/feature.ma.RDS"))

        
  # Plot all feature ranges ####
        pdf(file = paste0(this.run,'/','feature.ranges.pdf'),   # The directory you want to save the file in
            width = dim, # The width of the plot in inches
            height = dim)
            feature.shift_range <- feature$position %>% apply(., 1, function(x) range(ppm[x],na.rm = T))
            subs <- ind2subR(1:length(feature.shift_range), m = nrow(feature.shift_range))
            plot(feature.shift_range, subs$cols, pch = ".", cex = .01)
            segments(feature.shift_range[1,], 1:ncol(feature.shift_range),
                     feature.shift_range[2,], 1:ncol(feature.shift_range),
                     lwd = 0.1)
        dev.off()
       
###############################################################################################################           
        
# The TINA part  ####
              
   # OPTICS-based #### 
    
    if (nrow(feature.ma$stack) > 1000){
      t1 <- Sys.time()
      results <- tina_combineFeatures_optics(feature.ma$stack,
                                             max.eps = 50,
                                             minPts = 2,
                                             eps.stepsize = .01,
                                             max.plots = 600,
                                             plot.loc = ".",
                                             plot.name = "feature_clusters.pdf",
                                             nfeats = 10000,
                                             dist.threads = parallel::detectCores() - 1) # pars$par$ncores
  
      # Label the "noise" points as individual clusters
        stray.labels <- seq_along(results$clusters[[1]]) + max(results$labels) 
        stray.feats <- results$clusters[[1]] %>% as.list
      
      # Update the results object
        results$clusters <- c(results$clusters[-1], stray.feats)
        clusters <- list(method = 'optics',
                         results = results,
                         cluster.labs = clusts2labs(results$clusters),
                         groups = results$clusters)
        
        print(Sys.time() - t1)
        # length(clusters$cluster.labs %>% unique)
    } else {
      clusters <- list(method = 'none',
                       results = NULL,
                       cluster.labs = 1:nrow(feature.ma$stack),
                       groups = as.list(1:nrow(feature.ma$stack)))
    }
            
    saveRDS(clusters, paste0(this.run, "/clusters.RDS"))
    
 ###############################################################################################################   
    
    # check each cluster's quality
    # - calculate all pairwise alignments + ls fits for cluster features
    #   - returned as pw rmse matrix, lags
    # - identify key feature
    # - align all to key
    # - rmse-based threshold for removing features from cluster
    
          t1 <- Sys.time()
            # *** Note: this is parallelized for ncores - 2
            # *** Note: currently not using rmse cutoff.
            clust.info <- checkClusters(clusters = clusters, feature = feature.ma, par.cores = pars$par$ncores)
          print(Sys.time() - t1)

          keys <- lapply(clust.info, function(ci){
            ci$key.feat
          }) %>% unlist
                
      saveRDS(clust.info, paste0(tmpdir, "/clust.info.RDS"))
      # clust.info <- readRDS(paste0(tmpdir, "/clust.info.RDS"))
      
  ############################################################################################################### 
    
    # Split clusters
    
      # ....
      
  ############ Plot the cleaned clusters ####
  if (pars$tina$plots$cleaned.clusters){
    # Produce plots ####
        message("Generating plots. Progress:")
        cluster.list <- clusters$groups
        
        everyNth <- every.nth(select = pars$storm$number.of.plots, 
                              from = cluster.list)
        
        plots <- pblapply(cluster.list[everyNth], function(x)
          {
                fs <- feature$stack[x, ,drop = FALSE] %>% trim.sides
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

      # Produced feature.ma object (max-aligned)
      # - feature stack 
      # - position stack 
      # - ss
      #   - matrix 
      # - region
      #   - this could be updated to hold the list of regions from
      # - driver.relative (driver position; column of [].stack)

  # FINAL Objects to export: ####
  # final feature object (everything should apply to max-aligned)
  # - stack
  # - position
  # - driver.relative
  # - sfe
  #   - list of feature-specific objects

      # Erase the feature profiles - not necessary
      features.specd <- lapply(features.specd, function(x) {
        x$feat$profile <- NULL
        x$feat$position <- NULL
        x$feat$corr <- NULL
        x
      })
      
      feature.final <- list(stack = feature.ma$stack,
                            position = feature.ma$position,
                            driver.relative = feature.ma$subset,
                            sfe = features.specd)
      saveRDS(feature.final, paste0(tmpdir, "/feature.final.RDS"))

      
  # final clusters object
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
      
      saveRDS(cluster.final, paste0(tmpdir, "/cluster.final.RDS"))
      
#         #### #####
  message('-------------------------------------------------------')
  message('-------------------  TINA Complete  -------------------')
  message('-------------------------------------------------------')
        
}