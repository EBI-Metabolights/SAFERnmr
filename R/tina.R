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
#' See
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
#' @importFrom plotly plot_ly
#' @importFrom tictoc tic toc
#' @importFrom gridExtra grid.arrange
#' @importFrom magrittr %>%
#'
tina <- function(pars){
  message('-------------------------------------------------------')
  message('-------------------      TINA       -------------------')
  message('-------------------------------------------------------')
  printTime()
  message('\n\n\n')
  
  
################ Read parameters file ##################
  
  tmpdir <- pars$dirs$temp

# Params ####

    fse.result <- readRDS(paste0(tmpdir,"/fse.result.RDS"))
      xmat <- fse.result$xmat
      ppm <- fse.result$ppm
      
      # Check fse.result for nulls
        fse.result %>% test_nullish # just double-checking
      
    bounds <- pars$tina$bounds        # only consider signatures within this region (ppm)
    min.subset <- pars$tina$min.subset          # don't keep features if their subsets are too small (basically does nothing)
    prom.ratio <- pars$tina$prom.ratio

########### Setup ###################################################################################

    # Run tina_setup to set up successful features

        feature <- tina_setup(fse.result$storm_features, xmat)
          feature %>% test_nullish
        feature %>% compress_features %>% debug_write("feature.RDS", pars)
        # feature <- readRDS(paste0(pars$dirs$temp, "/debug_extra.outputs", "/feature.RDS")) %>% expand_features
          
        message('\n\t---- Feature Filtering ----')
        
    # Run feature filter function to get the filter:
    #   - in defined ppm range?
    #   - long enough runs?
    #   - large enough subset?
    #   - has at least 1 true peak > 0.3 of the range of intensities?


    # Build feature filter
          printTime()
          filts  <- filterFeatures(feature, ppm = ppm,
                                    ppm.range = bounds, min.runlength = fse.result$noisewidth*2,
                                    min.subset = min.subset, prom.ratio = prom.ratio, give = "filter",
                                    max.features = nrow(feature$stack))

          filts %>% debug_write("filts.RDS", pars)
          # filts <- readRDS(paste0(pars$dirs$temp, "/debug_extra.outputs", "/filts.RDS"))
          
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
          feature %>% test_nullish('feature')
        message('\n\tFeature filtering complete. Saving results...')
        
        feature %>% compress_features %>% debug_write("feature.filtered.RDS", pars)
        # feature <- readRDS(paste0(pars$dirs$temp, "/debug_extra.outputs", "/feature.filtered.RDS")) %>% expand_features
        
################ Plotting Filtered features #######################
  tryCatch(expr = { 
  # Plot the good features that passed: ####
    if (pars$tina$plots$filtered.features){
        message('printing good features to pdf...')

        everyNth <- every_nth(select = pars$storm$number.of.plots,
                              from = sum(filt))

        plot_stormRefRegions_grid(xmat = NULL, ppm,
                                  fse.result$storm_features %>% .[filt] %>% .[everyNth], # if not doing a small region
                                  plotLoc = paste0(tmpdir,'/plots/'),
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

              everyNth <- every_nth(select = pars$storm$number.of.plots,
                                    from = sum(bad.ones))
              message('Plotting every ', sum(everyNth), ' / ', length(everyNth), ' features...')
              plot_stormRefRegions_grid(xmat = NULL, ppm,
                                        fse.result$storm_features %>% .[bad.ones] %>% .[everyNth], # if not doing a small region
                                        plotLoc = paste0(tmpdir,'/plots/'),
                                        filename = paste0('failed_filter_', filt.name, '.pdf'),
                                        calcStocsy = FALSE, n_xticks = 4)
        }
      })
      rm(res)
    }
      
      

  }, error = function(cond){message('Plotting filterd features failed in TINA. May want to look into this.')}
  )
########### SFE  ################################################################################## 

    # Force garbage collection to slim down workspace
        # Remove everything except
        # - feature
        # - xmat
        # - ppm
        # - pars
        rm(filts)
        rm(fse.result)
        gc()
        
    # Do spec-feature extraction for all features
        message('\n\n\t---- Spec-Feature Extraction (SFE) ----')
        printTime()
        t1 <- Sys.time()
            features.specd <- parallel::mclapply(1:nrow(feature$stack),
                                                 FUN = function(i){
                tryCatch(
                    expr = {
                            res <- 
                             sfe2(feature, i,
                                 xmat,
                                 ppm,
                                 r.thresh = pars$storm$correlation.r.cutoff)
                              res %>% test_nullish()
                              res
                    },
                    error = function(cond){
                             NA
                    }
                )                                   
              }, mc.cores = pars$par$ncores)

            message('\n\tParallel sfe done on ', length(features.specd), ' features.')
        print(Sys.time() -t1)

        # Nullish is not acceptable
        features.specd %>% test_nullish('features.specd')
        
        features.specd %>% debug_write("features.specd.RDS", pars)
        # features.specd <- readRDS(paste0(pars$dirs$temp, 
        #                                  "/debug_extra.outputs", 
        #                                  "/features.specd.RDS"))

        
########### Adjust feature object with sfe results  ############################################################
    # Apply/add new feature info to old features #####
      # * any features that failed sfe need to be removed
      # profile stack must be updated
      # - needs to allow for resizing?
      # - sfe will return trimmed profiles and undo any alignment of features
      
        message('\n\tapplying sfe results to feature stack...')
        
        passed.sfe <- !(features.specd %>% is.na)
        
        if (!any(passed.sfe)){stop('No features passed sfe.')} 
        
          # Go through feature object and apply filter
            feature$stack <- feature$stack[passed.sfe, ,drop = F]
            feature$position <- feature$position[passed.sfe, ,drop = F]
            feature$subset$ss.all <- feature$subset$ss.all[passed.sfe, ,drop = F]
            feature$subset$sizes <- feature$subset$sizes[passed.sfe]
            feature$driver.relative <- feature$driver.relative[passed.sfe]
            features.specd <- features.specd[passed.sfe]
            
        # Get the width of the new feature mat ####
          n.cols <- lapply(features.specd, function(res){
            res$feat$profile %>% length
          }) %>% unlist %>% max(na.rm = T)
            
          # Positions match stack

        # Make new feature mat using sfe-derived profile (left-justified) ####
          feature$stack <- lapply(features.specd, function(res){
            
            needed <- n.cols - length(res$feat$profile)
            
            return(  c(res$feat$profile, rep(NA, needed))  )
            
          }) %>% do.call(rbind, .)
            
          # *** note: at this point, positions do NOT match stack ***
          
          feature$position <- lapply(features.specd, function(res){

            needed <- n.cols - length(res$feat$position)

            return(  c(res$feat$position, rep(NA, needed))  )

          }) %>% do.call(rbind, .)
          
          feature$driver.relative <- lapply(features.specd, function(res){
            # driver.relative was somewhere in feature matrix columns. Now that features
            # are all left-justified, this will change. 
            return(  res$feat$driver.relative  )

          }) %>% do.call(rbind, .)
            
            feature %>% test_nullish('feature')

    # Align to feature maximum ####

        message('\n\taligning features to max peak...')
        feature.ma <- align_max(feature, scaling = FALSE)
           # if max align fails on a feature, remove it
           features.specd <- features.specd[feature.ma$alignment.success]
           feature.ma %>% test_nullish('feature.ma')
        
        feature.ma %>% compress_features %>% debug_write("feature.ma.RDS", pars)
        # feature.ma <- readRDS(paste0(pars$dirs$temp,
        #                              "/debug_extra.outputs",
        #                              "/feature.ma.RDS")) %>% expand_features
        
        rm(feature)

##################################################################################################################                
    # Plot all feature ranges ####
  tryCatch(expr = { 
          pdf(file = paste0(tmpdir,'/plots/','feature.ranges.pdf'))
        
              feature.shift_range <- feature.ma$position %>% apply(., 1, function(x) range(ppm[x],na.rm = T))
              subs <- ind2subR(1:length(feature.shift_range), m = nrow(feature.shift_range))
              plot(feature.shift_range, subs$cols, pch = ".", cex = .01)
              segments(feature.shift_range[1,], 1:ncol(feature.shift_range),
                       feature.shift_range[2,], 1:ncol(feature.shift_range),
                       lwd = 0.1)
          dev.off()
    }, error = function(cond){message('Plotting feature ranges failed in TINA. May want to look into this.')}
  )
            
          ######################################################################################################## 
                    
    # FINAL feature object to export: ####
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
        
      # Double-check for null features after sfe. Remove from feature.ma (this is much lighter than feature.final,
      # and will be used for downstream work in this file. 
        nullfeats <- feature.ma$stack %>% apply(1, function(x) all(is.na(x))) %>% unlist
          if (any(nullfeats)){warning('TINA: there were ', sum(nullfeats),' null features after sfe. Removing...')
            feature.ma$stack <- feature.ma$stack[!nullfeats, ,drop = F]
            feature.ma$position <- feature.ma$position[!nullfeats, ,drop = F]
            feature.ma$subset$ss.all <- feature.ma$subset$ss.all[!nullfeats, ,drop = F]
            feature.ma$subset$sizes <- feature.ma$subset$sizes[!nullfeats]
            feature.ma$region$sizes <- feature.ma$region$sizes[!nullfeats]
            feature.ma$driver.relative <- feature.ma$driver.relative[!nullfeats]
          }
            
        message('\n\twriting post-sfe features to file...')
        feature.final <- list(stack = feature.ma$stack,
                              position = feature.ma$position,
                              driver.relative = feature.ma$subset,
                              sfe = features.specd)
        
        feature.final %>% test_nullish('feature.final')
         
        saveRDS(feature.final %>% compress_features, paste0(tmpdir, "/feature.final.RDS"))
        # feature.final <- readRDS(paste0(tmpdir, "/feature.final.RDS")) %>% expand_features
        
        rm(features.specd)
        gc()
        
###############################################################################################################           
      do.clustering <- tryCatch(
        {
          pars$tina$do.clustering
        }, error = function(cond){FALSE}
      )
        
      tryCatch(
        {
          # This doesn't return clusters - just writes the data to file
            cluster_features(pars, feature.final, min.features = 1000, do.clustering)
            # cluster.final <- readRDS('/Users/mjudge/Documents/current_run/cluster.final.RDS')
        },
        error = function(cond)
          {
            message('\n\tClustering failed, trying again without')
            cluster_features(pars, feature.final, min.features = 1000, do.clustering = FALSE)
          }
      )
      
        
#         #### #####
  message('-------------------------------------------------------')
  message('-------------------  TINA Complete  -------------------')
  message('-------------------------------------------------------')
        
}


