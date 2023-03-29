#' Tina Function. Tina is not annotation.
#'
#' Function for running TINA stage of pipe pipeline.
#'
#' @param pars A list of parameters for the TINA pipeline.
#'
#' @return A list containing the results and cluster labels from the TINA pipeline.
#'
#' @export
#' @importFrom umap umap
#' @importFrom apcluster apcluster
#' @importFrom plotly plot_ly
#' @importFrom shiny Shiny
#' @importFrom ggridges ggridges
#' @importFrom tictoc tic toc
#' @importFrom gridExtra grid.arrange
#' @importFrom coop pcor covar
tina <- function(pars){
  message('-------------------------------------------------------')
  message('-------------------      TINA       -------------------')
  message('-------------------------------------------------------')
  message('\n\n\n')
  
################ Read parameters file ##################
  
  tmpdir <- pars$dirs$temp
  this.run <- paste0(tmpdir)


# 
#   # getPackages.par()
#   packages <- c("gridExtra","umap","tictoc","apcluster","plotly","shiny", "ggridges")
#       new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
#       if(length(new.packages) > 1){install.packages(new.packages)} 
#       pacman::p_load(packages, character.only = TRUE)

# Params ####

    fse.result <- readRDS(paste0(this.run,"/fse.result.RDS"))
    bounds <- pars$tina$bounds        # only consider signatures within this region (ppm)
    # rcutoff <- 0.99           # can also be set with "elbow" method
    # pcutoff <- 0.001          # cutoff for afc correlations (basically does nothing)
    # overlap.fraction <- 0.5   # don't compare features if they don't overlap by at least half their points
    #                           # additionally, when calculating the combined feature shape, only use points
    #                           # for which this % of features have data
    min.subset <- pars$tina$min.subset          # don't keep features if their subsets are too small (basically does nothing)
    prom.ratio <- pars$tina$prom.ratio

########### Setup #############
        
  
    # Run tina_setup to set up successful features
        
        feature <- tina_setup(fse.result$storm_features, fse.result$xmat)

    # Run feature filter function to get the filter:
    #   - in defined ppm range?
    #   - long enough runs? 
    #   - large enough subset?
    #   - has at least 1 true peak > 0.3 of the range of intensities?
        
        filt  <- filterFeatures(feature, ppm = fse.result$ppm, 
                                  ppm.range = bounds, min.runlength = fse.result$noisewidth*3,
                                  min.subset = min.subset, prom.ratio = prom.ratio, give = "filter")
        
    # Re-run tina_setup on the filtered features
        feature <- tina_setup(fse.result$storm_features[filt], fse.result$xmat)
        
    # Do spec-feature extraction for features
    
    # Align to feature maximum
        feature.ma <- align.max(feature, scaling = FALSE)
        
    # 
        saveRDS(feature.ma, paste0(tmpdir, "/feature.RDS"))

        
  # Plot all feature ranges ####
            # feature.shift_range <- feature$position %>% apply(., 1, function(x) range(ppm[x],na.rm = T))
              # subs <- ind2subR(1:length(feature.shift_range), m = nrow(feature.shift_range))
            # plot(feature.shift_range, subs$cols, pch = ".", cex = .01)
            # segments(feature.shift_range[1,], 1:ncol(feature.shift_range),
            #          feature.shift_range[2,], 1:ncol(feature.shift_range),
            #          lwd = 0.1)
       
        
        
# The TINA part (skip for the moment) ####
# ##### UMAP -> APcluster ####

    results <- tina_combineFeatures(feature.ma$stack,
                                     doUMAP = TRUE,
                                     umap.n_neighbors = pars$tina$umap$n_neighbors, #30,
                                     umap.min_dist = pars$tina$umap$min_dist, #0.8,
                                     ap.rcutoff = pars$tina$affinity.propagation$rcutoff,
                                     ap.q = pars$tina$affinity.propagation$q, # if not converging, try reducing a little.
                                     ap.lam = pars$tina$affinity.propagation$lambda, # oscillations possible. dampening factor
                                     plot.loc = paste0(this.run,"/"),
                                     plot.name = paste0("clusters_umap_ap_r",pars$tina$affinity.propagation$rcutoff,".pdf"))
        
    clusters <- list(results = results,
                     cluster.labs = clusts2labs(results$clusters))
    
    saveRDS(clusters, paste0(this.run, "/clusters.RDS"))

  message('-------------------------------------------------------')
  message('-------------------  TINA Complete  -------------------')
  message('-------------------------------------------------------')
        
}