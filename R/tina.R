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
#' @importFrom umap umap
#' @importFrom apcluster apcluster
#' @importFrom plotly plot_ly
#' @importFrom shiny Shiny
#' @importFrom ggridges ggridges
#' @importFrom tictoc tic toc
#' @importFrom gridExtra grid.arrange
#' @importFrom coop pcor covar
tina <- function(pars) {
  message("-------------------------------------------------------")
  message("-------------------      TINA       -------------------")
  message("-------------------------------------------------------")
  message("\n\n\n")

  ################ Read parameters file ##################

  tmpdir <- pars$dirs$temp
  this.run <- paste0(tmpdir)

  # Params ####

  fse.result <- readRDS(paste0(this.run, "/fse.result.RDS"))
  bounds <- pars$tina$bounds # only consider signatures within this region (ppm)
  #                           # additionally, when calculating the combined feature shape, only use points
  #                           # for which this % of features have data
  min.subset <- pars$tina$min.subset # don't keep features if their subsets are too small (basically does nothing)
  prom.ratio <- pars$tina$prom.ratio

  ########### Setup #############


  # Run tina_setup to set up successful features

  feature <- tina_setup(fse.result$storm_features, fse.result$xmat)

  # Run feature filter function to get the filter:
  #   - in defined ppm range?
  #   - long enough runs?
  #   - large enough subset?
  #   - has at least 1 true peak > 0.3 of the range of intensities?

  filt <- filterFeatures(feature,
    ppm = fse.result$ppm,
    ppm.range = bounds, min.runlength = fse.result$noisewidth * 3,
    min.subset = min.subset, prom.ratio = prom.ratio, give = "filter"
  )

  # Re-run tina_setup on the filtered features
  feature <- tina_setup(fse.result$storm_features[filt], fse.result$xmat)

  # Do spec-feature extraction for features

  # Align to feature maximum
  feature.ma <- align.max(feature, scaling = FALSE)

  #
  saveRDS(feature.ma, paste0(tmpdir, "/feature.RDS"))

  # The TINA part (skip for the moment) ####
  # ##### UMAP -> APcluster ####

  results <- tina_combineFeatures(feature.ma$stack,
    doUMAP = TRUE,
    umap.n_neighbors = pars$tina$umap$n_neighbors, # 30,
    umap.min_dist = pars$tina$umap$min_dist, # 0.8,
    ap.rcutoff = pars$tina$affinity.propagation$rcutoff,
    ap.q = pars$tina$affinity.propagation$q, # if not converging, try reducing a little.
    ap.lam = pars$tina$affinity.propagation$lambda, # oscillations possible. dampening factor
    plot.loc = paste0(this.run, "/"),
    plot.name = paste0("clusters_umap_ap_r", pars$tina$affinity.propagation$rcutoff, ".pdf")
  )

  clusters <- list(
    results = results,
    cluster.labs = clusts2labs(results$clusters)
  )

  saveRDS(clusters, paste0(this.run, "/clusters.RDS"))

  message("-------------------------------------------------------")
  message("-------------------  TINA Complete  -------------------")
  message("-------------------------------------------------------")
}
