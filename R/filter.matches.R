#' Match Filtering
#'
#' Filter and process matched peaks based on user-specified criteria, such as singlet removal and ppm distance filtering.
#' Also calculate backfits of ref-feats (reference subsignatures fished out by feature matches to reference spectra) to each individual dataset spectrum.
#' Note: ref-feat fit is obtained using the original feature, then backfit feasibility scores are calculated. BFF scores indicate the extent to which ref-feat resonances do not exceed actual spectral signal. Specifically, each single ref-feat resonance positive residual is evaluated as a fraction of the total ref-feat height. Because the absence of a feasible fit for any single reference resonance invalidates a match, the worst-violating resonance gives the score for the whole ref-feat in a particular spectrum. Note that a ref-feat receives a BFF score for each spectrum it is fit to.
#'
#' @param pars A list of input parameters.
#' @return A list of filtered and processed matched peak information, including back-fits to the original spectra.
#' @import yaml, magrittr, pbapply
#' @export filter.matches
filter.matches <- function(pars) {
  message("--------------------------------------------------------------")
  message("-------------------     Match Filtering    -------------------")
  message("--------------------------------------------------------------")
  message("\n\n\n")

  ################ Read parameters file ##################


  tmpdir <- pars$dirs$temp
  this.run <- paste0(tmpdir)

  ##################################################################################################################
  # Read data and set up ####

  message("Loading data from files...\n\n\n")

  fse.result <- readRDS(paste0(this.run, "/fse.result.RDS"))
  feature <- readRDS(paste0(this.run, "/feature.RDS"))
  matches <- readRDS(paste0(this.run, "/matches.RDS"))

  # Format matches ####

  matches.split <- split(matches, names(matches))
  rm(matches)
  match.info <- do.call(rbind, matches.split$matches)
  rownames(match.info) <- NULL

  # peak.qualities aligns with fit and match.info now, but feat number does not index it (there are missing features)!
  pq.featureNumbers <- unique(match.info[, "feat"]) # this does not sort (just for good measure)

  peak.qualities <- matches.split$peak.quality
  fits.feature <- matches.split$fits %>% unlist(recursive = F)
  rm(matches.split)


  # Do the filtering (functionalized)
  res <- filter.matches_singlets(
    match.info, fits.feature,
    peak.qualities, pq.featureNumbers,
    pars$matching$filtering$res.area.threshold
  )
  match.info <- res$match.info
  fits.feature <- res$fits.feature

  res <- filter.matches_shiftDelta(match.info, feature,
    ppm = fse.result$ppm, fits.feature,
    ppm.tol = pars$matching$filtering$ppm.tol
  )
  match.info <- res$match.info
  fits.feature <- res$fits.feature

  ######################### Back-fit reference to spectra  #############################

  message("Back-fitting ref-feats to each spectrum in the relevant subset...\n\n")
  xmat <- fse.result$xmat
  ppm <- fse.result$ppm

  # Back-fit each matched reference region to the subset spectra

  m.inds <- 1:nrow(match.info)
  t1 <- Sys.time()
  backfits <- backfit_ref.feats.2.subset.specs(m.inds, fits.feature, match.info,
    feature,
    xmat, ppm,
    plots = F
  ) # plots are heavy and time-expensive!
  Sys.time() - t1
  message("Saving backfits...\n\n\n")
  saveRDS(backfits, paste0(this.run, "/backfits.RDS"))

  # ########## save filtered data ########################################################################

  message("Saving split and filtered match data...\n\n")

  saveRDS(match.info, paste0(this.run, "/match.info.RDS"))

  saveRDS(fits.feature, paste0(this.run, "/fits.RDS"))

  saveRDS(peak.qualities, paste0(this.run, "/peak.qualities.RDS"))

  message("-----------------------------------------------------------------")
  message("-----------------  Matching Filtering Complete ------------------")
  message("-----------------------------------------------------------------")
}
