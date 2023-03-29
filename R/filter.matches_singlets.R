#' Filter matches with singlets
#'
#' This function removes matches that have only one peak in their respective feature, reference or feature-not-never-fit regions.
#'
#' @param match.info The match information data.frame obtained within filter.matches.
#' @param fits.feature A list of fitted features.
#' @param peak.qualities A list of peak quality vectors (~ feature points' relevance in the reference database)
#' @param pq.featureNumbers A vector of feature numbers to index the peak.qualities list. This is necessary when only a subset of features are matched and feature fits are filtered.
#' @param res.area.thresh The minimum reference spectrum resonance area that must be accounted for by the fit feature in order to consider it matched.
#'
#' @return A list containing the filtered match information and the filtered fitted features.
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlist list.unlist
#'
#' @export
filter.matches_singlets <- function(match.info, fits.feature, peak.qualities, pq.featureNumbers, res.area.thresh) {


  ########################################################################################################################
  ## Removing singlets ####

  # # Dependencies ####
  #   source('./../extractPeaks_corr.R')
  #   source('./../pk.maxs.R')
  #   source('./../pk.bounds.R')
  #   source('./../corr_expand.R')
  #   source('./../ind2subR.R')
  #   source('./../stackplot.R')


  ########### singlet filter for fit features ################

  message("filtering out matches for fit features with 1 or fewer resonances...")
  match.info[, "numpeaks.feat"] <- lapply(fits.feature, function(ff) {
    return(pk.maxs(ff$feat.fit, mask = !is.na(ff$residuals)) %>% length())
  }) %>% unlist()
  fits.feature <- fits.feature[match.info$numpeaks.feat > 1]
  match.info <- dplyr::filter(match.info, numpeaks.feat > 1)
  message(nrow(match.info), " matches remaining\n\n")

  ########### singlet filter for fit ref regions ################

  message("filtering out matches involving ref features with 1 or fewer resonances...")
  match.info[, "numpeaks.ref"] <- lapply(fits.feature, function(ff) {
    return(pk.maxs(ff$spec.fit, mask = !is.na(ff$residuals)) %>% length())
  }) %>% unlist()
  fits.feature <- fits.feature[match.info$numpeaks.ref > 1]
  match.info <- dplyr::filter(match.info, numpeaks.ref > 1)
  message(nrow(match.info), " matches remaining\n\n")

  ########### singlet filter for feature-not-never-fit regions ################

  message("filtering out matches involving not-never-fit feature regions with 1 or fewer resonances...")
  match.info[, "numpeaks.feat.nnf"] <- lapply(1:nrow(match.info), function(m) {
    ff <- fits.feature[[m]]
    fnum <- match.info[m, "feat"]

    peak.quality <- peak.qualities[[which(fnum == pq.featureNumbers)]]
    f.adj <- ff$feat.fit - ff$feat.fit * peak.quality
    f.adj <- f.adj - min(f.adj, na.rm = T)
    return(pk.maxs(f.adj, mask = !is.na(f.adj)) %>% length())
  }) %>% unlist()
  fits.feature <- fits.feature[match.info$numpeaks.feat.nnf > 1]
  match.info <- dplyr::filter(match.info, numpeaks.feat.nnf > 1)
  message(nrow(match.info), " matches remaining\n\n")

  ########### singlet filter for refpeaks.matched ################

  message(
    "filtering out matches for which 1 or fewer resonances had at least ", pars$matching$filtering$res.area.threshold,
    " of their area explained by the matching feature..."
  )
  res.area.threshold <- res.area.thresh
  match.info[, "refpeaks.matched"] <- lapply(fits.feature, function(ff) {
    # ff <- fits.feature[[m]]
    pks <- extractPeaks_corr(ff$spec.fit, plots = F)
    pk.coverage <- lapply(pks$bounds[pks$truePeak], function(p) {
      # p <- pks$bounds[[2]]
      pkinds <- p %>%
        unlist() %>%
        fillbetween()
      return(sum(ff$feat.fit[pkinds], na.rm = T) / sum(ff$spec.fit[pkinds], na.rm = T))
    })
    return(sum(pk.coverage > (1 - res.area.threshold) & pk.coverage < (1 + res.area.threshold)))
  }) %>% unlist()
  fits.feature <- fits.feature[match.info$refpeaks.matched > 1]
  match.info <- dplyr::filter(match.info, refpeaks.matched > 1)
  message(nrow(match.info), " matches remaining\n\n")

  return(list(
    match.info = match.info,
    fits.feature = fits.feature
  ))
}
