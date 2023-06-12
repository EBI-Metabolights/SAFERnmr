#' Calculate how often each feature point is fit across the reference database
#'
#' Get weights for feature profile points according to how often they fit against reference DB peaks.
#' mean of the % of the feature signal that had no corresponding ref signal (extra peak in feature)
#' useful to scale this by some exponential (e.g. 2) to select extremely high percentages
#' lapply(allmatches.fits, function(x) x$overshoot.pct.feat) %>% do.call(rbind,.)
#'
#' @param pos.res.pct.feat matrix of positive residuals for a feature to all references fit as a percent of the feature resonance height
#' @param scale.exp a scaling exponent for the weights to select for values very close to 1 (default is 2)
#'
#' @importFrom magrittr %>%
#'
#' @return A vector of weights, one for each point in the feature profile,
#' indicating how rarely the feature point was fit in the ref DB matches the
#' which can be (e.g.) multiplied by and then subtracted from feature
#' to squash not-never-fit peaks in the feature, or (1-result) * rmse values to
#' downweight the effect of not-never-fit peaks (those not represented in the ref DB,
#' and more likely to be false positives in feature extraction).
#'
#'
#'
#' @export
bad_peaks <- function(pos.res.pct.feat, scale.exp = 2){
        mat <- pos.res.pct.feat #res.mat
        # mat <- mat * allmatches.feat[, 'rval']
        # mat <- mat * (1/allmatches.feat[, 'rmse'])
        # mat[mat<0] <- 0
        mat[is.na(mat)] <- 0
        
        # heatmap(mat, Colv = NA, scale = 'none')

      bad.peaks <- colMeans(mat, na.rm = T)^scale.exp #* fit$feat.fit
      # bad.peaks[is.na(bad.peaks)] <- 0
      # simplePlot(bad.peaks^10, xdir = "n")
      # simplePlot(rbind(fit$feat.fit, fit$feat.fit * bad.peaks))
      # simplePlot(rbind(fit$feat.fit - fit$feat.fit * bad.peaks, fit$feat.fit))
    # Just record the probability that it's a bad peak
      # bad.peaks %>% plot
    return(bad.peaks)
}