#' Calculate the probability that a peak is "bad"
#'
#' Detect likelihood of each point being a bad peak
#' mean of the % of the feature signal that had no
#' corresponding ref signal (extra peak in feature)
#' useful to scale this by some exponential (e.g. 2) to select extremely high percentages
#' (higher values does less squashing of the bad peaks)
#' opf.mat <- lapply(allmatches.fits, function(x) x$overshoot.pct.feat) %>% do.call(rbind,.)
#' output is a vector which can be (e.g.) multiplied by and then subtracted from fit$fit.feat
#' @param pos.res.pct.feat matrix of positional residue differences between a 
#'   reference and a feature spectrum, as output by \code{positionalResidue}.
#' @param scale.exp a scaling exponent for the probability values (default is 2)
#'
#' @return A vector of probabilities, one for each peak in the feature spectrum,
#'   indicating the probability that the peak is "bad".
#'
#'
#' @export
#' @importFrom stats colMeans
#' @importFrom graphics heatmap
#' @importFrom magrittr %>%
bad.peaks <- function(pos.res.pct.feat, scale.exp = 2){
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