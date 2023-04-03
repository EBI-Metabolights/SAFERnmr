#' Match a feature against a reference using FFT-based convolution and correlation
#' picks the top max.hits peaks in the convolution, then calculates correlation and
#' pvalue of the PCC at those points in the ref.
#' Returns a dataframe with the xcorr information.
#'
#' @param f.num numeric: Feature number
#' @param r.num numeric: Reference number
#' @param feat numeric: Vector of feature intensities
#' @param ref numeric: Vector of reference intensities
#' @param feat.ft.c numeric: FFT of feature intensities, conjugated
#' @param ref.ft numeric: FFT of reference intensities
#' @param pad.size numeric: Size of the padding applied during FFT-based convolution
#' @param max.hits numeric: Maximum number of candidate lags to consider
#' @param r.thresh numeric: Correlation threshold for considering a match
#' @param p.thresh numeric: P-value threshold for considering a match
#'
#' @return A data frame with the following columns:
#'   - feat: Feature number or identifier
#'   - ref: Reference number or identifier
#'   - lag: The lag that produces the highest correlation
#'   - rval: The correlation coefficient
#'   - pval: The p-value of the correlation
#'   - pts.matched: Number of data points used in the correlation
#'   - pts.feat: Number of data points in the feature vector
#'   - feat.start: The starting index of the feature vector used in the correlation
#'   - feat.end: The ending index of the feature vector used in the correlation
#'   - ref.start: The starting index of the reference vector used in the correlation
#'   - ref.end: The ending index of the reference vector used in the correlation
#'
#'
#' @export
feature_match2ref_slim <- function(f.num, r.num, feat, ref,
                                   feat.ft.c, ref.ft,
                                   pad.size,
                                   max.hits = 5,
                                   r.thresh = 0.8, p.thresh = 0.01) {

  # Do the FFT-based conv/xcorr ####

  r <- (feat.ft.c * ref.ft) %>%
    fftw::FFT(., inverse = TRUE) %>%
    Re() %>%
    c()

  # Get maxima (candidate lags) ####

  lmxs <- localMaxima(r)

  # Sort maxima
  lags <- lmxs[order(r[lmxs], decreasing = T)] # sort by xcorr peak height

  # Restrict lags to those not inside padding (padding is really just for end
  # effects in the FT, not for actual comparison. Perhaps it's necessary to
  # pad the matrix twice?). Also only take the top n hits (to keep calculations
  # reasonable).

  lags <- lags[lags >= pad.size] %>% .[1:max.hits]

  # Loop though candidate lags and evaluate fit at each one ####

  feat <- t(c(feat))
  ref <- t(c(ref))
  inds.trim.feat <- trim.sides(feat, out = "inds")


  fits <- lapply(lags, function(lag) {
    # Calculate corr at each shift
    ref.pos <- lag - pad.size + inds.trim.feat
    use <- !is.na(feat[inds.trim.feat] + ref[ref.pos])

    # Make sure there are enough points to do a correlation:
    if (sum(use) < 3) {
      return(NULL)
    }

    r <- suppressWarnings(
      cor(feat[inds.trim.feat[use]],
        ref[ref.pos[use]],
        use = "pairwise.complete.obs",
        method = "pearson"
      )
    )
    return(data.frame(
      ref.start = min(ref.pos),
      ref.end = max(ref.pos),
      pts.matched = sum(use),
      rval = r
    ))
  }) %>% do.call(rbind, .)

  # In case of infinite rvals, set to zero:
  fits$rval[is.infinite(fits$rval)] <- 0
  r.pass <- fits$rval > r.thresh

  # Calculate pvals using t-distribution
  a <- -abs(fits$rval * sqrt((fits$pts.matched - 2) / (1 - fits$rval^2)))
  pvals <- 2 * pt(a, (fits$pts.matched - 2))
  p.pass <- pvals < p.thresh

  # average fit intensity as a fraction of the feature signal (want to fit parts that are dominant)

  r.p.pass <- which(r.pass & p.pass)

  if (!any(r.p.pass)) {
    return(NULL)
  }


  matches.ranked <- order(pvals[r.p.pass]) %>% r.p.pass[.]

  matches <- data.frame(
    feat = f.num,
    ref = r.num,
    lag = lags[matches.ranked],
    rval = fits$rval[matches.ranked],
    pval = pvals[matches.ranked],
    pts.matched = fits$pts.matched[matches.ranked],
    pts.feat = length(inds.trim.feat),
    feat.start = min(inds.trim.feat),
    feat.end = max(inds.trim.feat),
    ref.start = fits$ref.start[matches.ranked],
    ref.end = fits$ref.end[matches.ranked]
  )

  # Record results
  return(matches)
}
