#' Statistical Decomposition for spectra.
#'
#' Given a spectral matrix and a list of picked peaks, loop
#' through each picked peak and perform a STOCSY using it as
#' a driver peak. Threshold the result at r >= 0.8, giving
#' pseudospectra (target list). Pick peak these using standard settings,
#' and store the peak clusters in ppeaks.
#'
#' @param peaks picked peaks for representative spectrum, numeric vector
#'  containing ppm indices of peaks
#' @param spec spectral matrix consisting of n columns (one for each ppm
#'  value). row  1 (named 'ppm') is the vector of ppm values. rows 2:m+1
#' (names undefined) are the m sample measurements for each ppm value.
#' @param params yaml object containg parameters for the pipeline.
#' @return list of:
#'   target.RDS: list with one element for each cluster, containing:
#'       [1] ppm vector
#'       [2] thresholded intensities for each ppm value
#'   ppeaks.RDS: list with picked peaks from each cluster.
#'   One element per cluster, containing:
#'       numeric vector
#'         col 1: (named "chemical-shift") ppm values
#'         col 2: (named "intensity") intensities (covariance from STOCSY)
#' @export
stat_decomp <- function(peaks, spec, params) {

    # Do the STOCSY deconvolution
    message("\nRunning STOCSY on all provided peaks...")
    target <- pbapply::pblapply(
        peaks, function(x) stocsydec(spec, x, params$sd_pars$cutoff)
    )

    # Peak pick the results
    message("\nPeak picking STOCSY results to get clusters...")
    ppeaks <- pbapply::pblapply(target, function(x) ppTargetSpec(x))

    # Save data
    message("\nSaving results. Please wait...")
    message("\nData written to target.RDS, ppeaks.RDS, and peaks.RDS .")
    message("\nStatistical Decomposition process completed.\n\n\n")
    return(list(
        target = target,
        ppeaks = ppeaks
    ))
}
