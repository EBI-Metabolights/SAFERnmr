#' Statistical Decomposition for spectra.
#'
#' Description
#' @param peaks rds peaks file loaded into memory.
#' @param spec rds spec file loaded into memory.
#' @param params yaml object containg parameters for the pipeline.
#'
#' @return list of stuff
#' @export
stat_decomp <- function(peaks, spec, params) {

    # Do the STOCSY deconvolution
    message("\nRunning STOCSY on all provided peaks...")
    target <- pblapply::pblapply(
        peaks, function(x) stocsydec(spec, x, params$sd_pars$cutoff)
    )

    # Peak pick the results
    message("\nPeak picking STOCSY results to get clusters...")
    ppeaks <- pblapply::pblapply(target, function(x) ppTargetSpec(x))

    # Save data
    message("\nSaving results. Please wait...")
    message("\nData written to target.RDS, ppeaks.RDS, and peaks.RDS .")
    message("\nStatistical Decomposition process completed.\n\n\n")
    return(list(
        target = target,
        ppeaks = ppeaks
    ))
}
