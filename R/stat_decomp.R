#' Statistical Decomposition for spectra.
#'
#' Description
#' @param peaks_loc rds peaks file loaded into memory.
#' @param spec_loc rds spec file loaded into memory.
#' @param params yaml object containg parameters for the pipeline.
#'
#' @return list of stuff
#' @export
stat_decomp <- function(peaks, spec, params) {
    library(stringr)
    # put speaq back in here.
    library(yaml)
    library(pbapply)
    library(devtools)
    library(MassSpecWavelet)
    # Read params

    # want to later must have the params object passed into the method,
    # rather than loading it.
    pars <- yaml.load_file(params)

    # Do the STOCSY deconvolution
    message("\nRunning STOCSY on all provided peaks...")
    target <- pblapply(
        peaks, function(x) stocsydec(spec, x, pars$sd_pars$cutoff)
    )

    # Peak pick the results
    message("\nPeak picking STOCSY results to get clusters...")
    ppeaks <- pblapply(target, function(x) ppTargetSpec(x))

    # Save data
    message("\nSaving results. Please wait...")
    saveRDS(target, "./data/target.RDS")
    saveRDS(ppeaks, "./data/ppeaks.RDS")
    saveRDS(peaks, "./data/peaks.RDS")

    message("\nData written to target.RDS, ppeaks.RDS, and peaks.RDS .")
    message("\nStatistical Decomposition process completed.\n\n\n")
    return(list(
        target = target,
        ppeaks = ppeaks
    ))
}
