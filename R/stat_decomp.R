#' Statistical Decomposition for spectra.
#'
#' Description
#' @param peaks_loc location of the rds peaks file.
#' @param spec_loc location of the rds spec file.
#' @param params location of the yaml params file.
#'
#' @return list of stuff
#' @export
stat_decomp <- function(peaks_loc, spec_loc, params) {
    library(stringr)
    # put speaq back in here.
    library(yaml)
    library(pbapply)
    library(devtools)
    library(MassSpecWavelet)
    # Read params
    pars <- yaml.load_file(params)

    # Load data
    if (missing(peaks_loc)) {
        peaks <- readRDS(str_c("./data/", "peaks.RDS"))
    } else {
        peaks <- readRDS(peaks_loc)
    }

    if (missing(spec_loc)) {
        spec <- readRDS(str_c("./data/", "spec.RDS"))
    } else {
        spec <- readRDS(spec_loc)
    }

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
    return(c(target, ppeaks, recursive = TRUE))
}
