#' statDecomp
#'
#' Michael T. Judge and Gonçalo Graça
#' v.0.1 19MAY2022
#' Description:
#'   Given a spectral matrix and a list of picked peaks, loop through each picked peak
#'   and perform a STOCSY using it as a driver peak. Threshold the result at r >= 0.8,  
#'   giving pseudospectra (target list). Pick peak these using standard settings, 
#'   and store the peak clusters in ppeaks.
#'
#' Input files:
#'   peaks.RDS: picked peaks for representative spectrum
#'     numeric vector containing ppm indices of peaks
#'   spec.RDS: spectral matrix consisting of n columns (one for each ppm value)
#'     row  1 (named 'ppm') is the vector of ppm values
#'     rows 2:m+1 (names undefined) are the m sample measurements for each ppm value
#' Output files:
#'   target.RDS: list with one element for each cluster, containing:
#'       [1] ppm vector
#'       [2] thresholded intensities for each ppm value
#'   ppeaks.RDS: list with picked peaks from each cluster. One element per cluster, containing:
#'       numeric vector
#'         col 1: (named "chemical-shift") ppm values
#'         col 2: (named "intensity") intensities (covariance from STOCSY)
#' Source dependencies:
#'   stocsydec.R
#'   ppTargetSpec.R
#' @param peaks_loc location of the peaks .RDS file.
#' @param spec_loc location of the spectra ,RDS file.
#' @param params workflow parameters as a .yaml file.
#' @return list object containing both picked peaks and target outputs.
#' @export
stat_decomp <- function(peaks_loc, spec_loc, params) {
  message(paste0("\n\n",
                "------------------------------------------------------------------\n",
                "---------------- Statistical Decomposition Module ----------------\n",
                "---------------- 	v1.0                      ----------------\n",
                "------------------------------------------------------------------\n\n"))
  library(stringr)
  library(speaq)
  library(yaml)
  library(pbapply)
  library(devtools)
  # Read params
    pars <- yaml.load_file(params)


  # we need some logic here that uses these presupplied .rds files as test cases,
  # but defers to user supplied .rds files if supplied.

  # Load data
    if (missing(peaks_loc)) {
      peaks <- readRDS(str_c("./data/","peaks.RDS"))
    } else {
      peaks <- readRDS(peaks_loc)
    }

    if (missing(spec_loc)) {
      spec <- readRDS(str_c("./data/","spec.RDS"))
    } else {
      spec <- readRDS(spec_loc)
    }

  ########################################################################################  

  # Do the STOCSY deconvolution 
    message("\nRunning STOCSY on all provided peaks...")
    target <- pblapply(peaks, function(x) stocsydec(spec, x, pars$sd_pars$cutoff))

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
