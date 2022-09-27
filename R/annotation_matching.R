#' Annotation Module
#'
#' Given picked peaks from statistical decomposition module, do a simple
#' comparison between each query chemical shift list (STOCSY-derived)
#' and each hmdb reference spectrum chemical shift list. Optimize the
#' matching using the Hungarian Algorithm for linear assignment, after
#' all pairs with distances exceeding the threshold are removed.
#' @param peaks initial peaks picked (driver for picked peaks and target)
#' @param target full resolution STOCSY result. One element per STOCSY cluster,
#' containing ppm vector and thresholded intensity for each ppm value
#' @param ppeaks peak picked STOCYS result. One element per STOCSY cluster,
#' containing a numeric vector: col1: (named "chemical shift") ppm values,
#' col2: (named "intensity") intensity (covariance from STOCSY).
#' @param spec spectral data (also ppm vector).#     File should contain:
#'  List containing compound entries from hmdb or other db, SORTED by hmdb ID:
#'     spec_ID
#'       > metadata
#'         [various fields]
#'       > spectrum_peaks
#'         double (n x 2 array = n peaks x (ppm + intensity in STOCSY))
#'       > biospecimen
#'         (char array; e.g.
#' "Blood" "Cerebrospinal Fluid (CSF)" "Feces" "Urine" )
#'       > dsource
#'         (char, data source e.g. 'spec_file' or 'metabocard')
#' @param params params.yaml file, see package desc
#' @param refdb reference database containing spectra to match against,
#' loaded into memory.
#' 
#'
#' @return nothing, export matches is called.
#' @export
annotation_matching <- function(
  peaks, target, ppeaks, spec, params, refdb) {


  # Reformat metadata (this may eventually be unnecessary):

  metadata <- plyr::rbind.fill(
    lapply(1:length(refdb), function(x) refdb[[x]]$metadata))

  # Run the matching function:
  message(
    paste0(
      "\n\nMatching STOCSY clusters using ",
       params$am_pars$refdb_file, " as the reference...\n"))

  matches <- pblapply::pblapply(
    1:length(target),
    function(x) {
      matchToMultiRef(
        target = ppeaks[[x]], # target (picked peaks from STOCSY decomp specs)
        driver_ppm = spec[1, peaks[x]], # driver_ppm
        references = refdb, # refs, driver_ppm - ?
        metadata = metadata,
        tol = params$am_pars$dist_thresh,
        matchMethod = params$am_pars$matchMethod
      )
    }
  )

  message(
    "\nData written to matchPairs.tsv, referenceList.tsv, and targetList.tsv .")
  message("\nAnnotation matching process completed.\n\n")

  return(matches)
}
