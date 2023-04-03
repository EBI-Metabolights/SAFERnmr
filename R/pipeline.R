#' Pipeline function to execute the entire ImperialNMRTool workflow
#'
#'  Feature shape based annotation pipeline
#'  /
#'  To run from here, you need the following files in ./data:
#'  params.yaml         # parameter file
#'  spectral.matrix.RDS # spectral matrix
#'  data.list.RDS       # library files
#'  lib.info.RDS        # index for library files
#'
#'  Running the setup script will accomplish this, as well as ensuring that the correct
#'  library is being used.
#'  This pipeline accepts a spectral matrix (preferably with minimal processing; e.g.
#'  no strong normalization, no scaling, no alignment) in the form of and RDS file
#'  with a matrix (row 1 = ppm vector, each additional row is a 1D NMR sample).
#'
#'
#' # MTJ 2023
#'
#' @param params_loc Path to a YAML file containing user-specified parameters
#' @param params_obj List of parameters
#' @return NULL
#' @importFrom yaml yaml.load_file
#' @importFrom base missing
#' @importFrom utils readRDS
#' @importFrom utils write.table
#' @importFrom utils saveRDS
#' @export
pipeline <- function(params_loc, params_obj) {
  if (isFALSE(missing(params_obj))) {
    run_params <- params_obj
  } else if (missing(params_loc)) {
    # load default params
    filepath <- system.file(
      "extdata", "default_params.yaml",
      package = "ImperialNMRTool"
    )
    run_params <- yaml::yaml.load_file(filepath, eval.expr = TRUE)
    default <- TRUE
  } else {
    # load supplied params
    run_params <- yaml::yaml.load_file(params_loc, eval.expr = TRUE)
  }


  ###################################################################################################################################
  ## Feature Shape Extraction

  fse(run_params)

  ###################################################################################################################################
  ## TINA / SAFARI
  # - filter out feature shapes which make no sense
  # - associate features whose shapes are highly similar

  tina(run_params)

  ###################################################################################################################################
  ## Match spectra to the database
  # - match using convolution-based cross-correlation
  # - parallelized

  match.features2refs.par(run_params)


  ###################################################################################################################################
  ## Match filtering
  # - filter out singlet matches
  # - filter for {pascalstriangle}(dirac)
  # - backfit ref subsignatures to individual dataset spectra

  filter.matches(run_params)


  ###################################################################################################################################
  ## Assess matches
  # - produce a score to summarize the quality of the matches:
  # - quantify the sum of the best evidence for each compound in each spectrum:
  #   - a ref is associated with a spectrum via multiple ref features
  #     (these are subsignatures)
  #   - a ref resonance can be covered by multiple ref features, each with a
  #     different backfit feasibility score ~[0,1]
  #   - for a given ref signature - dataset spectrum association, we want to add
  #     up the best evidence (i.e. use the evidence from the best-fitting ref
  #     features) at each point in the reference signature
  #   - Fraction of ref signature explained [0,1], scaled by the best backfit
  #     feasibility score at each point in the ref signature
  #   - calls a script "pair.score.summation.R" which has a parallelized section
  #   *** Format as MAF file and print match plots on request ***

  score.matches(run_params)

  matches <- readRDS(paste0(run_params$dirs$temp, "/matches_scored_named.RDS"))
  write.table(matches, sep = "\t", file = "matches_scored_named.txt", row.names = F, col.names = T)
  saveRDS(sessionInfo(), paste0(run_params$dirs$temp, "/session.info.RDS"))
}
