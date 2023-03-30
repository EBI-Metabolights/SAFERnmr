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
#'  FSE: the spectral matrix is decomposed into compound features using feature
#'  shape extraction. First, a local STOCSY is performed at every point along the
#'  spectrum (within a sliding window of ~ 100 points; enough to capture multiple
#'  resonances within any multiplet). For each of these STOCSYs, the central peak
#'  in the correlation profile (correlation pocket; corrpocket) typically captures
#'  a resonance, as the correlation is 1 by definition at the STOCSY'd point, and
#'  typically falls off as you approach the boundaries of the resonance. This is
#'  taken, with the next highest correlation peak within the window, to form a
#'  rough statistical description of two resonances which have an correlation in
#'  intensity across samples, albeit separated by chemical shift. We term this a
#'  'protofeature'. Importantly, each point and its associated window will capture
#'  the dominant protofeature most associated with that point. This changes for
#'  adjacent points, and many protofeatures will be duplicated multiple times. A
#'  protofeature should be considered as a rough hypothesis about a statistical
#'  association between two resonances, which happen to be sufficiently aligned
#'  so that they produce a coherent signal.
#'  If the windows are all aligned, we can plot the % of central corrpeaks containing
#'  each window point. From this distribution, it is clear that nearly all
#'  protofeatures, including those from noise peaks and real peaks, include the
#'  most central window points. As such, these cannot reliably be used to identify
#'  noise. However, the correlation peak about noise tends to be much smaller, and
#'  characteristically so. As such, a noisewidth can be estimated from this
#'  distribution. This is the origin of the noise.percentile cutoff, which is applied
#'  like so: "given a noise.percentile = 0.99, consider only those protofeatures
#'  for which both peaks have a width > the smallest 1% of peaks". Reducing this
#'  number therefore gives a more selective cut. Protofeatures are also filtered so
#'  that the
#'
#'  STORM: Joram Posma's STORM has been adapted and optimized to accept these
#'  protofeatures in the following ways:
#'  - first, since many of the protofeatures are noise,
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
    run_params <- yaml::yaml.load_file(filepath, eval.expr=TRUE)
    default <- TRUE
  } else {
    # load supplied params
    run_params <- yaml::yaml.load_file(params_loc, eval.expr=TRUE)
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
