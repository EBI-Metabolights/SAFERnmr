###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
## Feature shape based annotation pipeline ####
## /
# To run from here, you need the following files in ./data:
  # params.yaml         # parameter file
  # spectral.matrix.RDS # spectral matrix
  # data.list.RDS       # library files
  # lib.info.RDS        # index for library files
  
# Running the setup script will accomplish this, as well as ensuring that the correct 
# library is being used. 
# This pipeline accepts a spectral matrix (preferably with minimal processing; e.g.
# no strong normalization, no scaling, no alignment) in the form of and RDS file 
# with a matrix (row 1 = ppm vector, each additional row is a 1D NMR sample). 
# 
# FSE: ####
# the spectral matrix is decomposed into compound features using feature
# shape extraction. First, a local STOCSY is performed at every point along the 
# spectrum (within a sliding window of ~ 100 points; enough to capture multiple
# resonances within any multiplet). For each of these STOCSYs, the central peak
# in the correlation profile (correlation pocket; corrpocket) typically captures 
# a resonance, as the correlation is 1 by definition at the STOCSY'd point, and 
# typically falls off as you approach the boundaries of the resonance. This is 
# taken, with the next highest correlation peak within the window, to form a 
# rough statistical description of two resonances which have an correlation in 
# intensity across samples, albeit separated by chemical shift. We term this a 
# 'protofeature'. Importantly, each point and its associated window will capture 
# the dominant protofeature most associated with that point. This changes for 
# adjacent points, and many protofeatures will be duplicated multiple times. A
# protofeature should be considered as a rough hypothesis about a statistical 
# association between two resonances, which happen to be sufficiently aligned 
# so that they produce a coherent signal. 
# If the windows are all aligned, we can plot the % of central corrpeaks containing
# each window point. From this distribution, it is clear that nearly all 
# protofeatures, including those from noise peaks and real peaks, include the 
# most central window points. As such, these cannot reliably be used to identify
# noise. However, the correlation peak about noise tends to be much smaller, and 
# characteristically so. As such, a noisewidth can be estimated from this 
# distribution. This is the origin of the noise.percentile cutoff, which is applied
# like so: "given a noise.percentile = 0.99, consider only those protofeatures 
# for which both peaks have a width > the smallest 1% of peaks". Reducing this 
# number therefore gives a more selective cut. Protofeatures are also filtered so
# that the 
# 
# STORM: Joram Posma's STORM has been adapted and optimized to accept these
# protofeatures in the following ways: 
# - first, since many of the protofeatures are noise, we provide failure modes 
#   and reporting for the following cases: 
#   "empty subset",          # empty subset (no spectrum contains signature)
#   "subset degenerated",    # 1-3 spectra in the subset (not enough spectra to 
#                              get a reliable correlation)
#   "reference degenerated", # signature degenerates to include < 3 points (not
#                              meaningful to correlate shapes)
#   "did not converge"       # subset continues to change after 24 iterations
#
# - additionally, the correlation r and p-value cutoff q are both used during 
#   both the subset selection and reference update steps. 
# - we also remove any regions of the reference for which there are fewer than 
#   minpeak values after r and p value thresholding. This helps avoid noise. 
# - finally, STORM uses the max covariance value of the ref shape in each iteration
#   as the driver. This can cause the algorithm to get lured away from small 
#   features, to nearby features with larger covariance. To mitigate this effect,
#   we require that each driver update step remain on the same covariance signal
#   peak throughout all iterations. Bear in mind that other nearby peaks will be
#   assessed if they are part of a correlation pocket. 
#   
# STORM extracts meaningful features using protofeatures to define the region of 
# interest and a rough sketch of the feature shape highly correlated with each 
# spectral point. In the future, HCA could be used to cluster potential starting
# feature shapes correlated with each driver, or the nonoptimal subset for each 
# point could be re-STORMed to detect any other feature shapes present. It is 
# also perfectly reasonable to combine feature shapes from different STORM runs 
# for a given dataset, as these comprise a list of somewhat independently tested
# feature shapes, and duplication is not an issue. 
#   
# TINA (TINA Is Not Alignment): ####
#
# We cannot eliminate all poor shapes, but there are a couple of useful heuristics
# which generally reduce unnecessary computation downstream. First, there are many 
# feature shapes which are either quite poor in quality, or do not contain 
# sufficient information to be useful for annotation. We keep features with:
#   - in defined ppm range (generally [-1, 11])
#   - long enough runs? (contains runs > noisewidth*3 adjacent points)
#   - large enough subset? (>= 5 spectra with feature, good correlation reliability)
#   - has at least 1 true peak > 0.3 of the range of intensities? (not just the
#     side of a broad peak; not monotonic)
# Since the same features will often be extracted multiple times (either in the
# same spectral region, or other regions; i.e. same peak, but misaligned), it is
# advantageous to reduce this redundancy by clustering feature shapes. We accomplish 
# this using a combination of UMAP projection and Affinity Propagation clustering.
# Before comparing feature shapes, we align them to their maximum intensity
# resonance. This is quick, and usually performs well enough. 
# UMAP uses euclidean distance as a default. In practice, UMAP is able to sort 
# feature shapes into tight clusters when low n_neighbors (e.g. 5) and min_dist 
# (e.g. 0.05) are used. This does not capture global relationships as well, but 
# we only use it to identify tight clusters.
# All pairwise euclidean distances are computed using 
# See 
# Ultimately, it doesn't matter much what clustering method is used, as this is 
# primarily a means of combining highly similar features to reduce the computational
# burden of pairwise comparisons to reference spectra. 
# 
# Matching: ####
# Matching is accomplished by cross-correlating all pairwise feature - reference 
# spectrum pairs. Since there are feature-level comparisons to make (i.e. a single
# feature across all matches), iteration over features is serial. Iteration over
# references for each feature is done in parallel. This comparison is not optimal,
# but works as a proof of concept and can be scaled up. Increasing numbers of 
# features results in a linear increase in computational time, but also more memory 
# usage for each core. Generally, it is advisable to leave 1 or 2 cores on one's 
# machine free for system operations and the main R instance. For each comparison 
# that is, for each feature - reference pair being tested, the top max.hits 
# (e.g., 5) convolution hits are assessed for pearson correlation coefficient 
# (r.thresh) and pvalue (p.thresh). For hits passing both thresholds, the feature
# is fit to the reference using least squares. RMSE is reported. 
# 
# Next, it's best to avoid penalizing a feature fit just because it contains false 
# positive points. Once a feature has been fit to all reference spectra (at 0-5  
# more places), resonances which were never fit are identified using peak poorness:
# 1) Line up all fits for that feature across the database. 
# 2) For each fit, take the positive residuals, divide by feature intensity. If 
#    values are close to 1, the feature was completely absent in the ref
#    signature at those point. This is a measure of ~ peak poorness for this database
# 3) Take the mean of those values across all fits for each point in the feature. 
# 4) Square peak-poorness to squash them down unless very high:
#     peak.quality = 1-(peak.poorness^2)
# 5) For each fit's non-missing values:
#     rmse.weighted = sum(residuals * peak.quality)/n.points
#     
# This metric appears to give more intuitive results. Missing reference resonances
# compared to the feature is preferable to the opposite, which will be measured 
# during backfitting.
# Matches are written to file. 
# Note: matching in parallel is much lighter (~50%) if run in a new R instance  
# outside of RStudio due to memory leaks/overhead. 
# 
#
# MTJ 2023 
# Pipeline wrapper function ####
#'  Feature shape based annotation pipeline ####
#'  /
#'  To run from here, you need paths to the following files:
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
#' 
#' 
#' @return NULL
#' 
#' 
#' @importFrom utils write.table
#' @import magrittr
#' @importFrom yaml yaml.load_file
#' 
#' @export
pipeline <- function(params_loc, params_obj) {
  
  # status <- NULL
  # status <- tryCatch({
                if (isFALSE(missing(params_obj))) {
                  run_params <- params_obj
                } else if (missing(params_loc)) {
                  # load default params
                  filepath <- base::system.file(
                    "extdata", "default_params.yaml",
                    package = "SAFER"
                  )
                  run_params <- yaml::yaml.load_file(filepath, eval.expr = TRUE)
                  default <- TRUE
                  
                } else {
                  # >>> This is the usual route <<<
                  # load supplied params
                  run_params <- yaml::yaml.load_file(params_loc, eval.expr = TRUE)
                  dir.create(run_params$dirs$temp, showWarnings = F)
                  file.copy(params_loc, paste0(run_params$dirs$temp,'/params.yaml'))
                }
              
                # if (run_params$galaxy$enabled == FALSE) {
                #   setup(run_params)
                # }
                
                pars <- run_params
                
                if (is.null(pars$dirs$temp)){
                  pars$dirs$temp <- '.'
                } else {
                  if(!dir.exists(pars$dirs$temp)){
                    dir.create(pars$dirs$temp, showWarnings = F)
                  }
                }
              
                # Don't use lib.info path from now on. Assume it's in tmpdir
                file.copy(pars$files$lib.info, paste0(pars$dirs$temp,'/lib.info.RDS'))
              
  #           }, error = function(cond){return('failed setup')})
  # 
  # if (!is.null(status)){return(status)}
  
################################################################################################################################### 
## Feature Shape Extraction

  # status <- tryCatch({
    fse(pars)
  #   }, error = function(cond){return('failed fse')})
  # 
  # if (!is.null(status)){return(status)}
  
################################################################################################################################### 
## TINA / SAFARI
# - filter out feature shapes which make no sense
# - associate features whose shapes are highly similar

  # status <- tryCatch({
    tina(pars)
    # }, error = function(cond){return('failed tina')})
  # 
  # if (!is.null(status)){return(status)}
  
################################################################################################################################### 
## Match spectra to the database 
# - match using convolution-based cross-correlation
# - parallelized
  
  # status <- tryCatch({
              match_features2refs_par_setup(pars)
  #           }, error = function(cond){return('failed match setup')})
  # 
  # if (!is.null(status)){return(status)}
  
  gc() # garbage collect before big parallel compute
  
  # status <- tryCatch({
              match_features2refs_par_explicit(pars)
  #           }, error = function(cond){return('failed matching')})
  # 
  # if (!is.null(status)){return(status)}
  
################################################################################################################################### 
## Match filtering
# - filter out singlet matches
# - filter for (deltappm)
# - backfit ref subsignatures to individual dataset spectra

  # status <- tryCatch({
              filter_matches(pars) 
  #           }, error = function(cond){return('failed match filtering')})
  #   
  # if (!is.null(status)){return(status)}
  
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
#   - calls a function "pair_score_summation.R" which has a parallelized section
#   *** Format as MAF file and print match plots on request ***

  # status <- tryCatch({
              score_matches(pars)
  #           }, error = function(cond){return('failed match scoring')})
  #   
  # if (!is.null(status)){return(status)}
  

  # matches <- readRDS(paste0(this.run, "/matches_scored_named.RDS"))
  # write.table(matches, sep = '\t', file = 'matches_scored_named.txt', row.names = F, col.names = T)
  message('Saving session info...')
  saveRDS(sessionInfo(), paste0(pars$dirs$temp, "/session.info.RDS"))
  status <- 'success'
  return(status)
}