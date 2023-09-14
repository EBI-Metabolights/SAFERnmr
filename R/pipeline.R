###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
## Feature shape based annotation pipeline ####
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
              
               pars <- run_params
                
                if (is.null(pars$dirs$temp)){
                  pars$dirs$temp <- '.'
                } else {
                  if(!dir.exists(pars$dirs$temp)){
                    dir.create(pars$dirs$temp, showWarnings = F)
                  }
                }
              
                # Don't use lib.info path from now on. Assume it's in tmpdir
                # file.copy(pars$files$lib.info, paste0(pars$dirs$temp,'/lib.info.RDS'))
              
  #           }, error = function(cond){return('failed setup')})
  # 
  # if (!is.null(status)){return(status)}
  info <- sessioninfo::package_info(pkgs = 'SAFER')
  
  message('------------------------------------------------------------------------------')
  message('-----------------------------       SAFER        -----------------------------')
  message('-----------------------------       v',info$loadedversion,'       -----------------------------')
  message('------------------------------------------------------------------------------')
  printTime()
  
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
  
  zip(files = paste0(pars$dirs$temp, '/debug_extra.outputs'), 
      zipfile = paste0(pars$dirs$temp, '/debug_extra.outputs.zip'))
  unlink(paste0(pars$dirs$temp, '/debug_extra.outputs'), recursive = TRUE)
  
  zip(files = paste0(pars$dirs$temp, '/plots'), 
      zipfile = paste0(pars$dirs$temp, '/plots.zip'))
  unlink(paste0(pars$dirs$temp, '/plots'),recursive = TRUE)
  
  zip(files = paste0(pars$dirs$temp, '/temp_data_matching'), 
      zipfile = paste0(pars$dirs$temp, '/temp_data_matching.zip'))
  unlink(paste0(pars$dirs$temp, '/temp_data_matching'), recursive = TRUE)
  
  message('Saving session info...')
  saveRDS(sessionInfo(), paste0(pars$dirs$temp, "/session.info.RDS"))
  status <- 'success'
  return(status)
}