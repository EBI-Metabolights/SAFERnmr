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
  
  start.time <- Sys.time()
  # status <- NULL
  # status <- tryCatch({
  # The goal of this section is to get:
  # - params_loc (an actual filepath of the params file)
  # - pars (the parameters used for the run)
  
                if (isFALSE(missing(params_obj))) {
                  pars <- params_obj
                  
                  # Write the params file to the temp dir
                  params_loc <- './params.yaml'
                  yaml::write_yaml(params_obj, file = params.loc)
                  
                } else if (missing(params_loc)) {
                  # load default params
                  filepath <- base::system.file(
                    "extdata", "default_params.yaml",
                    package = "SAFER"
                  )
                  
                  params_loc <- filepath
                  pars <- yaml::yaml.load_file(params_loc, eval.expr = TRUE)
                  default <- TRUE
                  
                } else {
                  # >>> This is the usual route <<<
                  # load supplied params
                  pars <- yaml::yaml.load_file(params_loc, eval.expr = TRUE)
                  
                }
               
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
  
               
  # Create a timestamped subdirectory under this tmp dir name and set tmpdir to that
    pars$dirs$temp <- paste0( pars$dirs$temp, '/', start.time %>% as.numeric %>% round )
    dir.create(pars$dirs$temp , showWarnings = F)  
    
    dir.create(pars$dirs$temp, showWarnings = F)
    file.copy(params_loc, paste0(pars$dirs$temp,'/params.yaml'), overwrite = TRUE)
               
  # Validate parameters before beginning
    
    pars.passed.checks <- valid_pars(pars)
    
    status <- 'completed setup'
    
    run.summary <- data.frame(local.id = start.time,
                              status = status,
                              pars.passed.checks = pars.passed.checks,
                              study = pars$study$id, 
                              spec.MHz = pars$study$spectrometer.frequency,
                              data = pars$files$spectral.matrix,
                              ref.library = pars$files$lib.data,
                              np = pars$corrpockets$noise.percentile,
                              protofeature.r = pars$corrpockets$rcutoff,
                              storm.r = pars$storm$correlation.r.cutoff,
                              storm.b = pars$storm$b,
                              min.subset = pars$tina$min.subset,
                              max.feats = pars$debug$throttle_features,
                              max.hits = pars$matching$max.hits,
                              match.r = pars$matching$r.thresh, # this will get updated by bf limit
                              ppm.tol = pars$matching$filtering$ppm.tol,
                              max.backfits = pars$matching$filtering$max.backfits,
                              par.ncores = pars$par$ncores,
                              galaxy = pars$galaxy$enabled,
                              verbose = pars$debug$all.outputs
                              )
    
  run.summary$write.time <- Sys.time()
  write.csv(x = run.summary, file = paste0(pars$dirs$temp, "/run.summary.csv"))
               
################################################################################################################################### 
## Feature Shape Extraction

  fse.summary <- tryCatch({
    
     run.sum <- fse(pars)

  }, error = function(cond){NULL})
  
  if (is.data.frame(fse.summary)){
    
    run.summary <- cbind(run.summary, fse.summary)
    run.summary$status <- 'completed fse'
    
  } else {
    
    run.summary$status <- 'failed fse' 
    
  }
  
  run.summary$write.time <- Sys.time()
  write.csv(x = run.summary, file = paste0(pars$dirs$temp, "/run.summary.csv"))
  
################################################################################################################################### 
## TINA / SAFARI
# - filter out feature shapes which make no sense
# - associate features whose shapes are highly similar

  tina.summary <- tryCatch({
    
     run.sum <- tina(pars)

  }, error = function(cond){NULL})
  
  if (is.data.frame(tina.summary)){
    
    run.summary <- cbind(run.summary, tina.summary)
    run.summary$status <- 'completed tina'
    
  } else {
    
    run.summary$status <- 'failed tina' 
    
  }
  
  run.summary$write.time <- Sys.time()
  write.csv(x = run.summary, file = paste0(pars$dirs$temp, "/run.summary.csv"))
  
  
################################################################################################################################### 
## Match spectra to the database 
# - match using convolution-based cross-correlation
# - parallelized
  
  # status <- tryCatch({
  
              
  match.setup.summary <- tryCatch({
    
     match_features2refs_par_setup(pars)

  }, error = function(cond){NULL})
  
  if (is.data.frame(match.setup.summary)){
    
    run.summary <- cbind(run.summary, match.setup.summary)
    run.summary$status <- 'completed match setup'
    
  } else {
    
    run.summary$status <- 'failed match setup' 
    
  }
  
  run.summary$write.time <- Sys.time()
  write.csv(x = run.summary, file = paste0(pars$dirs$temp, "/run.summary.csv"))
  
  gc() # garbage collect before big parallel compute
  
  
  match.summary <- tryCatch({
    
     match_features2refs_par_explicit(pars)

  }, error = function(cond){NULL})
  
  if (is.data.frame(match.summary)){
    
    run.summary <- cbind(run.summary, match.summary)
    run.summary$status <- 'completed matching'
    
  } else {
    
    run.summary$status <- 'failed matching'
    
  }
  
  run.summary$write.time <- Sys.time()
  write.csv(x = run.summary, file = paste0(pars$dirs$temp, "/run.summary.csv"))
  
################################################################################################################################### 
## Match filtering
# - filter out singlet matches
# - filter for (deltappm)
# - backfit ref subsignatures to individual dataset spectra


  match.filt.summary <- tryCatch({
    
     filter_matches(pars)

  }, error = function(cond){NULL})
  
  if (is.data.frame(match.filt.summary)){
    
    run.summary <- cbind(run.summary, match.filt.summary)
    run.summary$status <- 'completed match filtering/backfitting'
    
  } else {
    
    run.summary$status <- 'failed match filtering/backfitting'
    
  }
  
  run.summary$write.time <- Sys.time()
  write.csv(x = run.summary, file = paste0(pars$dirs$temp, "/run.summary.csv"))
  
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

  
  match.score.summary <- tryCatch({
    
     score_matches(pars)

  }, error = function(cond){NULL})
  
  if (is.data.frame(match.score.summary)){
    
    run.summary <- cbind(run.summary, match.score.summary)
    run.summary$status <- 'completed match scoring'
    
  } else {
    
    run.summary$status <- 'failed match scoring'
    
  }
  
  run.summary$write.time <- Sys.time()
  write.csv(x = run.summary, file = paste0(pars$dirs$temp, "/run.summary.csv"))
  
################################################################################################################################### 
## Summarize Results
# - produce a score to summarize the quality of the matches:

  run.summary$total.time <- as.character((Sys.time()-start.time) %>% round(1))   
  
  si <- Sys.info() %>% as.list()
  
  run.summary$safer.version <- packageVersion("SAFER") %>% as.character
  run.summary$r.version <- R.version$version.string
  run.summary$r.platform <- R.version$platform
  run.summary$system.version <- si$version
  run.summary$run.id <- start.time %>% as.numeric %>% round
  
  run.summary$write.time <- Sys.time()
  names(run.summary) <- names(run.summary) %>% stringr::str_replace_all('\\.', '_') # change dots for _ to make shell stuff easier on FTP
  print(t(run.summary))     
  write.csv(x = run.summary, file = paste0(pars$dirs$temp, "/run.summary.csv"))          
  # run.summary <- read.csv(paste0(pars$dirs$temp, "/run.summary.csv"))
  
################################################################################################################################### 
## Clean up 
# - produce a score to summarize the quality of the matches:
              
  paste0(pars$dirs$temp, '/debug_extra.outputs')
  if (dir.exists(paste0(pars$dirs$temp, '/debug_extra.outputs'))){
    zip(files = paste0(pars$dirs$temp, '/debug_extra.outputs'), 
        zipfile = paste0(pars$dirs$temp, '/debug_extra.outputs.zip'))
    unlink(paste0(pars$dirs$temp, '/debug_extra.outputs'), recursive = TRUE)
  }
  
  if (dir.exists(paste0(pars$dirs$temp, '/plots'))){
    zip(files = paste0(pars$dirs$temp, '/plots'), 
        zipfile = paste0(pars$dirs$temp, '/plots.zip'))
    unlink(paste0(pars$dirs$temp, '/plots'),recursive = TRUE)
  }
  
  if (dir.exists(paste0(pars$dirs$temp, '/temp_data_matching'))){
    zip(files = paste0(pars$dirs$temp, '/temp_data_matching'), 
        zipfile = paste0(pars$dirs$temp, '/temp_data_matching.zip'))
    unlink(paste0(pars$dirs$temp, '/temp_data_matching'), recursive = TRUE)
  }
  
  message('Saving session info...')
  saveRDS(sessionInfo(), paste0(pars$dirs$temp, "/session.info.RDS"))
  status <- 'success'
  return(status)
}