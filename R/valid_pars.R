#' Parameter validation (using params object)
#'
#' Correct any obviously incorrect params. Replace any empties with default params.
#' 
#' Then go through a bunch of checks to see if params are reasonable. 
#'
#' @param pars A list of input parameters.
#' @return Adjusted parameters
#' @importFrom magrittr %>%
#' 
#' @export
valid_pars <- function(pars){
  # Accessory functions: ####
    detect.convert.scientificNotation <- function(x){
      
      matches <- grep("\\b[+-]?[0-9]+(\\.[0-9]+)?[eE][+-]?[0-9]+\\b", 
           x, perl = TRUE, value = TRUE)
      
      # If string contains one instance (expected)
      if (length(matches) == 1 && length(matches[[1]]) > 0) {
        
        return(as.numeric(matches[[1]]))
        
      } 
      
      else {
        
        return(x)
        
      }
      
    }
    convert_sciNotation_pars <- function(pars){
      
      pars <- tryCatch(
        {
          # rrapply::rrapply(pars,
          #                  condition = \(x) is.character(x),
          #                  f = \(x) detect.convert.scientificNotation(x),
          #                  how = "replace")
          rapply(pars,
                 class = c("character"),
                 function(x) detect.convert.scientificNotation(x),
                 how = 'replace')
          
        }, error = function(cond){
          
          return(pars)
          
        }
      )
      
      return(pars)
    }
    in_range <- function(a, b, c){
      bc <- sort(c(b,c))
      bc[1] <= a & a <= bc[2]
    }
    check_each_par_value <- function(pars){
  
      error_messages <- NULL

      # Validate dirs ####
      
        dirs <- pars$dirs
        if(!is.null(dirs)){
          if(!is.null(dirs$temp)){
            if(!is.character(dirs$temp)){
              error_messages <- c(error_messages, paste("'dirs$temp' must be a character string, but it is of type", class(dirs$temp)))
            }
            if(!file.exists(dirs$temp)){
              error_messages <- c(error_messages, paste("'dirs$temp' does not exist as a file path:", dirs$temp))
            }
          }
          
          if(!is.null(dirs$lib)){
            if(!is.character(dirs$lib)){
              error_messages <- c(error_messages, paste("'dirs$lib' must be a character string, but it is of type", class(dirs$lib)))
            }
            if(!file.exists(dirs$lib)){
              error_messages <- c(error_messages, paste("'dirs$lib' does not exist as a file path:", dirs$lib))
            }
          }
        }
        
      
      # Validate files ####
      
        files <- pars$files
        if(!is.null(files)){
          if(!is.null(files$spectral.matrix)){
            if(!is.character(files$spectral.matrix)){
              error_messages <- c(error_messages, 
                                  paste("'files$spectral.matrix' must be a character string, but it is of type", class(files$spectral.matrix)))
            }
            if(!tools::file_ext(files$spectral.matrix) == 'RDS'){
              error_messages <- c(error_messages, 
                                  paste("'files$spectral.matrix' should be the full path to an .RDS file containing a spectral matrix with ppm values in the first row, and corresponding intensity values for each spectrum in the following rows."))
            }
            if(!file.exists(files$spectral.matrix)){
              error_messages <- c(error_messages, 
                                  paste("'files$spectral.matrix' does not exist as a file path:", files$spectral.matrix))
            }
          }
          
          if(!is.null(files$lib.data)){
            if(!is.character(files$lib.data)){
              error_messages <- c(error_messages, 
                                  paste("'files$lib.data' must be a character string, but it is of type", class(files$lib.data)))
            }
            if(!tools::file_ext(files$lib.data) == 'RDS'){
              error_messages <- c(error_messages, 
                                  paste("'files$lib.data' should be the full path to .RDS file containing a library data in the SAFER format."))
            }
            if(!file.exists(files$lib.data)){
              error_messages <- c(error_messages, 
                                  paste("'files$lib.data' does not exist as a file path:", files$lib.data))
            }
          }
          
          ## lib.info is no longer used
          # if(!is.null(files$lib.info)){
          #   if(!is.character(files$lib.info)){
          #     error_messages <- c(error_messages, 
          #                         paste("'files$lib.info' must be a character string, but it is of type", class(files$lib.info)))
          #   }
          #   if(!file.exists(files$lib.info)){
          #     error_messages <- c(error_messages, 
          #                         paste("'files$lib.info' does not exist as a file path:", files$lib.info))
          #   }
          # }
          
        }
        
      
      # Validate corrpockets ####
        
        corrpockets <- pars$corrpockets
        if(!is.null(corrpockets)){
          if(!is.null(corrpockets$half.window)){
            if (!is.numeric(corrpockets$half.window)) {
              error_messages <- c(error_messages, paste("corrpockets$half.window should be numeric. Current value:", 
                                                        corrpockets$half.window))
            } else {
              if(!in_range(corrpockets$half.window, 0.02, 0.1)){
                  c(error_messages, 
                    paste0("corrpockets$half.window should be a numeric representing the expected upper bound of first-order coupling constants in the spectrum (in ppm). In other words: in ppm units, how far apart do you expect the resonances in a doublet or triplet to be? This is normally 0.02 - 0.08 ppm. ",
                           " Current value: ", corrpockets$half.window, " ppm.")
                    )
                }
            }
          }
          
          if(!is.null(corrpockets$noise.percentile)){
            if (!is.numeric(corrpockets$noise.percentile)) {
              error_messages <- c(error_messages, paste("corrpockets$noise.percentile should be numeric. Current value:", 
                                                        corrpockets$noise.percentile))
            } else {
                if (!in_range(corrpockets$noise.percentile, 0.9, 0.9999)){
                  error_messages <- c(error_messages, 
                                      paste0("corrpockets$noise.percentile normally sits between 0.9 and 0.9999.", 
                                             " Current value: ", corrpockets$half.window))
                }
            }
          }
          
          if(!is.null(corrpockets$only.region.between)){
            if (!is.null(corrpockets$only.region.between)) {
              error_messages <- c(error_messages, 
                                  paste("corrpockets$only.region.between should be 2 numerics provided as e.g. c(0, 10) or empty(NULL). ",
                                        "Current value:",corrpockets$only.region.between))
            } 
          }
          
          if(!is.null(corrpockets$rcutoff)){
            if (!is.numeric(corrpockets$rcutoff)) {
              error_messages <- c(error_messages, paste("corrpockets$rcutoff should be numeric. Current value:", corrpockets$rcutoff))
            } else {
                if (!in_range(corrpockets$rcutoff, 0.5, 0.9999)){
                  error_messages <- c(error_messages, 
                                      paste0("corrpockets$rcutoff normally sits between 0.5 and 0.9. More permissive values are recommended, since STORM threshold will have a much greater effect on feature shape.", 
                                             " Current value: ", corrpockets$rcutoff))
                }
            }
          }
        }
      
      
      # Validate storm ####
        
        storm <- pars$storm
        if(!is.null(storm)){
          if(!is.null(storm$correlation.r.cutoff)){
            if (!is.numeric(storm$correlation.r.cutoff)) {
              error_messages <- c(error_messages, paste("storm$correlation.r.cutoff should be numeric. Current value:", storm$correlation.r.cutoff))
            } else {
                if (!in_range(storm$correlation.r.cutoff, 0.7, 0.9)){
                  error_messages <- c(error_messages, 
                                      paste0("storm$correlation.r.cutoff normally sits between 0.7 and 0.9. More permissive values are not recommended, since STORM is intended to pull out highly similar feature shapes, and this will compromise them and lead to much greater numbers of less useful matches. However, if the cutoff is too high, it will lead to far fewer features, and may result in insufficient spectral coverage.", 
                                             " Current value: ", storm$correlation.r.cutoff))
                }
            }
          }
          
          if(!is.null(storm$q)){
            if (!is.numeric(storm$q)) {
              error_messages <- c(error_messages, paste("storm$q should be numeric. Current value:", storm$q))
            }
          }
          
          if(!is.null(storm$b)){
            if (!is.numeric(storm$b)) {
              error_messages <- c(error_messages, paste("storm$b should be numeric. Current value:", storm$b))
            } else {
                if (!in_range(storm$b, 1, 3)){
                  error_messages <- c(error_messages, 
                                      paste0("storm$b is a window size for STORM window expansion, as a multiple of the distance between the primary and secondary correlation peaks. Normally between 1 and 3.",
                                             " Current value: ", storm$b))
                }
            }
          }
          
          if(!is.null(storm$number.of.plots)){
            if (!is.numeric(storm$number.of.plots)) {
              error_messages <- c(error_messages, paste("storm$number.of.plots should be numeric. Current value:", storm$number.of.plots))
            }
          }
        }
      
      
      # Validate matching ####
      
        matching <- pars$matching
        if(!is.null(matching)){
        if(!is.null(matching$cluster.profile)){
          if (!is.character(matching$cluster.profile)) {
            error_messages <- c(error_messages, paste("matching$cluster.profile should be a character. Current value:", matching$cluster.profile))
          } else {
              if (!any(c("representative.feature") %in% matching$cluster.profile)){
                error_messages <- c(error_messages,
                                    paste0("matching$cluster.profile is currently only allowed to be 'representative.feature'",
                                           " Current value: ", storm$cluster.profile))
              }
          } 
        }
        
        if(!is.null(matching$ref.sig.SD.cutoff)){
          if (!is.numeric(matching$ref.sig.SD.cutoff)) {
            error_messages <- c(error_messages, paste("matching$ref.sig.SD.cutoff should be numeric. Current value:", matching$ref.sig.SD.cutoff))
          }
        }
        
        if(!is.null(matching$max.hits)){
          if (is.numeric(matching$max.hits)){
            matching$max.hits <- matching$max.hits %>% round %>% as.integer
          }
          if (!is.integer(matching$max.hits)) {
            error_messages <- c(error_messages, paste("matching$max.hits should be an integer > 0. Current value:", matching$max.hits))
          } else {
              if (!in_range(matching$max.hits, 1, 10)){
                error_messages <- c(error_messages,
                                    paste0("matching$max.hits should be an integer > 0 to indicate the number of times a feature can be matched to the same PCRS/spectrum. Typical values are 3-5. This allows similar hits to be evaluated in more detail at the fitting step, instead of simply choosing the one with the largest convolution value (not always the best indicator of a high-quality match).",
                                           " Current value: ", matching$max.hits))
              }
          } 
        }
        
        if(!is.null(matching$r.thresh)){
          if (!is.numeric(matching$r.thresh)) {
            error_messages <- c(error_messages, paste("matching$r.thresh should be numeric. Current value:", matching$r.thresh))
          } else {
              if (!in_range(matching$r.thresh, 0.7, 1)){
                error_messages <- c(error_messages,
                                    paste0("matching$r.thresh is typically no lower than 0.7, and cannot exceed 1.",
                                           " Current value: ", matching$r.thresh))
              }
          } 
        }
        
        if(!is.null(matching$p.thresh)){
          if (!is.numeric(matching$p.thresh)) {
            error_messages <- c(error_messages, paste("matching$p.thresh should be numeric. Current value:", matching$p.thresh))
          }
        }
        
        if(!is.null(matching$filtering$res.area.threshold)){
          if (!is.numeric(matching$filtering$res.area.threshold)) {
            error_messages <- c(error_messages, paste("matching$filtering$res.area.threshold should be numeric. Current value:", matching$filtering$res.area.threshold))
          } else {
              if (!in_range(matching$filtering$res.area.threshold, .1, .5)){
                error_messages <- c(error_messages,
                                    paste0("matching$filtering$res.area.threshold is typically between 0.1 and 0.5. It represents the resonance area accounted for necessary to consider a it 'matched' for the purpose of detecting singlet matches.",
                                           " Current value: ", matching$filtering$res.area.threshold))
              }
          } 
        }
        
        if(!is.null(matching$filtering$ppm.tol)){
          if (!is.numeric(matching$filtering$ppm.tol)) {
            error_messages <- c(error_messages, paste("matching$filtering$ppm.tol should be numeric. Current value:", matching$filtering$ppm.tol))
          } else {
              if (!in_range(matching$filtering$ppm.tol, 0, .5)){
                error_messages <- c(error_messages,
                                    paste0("matching$filtering$ppm.tol is typically between 0.001 and 0.1. Increasing this greatly increases the number of potential matches for each feature.",
                                           " Current value: ", matching$filtering$ppm.tol))
              }
          } 
        }
        
        if(!is.null(matching$filtering$max.backfits)){
          if (!is.numeric(matching$filtering$max.backfits)) {
            error_messages <- c(error_messages, paste("matching$filtering$max.backfits should be a numeric Current value:", matching$filtering$max.backfits))
          } else {
              if (!in_range(matching$filtering$max.backfits, 1E3, 1E7)){
                error_messages <- c(error_messages,
                                    paste0("matching$filtering$max.backfits a key parameter for preventing the pipeline from reaching astronomical compute times on a given platform. For 48 cores and 500 Gb of RAM, 1E7 seems to be a good upper limit. 1E8 still won't finish after a few days. We are exploring the degree to which this can be reduced without dramatically changing the quality of results.",
                                           " Current value: ", matching$filtering$max.backfits))
              }
          } 
        }
      }
      
      
      # Validate par ####
        
        par <- pars$par
        if(!is.null(par)){
          if(!is.null(par$ncores)){
            if(is.numeric(par$ncores)){
              par$ncores <- as.integer(par$ncores %>% round)
            }
            if (!is.integer(par$ncores)) {
              error_messages <- c(error_messages, paste("par$ncores should be an integer. Current value:", par$ncores))
            } else {
                if (!in_range(par$ncores, 1, parallel::detectCores()-1)){
                  error_messages <- c(error_messages,
                                      paste0("par$ncores sets the number of cores available for parallel operations in the pipeline. Ensure there are enough cores available (checked here using the range: c(1, parallel::detectCores()-1)), which evaluates to ", c(1, parallel::detectCores()-1),", and that there is enough RAM to allocate to each one. The average dataset will require up to 10-15Gb per core, depending on several factors such as the reference library size, the number of features, the size of the spectral matrix (number of samples and number of points), the average subset size, etc.. ",
                                             " Current value: ", par$ncores))
                }
            } 
          }
          
          if(!is.null(par$type)){
            if (!is.character(par$type)) {
              error_messages <- c(error_messages, paste("par$type should be either FORK or PSOCK, and sets the cluster type for the matching loop. All other parallel operations should be using mclapply; see documentation for that if more detail is required. Current value:", par$type))
            }
          }
        }
      
      
      # Validate galaxy ####
        
        galaxy <- pars$galaxy
        if(!is.null(galaxy)){
          if(!is.null(galaxy$enabled)){
            if (!is.logical(galaxy$enabled)) {
              error_messages <- c(error_messages, paste("galaxy$enabled should be logical. Current value:", galaxy$enabled))
            }
          }
        }
      
      
      # Validate debug ####
      
        debug <- pars$debug
        if(!is.null(debug)){
        if(!is.null(debug$enabled)){
          if (!is.logical(debug$enabled)) {
            error_messages <- c(error_messages, paste("debug$enabled should be logical. Current value:", debug$enabled))
          }
        }
        
        if(!is.null(debug$throttle_features)){
          if (!is.numeric(debug$throttle_features)) {
            error_messages <- c(error_messages, paste("debug$throttle_features should be numeric. Current value:", debug$throttle_features))
          }
        }
        
        if(!is.null(debug$all.outputs)){
          if (!is.logical(debug$all.outputs)) {
            error_messages <- c(error_messages, paste("debug$all.outputs should be logical. Current value:", debug$all.outputs))
          }
        }
      }
      
      # ####
      
      return(error_messages)
    }
  
  # Blanket checks: ####
  
    # Correct any scientific notations that didn't get converted by yaml
      pars <- convert_sciNotation_pars(pars)
    
    # Convert any empties or missings to defaults
    
      
      
  # Field-specific Validators: ####
  
    warnings <- check_each_par_value(pars)
    
    # Print all error messages
    
    validation.pass <- length(warnings) == 0
    
  # Report ####
  
    if(!validation.pass) {
      
      cat("Parameter validation warnings:\n")
      i <- 0
      for (msg in warnings) {
        cat(i, ' - ', msg, "\n\n")
        i <- i + 1
      }
      
    } else {
      
      cat("Validation successful. No errors found.\n")

    }
    
    return(validation.pass)
    
}

 



