#' Parameter validation (using params object)
#'
#' Go through a bunch of checks to see if params are reasonable
#'
#' @param pars A list of input parameters.
#' @return Adjusted parameters
#' @importFrom magrittr %>%
#' @importFrom rrapply rrapply
#' 
#' @import pbapply
#' 
#' @export
validate_pars <- function(pars){
  
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
          rrapply::rrapply(pars,
                           condition = \(x) is.character(x),
                           f = \(x) detect.convert.scientificNotation(x),
                           how = "replace")
          
        }, error = function(cond){
          
          return(pars)
          
        }
      )
      
      return(pars)
    }
  
  # Blanket checks: ####
  
    # Correct any scientific notations that didn't get converted by yaml
      pars <- convert_sciNotation_pars(pars)
    
  
  # Field-specific Validators: ####
  
  
  # Detect any parameters outside of normal ranges
  
    # Validate dirs
    if(!is.null(dirs)){
      if(!is.null(dirs$temp)){
        stopifnot(is.character(dirs$temp))
        stopifnot(file.exists(dirs$temp))
      }
      
      if(!is.null(dirs$lib)){
        stopifnot(is.character(dirs$lib))
        stopifnot(file.exists(dirs$lib))
      }
    }
    
    # Validate study
    if(!is.null(study)){
      if(!is.null(study$id)){
        stopifnot(is.character(study$id))
      }
      
      if(!is.null(study$spectrometer.frequency)){
        stopifnot(is.numeric(study$spectrometer.frequency))
      }
    }
    
    # Validate files
    if(!is.null(files)){
      if(!is.null(files$spectral.matrix)){
        stopifnot(is.character(files$spectral.matrix))
        stopifnot(file.exists(files$spectral.matrix))
      }
      
      if(!is.null(files$lib.data)){
        stopifnot(is.character(files$lib.data))
        stopifnot(file.exists(files$lib.data))
      }
      
      if(!is.null(files$lib.info)){
        stopifnot(is.character(files$lib.info))
        stopifnot(file.exists(files$lib.info))
      }
    }
    
    # Validate corrpockets
    if(!is.null(corrpockets)){
      if(!is.null(corrpockets$half.window)){
        stopifnot(is.numeric(corrpockets$half.window))
      }
      
      if(!is.null(corrpockets$noise.percentile)){
        stopifnot(is.numeric(corrpockets$noise.percentile))
      }
      
      if(!is.null(corrpockets$only.region.between)){
        stopifnot(is.null(corrpockets$only.region.between))
      }
      
      if(!is.null(corrpockets$rcutoff)){
        stopifnot(is.numeric(corrpockets$rcutoff))
      }
    }
    
    # Validate storm
    if(!is.null(storm)){
      if(!is.null(storm$correlation.r.cutoff)){
        stopifnot(is.numeric(storm$correlation.r.cutoff))
      }
      
      if(!is.null(storm$q)){
        stopifnot(is.numeric(storm$q))
      }
      
      if(!is.null(storm$b)){
        stopifnot(is.numeric(storm$b))
      }
      
      if(!is.null(storm$number.of.plots)){
        stopifnot(is.numeric(storm$number.of.plots))
      }
    }
    
    # Validate tina
    if(!is.null(tina)){
      if(!is.null(tina$bounds)){
        stopifnot(is.numeric(tina$bounds))
      }
      
      if(!is.null(tina$min.subset)){
        stopifnot(is.numeric(tina$min.subset))
      }
      
      if(!is.null(tina$prom.ratio)){
        stopifnot(is.numeric(tina$prom.ratio))
      }
      
      if(!is.null(tina$do.clustering)){
        stopifnot(is.logical(tina$do.clustering))
      }
      
      if(!is.null(tina$clustering$max.eps)){
        stopifnot(is.numeric(tina$clustering$max.eps))
      }
      
      if(!is.null(tina$clustering$minPts)){
        stopifnot(is.numeric(tina$clustering$minPts))
      }
      
      if(!is.null(tina$clustering$eps.stepsize)){
        stopifnot(is.numeric(tina$clustering$eps.stepsize))
      }
      
      if(!is.null(tina$nfeats)){
        stopifnot(is.numeric(tina$nfeats))
      }
      
      if(!is.null(tina$plots$max.plots)){
        stopifnot(is.numeric(tina$plots$max.plots))
      }
      
      if(!is.null(tina$plots$filtered.out)){
        stopifnot(is.logical(tina$plots$filtered.out))
      }
      
      if(!is.null(tina$plots$filtered.features)){
        stopifnot(is.logical(tina$plots$filtered.features))
      }
      
      if(!is.null(tina$plots$cleaned.clusters)){
        stopifnot(is.logical(tina$plots$cleaned.clusters))
      }
    }
    
    # Validate matching
    if(!is.null(matching)){
      if(!is.null(matching$cluster.profile)){
        stopifnot(is.character(matching$cluster.profile))
      }
      
      if(!is.null(matching$ref.sig.SD.cutoff)){
        stopifnot(is.numeric(matching$ref.sig.SD.cutoff))
      }
      
      if(!is.null(matching$max.hits)){
        stopifnot(is.numeric(matching$max.hits))
      }
      
      if(!is.null(matching$r.thresh)){
        stopifnot(is.numeric(matching$r.thresh))
      }
      
      if(!is.null(matching$p.thresh)){
        stopifnot(is.numeric(matching$p.thresh))
      }
      
      if(!is.null(matching$filtering$res.area.threshold)){
        stopifnot(is.numeric(matching$filtering$res.area.threshold))
      }
      
      if(!is.null(matching$filtering$ppm.tol)){
        stopifnot(is.numeric(matching$filtering$ppm.tol))
      }
      
      if(!is.null(matching$filtering$max.backfits)){
        stopifnot(is.character(matching$filtering$max.backfits))
      }
    }
    
    # Validate par
    if(!is.null(par)){
      if(!is.null(par$ncores)){
        stopifnot(is.numeric(par$ncores))
      }
      
      if(!is.null(par$type)){
        stopifnot(is.character(par$type))
      }
    }
    
    # Validate galaxy
    if(!is.null(galaxy)){
      if(!is.null(galaxy$enabled)){
        stopifnot(is.logical(galaxy$enabled))
      }
    }
    
    # Validate debug
    if(!is.null(debug)){
      if(!is.null(debug$enabled)){
        stopifnot(is.logical(debug$enabled))
      }
      
      if(!is.null(debug$throttle_features)){
        stopifnot(is.numeric(debug$throttle_features))
      }
      
      if(!is.null(debug$all.outputs)){
        stopifnot(is.logical(debug$all.outputs))
      }
    }

  return(pars)
}

 



