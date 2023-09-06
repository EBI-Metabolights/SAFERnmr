#' debug_write #########################################################################################################################
#'
#' If saving all outputs to the pipeline, then save obj under f.name in  directory
#' Uses pars$debug$enabled to check if debugging. If so, check pars$debug$all.outputs 
#' to see whether or not to save. 
#' 
#' Saving location is pars$dirs$temp /debug_extra.outputs 
#'
#' @param obj some object you want to save if debugging
#' @param f.name RDS filename the object should be saved under (no directory info, just name with extension)
#' @param pars parameters from the run; contains info on whether of not to save
#'
#' @export
debug_write <- function(obj, f.name, pars){
  tryCatch(
    {
        if (pars$debug$enabled){
          if (pars$debug$all.outputs){
            dir.create(paste0(pars$dirs$temp, "/debug_extra.outputs"), showWarnings = F)
            saveRDS(obj, paste0(pars$dirs$temp, "/debug_extra.outputs/", f.name))
          }
        }
      return(paste0('wrote "', f.name, '" to debug_extra.outputs...'))
    },
    error = function(cond){
      return(paste0('warning: debug_write(): "', f.name, '" could not be saved.'))
    }
  )
}
