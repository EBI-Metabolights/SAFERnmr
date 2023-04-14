#' Setup function for initializing the analysis environment
#'
#' This function performs various setup tasks, such as copying library data and metadata, and finding the best fit for the spectrometer frequency.
#'
#' @param params_obj An object containing various input parameters, such as the study and library directories.
#'
#' @return None
#'
#' @import yaml
#' @import magrittr
#' @importFrom base message
#' @importFrom base file.copy
#' @importFrom base stop
#' @importFrom utils readRDS
#'
#' @export
#'
#' @examples
#' setup(params_obj)
#'
#' @keywords setup, initialization, data
#' @rdname setup
setup <- function(params_obj) {


  # ###################################################################################################################################

  message("-------------------------------------------------------")
  message("-------------------      Setup      -------------------")
  message("-------------------------------------------------------")
  message("\n\n\n")


  library(yaml)
  library(magrittr)


  ##################################################################################################################
  ## Library data copy ####

  # Which spectrometer frequency out of those available is closest to the dataset?
  # - get spec freq from data
  sf <- params_obj$study$spectrometer.frequency

  # - get spec freqs for available libraries
  field.strengths <- readRDS(paste0(params_obj$dirs$lib, "/available.field.strengths_MHz.RDS"))

  # - find best fit
  spec.freq <- (field.strengths - sf) %>%
    abs() %>%
    which.min() %>%
    field.strengths[.]
  lib.file <- paste0(params_obj$dirs$lib, "/data.list_", spec.freq, "MHz.RDS")

  # Library info file:
  message("Copying library info from ", params_obj$dirs$lib, "/gissmo_lib.info.RDS")
  message(" to ", paste0(this.run, "/lib.info.RDS"), "\n\n\n")

  copied <- file.copy(
    from = paste0(pars$dirs$lib, "/gissmo_lib.info.RDS"),
    to = paste0(this.run, "/lib.info.RDS"),
    overwrite = T
  )
  if (!copied) {
    stop(
      "Library info was not copied. Ensure that both file (",
      paste0(pars$dirs$lib, "/gissmo_lib.info.RDS"),
      ") and destination (",
      this.run, ") exist."
    )
  }

  # Library data:
  message("Copying library from ", lib.file)
  message(" to ", paste0(this.run, "/lib.data.RDS"), "\n\n\n")

  copied <- file.copy(
    from = lib.file,
    to = paste0(this.run, "/lib.data.RDS"),
    overwrite = T
  )
  if (!copied) {
    stop(
      "Library was not copied. Ensure that both file (",
      lib.file,
      ") and destination (",
      this.run, ") exist."
    )
  }

  ################################################################################################

  message("-------------------------------------------------------")
  message("--------------------  Setup Complete.  ----------------")
  message("-------------------------------------------------------")
}
