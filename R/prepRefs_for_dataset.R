#' Prepare reference spectra for matching to a dataset
#' To match features from a dataset to reference spectra in fine detail, they need
#' to have the same spectral axis, more or less. This function interpolates each
#' reference spectrum in the list to the dataset ppm axis, and then filters its
#' signal to avoid inclusion of noise/zeros, which are replaced with NAs so that
#' cross-correlations, etc. are minimized. In future versions, a solvent region
#' could be provided for each (or all) that will be replaced with NAs to avoid
#' matching that peak.
#'
#'  Takes gissmo (or other) reference spectra as a list with
#'   - tag (short name, e.g. "bmse000005" - not necessarily unique)
#'   - ref.name (should be unique with field strength etc.)
#'   - ppm axis
#'   - ref spectrum data
#'
#' @param data.list a list of reference spectra, where each element of the list is a list containing the fields "tag", "ref.name", "compound.name", "ppm", and "data"
#' @param ppm.dataset the ppm axis of the dataset to which the reference spectra will be matched
#' @param ref.sig.SD.cutoff the cutoff for reference spectra signal, expressed as a fraction of the standard deviation
#'
#' @return a list of processed reference spectra, where each element of the list is a list containing the fields "tag", "ref.name", "compound.name", "ppm", "data", and "mapped". The "mapped" field contains a list with the fields "ppm", "data", "sig.cutoff", and "sd.cutoff"
#'
#' @examples
#' data.list <- list(ref1, ref2, ref3)
#' ppm.dataset <- c(0, 1, 2, 3, 4)
#' ref.sig.SD.cutoff <- 0.01
#' prepRefs_for_dataset(data.list, ppm.dataset, ref.sig.SD.cutoff)
#'
#' @export
prepRefs_for_dataset <- function(data.list,                 # list of gissmo (or other) ref spectra and ppm axes
                                 ppm.dataset,                       # ppm axis of dataset that will be matched
                                 ref.sig.SD.cutoff = 0.01,  # cutoff for ref spectra signal (fraction of SD)
                                 n.cores = 1)               # this can be parallelized
  {
  
    
    # For each ref, process for this matching ####
      ppm <- ppm.dataset
      ref.sig.SD.cutoff <- 0.01
      
      data.list.processed <- 
          # pblapply(data.list, function(ref.spec)
            # ref.spec <- data.list[[2]]
          mclapply(data.list, function(ref.spec)
            {
              
              # Linear interpolation to dataset ppm axis ####
                # ref.spec <- data.list[[1]]
                #  simplePlot(ref.spec$data, xvect = ref.spec$ppm)

                ref.int <- approx(ref.spec$ppm, ref.spec$data, 
                                   ppm, # get values for these ppms (can be selective, or if goes oob, then NA)
                                   method="linear", rule=1, f=0)
                ref.ppm <- ref.int$x
                ref.data <- ref.int$y - min(ref.int$y, na.rm = TRUE)
                # simplePlot(ref.data, xvect = ref.ppm)
                
              # # If not in the original ppm range, interpolation makes no sense. 
              # # Replace with zeros.
              #   ppm.in.ref <- ref.spec$ppm %>% range %>% vectInds(., ref.ppm) %>% fillbetween
              #   ref.data[-ppm.in.ref] <- 0
                
              # Set anything below threshold to NA to avoid extraneous comparisons ####
                ncutoff <- median(ref.data, na.rm = TRUE) + sd(ref.data, na.rm = TRUE) * ref.sig.SD.cutoff
                ref.data[ref.data < ncutoff] <- NA
                ref.data <- ref.data - min(ref.data, na.rm = TRUE)
                ref.data <- ref.data / sum(ref.data, na.rm = T)
                
                data.compressed <- tryCatch(
                  {
                    co_compress(stack.list = list(data = ref.data, ppm = ref.ppm), 
                                sparse.val = NA, key = 'data')
                  }, 
                  error = function(cond){
                    NA
                  }
                )
                
                if (length(data.compressed) == 1){stop('compression failed for ', ref.spec$tag)}
                
                
              # Return a copy of the list element with mapped data added ####
                return( list(tag = ref.spec$tag,
                             ref.name = ref.spec$ref.name,
                             compound.name = ref.spec$compound.name,
                             id = ref.spec$id,
                             # ppm = ref.spec$ppm,
                             # data = ref.spec$data,
                             mapped = list(ppm = NULL,
                                           data = NULL,
                                           data.compressed = data.compressed,
                                           sig.cutoff = ncutoff,
                                           sd.cutoff = ref.sig.SD.cutoff),
                             info = ref.spec$info
                            )
                       )
            }, mc.cores = n.cores
          )
      
      # ref.mat <- lapply(data.list.processed, function(x) x$mapped$data) %>% do.call(rbind, .)
        # simplePlot(ref.mat)
      # lapply(data.list.processed, function(x) length(x$mapped$data.compressed) == 1) %>% unlist %>% any
  
  return(data.list.processed)
  
}