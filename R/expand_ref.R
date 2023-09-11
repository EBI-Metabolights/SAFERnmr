#' expand.ref
#' 
#' Expand library data processed
#' 
#' For a given ref object, expand the compressed data and ppm vector. The ppm
#' vector will only expand for non-na values; thus, just replace with ppm for 
#' the dataset (it's just a copy, anyways). 
#' 
#' @param ref list of dataset-processed library data for a compound
#' @param ppm dataset ppm vector the data were interpolated to
#' 
#' @return ref with mapped$data and mapped$ppm expanded, ready for use
#' @export
#' 
expand_ref <- function(ref, ppm){
  dat <- ref$mapped$data.compressed %>% expand_stacklist
    ref$mapped$data <- dat$data
    ref$mapped$ppm <- ppm
  return(ref)
}