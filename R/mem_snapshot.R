#' Write a file with current R memory usage
#'
#' 
#' @param fname filename to write to (usually timestamped)
#' @return writes a file with current R memory usage
#'
#' @export
mem_snapshot <- function(fname){
  fileConn<-file(fname)
  writeLines(paste0(as.character(round(pryr::mem_used() / 1000000000, 2)), ' Gb'), fileConn)
  close(fileConn)
}