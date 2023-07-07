#' Accessory function to print timestamp to message throughout pipeline
#'
#' @return timestamp (message)
#' 
#' @export
printTime <- function(){
  message('\n\ttimestamp : ', Sys.time(), '\n')
}