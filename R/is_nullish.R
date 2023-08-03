#' is_nullish: contaminated with NaN, NULL, Inf, or empty values?
#' 
#' Check sure if elements are actually usable values (i.e.
#' not: empty, NaN, Inf, or NULL). For vectors and matrices, this 
#' behaves like a composite is.na() like function.
#' 
#' Works on basic types: numeric, integer, complex, character, logical
#' And data structures: vector, matrix, list, data.frame, table, nested list,
#' factor
#' 
#' If x is a list, this will recursively check each element and 
#' return TRUE for each one if all elements are usable values. Top level
#' is preserved. *If error message from lapply, returns TRUE.
#' 
#' Examples:
#'    is_nullish(c())
#'    
#'    obj <- c(1.2i,4, NaN)
#'    is_nullish(obj)
#'    
#'    obj <- c(1.2,4, NaN)
#'    is_nullish(obj) 
#'    
#'    obj <- c(2,4, NAN)
#'    is_nullish(obj)
#'    
#'    obj <- c(T,F,T)
#'    is_nullish(obj)
#'    
#'    obj <- list(a = 1, b = NULL, c = NA, d = list(a =  NULL, b = 1))
#'    is_nullish(obj)
#'    
#'    tab <- matrix(c(7, 5, 14, NA, 3, 2, NaN, 6, Inf), ncol=3, byrow=TRUE)
#'    colnames(tab) <- c('colName1','colName2','colName3')
#'    rownames(tab) <- c('rowName1','rowName2','rowName3')
#'    is_nullish(tab)
#'    is_nullish(as.table(tab))
#'    
#'    obj <- data.frame(filed1 = NaN,
#'                      f2 = 2)
#'    is_nullish(obj)
#'    
#'    obj <- factor(c("yes", "no", "no", "yes", "yes"))
#'    is_nullish(obj)
#'
#'#' @param x an object of (any?) type
#'
#' @return whether the object has NaN, NULL, empty, or Inf values
#'
#' @export
is_nullish <- function(x){
  
  # Check if is EMPTY or NULL
    if (length(x) == 0){
      
      return(TRUE)
      
    } else {
      
      
        if (is.list(x)){
          # If list, then recursively check:
          
            return(lapply(x, function(x.elem) is_nullish(x.elem) %>% as.logical %>% any) %>% unlist)
          
        } 
      
        if (is.table(x)){
          
          # If table, then coerce to vect
            
            return(x %>% c %>% is_nullish %>% any)
          
        } else {
          # If not a list, then check if elements are Inf or NaN
            return(is.nan(x) | is.infinite(x))
          
        }
      
    }
  
}

test_nullish <- function(obj, msg=NULL){
  
  if (any(is_nullish(obj))){
    call.txt <- deparse(sys.call(-1))
    #deparse(sys.calls()[[sys.nframe()-1]]) # function name
    if (is.null(msg)){
      return(stop('"',call.txt, '" ...object is nullish.', call. = F))
    }
  } else {
      return(stop(msg, call. = F))
  }
}
