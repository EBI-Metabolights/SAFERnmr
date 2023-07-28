#' Sparse (NA-gapped) matrix compression
#'
#' Compresses matrix to its non-NA values and records their positions
#' to enable perfect reconstruction on-the-fly. Intended for sparse matrices
#' and parallel operations as well as general memory management.
#' Example: c.stack <- compress_stack(rbind(c(NA, NA, 1, 3, 5, 3, NA, 1, NA), c(NA, 1, 3, 5, NA, 3, NA, 1, 3), c(3, 5, NA, NA, 1, NA, NA, 1, 3)))
#'
#' @param stack NA-laden matrix with real values you want to retain
#' @return compressed matrix (list format)
#'          - pos : linear index within stack
#'          - vals: non-NA vals within stack
#'          - m : number of rows in stack
#'          - n : number of columns in stack
#' 
#'          
#' @importFrom magrittr %>%
#'
#' @export
compress_stack <- function(stack){
  
  # stack <- ref.mat
  keep <- which(!is.na(stack))
  # return cstack
  return(list(pos = keep,
              vals = stack[keep],
              m = nrow(stack),
              n = ncol(stack))
         )
}

##########################################################################################################################
#' Select rows from compressed matrix (c.stack)
#' selection contains the necessary information to reconstruct selected rows,
#' but nothing else. Returns c.stack object. NOTE: positions are still with 
#' regard to the FULL stack; however, decompressing with cstack_expandRows will
#' only return the selection (row order preserved).
#'
#' 
#' @param c.stack compressed matrix object
#' @return compressed matrix with only the selected row data
#'          - pos : linear index within FULL stack
#'          - vals: non-NA vals from FULL stack 
#'          - m : number of rows in stack
#'          - n : number of columns in stack
#' Example selection <-  cstack_selectRows(c.stack, row.nums = c(1,2))
#'
#' @export
cstack_selectRows <- function(c.stack, row.nums){
  
  stack.coords <- ind2subR(c.stack$pos, c.stack$m)
  keep <- which(stack.coords$rows %in% row.nums)
  c.stack$pos <- c.stack$pos[keep]
  c.stack$vals <- c.stack$vals[keep]
  
  return(c.stack)
}


##########################################################################################################################

#' Decompress c.stack object. Works on full stack and selection. Will return a 
#' matrix with the number of rows in the c.stack obj, in their original order.
#'
#' @param c.stack compressed matrix object
#' @return compressed matrix with only the selected row data
#'          - pos : linear index within stack
#'          - vals: non-NA vals within stack
#'          - m : number of rows in stack
#'          - n : number of columns in stack
#' @examples rows.expanded <- cstack_expandRows(c.stack)
#'           rows.expanded <- cstack_expandRows(selection)
#'
#' @export

cstack_expandRows <- function(cstack){
  # Convert the stack positions for the selection back to stack coords
    stack.pos <- ind2subR(cstack$pos, cstack$m)
    stack.pos$rows <- stack.pos$rows %>% dplyr::dense_rank()
    
  # Build and fill a matrix with the selected rows expanded
  # * cstack$vals are what was compressed
  # * stack.pos is cstack$pos -> coordinates -> renamed rows to 1:length(unique)
    rows.expanded <- matrix(NA, max(stack.pos$rows), cstack$n)
    rows.expanded[sub2indR(stack.pos$rows,
                           stack.pos$cols,
                           nrow(rows.expanded))] <- cstack$vals

  return(rows.expanded)
}


