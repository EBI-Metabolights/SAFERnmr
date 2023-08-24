# compress_stack #########################################################################################################################
#' Sparse (NA-gapped) matrix compression
#'
#' Compresses matrix to its non-NA values and records their positions
#' to enable perfect reconstruction on-the-fly. Intended for sparse matrices
#' and parallel operations as well as general memory management.
#' 
#' In the case of feature objects, the profile and position mats can be 
#' compressed together. 
#' 
#' Example: c.stack <- compress_stack(rbind(c(NA, NA, 1, 3, 5, 3, NA, 1, NA), c(NA, 1, 3, 5, NA, 3, NA, 1, 3), c(3, 5, NA, NA, 1, NA, NA, 1, 3)))
#'
#' @param stack sparse matrix with real values you want to retain
#' @param sparse.val (optional) value taking up too much space (i.e. NA default, or 0, or other)
#' 
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
compress_stack <- function(stack, sparse.val = NA){
  
  # Handle single row cases
    if (is.vector(stack)){stack <- matrix(stack, nrow = 1)}
  
  # If using some other value, convert to NA
    if (!is.na(sparse.val)){
      stack[stack == sparse.val] <- NA
    }
  
  # stack <- ref.mat
  keep <- which(!is.na(stack))
  
  # return cstack
  return(list(pos = keep,
              vals = stack[keep],
              m = nrow(stack),
              n = ncol(stack),
              sparse.val = sparse.val)
         )
}

# cstack_selectRows #########################################################################################################################
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
  # c.stack$m <- length(unique(stack.coords$rows[keep])) # reset m for subset matrixrows? NO: pos is in terms of m
  
  return(c.stack)
}


# cstack_expandRows #########################################################################################################################

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
# compress_features #####################################################################################################
#' Sparse (NA-gapped) matrix compression applied to features object.
#'
#' In the case of feature objects, the profile and position mats can be 
#' compressed together. 
#' 
#' Example: c.stack <- compress_stack(rbind(c(NA, NA, 1, 3, 5, 3, NA, 1, NA), c(NA, 1, 3, 5, NA, 3, NA, 1, 3), c(3, 5, NA, NA, 1, NA, NA, 1, 3)))
#'
#' @param feature feature object with 
#'          - stack (profiles on rows)
#'          - position (cols in xmat on rows), same size as stack
#'          - driver.relative and sfe optional
#'          
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
compress_features <- function(feature){
  # How to store the stack
    c.stack <- compress_stack(feature$stack)
    
  # Position stack only requires the first element if profiles are NA-gapped
    min.col.pos <- Rfast::rowMins(feature$position)
    min.col.x <- feature$position[sub2indR(rows = 1:nrow(feature$position),
                                           cols = min.col.pos, 
                                           m = nrow(feature$position))]
    
    first.col.in.x <- min.col.x - (min.col.pos - 1)
    
    # The first value in each row is stored. To expand for a given row, simply 
    # construct an NA vector of size n (number of columns for feat obj matrices), 
    # then fill elements 
    
  # sfe and driver.relative
    # not sure if this can be compressed
    # if (is.null(feature$driver.relative)){feature$driver.relative <- NA}
    # if (is.null(feature$sfe)){feature$sfe <- NA}
    
  compressed.feature <- list(stack = c.stack,
                             position = data.frame(xcol = first.col.in.x,
                                                   row = 1:nrow(feature$stack))
                             # driver.relative = feature$driver.relative,
                             # sfe = feature$sfe
                             )
  return(compressed.feature)
}

# select_features #####################################################################################################
#' Compressed features object subsetting
#'
#' Pull selected features from compressed obj.
#' 
#' just.rows can be used to extract a subset of features or all, just in a different order.
#' c.stack <- compress_stack()
#' Example: 
#' # Make example feature obj
#'  feature <- list(profile = rbind(c(NA, NA, 1, 3, 5, 3, NA, 1, NA), c(NA, 1, 3, 5, NA, 3, NA, 1, 3), c(3, 5, NA, NA, 1, NA, NA, 1, 3)))
#'  feature$position <- outer(1:3, 1:9, "+")
#'  feature$position[is.na(feature$profile)] <- NA
#' # Compress
#'  compressed.feature <- compress_features(feature)
#' # Decompress all
#'  feature.exp <- expand_features(compressed.feature)
#' # Decompress rows 1 & 2
#'  feature.exp <- expand_features(compressed.feature, c(1,2))
#'  
#' @param compressed.feature compressed feature obj
#' @param just.rows row numbers of feature selection (or reordered row numbers)
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
select_features <- function(compressed.feature, row.nums=NULL){
  
  # # Get easy access to m and n
  #   m <- compressed.feature$stack$m
  #   n <- compressed.feature$stack$n
   
  # Do different things if selecting rows
  
    # Select rows from stack 
      cstack.selection <- cstack_selectRows(c.stack = compressed.feature$stack, 
                                            row.nums = row.nums)
      
      compressed.feature$stack <- cstack.selection
      
    # Select rows from pos 
    
      # Which position elements correspond to the row numbers being pulled?
        rel.rows <- compressed.feature$position$row %in% row.nums
        compressed.feature$position <- compressed.feature$position[rel.rows,] # incorrect; row.nums does not index pos.
        
      # compressed.feature$driver.relative$ss.all <- compressed.feature$driver.relative$ss.all[row.nums, ,drop=F]
      # compressed.feature$driver.relative$sizes <- compressed.feature$driver.relative$sizes[row.nums]
      # compressed.feature$sfe <- compressed.feature$sfe[row.nums]
      
  return(compressed.feature)
}


# expand_features #####################################################################################################
#' Compressed features object decompression.
#'
#' In the case of feature objects, the profile and position mats can be 
#' compressed together to get further compression gains.
#' 
#' just.rows can be used to extract a subset of features or all, just in a different order.
#' c.stack <- compress_stack()
#' Example: 
#' # Make example feature obj
#'  feature <- list(profile = rbind(c(NA, NA, 1, 3, 5, 3, NA, 1, NA), c(NA, 1, 3, 5, NA, 3, NA, 1, 3), c(3, 5, NA, NA, 1, NA, NA, 1, 3)))
#'  feature$position <- outer(1:3, 1:9, "+")
#'  feature$position[is.na(feature$profile)] <- NA
#' # Compress
#'  compressed.feature <- compress_features(feature)
#' # Decompress all
#'  feature.exp <- expand_features(compressed.feature)
#' # Decompress rows 1 & 2
#'  feature.exp <- expand_features(compressed.feature, c(1,2))
#'  
#' @param compressed.feature compressed feature obj
#' @param just.rows row numbers of feature selection (or reordered row numbers)
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
expand_features <- function(compressed.feature, row.nums=NULL){
  # Get easy access to m and n
    m <- compressed.feature$stack$m
    n <- compressed.feature$stack$n
    
  # Do different things if selecting rows
  if (!is.null(row.nums)){
    # Select rows, then expand stack
      compressed.feature <- select_features(compressed.feature, row.nums = row.nums)
  }
    
  # Full stack expansion of whatever is in compressed.feature at this point
    # Profiles
      stack <- cstack_expandRows(compressed.feature$stack)

    # Positions 
      c.pos <- compressed.feature$position$xcol

    # Re-create Positions matrix
      # The first value in each row is stored. To expand for given row(s), simply 
      # construct a linear count from c.pos, then replace NAs according to stack.
        
        cols.to.add <- n - 1
        
        pos.stack <- lapply(1:nrow(stack), function(row.m){
          return(c.pos[row.m]:(c.pos[row.m] + cols.to.add))
        }) %>% do.call(rbind,.)
      
        pos.stack[is.na(stack)] <- NA  
  

  # Build output 
    feature <- list(profile = stack,
                    position = pos.stack
                    # driver.relative = compressed.feature$driver.relative,
                    # sfe = compressed.feature$sfe
                    )
  return(feature)
}

# co_compress #####################################################################################################
#' Sparse (NA-gapped) matrix compression from one matrix applied to others
#'
#' Generalized case of stack.lists matrix 
#' 
#' Example: c.stack <- compress_stack(rbind(c(NA, NA, 1, 3, 5, 3, NA, 1, NA), c(NA, 1, 3, 5, NA, 3, NA, 1, 3), c(3, 5, NA, NA, 1, NA, NA, 1, 3)))
#'
#' @param stack.list stack.list object with 
#'          - stack (profiles on rows)
#'          - position (cols in xmat on rows), same size as stack
#'          - driver.relative and sfe optional
#'          
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
co_compress <- function(stack.list, sparse.val = 0, key = 1, apply.to = NULL){
  # stack.list <- spec %>% lapply(as.numeric)
  # key <- 2
  
  if (is.character(key)){
    # Assume field name provided
      key <- which(names(stack.list) %in% key)
      
  }
  
  if (is.null(apply.to)){apply.to <- 1:length(stack.list[-key])}
  
  key.stack <- stack.list[[key]]
    
  # How to store the stack
    c.stack <- compress_stack(key.stack, sparse.val = sparse.val)
    stack.list[[key]] <- c.stack$vals
    
  # Loop through each stack in the list and apply the same positions
    
    stack.list[apply.to] <- lapply(stack.list[apply.to], function(x) x[c.stack$pos])
               
  # Replace vals with list of vals
    c.stack$vals <- stack.list
    
                             
  return(c.stack)
}

