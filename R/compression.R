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
#' Note: for spectral matrices, it is rare that all columns have values. As such,
#' spectral matrices MUST not be transposed; i.e. each column should be
#' a ppm value and each row is a spectrum. Otherwise ppms will be lost upon 
#' reconstruction!
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
  
  # Check if there's anything in there. 
    if (is.null(stack)){stop('compress_stack() does not work on NULL')}
  
  # Handle single row cases
    if (is.vector(stack)){stack <- matrix(stack, nrow = 1)}
  
  # If using some other value, convert to NA
    if (!is.na(sparse.val)){
      stack[stack == sparse.val] <- NA
    }
  
  is.val <- !is.na(stack)
  if (!any(is.val)){stop('undefined when all values in stack == sparse.val')}
  
  # stack <- ref.mat
  keep <- which(is.val)
  
  # compress positions using range expansion
    
    ranges <- compress_pos(keep)
    if (!all.equal(expand_pos(ranges), keep)){stop('Error in compress_stack(): compress_pos/expand_pos failed')}
    
  # return cstack
  return(list(pos = ranges,
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
  pos <- expand_pos(c.stack$pos)
  stack.coords <- ind2subR(pos, c.stack$m)
  keep <- which(stack.coords$rows %in% row.nums) 
  c.stack$pos <- pos[keep] %>% compress_pos # subset the vector and range-compress
  c.stack$vals <- c.stack$vals[keep]
  # reset m for subset matrixrows? NO: pos is in terms of m, and must remain that way
  # lest we need to recompute them all. 
  
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
    pos <- expand_pos(cstack$pos)
    stack.pos <- ind2subR(pos, cstack$m)
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
    
  # Feature positions are a special case;
  # they only require the first element of each row, since the corresponding
  # profiles are NA-gapped and retain the gap information:
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
#' @param stack.list list of matrices with matching positions (same sizes)
#' @param sparse.val passthrough to compress_stack()
#' @param key which matrix (name or list index) to use as the compression model (default = 1). If name, will be converted to list index.
#' @param apply.to which matrices (other than key) to apply compression scheme to
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
co_compress <- function(stack.list, sparse.val = NA, key = 1, apply.to = NULL, stack.names = NULL){
  # stack.list <- list(ref.data = ref.data, ref.ppm = ref.ppm)
  # key <- 2
  
  if (is.character(key)){
    # Assume field name provided
      key <- which(names(stack.list) %in% key)
      
  }
  
  if (is.null(apply.to)){apply.to <- 1:length(stack.list) %>% .[-key]}
  
  key.stack <- stack.list[[key]]
    
  # How to store the stack
    c.stack <- compress_stack(key.stack, sparse.val = sparse.val)
    stack.list[[key]] <- c.stack$vals
    
  # Loop through each stack in the list and apply the same positions
    
    pos <- expand_pos(c.stack$pos) # decompress position filter
    # loop through stacks and apply pos filter:
    # (names actually get preserved in the assignment)
    stack.list[apply.to] <- lapply(stack.list[apply.to], function(x) x[pos])
               
  # Replace vals with list of vals
    c.stack$vals <- stack.list

  return(c.stack)
}

# compress_pos #####################################################################################################
#' Since positions (e.g. "keep" in compress_stack()) are just runs of increasing inds,
#' there isn't really a need to keep all of them. Rather, just store the start and end
#' points of the runs. 
#' 
#' Example: all.equal(expand_pos(compress_pos(keep)), keep)
#'
#' @param keep vector of non-sparse position inds
#' @return ranges of the incremental runs in keep (can be directly expanded to whatever was in keep)
#' 
#'          
#' @importFrom magrittr %>%
#'
#' @export
compress_pos <- function(keep){
    # compress positions ("keep" in compress_stack()) using range expansion
    # Position will increase by 1 in a range, but > 1 between ranges.
    # Get the range starts and ends (within pos). A single point will
    # be both start and end. 
    
      jumps <- which(diff(keep) > 1)
      rng.end <- c(jumps, length(keep)) # the last index is always an end
      rng.start <- c(1, jumps+1) # the first index is always a start
      ranges <- data.frame(start = rng.start %>% keep[.],
                           end = rng.end %>% keep[.])
      return(ranges)
}

# expand_pos #####################################################################################################
#' Expand ranges to a vector
#'
#' @param keep vector of non-sparse position inds
#' @return ranges of the incremental runs in keep
#' 
#' Example: all.equal(expand_pos(compress_pos(keep)), keep)
#'          
#' @importFrom magrittr %>%
#'
#' @export
expand_pos <- function(ranges){
    
  # expand each range, then concatenate
  pos <- apply(ranges, 1, fillbetween) %>% unlist %>% c
  
  return(pos)
}

# expand_stacklist #####################################################################################################
#' Get each stack(s) from co-compressed stacklist
#'
#' @param c.stack compressed stack where vals is a list of co-compressed stacks
#' @param which.stacks inds or names of stacks to return from vals 
#' @return ranges of the incremental runs in keep
#' 
#' Example: all.equal(expand_pos(compress_pos(keep)), keep)
#'          
#' @importFrom magrittr %>%
#'
#' @export
expand_stacklist <- function(c.stack, which.stacks = NULL){
  # c.stack <- x$mapped$data.compressed

  if (!is.list(c.stack$vals)){warning('c.stack$vals is not a list'); return(NULL)}
  if (is.null(which.stacks)){which.stacks <- 1:length(c.stack$vals)}
  
  # Convert the stack positions for the selection back to stack coords
    pos <- expand_pos(c.stack$pos)
    stack.pos <- ind2subR(pos, c.stack$m)
    stack.pos$rows <- stack.pos$rows %>% dplyr::dense_rank()
    
  # Build and fill a matrix with the selected rows expanded
  # * cstack$vals are what was compressed
  # * stack.pos is cstack$pos -> coordinates -> renamed rows to 1:length(unique)
    rows.expanded <- matrix(NA, max(stack.pos$rows), c.stack$n)
    inds <- sub2indR(stack.pos$rows,
                             stack.pos$cols,
                             nrow(rows.expanded))
    
    # Loop through and put expanded stacks back in place (to preserve names)
    c.stack$vals[which.stacks] <- lapply(c.stack$vals[which.stacks], function(stack.vals){
      rows.expanded[inds] <- stack.vals
      return(rows.expanded)
    })
    
  # Return the requested matrices 
    return(c.stack$vals[which.stacks])
    
}
