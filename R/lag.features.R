#' From a feature stack and lag table, align features to a given feature. Works
#' on pairs as well. Returns the matrix containing aligned profiles (e.g., for plotting).
#'
#' lag.table example:
#'   f1 f2 lag.in.f2 pos.big          val
#'   1  1  2         0    9172 3.705236e+18
#'   2  2  1         0    3057 3.705236e+18
#' (last 2 columns unused here, but are in the output from feat.align())
#' 
#' lag.table f1 and f2 values need to index rows of featureStack.
#' lag.table is the key to extract rows from featureStack and align them.
#'
#' @param featureStack matrix containing feature profiles
#' @param lag.table lag information table from feat.align or similar
#' @param to feature number to align to (must be in f1 or f2) 
#'
#' @return matrix with aligned feature profiles. 
#'
#'
#' @export
lag.features <- function(featureStack,  # each row is a feature
                         lag.table,     # has fields f1, f2, lag.in.f2
                         to = NULL)        # which feature to align to
{       
  # featureStack <- feature.stack
    # lag.table <- pair
    if (is.null(to)){to <- min(lag.table$f1)}
  
  # Use the lag.table to adjust the relative alignments of the relevant rows of featureStack
    
    # Define relevant rows 
      # Remove duplicate matches ####
      
        lag.table <- cbind(lag.table$f1,lag.table$f2) %>% t %>% sortPairs %>% t %>% duplicated %>% "!"(.) %>% lag.table[.,]
      
      # Pull out all matches involving the key feature ("to") ####
      
        hits <- pw.lags.relative.to(lag.table, feat.num = to)

      # Order the rows ####
        
        key.in.f1 <- hits$f1 == to
        key.in.f2 <- hits$f2 == to
        # Pull out the features that weren't the feature (can't unique - need to keep in order; 'to' feature is first)
        rows.used <- c(to,  cbind(hits$f1, hits$f2) %>% .[!cbind(key.in.f1, key.in.f2)])
        
        fs <- featureStack[rows.used, ]
          # simplePlot(trim.sides(fs))
  
        lags <- -c(0, hits$lag.in.f2) # 0 is for the key feature
        
        
      # Calculate matrix of column inds shifted by lags
        indsmat <- outer(  1:ncol(fs)  , lags, "-") %>% t
        indsmat <- indsmat - min(indsmat) + 1
        linds <- sub2indR(rows = 1:nrow(indsmat), cols = indsmat, m = nrow(indsmat))
        tmpmat <- matrix(NA, nrow(indsmat), max(indsmat))
        tmpmat[linds] <- fs
        tmpmat <- tmpmat[order(rows.used), ] # re-sort the rows so they match the relative order passed in
  
  return(tmpmat)
    # simplePlot(trim.sides(tmpmat))
}