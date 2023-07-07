#' Since pw feature alignments can have the feauture of interest in the f1 or f2
#' column, provide a function which formats lags relative to the specified feature.
#' Note: Reversing lags also means negating them. 
#'
#' @param lag.table lag table from feat_align, etc. 
#' @param feat.num feature number to arrange the lag table relative to
#'
#' @return subsetted and rearranged lag table
#' @export
#'
#'
#'
#' @importFrom magrittr %>%
pw_lags_relative_to <- function(lag.table, feat.num){
  # lag.table <- ci$lag.table
  
  key.in.f1 <- lag.table$f1 == feat.num
  key.in.f2 <- lag.table$f2 == feat.num
  
  # First, if the key feature is in f1, it was aligned to the corresponding f2, 
  # and the lag must be negated. Easier to do before subsetting.
  
    lag.table$f1[key.in.f1] <- lag.table$f2[key.in.f1]
    lag.table$f2[key.in.f1] <- feat.num
    lag.table$lag.in.f2[key.in.f1] <- -lag.table$lag.in.f2[key.in.f1]

  # Now subset:
    use <- key.in.f1 | key.in.f2
    lag.table <- lag.table[use, ]
    key.in.f1 <- key.in.f1[use]
    key.in.f2 <- key.in.f2[use]

    return(lag.table)
}