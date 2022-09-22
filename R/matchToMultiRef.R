#' matchToMultiRef
#'
#' Matches peaks from one pseudo-spectrum corresponding to one driver peak
#' to all database peaks.
#' Returns a data frame with the results ranked by score.
#' @param references reference spectra as a db file.
#' @param metadata metadata that is included in the final merged dataframe.
#' @param target target driver peak.
#' @param driver_ppm parts per million of driver peak.
#' @param match_method matching methodology as a string.
#' @param tol tolerance in ppm.
#' @param Itol tbc.
#' @param intensity flag indicating whether to figure intensity into matching.
#' @return Dataframe of results ranked by score.
#' @export
matchToMultiRef <- function(references, metadata, target, driver_ppm,
                            matchMethod = "", tol = 0.02, Itol = 20,
                            intensity = FALSE) {
  library(flexclust)
  # Run either the old or new matching function (v2 has matchMethod options, 'basic' method = old)
  if (!matchMethod == "") {

    # note from callum - does assigning intensity to FALSE as is done below not mean it will always be FALSE,
    # even if the function is called with this flag set to TRUE?
    res <- lapply(
      1:length(references),
      function(x) {
        matchPeaksToRef2(
          target = target,
          driver_ppm = driver_ppm,
          reference = references[[x]],
          tol = tol,
          Itol = Itol,
          matchMethod = matchMethod,
          intensity = FALSE
        )
      }
    )
  } else {
    res <- lapply(
      1:length(references),
      function(x) {
        matchPeaksToRef(target,
          driver_ppm,
          references[[x]],
          tol = tol,
          Itol = Itol,
          intensity = FALSE
        )
      }
    )
  }


  # The ppms of matched peaks are removed from further progress for now
  # as it causes problems when converting the result to data frame format
  # MTJ: we'll add back all this data and more later, but for now just work with
  # the subset so we can make minimal modifications to GG's code:

  includeFields <- c(
    "hmdb_id", "spectrum_id", "matches", "reference_peaks", "score")
  clean_res <- lapply(1:length(res), function(x) res[[x]][includeFields])


  # Convert the list to data frame and rename the column names

  clean_res <- data.frame(
    matrix(unlist(clean_res), nrow = length(clean_res), byrow = TRUE))

  colnames(clean_res) <- includeFields

  # Combine based on common colnames ("hmdb_id","spectrum_id")

  final_res <- merge.data.frame(
    clean_res,
    metadata
  )


  # Rank the results by score
  sorted_res <- final_res[order(final_res$score, decreasing = TRUE), ]
  ranked_res <- sorted_res
  ranked_res$score <- as.numeric(ranked_res$score)

  # Trim all non-matches
  ranked_res <- ranked_res[which(!ranked_res$score == 0), ]

  # "dense ranking"
  ranked_res$rank <- as.integer(factor(-ranked_res$score))

  # Handle the (hopefully rare) case of no matches
  if (nrow(ranked_res) > 0) {

    # Add driver peak column
    ranked_res[, "driver_peak_ppm"] <- driver_ppm

    # Using the filtered results, pull out the list elements for
    # potential matches (contains more data)

    # Match up the HMDB id and spec IDs
    filtOn <- c("hmdb_id", "spectrum_id")
    filtIDs <- ranked_res$hmdb_id
    filtSpec <- ranked_res$spectrum_id

    res <- res[as.numeric(rownames(ranked_res))]
  } else {
    ranked_res <- NULL
    res <- NULL
  }

  # Report the results

  # Stick the list and df results into a list so we can pass them both out
  out <- list(scores_matrix = ranked_res, filteredResults = res)

  # Return them

  return(out)
}
