#' Export Matches to filesystem.
#'
#' Export matches to file for visualization and storage
#' Ideally, this function will allow:
#' cleanMatches(matches,ranklimit = 1,annot_field = "hmdb_id")
#' Export matches
#' matches object is a list with two objects for each STOCSY cluster matched
#' (including those with null matches):
#'   - scores_matrix
#'     - data frame containing 1 row for each non-zero scored match, and columns
#'       with match info/metadata such as hmdb_id, ref spectrum id, (Jaccard) scores,
#'      number of matches, number of possible matches in ref, compound name, instrument
#'       properties, experiment properties, etc.
#'   - filteredResults
#'     - list of lists, each list contains additional info for each row in
#'       scores_matrix which could not fit into a data frame. Some info is
#'       duplicated in the scores_matrix, but the important data here is:
#'       - peaklists: reference and target peaklists (matched lists)
#'       - matchlist: detailed match info for the matched peaks only (both ref and target):
#'         - "ppms" of the peaks
#'         - relative "intensities" of the reference peaks (if known)
#'         - "inds" in the starting ref or target list, respectively
#' Extract data and put in tall matrix (all matched peaks).
#' Define Format: In order to plot the matches as two lists of scatter points and
#' lines connecting the matched target and ref points, we need:
#'   - all ref peak positions and intensities
#'   - all target peaks and intensities
#'   - the connection from each matched ref to its target
#'     - peaks with no match could be coded as NaN
#'   - indicate the limiting peak (intensity-wise) for a match
#' @param matches matches object, see above description for further information.
#' @param X tbc.
#' @param rankLimit rank limit, defaults to 1.
#' @return N/A, matches are written to file.
#' @export
exportMatches <- function(matches,
                          X = NULL,
                          rankLimit = 1,
                          output_dir) {


  # Make flatfiles with all data:

  # Initialize dataframes
  matchPairs <- NULL
  referenceList <- NULL
  targetList <- NULL
  # length(matches)

  # Figure out which matches are empty and exclude them
  hasData <- unlist(lapply(
    1:length(matches),
    function(x) any(matches[[x]]$scores_matrix$rank <= rankLimit)
  ))

  for (i in which(hasData)) {

    # Each iteration is for another STOCSY cluster
    clustNumber <- i

    # Pull out the rows/lists that pass the rank threshold:
    included <- matches[[i]]$scores_matrix$rank <= rankLimit
    tmp <- matches[[i]]$filteredResults[included]

    for (j in 1:length(tmp)) {

      # Each iteration is for another matched compound for this cluster
      thispair <- tmp[[j]]

      # Record the new match pair:
      mp <- data.frame(
        cluster = clustNumber,
        spec_id = thispair$spectrum_id,
        ref_id = thispair$hmdb_id,
        ref_name = thispair$hmdb_name,
        jaccard_score = thispair$score
      )

      # Cat to the big df
      matchPairs <- rbind(matchPairs, mp)

      # Record the reference peaks:
      # Initialize tmp df with ID
      rp <- data.frame(matchPair_ID = nrow(matchPairs))
      # NOTE: this only works because we already added the new match pair;
      #       computing would probably be better.

      # Add the peak info (fills in IDs to match)
      rp <- cbind(rp, thispair$peaklists$reference)

      # Add the indices:
      rp$target <- NA
      rp$target[thispair$matchlist$ref$inds] <- thispair$matchlist$target$inds
      # Note: if there are repeated indices, results will be correct, but not predictable

      # Record the row numbers
      rp$peakNumber <- 1:nrow(rp)

      # Add the rows to referenceList
      referenceList <- rbind(referenceList, rp)

      # Record the target peaks:
      # Initialize tmp df with ID
      tp <- data.frame(matchPair_ID = nrow(matchPairs))


      # Add the peak info (fills in IDs to match)
      tp <- cbind(tp, thispair$peaklists$target)

      # Add the indices:
      tp$ref <- NA
      tp$ref[thispair$matchlist$target$inds] <- thispair$matchlist$ref$inds
      # Note: if there are repeated indices in ref, results will be correct but duplicated

      # Record the row numbers
      tp$peakNumber <- 1:nrow(tp)

      # Add the rows to referenceList
      targetList <- rbind(targetList, tp)
    }
  }

  # Reorder cols as we see fit:
  matchPairs$ID <- 1:nrow(matchPairs)
  matchPairs <- matchPairs[, c("ID", "cluster", "spec_id", "ref_id", "ref_name", "jaccard_score")]
  referenceList <- referenceList[, c("matchPair_ID", "peakNumber", "chemical-shift", "intensity", "target")]
  colnames(referenceList) <- c("matchPair_ID", "peakNumber", "ppm", "intensity", "target")
  # ^ Matlab will pitch a fit about hyphens in colnames

  targetList <- targetList[, c("matchPair_ID", "peakNumber", "chemical-shift", "intensity", "ref")]
  colnames(targetList) <- c("matchPair_ID", "peakNumber", "ppm", "intensity", "ref")
  # ^ Matlab will pitch a fit about hyphens in colnames

  # Write the matrices to tsv files (avoid issues with commas in cmpd names)

  write.table(matchPairs,
    file = stringr::str_c(output_dir, "matchPairs.tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )

  write.table(referenceList,
    file = stringr::str_c(output_dir, "referenceList.tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )

  write.table(targetList,
    file = stringr::str_c(output_dir, "targetList.tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )

}
