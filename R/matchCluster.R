#' @name matchClusterR
#'
#' @description MTJ 28FEB2022
#' Requires flexclust, pracma libraries
#' Calculate all pairwise distances between all target cluster peaks and reference cluster peaks and filter for tolerance. Match each peak with its closest ref peak:
#' A)
#' start with lowest overall distance, assign target[i] to ref[r] peak, record ref[r] as unavailable (remove its matches in the distance matrix)
#' find next lowest distance match, and repeat until no more matches
#' – easier to implement right now, necessary update if staying with GG code
#' B)
#' find the best set of matches which satisfy the following:
#' ppm ranks are preserved (no criss-crossing of match lines)
#' ref peaks can be assigned to at most one target peak
#' best match set determined by sum of distances
#' –> alternative, potentially better qualities and potentially faster (no loop/complex decision tree)
#' install.packages(here(RcppHungarian), repos=NULL, type="source")
#' @param t target peaks.
#' @param r reference peaks.
#' @param tol tolerance in ppm.
#' @param method matching methodology as a string, defaults to 'itmin'.
#' Can also be substructure based 'substruct'. optimal 'opt' not implemented.
#' @return List of matches.
#' @export
matchCluster <- function(t, r, tol = 0.02, method = "itmin") {

  # Initialize reporting vectors

  matched_targets <- NA
  length(matched_targets) <- length(t)
  target_inds <- matched_targets

  matched_refs <- NA
  length(matched_refs) <- length(r)
  ref_inds <- matched_refs

  distances <- matched_targets

  # Match using chemical shift only

  ## New code:

  # Calculate distance matrix; threshold on tol
  d <- flexclust::dist2(t, r, method = "euclidean") # applies abs() - what if it didn't?
  mask <- d <= tol

  # Initialize common format for export:
  matched_targets <- NA
  target_inds <- NA
  matched_refs <- NA
  ref_inds <- NA
  distances <- NA


  if (method == "itmin") {
    # Iterated minimum method
    # Make the dist matrix only with reasonable pairs:
    d[!mask] <- NA

    tlen <- nrow(d)
    rlen <- ncol(d)
    possibleMatches <- sum(rowSums(!is.na(d)) >= 1)

    # While matches are still available (max # matches (non-NA) = # of possible targets)
    i <- 1L
    while (any(!is.na(d)) & (length(na.omit(matched_targets)) <= possibleMatches)) {

      # Which is the lowest available distance?

      subs <- ind2subR(which.min(d), tlen)

      # Check if there was even a min
      if (length(subs[[1]]) == 1) {
        # Make the ref peak assignment

        matched_targets[i] <- t[subs$rows]
        target_inds[i] <- subs$rows
        matched_refs[i] <- r[subs$cols]
        ref_inds[i] <- subs$cols
        distances[i] <- d[subs[[1]], subs[[2]]]

        # Mark the ref peak unavailable in dist matrix

        d[, subs$cols] <- NA
        d[subs$rows, ] <- NA
        i <- i + 1
      } else {
        break
      } # If there is nothing left in d
    }
  } else if (method == "substruct") {

    # Pull out substructures from ref


    # Match substructs instead of ref peaks using convolution

    #
  } else if (method == "hca") {

    # Pull out substructures from ref

    # Match substructs instead of ref peaks using convolution

    #
  } else if (method == "dtw") {

    # Pull out substructures from ref


    # Match substructs instead of ref peaks using convolution

    #
  } else if (method == "hungarian" || method == "hungarian_scaled") {
    # Note: mask computed earlier
    # hungarian algorithm RcppHungarian
    # Note: option to sqrt scale the out-of-tolerance resonances to discourage them

    # Remove all peaks which have no matches within tolerance:

    dlabels <- d
    dlabels[1:numel(d)] <- 1:numel(d)

    # * New indices for reduced cost matrix * (tab indicates indexing level)
    redRows <- which(rowSums(mask) != 0)
    redCols <- which(rowSums(t(mask)) != 0)
    redDlabs <- dlabels[redRows, redCols]
    redD <- matrix(d[redRows, redCols], nrow = length(redRows))

    # Check to make sure there's even anything in the reduced matrix:

    if (sum(size(redD)) >= 2) {
      # Exponential scaling of out-of-tolerance resonances (optional)
      if (method == "hungarian_scaled") {
        oot <- redD > tol & redD < 1 & redD > 0
        redD[oot] <- sqrt(redD[oot])
      }

      # Calculate optimal matchings (d cannot have NAs). New indices (matchpa)

      matchPairs <- HungarianSolver(redD)

      # Which have matches - beware of R trying to turn 1xn into vector (set nrow to force it to not transpose)

      keepRows <- which(rowSums(matchPairs$pairs != 0) == 2)
      matchedTol <- matrix(matchPairs$pairs[keepRows, ], nrow = length(keepRows))

      # Convert match inds to rd inds
      match_rdinds <- sub2indR(matchedTol[, 1], matchedTol[, 2], nrow(redD))

      # < convert back to d inds
      match_dInds <- redDlabs[match_rdinds]
      match_dInds <- match_dInds[mask[match_dInds]]
      finalPairs <- ind2subR(match_dInds, nrow(d))

      # Put in common format for export:
      matched_targets <- t[finalPairs$rows]
      target_inds <- finalPairs$rows
      matched_refs <- r[finalPairs$cols]
      ref_inds <- finalPairs$cols
      distances <- matched_targets - matched_refs
    }
    # If not, defaults take over
  } else {
    return(NULL)
  } # bail out

  # # Fair algorithm?
  #   # Rank targets based on number of viable refs

  # Remove NAs from the list and return
  # Add the indices
  return(list(
    matchedTargets = na.omit(matched_targets),
    matchedTargetInds = na.omit(target_inds),
    matchedRefs = na.omit(matched_refs),
    matchedRefInds = na.omit(ref_inds),
    distances = na.omit(distances),
    method = method
  ))
}
