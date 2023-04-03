#' Cross-correlation with local alignment
#'
#' This function calculates the cross-correlation between a reference vector and a matrix of spectra, where each
#' spectrum is a row of the matrix. The function aligns the reference vector with each spectrum by shifting
#' the spectrum and then calculating the cross-correlation between the aligned spectrum and the reference vector.
#' The function returns the best shift for each spectrum, along with the corresponding correlation score, p-value,
#' and overlap size. The function can also generate a plot showing the aligned spectra and reference vector.
#'
#' @param x A matrix where each row is a spectrum.
#' @param ppm A vector of column names for the matrix x.
#' @param signature.vals A vector of values matching the signature index window.
#' @param signature.idx.wind The reference indices (across all rows of currentInds).
#' @param currentInds The window(s) in which the reference sits in x (often larger than the reference).
#' @param min.overlap The minimum number of points on which to base the correlation.
#' @param slide If "internal", the function just slides the reference vector until it hits the bounds of currentInds;
#' if "external", the function pads currentInds to ensure that all reference vector points cross the bounds of currentInds;
#' if "limited", the user provides the number of lags (duplicated for both directions).
#' @param lag.limit The limit on the number of lags to consider (for "limited" slide option).
#' @param plots A logical value indicating whether to generate a plot showing the aligned spectra and reference vector.
#'
#' @return A list with the following elements:
#' \item{best.shifts}{A vector of the best shift for each spectrum.}
#' \item{best.rvals}{A vector of the correlation scores for each spectrum.}
#' \item{best.pvals}{A vector of the p-values for each spectrum.}
#' \item{shiftedInds}{A matrix of the indices in x for the aligned spectra.}
#' \item{shiftedSpecs}{A matrix of the aligned spectra.}
#' \item{ref.aligned}{A vector of the reference vector values (aligned with the spectra).}
#' \item{overlap.sizes}{A vector of the overlap sizes for each spectrum.}
#' \item{g}{A plot showing the aligned spectra and reference vector (if \code{plots = TRUE}).}
#'
#' @examples
#' x <- matrix(rnorm(100), ncol = 10)
#' ppm <- seq(100)
#' signature.vals <- rnorm(10)
#' signature.idx.wind <- seq(5, 14)
#' currentInds <- matrix(seq(1, 100, by = 10), ncol = 2)
#' xcorr_localign(x, ppm, signature.vals, signature.idx.wind, currentInds, min.overlap = 3, slide = "internal", lag.limit = NA, plots = FALSE)
#'
#' @importFrom pracma Reshape
#' @importFrom tidyverse *
#' @importFrom stats xcorr
#' @importFrom matrixStats colVars
#' @importFrom cowplot plot_grid
#' @export
xcorr_localign <- function(x, ppm, # full matrix and colnames
                           signature.vals, # vals matching signature.idx.wind
                           signature.idx.wind, # ref indices (across all rows of currentInds)
                           currentInds, # the window(s) in which the ref sits in x (often larger than ref)
                           min.overlap = 3, # minimum amount of points on which to base the correlation
                           slide = "internal", # if internal, just slide ref till it hits currentInds bounds
                           # if external, pad currentInds to ensure all ref points cross bounds
                           # if "limited", provide number of lags (duplicated for both directions)
                           lag.limit = NA,
                           plots = F) {

  # currentInds needs to be the nrow matrix x ncol (specData). It's the inds of x.
  #
  # x-corr based on current ref within region -/+ span
  # Extract out the matrix for this feature (mapped for each row)

  specData <- shiftedInds_toVals(x, currentInds)
  spec.idx <- 1:ncol(specData)
  ref.idx <- signature.idx.wind

  # Apply to each spectrum and collect scores grid
  # These lags indicate the shift (# of elements starting from the leftmost
  # element) that the column slices of specData where ncol(slice) = length(ref),
  # will start from. These slices will be stacked and compared to the ref values.

  if (is.na(lag.limit)) {
    lag.limit <- length(lags)
  }

  # Internal (if already padded, = complete)
  if (slide == "internal") {
    start.left <- min(spec.idx) - min(ref.idx)
    start.right <- max(spec.idx) - max(ref.idx)
    lags <- start.left:start.right
  }
  if (slide == "limited") {
    start.left <- max(c(min(ref.idx) - lag.limit, 1))
    start.right <- min(c(min(ref.idx) + lag.limit, ncol(specData) - span(ref.idx)))
    lags <- start.left:start.right
  }
  if (slide == "external") {
    # Will not work with plotting at the moment
    # Complete (external, pad data)
    lags <- (min(spec.idx) - max(ref.idx)):(max(spec.idx) - min(ref.idx))
    pad.by <- ncol(specData) - 1
    ref.idx <- ref.idx + pad.by
    specData <- padmat(specData, use = NA, col.by = pad.by)
    currentInds <- apply(currentInds, 1, function(x) {
      x %>%
        range() %>%
        "+"(., c(-pad.by, pad.by)) %>%
        fillbetween()
    }) %>% t()
  }

  # Calculate xcorr for each pair:

  results <- lapply(
    1:nrow(specData),
    function(specrow) {
      # specrow <- 10
      # print(specrow)
      xcorr(specData[specrow, ], # %>% simplePlot
        signature.vals,
        ref.idx,
        lags,
        min.overlap = min.overlap
      )
    }
  )

  lags <- lapply(1:length(results), function(x) results[[x]]$lags) %>% do.call(rbind, .)


  # Get all the rvals and pvals so we don't hardcode a filter. Still provide the filter, but don't apply

  rvals <- lapply(1:length(results), function(x) results[[x]]$rvals) %>% do.call(rbind, .)
  pvals <- lapply(1:length(results), function(x) results[[x]]$pvals) %>% do.call(rbind, .)
  overlaps <- lapply(1:length(results), function(x) results[[x]]$overlaps) %>% do.call(rbind, .)


  # # Find the best shift for each spectrum
  #   # Used to use a > 1SD of the rvals filter, which worked okay, but did't help
  #   with variable overlapping obeservations (e.g. signatures with NAs in them, or
  #   calculations on padding). rvals aren't reliable for that, but pvals are.

  # Just use minimum p value for positive non-NA rvals on each row:

  pvals[rvals < 0] <- NA # only positive rvals
  bestFit.idx <- rep(min(ref.idx), nrow(specData)) # default is no shift


  # Find min pval for positive rvals in each row
  # Don't modify the default value for the rows which have only NAs (no positive rvals)

  hasMin <- apply(pvals, 1, function(x) !all(is.na(x))) %>% which()
  if (length(hasMin) == 0) {
    return(errorOut())
  }
  bestFit.idx[hasMin] <- apply(pvals[hasMin, , drop = FALSE], 1, which.min)

  # Using the column inds bestFit.idx and the row numbers, get linear inds for
  # the pvals and rvals matrices.
  best <- sub2indR(rows = 1:nrow(rvals), cols = bestFit.idx, nrow(rvals))


  # Extract the rvals and pvals and lags for each spectrum
  bestFit.score <- rvals[best]
  bestFit.pval <- pvals[best]
  bestFit.overlap <- overlaps[best]
  best.shifts <- lags[best] # Convert idx in shift matrix back to lag
  best.shifts[best.shifts %>% is.na()] <- 0

  # Build the aligned matrix, Line up the spectra (calculate inds):

  # The x matrix inds for each row are given by the currentInds matrix, which
  # spans specData (window). ref.idx indexes within that window, i.e. can be
  # used to pull the columns of currentInds, in turn giving the actual x matrix
  # locations of ref.idx in x.

  # Which part of specData aligns with ref.idx?
  # * There will be a different lag for each row.

  shiftedInds <- matrix(NA, nrow = nrow(currentInds), ncol = ref.idx %>% span())

  # shiftedInds <- outer(ref.idx, best.shifts, "+") %>% apply(., 1, (range %>% fillbetween))
  # shiftedInds[] <- currentInds[shiftedInds]

  for (i in 1:nrow(currentInds)) {
    shiftedInds[i, ] <- currentInds[i, (ref.idx + best.shifts[i]) %>% range() %>% fillbetween(),
      drop = FALSE
    ] # shouldn't be NAs
  }

  # Linear inds to access vals in x and put in shiftedSpecs

  lininds.shifted <- sub2indR(
    rows = rep(1:nrow(specData), ncol(shiftedInds)),
    cols = shiftedInds,
    m = nrow(specData)
  ) %>% pracma::Reshape(., nrow(specData), ncol(shiftedInds))

  # Create shiftedSpecs (matrix in which specData region matching ref has been moved
  # to the ref location)
  # Pluck the data out of x using the shifted inds

  shiftedSpecs <- matrix(NA, nrow = nrow(specData), ncol = ncol(specData))

  shiftedSpecs[, ref.idx %>% range() %>% fillbetween()] <- x[lininds.shifted]

  # Plot the ref signature followed by the shiftedSpecs
  #   What I want to see is the rows of shiftedSpecs, shifted.

  # Create ref vector (get spacing with NAs for disjoint refs)
  refshape <- matrix(NA, nrow = 1, ncol = ncol(specData))
  refshape[ref.idx] <- signature.vals

  # Make the plot
  g <- NULL
  if (plots) {
    # Calculate plot bounds (ref bounds)
    bounds <- 1:ncol(shiftedSpecs)

    # Combine ref with shapes from specData (aligned) to compare
    compSpec <- rbind(refshape[, bounds], shiftedSpecs[, bounds])

    # Build a "ref" for use in plot_addRef
    # Get indices of ref vals in the smaller matrix (bounds)
    idx <- compSpec[1, ] %>%
      is.na() %>%
      "!"(.) %>%
      which()

    # Build the list
    ref <- list(
      signature = list(
        idx.wind = idx,
        vals = compSpec[1, idx]
      ),
      wind = list(inds = compSpec[1, ] %>% seq_along())
    )

    # Make the plot
    g <-
      simplePlot(compSpec[-1, , drop = FALSE], xvect = ppm[bounds]) %>%
      plot_addRef(compSpec[-1, , drop = FALSE], ref, ppm[bounds], scaling = TRUE)
    plot(g)
    # stackplot(compSpec,vshift = .05,hshift = 0)
  }

  # Return results

  return(list(
    best.shifts = best.shifts,
    best.rvals = bestFit.score,
    best.pvals = bestFit.pval,
    shiftedInds = shiftedInds,
    shiftedSpecs = shiftedSpecs,
    ref.aligned = refshape,
    overlap.sizes = bestFit.overlap,
    g = g
  ))
}

#' Create an error output object with default values.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{best.shifts}{A numeric vector of length 1, set to 0.}
#'   \item{best.rvals}{A numeric vector of length 1, set to 0.}
#'   \item{best.pvals}{A numeric vector of length 1, set to 1.}
#'   \item{shiftedInds}{A NULL object.}
#'   \item{shiftedSpecs}{A NULL object.}
#'   \item{ref.aligned}{A NULL object.}
#'   \item{overlap.sizes}{A numeric vector of length 1, set to 0.}
#'   \item{g}{A NULL object.}
#' }
#'
#' @examples
#' errorOut()
#'
#' @export
errorOut <- function() {
  return(list(
    best.shifts = 0,
    best.rvals = 0,
    best.pvals = 1,
    shiftedInds = NULL,
    shiftedSpecs = NULL,
    ref.aligned = NULL,
    overlap.sizes = 0,
    g = NULL
  ))
}
