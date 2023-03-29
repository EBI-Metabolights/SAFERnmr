#' Filter matches based on ppm shift difference
#'
#' This function filters the matches in match.info based on the difference in 
#' ppm shift between the feature and the match. Only matches with a ppm shift 
#' difference less than or equal to ppm.tol are retained.
#'
#' @param match.info A data frame containing the match information
#' @param feature A feature object with position information
#' @param ppm A numeric vector of ppm values
#' @param fits.feature A logical vector indicating whether each match fits the feature
#' @param ppm.tol The maximum allowed ppm shift difference between the feature and the match
#'
#' @return A list with two elements: the filtered `match.info` data frame, and the filtered `fits.feature` logical vector
#'
#' @export
filter.matches_shiftDelta <- function(match.info,
                                      feature,
                                      ppm,
                                      fits.feature,
                                      ppm.tol = 0.5) {


  match.info[, "ppm.difference"] <- pblapply(1:nrow(match.info), function(m) {
    fnum <- match.info$feat[m]

    f.inds.trim <- match.info[m, c("feat.start", "feat.end")] %>% as.numeric()
    feat.range <- f.inds.trim %>%
      feature$position[fnum, .] %>%
      range(na.rm = T) %>%
      ppm[.]
    ref.range <- match.info[m, c("ref.start", "ref.end")] %>%
      as.numeric() %>%
      ppm[.]
    return(mean(feat.range) - mean(ref.range))
    # match.info$feat.end %>% sort %>% plot
  }) %>% unlist()

  # scattermore::scattermoreplot(sort(match.info[,'ppm.difference']),
  #                              1:nrow(match.info), xlab = "ppm difference", ylab = "match")
  #
  # scattermore::scattermoreplot(sort(abs(match.info[,'ppm.difference'])),
  #                              1:nrow(match.info), xlab = "ppm difference", ylab = "match")
  # later, potentially: as fraction of known shift range?


  keep <- abs(match.info$ppm.difference) <= ppm.tol
  match.info <- match.info[keep, ]
  fits.feature <- fits.feature[keep]

  #####


  return(list(
    match.info = match.info,
    fits.feature = fits.feature
  ))
}
