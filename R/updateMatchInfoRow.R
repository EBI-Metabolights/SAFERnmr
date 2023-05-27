#' In match propagation, take a copy of a match.info row and update the relevant 
#' fields. Return the updated row. Accessory fcn, called by filter.matches.R.
#' 
#' Within a feature cluster, only one representative (or the mean profile) was matched
#' to all the refs. Using the matches to that, and the lags of each cluster member
#' to the key feature/mean profile, update the fit data for each cluster member. 
#' It's annoying to do these updates in the filter.matches function once we 
#' recalculate the fits, so I do it in this one and return the row to add to 
#' match.info. 
#'
#' @param mi.row a row of match.info data frame 
#' @param fit pre-computed fit for intracluster lag-adjusted feature (the other cluster member) to the ref feat
#'
#' @return row of match.info, updated for other cluster members. 
#'
#' @export
updateMatchInfoRow <- function(mi.row, fit)
{

                        feat <- fit$feat.fit
                        spec <- fit$spec.fit
                        
                      # Update rval
                        use <- !is.na(feat + spec)
                        mi.row$rval <- cor(feat[use], spec[use])
                      
                      # Update pval using t-distribution
                      
                        a <- -abs(mi.row$rval * sqrt( (fit$pts.matched-2) /(1-mi.row$rval^2)))
                        mi.row$pval <- 2*pt(a,(fit$pts.matched-2))

                      # Update pts.matched
                      
                        mi.row$pts.matched <- fit$pts.matched
                        mi.row$pts.feat <- sum(!is.na(feat))
                      
                      # Update fit coefficients
                      
                        mi.row$fit.intercept <- fit$fit[1]
                        mi.row$fit.scale <- fit$fit[2]
                      
                      # Update scores
                      
                        mi.row$wasserstein.score <- fit$wasserstein.score <- score.wasserstein(feat, spec)
                        mi.row$sum.residuals <- fit$sum.residuals
                        mi.row$rmse <- fit$rmse
                        mi.row$rmse.weighted <- NA # can't get, keep previous?

                      # Update peak-level info (don't update, just keep)
                      
                        # mi.row$numpeaks.feat <- 
                        # mi.row$numpeaks.ref <- 
                        # mi.row$numpeaks.feat.nnf <- NA # can't get, keep previous?
                        # mi.row$refpeaks.matched <- 
  return(mi.row)
}