#' Calculate sliding window correlation matrix, and label correlation peaks (pockets)
#'
#'
#' @param x Input matrix.
#' @param ws Size of sliding window.
#' @param extractPockets If \code{TRUE}, calculate correlation pockets for each column of \code{x} within each sliding window.
#' @param plotting If \code{TRUE}, produce a stackplot of the correlation matrix with pockets highlighted.
#' @param vshift Vertical shift in plot.
#' @param ppm X-axis values for plot. If not specified, uses column numbers of \code{x}.
#'
#' @return A list containing the following items:
#' \item{corr_compact}{The compact form of the sliding correlation matrix.}
#' \item{indsmat}{The offset matrix used to calculate the sliding correlation matrix.}
#' \item{cov_compact}{The compact form of the sliding covariance matrix.}
#' \item{isPocket}{A logical matrix indicating whether each point in the sliding window is within a correlation pocket.}
#' \item{plot}{If \code{plotting} is \code{TRUE}, a stackplot of the correlation matrix. If extractPockets = T, show only the pockets.}
#' \item{window}{The relative inds in the sliding window.}
#' \item{center}{The center position of the sliding window.}
#'
#' @examples
#' # Generate a test matrix
#' x <- matrix(rnorm(2000), nrow = 100)
#'
#' # Calculate sliding correlation matrix with pockets and plot
#' slidingCorr(x, ws = 20, extractPockets = TRUE, plotting = TRUE)
#'
#' @importFrom pracma Reshape
#' @importFrom magrittr %>%
#'
#' @export
slidingCorr <- function(x,ws, extractPockets = FALSE, plotting = TRUE, vshift = 20, ppm = NULL){
  #ws <- 250
  
  if (is.null(ppm)){ppm <- 1:ncol(x)}
  
  wind <- -ws:ws
  os <- ws+1
  
  corrmat <- matrix(data = NA, nrow = 2*ws+1, ncol = ncol(x))
  covmat <- corrmat
  indsmat <- outer(wind, 1:ncol(x), "+")
    oob <- indsmat < 1 | indsmat > ncol(x)
    indsmat[oob] <- NA
    
  in.bounds <- (indsmat %>% is.na %>% "!"(.)) %>% pracma::Reshape(., nrow(indsmat), ncol(indsmat))
    
    
    for (j in 1:ncol(x)){
      use <- in.bounds[,j] %>% which
      corrmat[use, j] <- a<-cor(x[, j], x[, indsmat[use,j]])
      covmat[use, j] <- a<-cov(x[, j], x[, indsmat[use,j]])
    }
  
  # Calc corr pocket for each
    if (extractPockets){
      pockets <- rep(FALSE, length(in.bounds)) %>% pracma::Reshape(., nrow(in.bounds), ncol(in.bounds))

      
      for (j in 1:ncol(corrmat)){
        use <- in.bounds[,j] %>% which
        bounds <- corr_expand(peak = (use %in% os) %>% which,
                              localMinima(corrmat[use,j]),
                              vRange = c(1,length(use))) %>% unlist %>% use[.]
        pockets[bounds[1]:bounds[2],j] <- TRUE
        
        
        
        # # Align within the center peak region
        #   
        #     corrmat[use,j]
        #     pockets[,j]
        #     
        #     # Get the data from the external shift regions of 
        #       x[, indsmat[use,j]]
        #       

        
      }
      
      
      
    } else {
      pockets <- NULL
    }
  
  # Calculate alignment for each 
    # if (extractPockets){
    #   pockets <- rep(FALSE, length(in.bounds)) %>% pracma::Reshape(., nrow(in.bounds), ncol(in.bounds))
    #   
    #   for (j in 1:ncol(corrmat)){
    #     use <- in.bounds[,j] %>% which
    #     bounds <- corr_expand(peak = (use %in% os) %>% which,
    #                           localMinima(corrmat[use,j]),
    #                           vRange = c(1,length(use))) %>% unlist %>% use[.]
    #     pockets[bounds[2]:bounds[1],j] <- TRUE
    #   }
    # } else {
    #   pockets <- NULL
    # }
  
    # in.bounds 2 is the correlation pockets
      
            # a[a < 0] <- 0
      
          # corrs <- a %>% t %>% rowSums(., na.rm = TRUE) %>% t
          # lengths <- filt2 %>% t %>% rowSums %>% t
          
          # Debug using b, a small subset of a
#######################################################################################################    
    g <- NULL
    if (plotting){
      if (extractPockets){
        filt2 <- pockets
      }else {
        filt2 <- in.bounds # default is plot the whole window
      }
      
      a <- corrmat
      a[!filt2] <- NA
      
        b <- a
        f <- filt2
        bi <- indsmat

        cmat <- matrix(NA, nrow = ncol(b), ncol = ncol(b))
    
        for (i in 1:ncol(b)){
          # Map the indices from the compact corrmat, a, -> expanded offset matrix, cmat
            colinds_c <- bi[f[,i], i]
            rowinds_b <- f[,i] %>% which
          
          # Pull the corrs into their spots on the nxn matrix
            cmat[i, colinds_c] <- b[rowinds_b, i]
        }
        # Make the plot using stackplot
        
          g <- stackplot(cmat, vshift = vshift, hshift = 0, xvect = ppm)
    }
  
    # How to define the corr pocket bounds and scores using this info???
    # scores <- 
      # stackStocsys?
      
  return(list(corr_compact = corrmat,
              indsmat = indsmat,
              cov_compact = covmat,
              isPocket = pockets,
              plot = g,
              window = wind,
              center = os))
}
