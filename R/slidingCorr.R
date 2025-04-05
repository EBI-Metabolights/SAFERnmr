#' Calculate sliding window correlation matrix, and label correlation peaks (pockets)
#'
#'
#' @param x Input matrix.
#' @param ws Size of sliding window.
#' @param extractPockets If \code{TRUE}, calculate correlation pockets for each column of \code{x} within each sliding window.
#' @param plotting If \code{TRUE}, produce a stackplot of the correlation matrix with pockets highlighted.
#' @param vshift Vertical shift in plot.
#' @param ppm X-axis values for plot - only used for plotting. If not specified, uses column numbers of \code{x}.
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
slidingCorr <- function(x,ws, extractPockets = FALSE, plotting = TRUE, vshift = 20, ppm = NULL, n.cores = 10){
  #ws <- 250
  
  ppm.inds <- 1:ncol(x)
  if (is.null(ppm)){ppm <- ppm.inds}
  
  wind <- -ws:ws
  driver_pos <- ws+1
  
  # Just use a vector as a template, will combine into a matrix later
  corr_template <- matrix(data = NA, nrow = 2*ws+1, ncol=1)

  # Map for keeping calculations in bounds of edges
  indsmat <- outer(wind, 1:ncol(x), "+")
    oob <- indsmat < 1 | indsmat > ncol(x)
    indsmat[oob] <- NA
    
  in.bounds <- (indsmat %>% is.na %>% "!"(.)) %>% pracma::Reshape(., nrow(indsmat), ncol(indsmat))
   
  # Set up multicore
    
    # split up the ppm vector into chunks
    # but first, randomize it so certain cores don't get stuck with no 
    ppm.rand <- sample(ppm.inds)
      unrand <- order(ppm.rand, decreasing = FALSE)
    chunk.size <- ceiling(length(ppm.rand) / n.cores)
    ppm.chunks <- split(ppm.rand, ceiling(seq_along(ppm.rand) / chunk.size))

    results <- mclapply(ppm.chunks, function(ppm.segment){
      # ppm.segment <- ppm.chunks[[1]]
      
      lapply(ppm.segment, function(j){
        # j <- ppm.segment[[1]]
        corrmat <- covmat <- corr_template
        use <- in.bounds[,j] %>% which
        corrmat[use] <- cor(x[, j], x[, indsmat[use,j]])
        covmat[use] <- cov(x[, j], x[, indsmat[use,j]])
        
        return(list(cors = corrmat,
                    covs = covmat,
                    j = j))
      })
      
    }, mc.cores = n.cores) %>% unlist(recursive = FALSE)
  
    cors <- lapply(results, function(result){
      list(cors = result$cors,
           j = result$j)
    })
    
    covs <- lapply(results, function(result){
      list(covs = result$covs,
           j = result$j)
    })
    
  # Now we have a list containing the local covariance and correlation calcultation for each spectral point.
    
  # Calculate the primary correlation peak for each
    if (extractPockets){
      
      is.pocket.template <- rep(FALSE, nrow(in.bounds))
      
      pockets.unsorted <- mclapply(cors, function(result){
        # result <- cors[[1]]
        j <- result$j
        
        use <- in.bounds[,j] %>% which
        bounds <- corr_expand(peak = (use %in% driver_pos) %>% which,
                              localMinima(result$cors[use]),
                              vRange = c(1,length(use))) %>% unlist %>% use[.]
        is.pocket <- is.pocket.template
        is.pocket[bounds[1]:bounds[2]] <- TRUE
        
        return(is.pocket)
      }, mc.cores = n.cores) %>% do.call(cbind,.)
    
      pockets <- pockets.unsorted[, unrand]

    } else {
      pockets <- NULL
    }
    
    # Unrandomize the results
    
    corrmat <- sapply(cors[unrand], `[[`, "cors")
    covmat  <- sapply(covs[unrand], `[[`, "covs")
    
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
      
  return(list(corr_compact = corrmat,
              indsmat = indsmat,
              cov_compact = covmat,
              isPocket = pockets,
              plot = g,
              window = wind,
              center = driver_pos))
}
