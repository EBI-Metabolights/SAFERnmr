#######################################################################
# 1. Gaussian smoothing
#######################################################################
smooth_gauss <- function(x, sigma){
  if (sigma <= 0) return(x)
  k <- ceiling(4*sigma)
  gx <- -k:k
  kern <- exp(-(gx^2)/(2*sigma^2))
  kern <- kern / sum(kern)

  y <- stats::filter(x, kern, sides = 2)

  # pad edges
  if (anyNA(y)){
    good <- which(!is.na(y))
    y[1:(min(good)-1)] <- y[min(good)]
    y[(max(good)+1):length(y)] <- y[max(good)]
  }

  as.numeric(y)
}

#######################################################################
# 2. Linear interpolation over NA runs
#######################################################################
fill_na_interp <- function(x){
  na <- is.na(x)
  if (!any(na)) return(x)

  good <- which(!na)
  bad  <- which(na)

  # ------------------------------------------------------------------
  # 1. Pad ends FIRST so that interpolation endpoints are well-defined
  # ------------------------------------------------------------------
  if (min(good) > 1){
    x[1:(min(good)-1)] <- x[min(good)]
  }
  if (max(good) < length(x)){
    x[(max(good)+1):length(x)] <- x[max(good)]
  }

  # After padding, recompute NA positions:
  na <- is.na(x)
  if (!any(na)) return(x)

  good <- which(!na)
  bad  <- which(na)

  # ------------------------------------------------------------------
  # 2. Now interpolate (and extrapolate if needed)
  # ------------------------------------------------------------------
  x[bad] <- approx(good, x[good], xout = bad, rule = 2)$y

  x
}


#######################################################################
# 3. Normalize to [0,1]
#######################################################################
normalize01 <- function(x){
  x <- x - min(x, na.rm=TRUE)
  r <- max(x, na.rm=TRUE)
  if (r == 0) return(rep(0, length(x)))
  x / r
}

#######################################################################
# 4. Simple local maximum counter
#######################################################################
count_local_maxima <- function(y){
  n <- length(y)
  if (n < 3) return(0)
  y <- fill_na_interp(y)

  sum(y[2:(n-1)] > y[1:(n-2)] & y[2:(n-1)] > y[3:n])
}

#######################################################################
# 5. SCALE-SPACE PEAK COUNTING + collapse-scale identification
#######################################################################
analyze_feature_collapse <- function(
    feat,
    scales = c(0, 1, 2, 4, 8, 16, 32),
    fill_na = TRUE,
    normalize = TRUE,
    plot = FALSE,
    peak_cex = .5,
    main.title = "Scale-Space Peak Count",
    sigma.labels = TRUE
){
  
  x <- feat
  if (fill_na) x <- fill_na_interp(x)
  if (normalize) x <- normalize01(x)

  # smooth at each scale
  smoothed <- lapply(scales, function(s) smooth_gauss(x, s))

  # count peaks
  peak_counts <- sapply(smoothed, count_local_maxima)

  # collapse scale = first σ with <= 1 peak
  collapse_scale <- NA
  idx <- which(peak_counts <= 1)
  if (length(idx) > 0) collapse_scale <- scales[idx[1]]

  ###################################################################
  # Optional plotting
  ###################################################################
  if (plot){
    offsets <- seq(0, length(scales)-1)*1.2

    plot(NULL,
         xlim = c(1, length(x)),
         ylim = c(-0.5, max(offsets)+1),
         xlab="Index",
         ylab="Scale-space",
         main=main.title)

    for (i in seq_along(scales)){
      y <- smoothed[[i]]
      off <- offsets[i]

      # curve
      lines(y + off, lwd = 1.5)

      # peaks
      idxs <- which(
        c(FALSE,
          y[2:(length(y)-1)] > y[1:(length(y)-2)] &
          y[2:(length(y)-1)] > y[3:length(y)],
          FALSE)
      )

      if (length(idxs) > 0){
        points(idxs, y[idxs] + off, pch=19, col="red", cex=peak_cex)
      }

      if (sigma.labels){
      text(5, off + 0.2,
           paste0("sigma=", scales[i], " (", peak_counts[i], " peaks)"),
           adj=0)
        
      }
    }

    if (!is.na(collapse_scale)){
      abline(h = offsets[which(scales == collapse_scale)],
             col="blue", lty=2)
      legend("topright",
             legend = paste("collapse sigma =", collapse_scale),
             col="blue", lty=2, bty="n")
    }
  }

  ###################################################################
  # Return results
  ###################################################################
  list(
    smoothed = smoothed,
    peak_counts = peak_counts,
    collapse_scale = collapse_scale
  )
}

#######################################################################
# DEMO
#######################################################################
# set.seed(1)
# n <- 600
# x <- seq(0, 10, length.out=n)
# feat <- 
#     0.4*dnorm(x, 4.5, 0.05) +
#     1.0*dnorm(x, 5.0, 0.08) +
#     1.0*dnorm(x, 5.2, 0.08) +
#     0.6*dnorm(x, 5.5, 0.06) +
#     rnorm(n, 0, 0.001)
# 
# feat[sample(1:n, 30)] <- NA
# 
# # Run
# res <- analyze_feature_collapse(feat, plot=TRUE, peak_cex=0.5)
# 
# print(res$peak_counts)
# print(res$collapse_scale)
# 
# feat[sample(1:n,30)] <- NA
# 
# 
# print(res$shape)
# print(res$collapse_scale)
# print(res$persistence)
