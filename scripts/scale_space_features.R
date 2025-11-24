# -------------------------------------------------------------
# Scale-space demo for a single feature vector: feat
# -------------------------------------------------------------

# Gaussian smoothing helper -----------------------------------
smooth_gauss <- function(x, sigma){
  # sigma = Gaussian SD in points
  if (sigma <= 0) return(x)

  k <- ceiling(4*sigma)
  gx <- -k:k
  kern <- exp(-(gx^2)/(2*sigma^2))
  kern <- kern / sum(kern)

  # convolution with padding at ends
  y <- stats::filter(x, kern, sides = 2)
  # Replace NA edges with nearest valid values
  if (anyNA(y)) {
    good <- which(!is.na(y))
    y[1:(min(good)-1)] <- y[min(good)]
    y[(max(good)+1):length(y)] <- y[max(good)]
  }
  as.numeric(y)
}

# Interpolate over NA gaps (recommended) -----------------------
fill_na_interp <- function(x){
  na <- is.na(x)
  if (!any(na)) return(x)

  good <- which(!na)
  bad  <- which(na)

  # pad edges
  if (min(good) > 1){
    x[1:(min(good)-1)] <- x[min(good)]
  }
  if (max(good) < length(x)){
    x[(max(good)+1):length(x)] <- x[max(good)]
  }

  # linear interpolation within
  x[bad] <- approx(good, x[good], xout = bad)$y
  x
}

# Normalization (optional) -------------------------------------
normalize01 <- function(x){
  x <- x - min(x, na.rm=TRUE)
  rng <- max(x, na.rm=TRUE)
  if (rng == 0) return(rep(0, length(x)))
  x / rng
}

# -------------------------------------------------------------
# Create stackplot of smoothing evolution
# -------------------------------------------------------------
plot_scale_space <- function(feat, scales = c(0, 1, 2, 4, 8, 12), 
                             normalize = TRUE, 
                             fill_na = TRUE){
  x <- feat
  if (fill_na) x <- fill_na_interp(x)
  if (normalize) x <- normalize01(x)

  # smooth at each scale
  smoothed <- lapply(scales, function(s) smooth_gauss(x, s))

  # y-offset for stacked plotting
  offsets <- seq(0, length(scales)-1) * 1.2  # vertical spacing

  # empty plot
  plot(NULL, xlim = c(1, length(x)), 
       ylim = c(-0.5, max(offsets)+1), 
       xlab = "Index", ylab = "Scale-space", 
       main = "Scale-Space Evolution of Feature")

  # draw each smoothed line
  for (i in seq_along(scales)){
    lines(
      smoothed[[i]] + offsets[i], 
      col = "black", lwd = 1.5
    )
    text(5, offsets[i] + 0.2, paste0("σ = ", scales[i]), adj = 0)
  }

  # return the smoothed series for further analysis
  invisible(smoothed)
}

# -------------------------------------------------------------
# Example usage
# -------------------------------------------------------------

# Load your feature (example)
# feat <- feature$stack[f, ]   # actual usage in your own code

# For demo: generate a fake multiplet
set.seed(1)
n <- 600
x <- seq(0, 10, length.out = n)
feat <- 
    0.4*dnorm(x, 4.5, 0.05) +
    1.0*dnorm(x, 5.0, 0.08) +
    1.0*dnorm(x, 5.2, 0.08) +
    0.6*dnorm(x, 5.5, 0.06) +
    rnorm(n, 0, 0.001)

# Introduce NA gaps to simulate real features
feat[sample(1:n, 25)] <- NA

# Plot scale-space evolution
plot_scale_space(feat)
