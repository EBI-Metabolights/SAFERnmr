

# -------------------------------------------------------------
# SCALE-SPACE STACKPLOT WITH PEAK TRACKING
# -------------------------------------------------------------
plot_scale_space_with_peaks <- function(feat,
                                        scales = c(0,1,2,4,8,12),
                                        fill_na = TRUE,
                                        normalize = TRUE,
                                        peak_cex = .5){
  x <- feat
  if (fill_na)     x <- fill_na_interp(x)
  if (normalize)   x <- normalize01(x)

  # smooth at each scale
  smoothed <- lapply(scales, function(s) smooth_gauss(x, s))

  # detect peaks
  peaks <- lapply(smoothed, local_maxima)

  # y-offsets for stacked visualization
  offsets <- seq(0, length(scales)-1) * 1.2

  # --- PLOT ---
  plot(NULL, xlim = c(1, length(x)),
       ylim = c(-0.5, max(offsets)+1),
       xlab = "Index",
       ylab = "Scale-space",
       main = "Scale-Space Evolution With Peak Tracking")

  # draw curves + peaks
  for (i in seq_along(scales)){
    y <- smoothed[[i]]
    off <- offsets[i]

    # curve
    lines(y + off, lwd = 1.5)

    # peaks
    pk <- peaks[[i]]
    if (length(pk) > 0){
      points(pk, y[pk] + off,
             pch = 19, col = "red",
             cex = peak_cex)
    }

    # annotate scale
    text(5, off + 0.2, paste0("σ = ", scales[i]), adj = 0)
  }

  invisible(list(smoothed = smoothed,
                 peaks = peaks))
}

# -------------------------------------------------------------
# DEMO EXAMPLE
# -------------------------------------------------------------
set.seed(1)
n <- 600
x <- seq(0, 10, length.out=n)

# synthetic multiplet: 3-peak cluster + noise
feat <- 
    0.4*dnorm(x, 4.5, 0.05) +
    1.0*dnorm(x, 5.0, 0.08) +
    1.0*dnorm(x, 5.2, 0.08) +
    0.6*dnorm(x, 5.5, 0.06) +
    rnorm(n, 0, 0.001)

# add NA gaps like real NMR features
feat[sample(1:n, 30)] <- NA

# call the scale-space plotter
plot_scale_space_with_peaks(feat)
