###############################################################################
# Paralellized (not vectorized) collapse-scale analyzer for multiple features
# 
# The point of this analysis is to determine the point at which smoothing produces
# <= 1 true peak in the feature profile. Singlet/shoulder peaks will quickly collapse,
# even if they have small bumps on them. Here we use Gaussian smoothing with 
# varying sigma levels, which double at each level. We simply need to catch the point 
# at which real features begin to emerge. 
# 
# -------------------------------------------------------------
# Gaussian smoothing
# -------------------------------------------------------------
smooth_gauss <- function(x, sigma){
  if (sigma <= 0) return(x)
  k <- ceiling(4*sigma)             # 
  gx <- -k:k                        # 99.99% of the distribution is between +/- 4*sigma
  kern <- exp(-(gx^2)/(2*sigma^2))  # gaussian formula relating intensities to x coords via sigma
  kern <- kern / sum(kern)          # normalize intensities
  y <- stats::filter(x, kern, sides = 2) # apply the filter to profile

  # fill NA edges by nearest valid value
  if (anyNA(y)){
    good <- which(!is.na(y))
    y[1:(min(good)-1)] <- y[min(good)]
    y[(max(good)+1):length(y)] <- y[max(good)]
  }
  as.numeric(y)
}

# -------------------------------------------------------------
# Linear interpolation over NA runs (just for filtering purposes)
# -------------------------------------------------------------
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

  # linear interpolation inside
  x[bad] <- approx(good, x[good], xout = bad)$y
  x
}

# -------------------------------------------------------------
# Normalize to [0,1]
# -------------------------------------------------------------
normalize01 <- function(x){
  x <- x - min(x, na.rm=TRUE)
  r <- max(x, na.rm=TRUE)
  if (r == 0) return(rep(0, length(x)))
  x / r
}

# -------------------------------------------------------------
# LOCAL MAXIMUM DETECTION (very simple, global rule)
# -------------------------------------------------------------
local_maxima <- function(y){
  # y[i] is a peak if it is greater than neighbors
  # boundaries excluded
  n <- length(y)
  if (n < 3) return(integer(0))

  # exclude NAs
  good <- which(!is.na(y))
  if (length(good) < n) {
    # small fix: smooth edges already filled but double check
    y <- fill_na_interp(y)
  }

  idx <- which(
    c(FALSE,
      y[2:(n-1)] > y[1:(n-2)] & y[2:(n-1)] > y[3:n],
      FALSE)
  )
  idx
}
###############################################################################
analyze_features_collapse <- function(
    feat_mat,                   # matrix: features on rows
    scales = c(0, 1, 2, 4, 8, 16, 36),
    fill_na = TRUE,
    normalize = TRUE,
    n_cores = 1,                # can parallelize if desired
    verbose = TRUE
){
  
  if (!is.matrix(feat_mat))
    feat_mat <- as.matrix(feat_mat)

  nfeat <- nrow(feat_mat)

  if (verbose)
    message("Analyzing ", nfeat, " features across ", length(scales), " scales…")

  # Use mclapply if >1 core
  FUN <- function(i){
    analyze_feature_collapse(
      feat = feat_mat[i, ],
      scales = scales,
      fill_na = fill_na,
      normalize = normalize,
      plot = FALSE
    )
  }

  if (n_cores > 1){
    require(parallel)
    out <- parallel::mclapply(1:nfeat, FUN, mc.cores=n_cores)
  } else {
    out <- lapply(1:nfeat, FUN)
    # i <- 5
  }

  # Extract collapse scales + peak counts into matrices
  collapse <- sapply(out, `[[`, "collapse_scale")
  peakcounts <- t(sapply(out, `[[`, "peak_counts"))

  colnames(peakcounts) <- paste0("sigma_", scales)

  list(
    collapse_scales = collapse,
    peak_counts = peakcounts,
    results = out
  )
}

###############################################################################
# MULTI-PANEL GRID PLOT OF SCALE-SPACE STACKPLOTS
###############################################################################
plot_feature_grid <- function(
    feat_mat,
    feat_ids = 1:16,
    scales = c(0,1,2,4,8,16,32),
    nrow = NULL, ncol = NULL,
    ...
){
  nf <- length(feat_ids)

  # automatic grid
  if (is.null(nrow) && is.null(ncol)){
    nrow <- ceiling(sqrt(nf))
    ncol <- ceiling(nf / nrow)
  } else if (is.null(nrow)){
    nrow <- ceiling(nf / ncol)
  } else if (is.null(ncol)){
    ncol <- ceiling(nf / nrow)
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow = c(nrow, ncol), mar=c(2.5,2.5,2.5,1))

  for (i in feat_ids){
    feat <- feat_mat[i, ]
    analyze_feature_collapse(
      feat,
      scales = scales,
      plot = TRUE,
      ...
    )
    title(paste("Feature", i), line = 1)
  }
}

# Suppose feature$stack is your matrix (features on rows)
feature.stack <- lapply(sats, feat_profile) %>% do.call(rbind,.)
feat_mat <- feature.stack

# Analyze everything (no plotting)
vec <- analyze_features_collapse(
  feat_mat,
  scales = c(0,1,2,4,8,16,32),
  n_cores = 4   # optional
)

# Peak counts
head(vec$peak_counts)

# Collapse scales
summary(vec$collapse_scales)

# Grid plot
plot_feature_grid(
  feat_mat,
  feat_ids = 1:25,
  scales = c(0,1,2,4,8,16,32), main.title = ""
)

###############################################################################
# PLOT N EXAMPLES FOR A GIVEN COLLAPSE SCALE AND SAVE TO PDF
###############################################################################
plot_examples_for_scale <- function(
    feat_mat,
    collapse_scales,
    target_scale,
    scales = c(0,1,2,4,8,16,32),
    n_examples = 25,
    ncol = 5,
    peak_cex = 0.5,
    pdf_file = NULL,
    seed = 1,
    sigma.labels = TRUE
){
  set.seed(seed)

  # Which features match that collapse scale?
  idx <- which(collapse_scales == target_scale)

  if (length(idx) == 0){
    warning("No features found with collapse scale = ", target_scale)
    return(invisible(NULL))
  }

  # Sample up to n_examples
  chosen <- sample(idx, min(n_examples, length(idx)))

  # Default PDF name if not supplied
  if (is.null(pdf_file)){
    pdf_file <- paste0("collapse_scale_", target_scale, "_examples.pdf")
  }

  pdf(pdf_file, width = 8, height = 10)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  # Grid layout: 2 x ceiling(n/2)
  nrow <- ceiling(length(chosen) / ncol)

  par(mfrow = c(nrow, ncol), mar = c(2,2,2,1))

  for (i in chosen){
    feat <- feat_mat[i, ]

    analyze_feature_collapse(
      feat,
      scales = scales,
      plot = TRUE,
      peak_cex = peak_cex,
      sigma.labels = sigma.labels,
      main.title = ""
    )

    title(main = paste("Feature", i),
          cex.main = 0.9)
  }

  dev.off()

  message("Saved ", length(chosen), 
          " examples for collapse scale ", target_scale,
          " → ", pdf_file)

  invisible(chosen)
}

plot_features_byScale <- function(
    feat_mat,
    collapse_scales,
    target_scale,
    n_examples = 25,
    ncol = 5,
    peak_cex = 0.5,
    pdf_file = NULL,
    seed = 1)
  {
  set.seed(seed)

  # Which features match that collapse scale?
  idx <- which(collapse_scales == target_scale)

  if (length(idx) == 0){
    warning("No features found with collapse scale = ", target_scale)
    return(invisible(NULL))
  }

  # Sample up to n_examples
  chosen <- sample(idx, min(n_examples, length(idx)))

  # Default PDF name if not supplied
  if (is.null(pdf_file)){
    pdf_file <- paste0(n_examples, "_features_collapse_scale_", target_scale, ".pdf")
  }

  pdf(pdf_file, width = 8, height = 10)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  # Grid layout: 2 x ceiling(n/2)
  nrow <- ceiling(length(chosen) / ncol)

  par(mfrow = c(nrow, ncol), mar = c(2,2,2,1))

  for (i in chosen){
    feat <- feat_mat[i, ]

    plot(feat, type='l', lwd=1)

    title(main = paste("Feature", i),
          cex.main = 0.5)
  }

  dev.off()

  message("Saved ", length(chosen), 
          " examples for collapse scale ", target_scale,
          " → ", pdf_file)

  invisible(chosen)
}

###############################################################################

# Usage: 
# Suppose feature$stack is your matrix (features on rows)
feature.stack <- lapply(sats, feat_profile) %>% do.call(rbind,.)
feat_mat <- feature.stack
scales <- c(0,1,2,4,8,12,16,32)

# Analyze everything (no plotting)
vec <- analyze_features_collapse(
  feat_mat,
  scales = scales,
  n_cores = 4   # optional
)

# Peak counts
head(vec$peak_counts)

# Collapse scales
summary(vec$collapse_scales)

# Plot degradation of each feature profile
res <-lapply(scales, function(x){
  plot_examples_for_scale(feat_mat,
      vec$collapse_scales,
      target_scale = x,
      scales = scales,
      n_examples = 25,
      ncol = 5,
      peak_cex = 0.5,
      pdf_file = NULL,
      seed = 1, 
      sigma.labels = FALSE)
})

# Plot with just feature profile
res <-lapply(scales, function(x){
  res <-plot_features_byScale(feat_mat,
      vec$collapse_scales,
      target_scale = x,
      n_examples = 25,
      ncol = 5,
      peak_cex = 0.5,
      pdf_file = NULL,
      seed = 1)
})

