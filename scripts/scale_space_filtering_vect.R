###############################################################################
# VECTORIZED collapse-scale analyzer for multiple features
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
    print(i)
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

# Grid plot of first 12 features
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

# Usage: 

plot_examples_for_scale(feat_mat,
    vec$collapse_scales,
    target_scale = 32,
    scales = c(0,1,2,4,8,16,32),
    n_examples = 25,
    ncol = 5,
    peak_cex = 0.5,
    pdf_file = NULL,
    seed = 1, 
    sigma.labels = FALSE)
