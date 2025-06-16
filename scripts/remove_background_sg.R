library(wavelets)

xmat=data$xmat
ppm=data$ppm
pexp <- expand_protofeature(p, xmat, ppm, half.window)
driver <- pexp$driver

  plot_protofeature(p,
            half.window = half.window, ppm = data$ppm,
            xmat,
            bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
            showPeaks = TRUE)
  p$driver
  p.width <- 
  abs(p$primary.upper - p$primary.lower) + 
    abs(p$secondary.upper - p$secondary.lower)
  
wind <- pexp$specRegion.inds

specRegion = pexp$specRegion
ppmRegion = pexp$ppmRegion

y <- specRegion[1,]
simplePlot(y, ppmRegion)

library(signal)

# Second derivative with Savitzky-Golay

d2y <- sgolayfilt(y, p = 3, n = p.width + (1-p.width%%2), m = 2)
# p = poly order
# n = window size (must be odd)
# m = derivative order

lapply()
# Threshold based on curvature
curvature_threshold <- quantile(abs(d2y), 0.75)
peak_mask <- abs(d2y) > curvature_threshold

# Create a filtered version preserving only "high curvature" points
y_filtered <- y
y_filtered[!peak_mask] <- NA  # Or replace with baseline estimate

# Plot
plot(y, type = "l", col = "gray", main = "Curvature-Based Peak Filtering")
lines(y_filtered, col = "blue", lwd = 2)


## Function

plot_SG_ROI_stack <- function(x, y, p, vspace=.1, quant_steps = seq(0.99, 0.01, by = -0.01), poly_deg = 3) {
  
  # 1️⃣ Define ROI window from p
  
  window_rng <- p$driver + range(c(p$primary.lower, p$secondary.lower, p$primary.upper, p$secondary.upper))
  window_min <- window_rng[1]
  window_max <- window_rng[2]
  
  # 2️⃣ Calculate ROI width in x units
  ROI_width <- window_max - window_min

  # 3️⃣ Determine SG window size (odd!)
  p.width <- ROI_width
  n_ROI <- p.width + (1 - p.width %% 2)  # Force odd window
  
  message("Using SG window n = ", n_ROI, ", covering ROI width of ", round(ROI_width, 3), " ", attr(x, "units", exact = TRUE))
  
  # 4️⃣ Compute second derivative
  d2y <- sgolayfilt(y, p = poly_deg, n = n_ROI, m = 2)
  
  # 5️⃣ Initialize matrix for filtered spectra
  filtered_mat <- matrix(NA, nrow = length(quant_steps), ncol = length(y))
  
  # 6️⃣ Base plot setup (stackplot effect by vertical offset)
  vshift <- max(y) * vspace  # vertical shift between lines
  plot(x, y, type = "n", ylim = c(min(y), max(y) + vshift * length(quant_steps)),
       xlab = "ppm", ylab = "Intensity",
       main = "SG ROI-Optimized Stackplot with Excluded Regions in Blue")
  
  # Add ROI vlines
  abline(v = ppm[c(window_min, window_max)], col = "red", lty = 2, lwd = 2)
  
  # 7️⃣ Loop over quantiles and plot
  for (i in seq_along(quant_steps)) {
    q <- quant_steps[i]
    
    # Threshold at this quantile
    curvature_threshold <- quantile(abs(d2y), q)
    peak_mask <- abs(d2y) > curvature_threshold
    
    # Create filtered version (masked-out regions = NA here)
    y_filtered <- y
    y_filtered[!peak_mask] <- NA
    
    # Store in matrix
    filtered_mat[i, ] <- y_filtered
    
    # Plot excluded regions (in blue)
    y_excluded <- y
    y_excluded[peak_mask] <- NA
    lines(x, y_excluded + vshift * (i - 1), col = "red", lwd = 2)
    
    # Plot included regions (high-curvature parts in black)
    lines(x, y_filtered + vshift * (i - 1), col = "black", lwd = 2)
  }
  
  # Optional legend
  legend("topright", legend = paste0(quant_steps * 100, "%"), 
         title = "Curvature Quantile",
         col = "black", lty = 1, lwd = 2, bg = "white")
  
  # Return the matrix for further use
  invisible(filtered_mat)
}

# Example call using your variables:
x <- ppmRegion
y <- specRegion[1,]

filtered_mat <- plot_SG_ROI_stack(x, y, p)

###########
# Plotly

library(plotly)

plot_SG_ROI_stack_plotly <- function(x, y, ppm, p, vspace = 0.001, quant_steps = seq(0.9, 0.1, by = -0.05), poly_deg = 3) {
  
  window_rng <- p$driver + range(c(p$primary.lower, p$secondary.lower, p$primary.upper, p$secondary.upper))
  window_min_idx <- window_rng[1]
  window_max_idx <- window_rng[2]
  
  # SG window size (still in index units)
  ROI_width <- window_max_idx - window_min_idx
  p.width <- ROI_width
  n_ROI <- p.width + (1 - p.width %% 2)
  
  message("Using SG window n = ", n_ROI, ", covering ROI width of ", ROI_width, " points.")
  
  # Second derivative
  d2y <- sgolayfilt(y, p = poly_deg, n = n_ROI, m = 2)
  
  # Prepare Plotly object
  plt <- plot_ly()
  
  # Vertical shift
  vshift <- max(y) * vspace
  
  # Precompute y limits
  y_axis_max <- max(y) + vshift * length(quant_steps)
  y_axis_min <- min(y)
  
  # Helper: split mask into continuous runs
  split_runs <- function(mask) {
    rle_mask <- rle(mask)
    runs <- vector("list", length(rle_mask$lengths))
    
    idx <- 1
    for (i in seq_along(rle_mask$lengths)) {
      len <- rle_mask$lengths[i]
      runs[[i]] <- idx:(idx + len - 1)
      idx <- idx + len
    }
    
    return(runs)
  }
  
  # Loop over quantiles
  for (i in seq_along(quant_steps)) {
    # i <- 0
    # i <- i + 1
    q <- quant_steps[i]
    
    curvature_threshold <- quantile(abs(d2y), q)
    peak_mask <- abs(d2y) > curvature_threshold
    
    y_offset <- vshift * (i - 1)
    
    # Split into runs
    runs <- split_runs(peak_mask)
    
    # Plot "included" runs (TRUE mask)
    for (r in runs[which(rle(peak_mask)$values)]) {
      # ri <- 0
      # ri <- ri + 1
      # r <- runs[which(rle(peak_mask)$values)] %>% .[[ri]]
      r <- unlist(r)
      plt <- plt %>%
        add_trace(x = x[r], y = y[r] + y_offset, type = 'scatter', mode = 'lines',
                  line = list(color = "red", width = 2),
                  name = paste0("Included Q", q * 100, "%"),
                  showlegend = FALSE)
    }
    
    # Plot "excluded" runs (FALSE mask)
    for (r in runs[which(!rle(peak_mask)$values)]) {
      # ri <- 0
      # ri <- ri + 1
      # r <- runs[which(!rle(peak_mask)$values)] %>% .[[ri]]
      r <- unlist(r)
        r <- c(min(r)-1, r, max(r)+1)
      plt <- plt %>%
        add_trace(x = x[r], y = y[r] + y_offset, type = 'scatter', mode = 'lines',
                  line = list(color = "black", width = 2),
                  name = paste0("Excluded Q", q * 100, "%"),
                  showlegend = FALSE)
    }
  }
  
  # Add ROI vlines — in ppm units:
  window_min_ppm <- ppm[window_min_idx]
  window_max_ppm <- ppm[window_max_idx]
  
  plt <- plt %>%
    add_segments(x = window_min_ppm, xend = window_min_ppm,
                 y = y_axis_min, yend = y_axis_max,
                 line = list(color = "red", width = 2, dash = "dash"), name = "ROI lower") %>%
    add_segments(x = window_max_ppm, xend = window_max_ppm,
                 y = y_axis_min, yend = y_axis_max,
                 line = list(color = "red", width = 2, dash = "dash"), name = "ROI upper")
  
  # Final layout
  plt <- plt %>% layout(title = "SG ROI-Optimized Stackplot (Interactive, Split Runs, PPM Correct)",
                        xaxis = list(title = "ppm"),
                        yaxis = list(title = "Intensity (offset)", range = c(y_axis_min, y_axis_max)),
                        showlegend = FALSE)
  
  return(plt)
}
plt <- plot_SG_ROI_stack_plotly(x, y, ppm, p)
plt

