simplePlot_x <- function(ymat = NULL, xvect = NULL, n_xticks = NULL, xdir = "reverse",
                       linecolor = "gray", opacity = 0.6, linewidth = 0.5){

  # Handle single vector case:
  if (is.vector(ymat))
    ymat <- ymat %>% c %>% as.matrix %>% t

  # Default x-axis
  if (is.null(xvect)) xvect <- 1:ncol(ymat)

  # Melt for ggplot
  df <- as.data.frame(t(ymat))
  colnames(df) <- 1:ncol(df)
  df$ppm <- xvect

  d <- reshape2::melt(df, id.vars = "ppm")
  colnames(d) <- c("ppm", "specNumber", "Spectral_Intensity")

  nlines <- nrow(ymat)

  # --------------------------
  # Handle vectorized colors
  # --------------------------
  if (length(linecolor) > 1 && length(linecolor) == nlines) {
    color_vals <- alpha(linecolor, opacity)
  } else {
    color_vals <- rep(alpha(linecolor, opacity), nlines)
  }

  # assign a color per spectrum
  d$color <- color_vals[as.numeric(d$specNumber)]

  # --------------------------
  # Handle vectorized linewidths
  # --------------------------
  if (length(linewidth) > 1 && length(linewidth) == nlines) {
    lw_vals <- linewidth
  } else {
    lw_vals <- rep(linewidth, nlines)
  }

  d$lw <- lw_vals[as.numeric(d$specNumber)]

  # --------------------------
  # Build plot
  # --------------------------
  g <- ggplot2::ggplot(d) +
    ggplot2::geom_line(
      aes(x = ppm,
          y = Spectral_Intensity,
          group = specNumber,
          color = specNumber,
          linewidth = lw),
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = color_vals) +
    ggplot2::scale_linewidth(range = range(lw_vals))

  # --------------------------
  # X-axis direction and ticks
  # --------------------------
  if (is.null(n_xticks)) {
    if (xdir == "reverse") {
      g <- g + ggplot2::scale_x_reverse(breaks = scales::breaks_pretty())
    } else {
      g <- g + ggplot2::scale_x_continuous(breaks = scales::breaks_pretty())
    }
  } else {
    if (xdir == "reverse") {
      g <- g + ggplot2::scale_x_reverse(breaks = scales::breaks_extended(n = n_xticks))
    } else {
      g <- g + ggplot2::scale_x_continuous(breaks = scales::breaks_extended(n = n_xticks))
    }
  }

  # --------------------------
  # Theme
  # --------------------------
  g <- g +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = element_text(colour = "black", size = 12),
                   legend.position = "none",
                   axis.text.y = element_blank(),
                   axis.title = element_blank(),
                   axis.ticks = element_blank(),
                   panel.border = element_blank(),
                   panel.grid = element_blank())

  return(g)
}
