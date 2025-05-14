library(wavelets)

xmat=data$xmat
ppm=data$ppm
pexp <- expand_protofeature(p, xmat, ppm, half.window)
driver <- pexp$driver

wind <- pexp$specRegion.inds

specRegion = pexp$specRegion
ppmRegion = pexp$ppmRegion

y <- specRegion[1,]
simplePlot(y, ppmRegion)
dwt_result <- dwt(y, filter = "haar", n.levels = 5)

# Plot approximation and detail coefficients
par(mfrow = c(5, 1), mar = c(2, 4, 2, 1))

# Approximation coefficients (lowest frequency content)
plot(dwt_result@V[[4]], type = "l", main = "Approximation Coefficients (Level 4)", ylab = "V4")

# Detail coefficients (high to low frequency content)
plot(dwt_result@W[[4]], type = "l", main = "Detail Coefficients (Level 4)", ylab = "W4")
plot(dwt_result@W[[3]], type = "l", main = "Detail Coefficients (Level 3)", ylab = "W3")
plot(dwt_result@W[[2]], type = "l", main = "Detail Coefficients (Level 2)", ylab = "W2")
plot(dwt_result@W[[1]], type = "l", main = "Detail Coefficients (Level 1)", ylab = "W1")



# Step 1: Decompose the signal
i <- 0
i <- i + 1

level <- i
dwt_result <- dwt(y, filter = "la8", n.levels = level)

# Step 2: Zero out low-frequency components (broad background)
# You can experiment with this! For example:
dwt_result@V[[level]] <- dwt_result@V[[level]] * 0  # zeros, same class and dimensions  # Remove approximation (baseline)
dwt_result@W[[level]] <- dwt_result@W[[level]] * 0  # Remove broadest detail (optional)
# Leave W1–W4 intact for fine structure

# Step 3: Reconstruct the signal
y_cleaned <- idwt(dwt_result)

# Step 4: Plot to compare
plot(y, type = "l", col = "gray", main = "Wavelet Background Removal")
lines(y_cleaned, col = "blue")
legend("topright", legend = c("Original", "Background Removed"), col = c("gray", "blue"), lty = 1)
