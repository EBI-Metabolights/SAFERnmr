xmat=data$xmat
ppm=data$ppm
pexp <- expand_protofeature(p, xmat, ppm, half.window)
driver <- pexp$driver

baseline <- predict(loess(y ~ seq_along(y), span = 0.8))
y_corrected <- y - baseline

plot(y, type = "l", col = "gray", main = "LOESS Background Removal")
lines(baseline, col = "blue")
legend("topright", legend = c("Original", "Background Removed"), col = c("gray", "blue"), lty = 1)


## Only run on the mins

