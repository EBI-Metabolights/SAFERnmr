#' ppTargetSpecR
#'
#' Function to peak-pick spectra decomposed from mixture spectra.
#' Requires function detectSpecPeaks from speaq package returns a
#' list of chemical shifts and intensities which can be optionally
#' normalised to the most intense peak for improved comparison
#' against HMDB spectra.
#' Goncalo Graca, 8 January 2021, g.gomes-da-graca-at-imperial.ac.uk
#' @param target target spectra to peak pick.
#' @param normalise flag indicating normalisation of output, default FALSE.
#' @return list of chemical shifts and intensities.
#' @export
ppTargetSpec <- function(target, normalise = FALSE) {
    tmp <- detect_spec_peaks(
        target, nDivRange = 128,
        scales = seq(1, 16, 2),
        baselineThresh = 1000,
        SNR.Th = -1,
        verbose = FALSE
        )
    s <- target[, tmp[[2]]]
    s <- t(s)
    colnames(s) <- c("chemical-shift", "intensity")
    s <- s[order(s[, 1]), ]
    if (normalise) {
        s[, 2] <- s[, 2] / max(s[, 2])
    }
    return(s)
}
