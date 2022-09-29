#' Stocsy Decomposition
#'
#' Function to run 'STOCSY decomposition' Goncalo Graca, 15 September
#'  2020, g.gomes-da-graca/at/imperial.ac.uk
#' @param spectra individual spectra file to run through decomp
#' @param idx index number of spextra files to be run through decomp
#' @param cthr cutoff threshold
#' @return Result of STOCS decomposition.
#' @export
stocsydec <- function(spectra, idx, cthr = 0.8) {
    ppm <- spectra[1, ]
    c <- cor(spectra[-1, ], spectra[-1, idx])
    cv <- cov(spectra[-1, ], spectra[-1, idx])
    cv[which(c < cthr)] <- 0
    cv[which(cv < 0)] <- 0
    result <- rbind(ppm, cv[, 1])
    return(result)
}
