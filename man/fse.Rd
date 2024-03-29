% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fse.R
\name{fse}
\alias{fse}
\title{FSE Function}
\usage{
fse(pars)
}
\arguments{
\item{pars}{A list of parameters for the function.}
}
\value{
A list of features extracted from the spectra.
}
\description{
Extracts feature shapes from a matrix of NMR spectra using the FSE algorithm.
# FSE: the spectral matrix is decomposed into compound features using feature
shape extraction. First, a local STOCSY is performed at every point along the
spectrum (within a sliding window of ~ 100 points; enough to capture multiple
resonances within any multiplet). For each of these STOCSYs, the central peak
in the correlation profile (correlation pocket; corrpocket) typically captures
a resonance, as the correlation is 1 by definition at the STOCSY'd point, and
typically falls off as you approach the boundaries of the resonance. This is
taken, with the next highest correlation peak within the window, to form a
rough statistical description of two resonances which have an correlation in
intensity across samples, albeit separated by chemical shift. We term this a
'protofeature'. Importantly, each point and its associated window will capture
the dominant protofeature most associated with that point. This changes for
adjacent points, and many protofeatures will be duplicated multiple times. A
protofeature should be considered as a rough hypothesis about a statistical
association between two resonances, which happen to be sufficiently aligned
so that they produce a coherent signal.
If the windows are all aligned, we can plot the % of central corrpeaks containing
each window point. From this distribution, it is clear that nearly all
protofeatures, including those from noise peaks and real peaks, include the
most central window points. As such, these cannot reliably be used to identify
noise. However, the correlation peak about noise tends to be much smaller, and
characteristically so. As such, a noisewidth can be estimated from this
distribution. This is the origin of the noise.percentile cutoff, which is applied
like so: "given a noise.percentile = 0.99, consider only those protofeatures
for which both peaks have a width > the smallest 1% of peaks". Reducing this
number therefore gives a more selective cut. Protofeatures are also filtered so
that they must be bidirectional (i.e. both peaks must indicate each other as
highest correlated), and must include runs > noise.width.
}
\details{
STORM: Joram Posma's STORM has been adapted and optimized to accept these
protofeatures in the following ways:
- first, since many of the protofeatures are noise, we provide failure modes
  and reporting for the following cases:
  "empty subset",          #' empty subset (no spectrum contains signature)
  "subset degenerated",    #' 1-3 spectra in the subset (not enough spectra to
                             get a reliable correlation)
  "reference degenerated", #' signature degenerates to include < 3 points (not
                             meaningful to correlate shapes)
  "did not converge"       #' subset continues to change after 24 iterations
- additionally, the correlation r and p-value cutoff q are both used during
  both the subset selection and reference update steps.
- we also remove any regions of the reference for which there are fewer than
  minpeak values after r and p value thresholding. This helps avoid noise.

STORM extracts meaningful features using protofeatures to define the region of
interest and a rough sketch of the feature shape highly correlated with each
spectral point. In the future, HCA could be used to cluster potential starting
feature shapes correlated with each driver, or the nonoptimal subset for each
point could be re-STORMed to detect any other feature shapes present. It is
also perfectly reasonable to combine feature shapes from different STORM runs
for a given dataset, as these comprise a list of somewhat independently tested
feature shapes, and duplication is not an issue.
}
