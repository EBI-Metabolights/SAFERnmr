## SAFER Reboot
## Local Pipeline 
## MTJ MAY2025
## devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')

pars <- setup_local()
  pars$corrpockets$only.region.between <- c(-0.5, 10)
data <- load_data_local(pars)
protofeatures <- compute_protofeatures(pars, data, noise.width.multiple = 2, top.n.peaks = 5, n.cores = 6)
# Protofeatures are seeds for STORM

# What we're looking for is multiplet shapes which recur across spectra.
sats <- compute_sats(pars, data, protofeatures)
