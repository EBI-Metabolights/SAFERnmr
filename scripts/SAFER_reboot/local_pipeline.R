## SAFER Reboot
## Local Pipeline 
## MTJ MAY2025

pars <- setup_local()
data <- load_data_local(pars)
protofeatures <- compute_protofeatures(pars, data, noise.width.multiple = 2, top.n.peaks = 5, n.cores = 6)

protofeatures$