setwd('/Users/mjudge/Documents/GitHub/icl_nmr_R')
devtools::document()

run_params <- '/Users/mjudge/Documents/current_run/params.yaml'
pars <- yaml::yaml.load_file(run_params, eval.expr = TRUE)


show.me.the.evidence(results.dir = '/Users/mjudge/Documents/current_run')

x <- readRDS('/Users/mjudge/Downloads/MTBLS430_1r_noesygppr1d.comp_spectralMatrix.RDS')
dim(x)
# looks like MTBLS1

x <- readRDS('/Users/mjudge/Downloads/MTBLS430_1r_noesygppr1d.comp_spectralMatrix.RDS')
dim(x)
# looks like MTBLS430
# possible data duplication?

x <- readRDS('/Users/mjudge/Downloads/MTBLS395_1r_cpmgpr1d.comp_spectralMatrix.RDS')
dim(x)
# big boi - this will test us