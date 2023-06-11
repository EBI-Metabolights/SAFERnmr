setwd('/Users/mjudge/Documents/GitHub/icl_nmr_R')
devtools::document()

setwd('/Users/mjudge/Documents/not_galaxy')
pipeline(params_loc = '/Users/mjudge/Documents/mtbls430/430.params.yaml')

run_params <- '/Users/mjudge/Documents/current_run/params.yaml'
pars <- yaml::yaml.load_file(run_params, eval.expr = TRUE)

devtools::document()
show.me.the.evidence(results.dir = '/Users/mjudge/Documents/current_run')
