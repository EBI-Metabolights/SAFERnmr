setwd('/Users/mjudge/Documents/GitHub/icl_nmr_R')
devtools::document()
run_params <- '/Users/mjudge/Documents/current_run/params.yaml'
pars <- yaml::yaml.load_file(run_params, eval.expr = TRUE)


show.me.the.evidence(results.dir = '/Users/mjudge/Documents/current_run')
