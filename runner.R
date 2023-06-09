
setwd('/Users/mjudge/Documents/GitHub/icl_nmr_R')
devtools::document()
run_params <- '/Users/mjudge/Documents/current_run/params.yaml'
pars <- yaml::yaml.load_file(run_params, eval.expr = TRUE)

ImperialNMRTool::filter.matches(pars)
ImperialNMRTool::score.matches(pars) 

devtools::document()
show.me.the.evidence(results.dir = '/Users/mjudge/Documents/current_run')

# pipeline('/Users/mjudge/Documents/current_run/params.yaml')


ImperialNMRTool::tina(pars)


devtools::document()
usethis::use_pipe()
devtools::document()

run_params <- '/Users/mjudge/Documents/current_run/params.yaml'
pars <- yaml::yaml.load_file(run_params, eval.expr = TRUE)
pipeline('/Users/mjudge/Documents/current_run/params.yaml')

run_params <- '/Users/mjudge/Documents/current_run/params.yaml'
pars <- yaml::yaml.load_file(run_params, eval.expr = TRUE)

ImperialNMRTool::tina(pars)


setwd('/Users/mjudge/Documents/GitHub/icl_nmr_R')
devtools::document()
usethis::use_pipe()
devtools::document()
# ImperialNMRTool::match.features2refs.par.setup(pars)
# ImperialNMRTool::match.features2refs.par.explicit(pars)

run_params <- '/Users/mjudge/Documents/current_run/params.yaml'
pars <- yaml::yaml.load_file(run_params, eval.expr = TRUE)
# ImperialNMRTool::match.features2refs.par.explicit(pars)
ImperialNMRTool::filter.matches(pars)
ImperialNMRTool::score.matches(pars) 

run_params <- '/Users/mjudge/Documents/current_run/params.yaml'
pars <- yaml::yaml.load_file(run_params, eval.expr = TRUE)

ImperialNMRTool::score.matches(pars) 

