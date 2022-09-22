#' hurricane
#'
#' Run the NMR pipeline start to finish:
#' 1] statistical decomposition
#' 2] annotation matching
#' 3] export results
#' 4] visualisation [in progress in upstream repository]
#' Also saves to RDS outputs from statistical decomposition
#' A default params.yaml file will be used if none is supplied.
#' @param params_loc location of the params.yaml file containing
#'  pipeline run parameters.
#' @return N/A but .RDS and .TSV files are saved, see description.
#' @export
hurricane <- function(params_loc) {

    if (missing(params_loc)) {
        # load default params
        run_params <- yaml::yaml.load_file("./extdata/default_params.yaml")
    } else {
        # load supplied params
        run_params <- yaml::yaml.load_file(params_loc)
    }
    # it may be that a user supplies their own params file but it
    # is not complete. What do we do then? Fill in any missing fields
    # with the corresponding values from default params?

    peaks <- readRDS(run_params$general_pars$peaks_location)
    spec <- readRDS(run_params$general_pars$spec_location)

    s_d_results <- stat_decomp(
        peaks = peaks,
        spec = spec,
        params = run_params)

    saveRDS(s_d_results$target, stringr::str_c(
        run_params$general_pars$output_dir, "target.RDS"))
    saveRDS(s_d_results$ppeaks, stringr::str_c(
        run_params$general_pars$output_dir, "ppeaks.RDS"))

    matches <- annotation_matching(
        peaks = peaks,
        target = s_d_results$target,
        ppeaks = s_d_results$ppeaks,
        spec = spec,
        params = run_params)

    exportMatches(
        matches = matches,
        X = spec,
        rankLimit = run_params$am_pars$rank_limit)

}
