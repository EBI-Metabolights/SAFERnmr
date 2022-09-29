# R Package for ICL NMR tool

A package 'wrapper' for the tools and functions written by michael judge and goncalo graca, as part of efforts to make the overall NMR pipeline functional in galaxy.


## Installation
You can install the package using the R CLI straight from this repository.

1 - run `install.package("devtools")` to install the devtools utility library which includes the `install_github` function.

2 - run `library(devtools)` to load the library into the session.

3 - run `install_github("EBI-Metabolights/icl_nmr_R")` to install the package.

4 - run `library(ImperialNMRTool)` to load the package into the session.

To quickly verify if the package is installed you can run `??hurricane` . This will bring up the manual page for this function, and if loaded means that the package itself was installed and loaded successfully. The manual pages for each function contain the information you will need to run them. 

## Running the workflow

If you are runnning this in a normal (non cluster / workflow manager) environment, it is easiest to do so in RStudio. If you have followed the installation steps above, the package is primed to be run. The entrypoint to the package, `ImperialNMRTool::hurricane()` takes one parameter, a params.yaml file. This file needs to be replete with the parameters you want to run the workflow with. This package has some example data that you can make use of, found in `/inst/extdata`.

### Parameters file
The parameters .yaml file for the workflow takes the following structure:


```
general_pars:
  peaks_location: peaks.RDS     # location of peaks.RDS file.
  spec_location: spec.RDS       # location of spec.RDS file.
  output_dir: ~/                # tells the workflow where to write the output files.

sd_pars:
 cutoff : 0.8                   # standard cutoff for STOCSY
 
am_pars:
 rank_limit: 5                  # rank limit for exporting matches IE take the top 5 matches discard the rest
 dist_thresh: 0.02              # no sense in allowing < precision of db peaks (0.01 ppm for hmdb)
 matchMethod: hungarian_scaled  # or basic, itmin, hungarian, hungarian_scaled
 refdb_file: hmdb_spectra_28FEB2022.RDS # location of reference spectra .RDS file.
 
vis_pars:
 matlab_root:  /Applications/MATLAB_R2021b.app/bin # deprecated, will be removed in future release
``` 
 The four parameters that need to be configured to your local setup in order for the workflow to run are peaks_location, spec_location, output_dor and refdb_file.

 Once this file is prepared, you can run the workflow from RStudio / the command line via `ImperialNMRTool::hurricane("/path/to/your/params.yaml")`
 

#### Help / contact
This repository is maintained by the [MetaboLights](https://www.ebi.ac.uk/metabolights/) team at the EMBL-EBI. You can send any requests or queries to metabolights-help@ebi.ac.uk , quoting this repository.