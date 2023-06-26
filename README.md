# icl_nmr_R

This is a beta repo for the SAFER approach to 1D NMR data annotation. The current branch is [v2.0.3](https://github.com/EBI-Metabolights/icl_nmr_R/tree/v2.0.3).

The guiding principles here are: 
- compound feature shapes, convey most of the information about chemical structure present in 1H1D NMR data
- feature shape + chemical shift confer far more specificity in spectral annotation than peak lists
- increasing the specificity of the information on which matches are based decreases the liklihood of false positives
- empirically/statistically derived relationships in the data are important in guiding annotation, and are more scalable than expert knowledge

SAFER (Spectral Annotation by Feature Extraction and Reference matching) consists of three major steps:
1) feature definition and extraction,
2) feature-based mapping between pure compound reference spectra (PCRSs) and
3) back-fitting of reference-extracted features to dataset spectra to gauge the believablity of the fits

In more detail: 
    1) Feature Shape Extraction
        - a modified version of SubseT Optimization by Reference Matching (STORM) is used to extract hypothetical feature shapes from a dataset of 1D NMR spectra
        - features are extracted and quantified in-place in each spectrum
        - singlet and other non-specific feature shapes are removed (leaving compound features)
        - feature shapes are clustered to simplify the set, yielding a set of shapes to be matched
    2) Matching to PCRSs
        - each feature shape is cross-correlated with each PCRS
        - the shape is least squares fit to the ref region
        - several metrics are recorded and cutoffs are applied for rvalue and pvalue
    3) Backfitting extracted ref-features to dataset spectra
        - feature shape is fit, and this fit is applied to the ref-feature
  
Each fit constitutes a potential association between a region in a reference spectrum and a region in a sample spectrum, as well as the fit values that match their intensities. There will typically be millions of these between the average dataset and the current 1300 PCRSs. These can be thought of as individual pieces of evidence for a given annotation for its region of a given sample spectrum. All the best evidence for each reference in each spectrum can be summed up and weighted by its quality to derive a metabolite-sample score, which is then linked back to each independent piece of peak-specific evidence. 

To use this package (still writing this):
   1) set up the params file
   2) ensure the 4 necessary files are present
   3) Run in R:
      devtools::document('replace_with_cloned_github_directory')
      pipeline('path_to_data_directory')

To Run the Results Viewer (R Shiny app): 
1) Ensure you have R and Rstudio installed - shiny needs these to run.
2) Clone or download the github repo (most functions won't be needed for this demo)
3) Download and expand the demo data here: https://drive.google.com/file/d/1cBw8ZyY703Z1S5httWJ_J7OStgMZPOnQ/view?usp=sharing  
   These are the results files from a recent run on an unaligned dataset, which can be 
   found here (in case you want that as well):
 
   https://www.ebi.ac.uk/metabolights/editor/MTBLS1, with the spectral matrix available here:
   http://ftp.ebi.ac.uk/pub/databases/metabolights/studies/mariana/spectral_matrices/MTBLS1_nmrML_missing_spectralMatrix.RDS 
   ^ This is an auto-generated spectral matrix file extracted from processed data from MTBLS1

In RStudio, run this to build the package locally (like using a library() call):

  devtools::document('replace_with_cloned_github_directory') # e.g. '/Users/mjudge/Documents/GitHub/icl_nmr_R'
  
Run this to start the app (replacing the filepath):

  show_me_the_evidence('replace_with_downloaded_data_directory') # e.g. '/Users/mjudge/Downloads/mtbls1_demo'

Instructions:
- use the heatmap or the search box to select a compound with high scores (or compound of interest)
- select a spetral region of the PCRS, as well as some samples (a subset of high-scoring samples is usually best for the stackplot)
- a stackplot will appear with evidence for the selected peak plotted in blue
- pan and zoom in the PCRS window to move around the spectrum, or use the vshift slide bar to adjust the spacing in the stackplot

Feel free to suggest improvements and report bugs on this repo!


