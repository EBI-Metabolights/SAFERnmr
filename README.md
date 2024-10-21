<img width="364" alt="image" src="https://github.com/user-attachments/assets/7a12cd62-04f7-4cd4-a21c-6e52c54f6686">
## SAFER
(Spectral Annotation by Feature Extraction and Reference matching)

[![DOI](https://zenodo.org/badge/531570333.svg)](https://zenodo.org/doi/10.5281/zenodo.10022482)


# Background
This is a beta repo for the SAFER approach to 1D NMR data annotation (v2.0.3).

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

# To use this package (still writing this):
   1) set up the params file `params.yaml`
   2) ensure the 4 necessary files are present in the locations given in `params.yaml`
   3) Run in R:
      ```
      devtools::document('replace_with_cloned_github_directory')
      pipeline('path_to_params.yaml')
      browse_evidence('path_to_data_directory')
      ```


Feel free to suggest improvements and report bugs on this repo!


