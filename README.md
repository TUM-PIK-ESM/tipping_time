# Code accompanying Ben-Yami et al. 2024 (Science Advances)
Important files:

TT_fingerprints/cdo_script.sh and TT_fingerprints/TT_get_fingerprints.py: The code for creating the AMOC fingerprints from the datasets.

DD23_method1.py : Older code to run method 1 from Ditlevesen & Ditlevsen 2023. In the final manuscript only the variance and autocorrelation timeseries are used from these results.

DD23_MLE.R : Adapted code from DD23 to apply their MLE method to the AMOC fingerprints in TT_fingerprints.

To get Figures 4 and 5, first run DD23_method1.py and DD23_MLE.R, and then use the generated files to run TT_plots.ipynb. 

TT_generalAlternativeModels/generalPlots.ipynb: The code for generating Figures and 1 and 2.

TT_AMOCAlternativeModels/AMOCPlots.ipynb: The code for generating Figure 3.

To download the datasets underlying the analyses of Figures 1-3, visit the following [zenodo folder](https://doi.org/10.5281/zenodo.12549739).
