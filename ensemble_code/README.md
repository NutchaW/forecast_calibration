# Folder structure

This folder contains R scripts for building ensemble forecasts using combination methods in the article and trained parameters (in the `trained_params` subfolder) used in the article to produce the application results. 

# Instructions

- The `cluster_scripts` folder contains R scripts for parameter estimation for BLP and beta mixture combination methods and shell scripts for running the analysis in a cluster computing environment. It also contains subfolders for parameters obtained from cross validation runs (`training params` for training runs and `final_train_pars` for validation runs)
- `combine_params.R` can be used to combine saved parameters in the `cluster_scripts` folder for convenience. The combined parameters are saved in the `final_train` folder.
- `cv_train.R` builds cross validated ensemble forecasts' binned predictive distributions from the BLP and beta mixture methods, and builds training and test ensemble forecasts' binned  from all combination methods. This script also calculates log scores and PIT values for all ensemble forecasts.  
- The `BLPwork_functions_app` folder contains `functions_app.R`, `ls_pit_reli_funs.R`, `tables_plots_functions.R` scripts. `functions_app.R` contains functions for parameter estimation and assembling ensemble PDFs and CDFs. `ls_pit_reli_funs.R` and `tables_plots_functions.R` contain functions for producing plots and tables in the result section of the article. `tables_plots_functions.R` also contains a function for calculating Cramer distance. 
