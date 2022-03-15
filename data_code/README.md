# How to process forecast data files 

To create data files for model building,

1) Copy the `model-forecast` and `scores` folders from https://github.com/FluSightNetwork/cdc-flusight-ensemble/tree/first-papers to the main directory of this repository.
2) Create a folder named `data_transformed` if it does not already exist.
3) Run the `get_alldisn.R` script to process forecast files in `model-forecast` and to save output files to the `data_transformed` folder
4) Run the `input_data.R` script to create additional necessary data files in `data_transformed` folder.