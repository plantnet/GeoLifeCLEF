# GeoLifeCLEF

This repository is related to the GeoLifeCLEF challenges. The details of each challenge, the data, and all other useful information are present on the challenge pages: 
* [GeoLifeCLEF 2022](https://www.kaggle.com/competitions/geolifeclef-2022-lifeclef-2022-fgvc9)
* [GeoLifeCLEF 2023](https://www.kaggle.com/competitions/geolifeclef-2023-lifeclef-2023-x-fgvc10)
* [GeoLifeCLEF 2024](https://www.kaggle.com/competitions/geolifeclef-2024)
* [GeoLifeCLEF 2025](https://www.kaggle.com/c/geolifeclef-2025)

## Code Base
In this repository, you will find dataloaders, sample_data, and examples to help using the challenge datasets.
- In ``data/sample_data/`` you will find a small sample of the dataset to try codes and loaders.
- ``example_patch_loading.ipynb`` and ``example_patch_loading.py`` give an example of pytorch dataset creation for CNN tensors taking into account different cases.
- ``example_time_series_loading.ipynb`` and ``example_time_series_loading.py`` give an example of pytorch dataset creation for time series tensors taking into account different cases.

## Environment
 We provide a conda environment containing the needed libraries to use this code.
 ```conda env create -f environment.yml```
 ```conda activate glc23```
