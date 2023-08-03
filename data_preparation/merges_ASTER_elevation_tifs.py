#!/usr/bin/python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          : Joaquim Estopinan
# email           : estopinan.joaquim@gmail.com
# python_version  : >=3.0
# license         : GPLv3
# description     : Combines "ASTER Global Digital Elevation Model V003" .tif files directly downloaded from
#                   https://search.earthdata.nasa.gov queried with the bounding box file "GLC23.geojson"
# ==============================================================================

# Imports
import os
import subprocess
import sys
from util.tif2cog import tif2cog
import numpy as np
from pathlib import Path
import glob
from datetime import datetime


# ### PARAMS ###
# In
tif_p   = Path("/gpfsscratch/rech/pcz/uzk84jj/ASTER_elevation_data/tifs")
pattern = "/*_dem.tif"

# COG Params
levels        = 6
num_threads   = "ALL_CPUS"
gdal_cachemax = 4096

# Out
out_p      = Path("/gpfsscratch/rech/pcz/uzk84jj/ASTER_elevation_data")
out_name   = "merged_DEM_WGS84"
# ##############



# Constructs target_path
print("****** target_path: ", tif_p)
# List all .tif from target_path
L_tifs = glob.glob(str(tif_p)+pattern)
print("****** L_tifs: ", L_tifs)
print("****** len(L_tifs): ", len(L_tifs))

# out name and folder
out_folder = out_p / out_name
print("****** out_folder: ", out_folder)
out_folder.mkdir(parents=True, exist_ok=True)


# List .txt and .tif files
list_txt = out_name + '.txt'
out_tif  = out_name + '.tif'

# Creates list .txt file
with open(out_folder / list_txt, 'w') as f:
    for line in L_tifs:
        f.write(f"{line}\n")

# --- Creates COG .tif file ---
check = tif2cog(input_file_list=str(out_folder/list_txt),
                output_file=str(out_folder/out_tif),
                levels=levels,
                num_threads=num_threads,
                cachemax=gdal_cachemax
                )
print("\n\t\t\t*** ***\n", check, "\n\t\t\t*** ***\n")


# Projects to WGS84, needs GDAL install
subprocess.call("gdalwarp -overwrite -multi -wo 'NUM_THREADS=ALL_CPUS' -wm 4096 -t_srs EPSG:4326 \
-tap -r near -tr 0.000208333 0.000208333 -co 'COMPRESS=DEFLATE' \
-co 'PREDICTOR=2' -co 'BIGTIFF=YES' merged_DEM.tif merged_DEM_WGS84.tif", shell=True)



# EOF
