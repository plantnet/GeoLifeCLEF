#!/usr/bin/python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Joaquim Estopinan
# email           :joaquim.estopinan@inria.fr
# ==============================================================================

# Imports
import os
import subprocess
import glob
from osgeo import gdal
from pyproj import CRS
os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd
from shapely.geometry import Polygon
import rasterio
import rasterio.mask
import os.path

"""
Use gdal functions to merge multiple hdf files contained in an input folder.

* See Answer 1 <https://gis.stackexchange.com/questions/220281/how-to-use-gdal-to-convert-hdf-multi-band-file-to-a-single-geotiff-file>

:param hdf_folder : folder containing the hdf files to merge
:param output_file: output file name.
"""
# *****************************
hdf_folder  = "/home/jestopin/PHD/data/GLC23/Land_cover_MODIS_Terra_Aqua_500m/"
output_file = "LandCover_MODIS_Terra-Aqua_500m.tif"
epsg_dest   = 4326
# Bounding box
min_x, min_y, max_x, max_y = -32.26344,  26.63842,  35.58677,  72.18392
# *****************************
output_file_new_espg         = hdf_folder + output_file[:-4] + "_epsg" + str(epsg_dest) + ".tif"
output_file_new_espg_cropped = output_file_new_espg[:-4] + "_cropped.tif"



hdf_files = glob.glob(hdf_folder+"*.hdf")
print("hdf_files\t\t:", hdf_files)
print("Number of hdf files:\t", len(hdf_files))

for i, hdf in enumerate(hdf_files):
    print("\n\ti\t:", i, "\n\thdf\t:", hdf)

    hdf_name = "hdf_"+str(i)+".tif"

    # Step 1: Converts each hdf subdataset into a temporary geotiff
    gdal_cmd = ["gdal_translate", "-sds", hdf, "tmp_outs.tif"]
    subprocess.call(" ".join(gdal_cmd), shell=True)

    # Step 2: Merges each subdataset .tif file into a multi-band geotiff (each subdataset in a distinct band)
    gdal_cmd = ["gdal_merge.py", "-separate", "-o", hdf_name, "tmp_outs*tif"]
    subprocess.call(" ".join(gdal_cmd), shell=True)

    # Step 3: Deletes temporary geotiffs
    gdal_cmd = ["rm", "tmp_outs*tif"]
    subprocess.call(" ".join(gdal_cmd), shell=True)


# Step 4: Merges each subdataset .tif file into a multi-band geotiff (each subdataset in a distinct band)
gdal_cmd = ["gdal_merge.py", "-o", hdf_folder+output_file, "hdf_*.tif"]
subprocess.call(" ".join(gdal_cmd), shell=True)

# Step 5: Deletes temporary hdf_*.tif files
gdal_cmd = ["rm", "hdf_*.tif"]
subprocess.call(" ".join(gdal_cmd), shell=True)


# Step 6: Reprojects geotiff in dest_src
f       = gdal.Open(hdf_folder+output_file)
warp    = gdal.Warp(output_file_new_espg, f, dstSRS=CRS.from_epsg(epsg_dest))


# Step 7: Crops reprojected raster in desired bounding box

# WGS84 Bounding box
lon_list = [min_x, max_x, max_x, min_x, min_x]
lat_list = [min_y, min_y, max_y, max_y, min_y]
polygon_geom = Polygon(zip(lon_list, lat_list))
polygon = gpd.GeoDataFrame(crs='epsg:4326', geometry=[polygon_geom])

# Reads src raster
with rasterio.open(output_file_new_espg) as src:
    # Metadata
    out_meta = src.meta
    # Adapts bounding box to src CRS
    polygon_CRS = polygon.to_crs(out_meta['crs'])
    # Only reads cropped data
    out_image, out_transform = rasterio.mask.mask(src, polygon_CRS['geometry'], crop=True)

# Prepares out metadata
out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform,
                 "compress": "deflate"})

# Writes out raster
with rasterio.open(output_file_new_espg_cropped, "w", **out_meta) as dest:
    dest.write(out_image)

# Verification
raster = rasterio.open(output_file_new_espg_cropped)
print("raster.profile\t:\n", raster.profile)

# End
