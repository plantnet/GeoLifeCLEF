​​import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, mapping
from tqdm import tqdm
import requests
import rasterio
import rasterio.mask
from io import BytesIO

df = pd.read_csv("/path/to/occurrences.csv", delimiter=';', low_memory=False)  # Read the CSV file into a pandas DataFrame
geometry = gpd.points_from_xy(df.decimalLongitude, df.decimalLatitude)  # Create a geometry column using the decimalLongitude and decimalLatitude columns from the DataFrame
geo_df = gpd.GeoDataFrame(df, geometry=geometry)  # Create a GeoDataFrame from the DataFrame with the added geometry column

bounds = geo_df.total_bounds  # Calculate the bounds of the GeoDataFrame
bounds[:2] -= 1  # Expand the bounds by 1 unit in all directions
bounds[2:] += 1  # Expand the bounds by 1 unit in all directions
min_x, min_y, max_x, max_y = bounds  # Separate the bounds into variables for convenience

lon_list = [min_x, max_x, max_x, min_x, min_x]  # Define the vertices of a polygon using the expanded bounds
lat_list = [min_y, min_y, max_y, max_y, min_y]  # Define the vertices of a polygon using the expanded bounds

polygon_geom = Polygon(zip(lon_list, lat_list))  # Create a polygon geometry from the lon_list and lat_list
polygon = gpd.GeoDataFrame(crs='epsg:4326', geometry=[polygon_geom])   # Create a GeoDataFrame with the polygon geometry and CRS (coordinate reference system)

base_url = 'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/monthly/'  # Define the base URL for the raster data

variables = ['pr', 'tasmax', 'tas', 'tasmin']  # Define the variables to iterate over
years = list(range(2000, 2020))  # Define the years to iterate over
months = [str(i).zfill(2) for i in range(1,13)]  # Define the months to iterate over

for variable in variables:  # Iterate over the variables
    for year in tqdm(years, desc=f'Variable {variable}'):  # Iterate over the years
        for month in months:  # Iterate over the months
            link = base_url + f'{variable}/CHELSA_{variable}_{month}_{year}_V.2.1.tif'  # Construct the link for the current variable, year, and month
            r = requests.get(link)  # Send an HTTP GET request to the link
            if r.status_code == 200:  # Check if the request was successful (status code 200)
                with rasterio.open(BytesIO(r.content)) as src:  # Open the received content as a raster dataset
                    out_image, out_transform = rasterio.mask.mask(src, polygon['geometry'], crop=True)  # Mask the raster dataset using the polygon geometry
                    out_meta = src.meta
                out_meta.update({"driver": "GTiff",  # Update the metadata for the output raster
                                 "height": out_image.shape[1],
                                 "width": out_image.shape[2],
                                 "transform": out_transform,
                                 "compress": "deflate"})
                with rasterio.open(f'Monthly_Variables/{month}_{year}_{variable}.tif', "w", **out_meta) as dest:
                    dest.write(out_image)  # Create a new raster file with the masked image and updated metadata
