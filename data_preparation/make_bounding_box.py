​​import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, mapping

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
