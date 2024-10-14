import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
from tqdm import tqdm
import requests
import rasterio
import rasterio.mask
from io import BytesIO
from concurrent.futures import ThreadPoolExecutor, as_completed


def download_and_process(variable, year, month, polygon, base_url, session):
    """
    Download and process a single raster file.

    This function downloads a raster file for a specific variable, year, and month.
    It masks the raster file using the given polygon and saves the masked raster as a new file.

    Args:
        variable (str): The climate variable to download (e.g., 'pr', 'tasmax').
        year (int): The year of the data.
        month (str): The month of the data, zero-padded (e.g., '01', '02').
        polygon (GeoDataFrame): A GeoDataFrame containing the polygon geometry for masking.
        base_url (str): The base URL for downloading the raster data.
        session (requests.Session): A requests session for making HTTP requests.

    Raises:
        Exception: If an error occurs during the download or file processing.
    """
    link = f"{base_url}{variable}/CHELSA_{variable}_{month}_{year}_V.2.1.tif"
    try:
        r = session.get(link, stream=True, timeout=10)
        if r.status_code == 200:
            with rasterio.open(BytesIO(r.content)) as src:
                out_image, out_transform = rasterio.mask.mask(src, polygon["geometry"], crop=True)
                out_meta = src.meta
            out_meta.update({
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform,
                "compress": "deflate"
            })
            with rasterio.open(f"./rasters{month}_{year}_{variable}.tif", "w", **out_meta) as dest:
                dest.write(out_image)
        else:
            print(f"Failed to download {link}, status code: {r.status_code}")
    except Exception as e:
        print(f"Error processing {link}: {e}")


def download_variable(variable, years, months, polygon, base_url):
    """
    Download and process all files for a specific variable over multiple years and months.

    This function uses a thread pool to download and process climate data for a given variable
    across all specified years and months. It applies a polygon mask to each downloaded file
    and saves the processed raster files.

    Args:
        variable (str): The climate variable to download (e.g., 'pr', 'tasmax').
        years (list): A list of years to download (e.g., [1979, 1980, ...]).
        months (list): A list of months to download (e.g., ['01', '02', ...]).
        polygon (GeoDataFrame): A GeoDataFrame containing the polygon geometry for masking.
        base_url (str): The base URL for downloading the raster data.
    """
    with requests.Session() as session:  # Reuse session for efficiency
        tasks = []
        with ThreadPoolExecutor(max_workers=10) as executor:
            with tqdm(total=len(years) * len(months), desc=f'Processing variable: {variable}') as pbar:
                for year in years:
                    for month in months:
                        task = executor.submit(download_and_process, variable, year, month, polygon, base_url, session)
                        tasks.append(task)

                for future in as_completed(tasks):
                    pbar.update(1)


def download_all(variables, years, months, polygon, base_url):
    """
    Download and process all climate variables over multiple years and months.

    This function iterates over all specified climate variables, downloading and processing
    raster files for each combination of variable, year, and month. It applies a polygon mask
    to each downloaded file and saves the processed raster files.

    Args:
        variables (list): A list of climate variables to download (e.g., ['pr', 'tasmax', 'tas']).
        years (list): A list of years to download (e.g., [1979, 1980, ...]).
        months (list): A list of months to download (e.g., ['01', '02', ...]).
        polygon (GeoDataFrame): A GeoDataFrame containing the polygon geometry for masking.
        base_url (str): The base URL for downloading the raster data.
    """
    for variable in variables:
        download_variable(variable, years, months, polygon, base_url)


if __name__ == "__main__":
    # Load and process the CSV into a GeoDataFrame
    df = pd.read_csv("/path/to/occurrences.csv", delimiter=";", low_memory=False)
    geometry = gpd.points_from_xy(df.lon,
                                  df.lat)  # Create geometry column (corrected lat/lon order)
    geo_df = gpd.GeoDataFrame(df, geometry=geometry)

    # Calculate bounds and create polygon geometry
    bounds = geo_df.total_bounds
    bounds[:2] -= 10  # Expand bounds by 1 degree
    bounds[2:] += 10
    min_x, min_y, max_x, max_y = bounds

    # Define the vertices of the expanded polygon
    lon_list = [min_x, max_x, max_x, min_x, min_x]
    lat_list = [min_y, min_y, max_y, max_y, min_y]
    polygon_geom = Polygon(zip(lon_list, lat_list))
    polygon = gpd.GeoDataFrame(crs="epsg:4326", geometry=[polygon_geom])

    # Define base URL, variables, years, and months
    base_url = "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/monthly/"
    variables = ["pr", "tasmax", "tas", "tasmin"]
    years = list(range(1979, 2000))
    months = [str(i).zfill(2) for i in range(1, 13)]

    download_all(variables, years, months, polygon, base_url)
