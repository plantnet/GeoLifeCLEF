import datetime
import os
import threading
import numpy as np
import pandas as pd
from PIL import Image
from raster_extractor_memory import Extractor
#from tqdm.auto import tqdm
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--band', dest='band', default='red')

# initialize timer
start_time = time.time()

## PARAMETERS 
cell_size = 100000
margin = 30
patch_size = 1
output = '/data/zenith/share/GeoLifeCLEF/2023/data/time_series/'

args = vars(parser.parse_args())
band = args['band']

def get_quarter_dates(year, quarter):
    if quarter == 1:
        return '{}.12.02..{}.03.20'.format(year-1, year)
    if quarter == 2:
        return '{}.03.21..{}.06.24'.format(year, year)
    if quarter == 3:
        return '{}.06.25..{}.09.12'.format(year, year)
    if quarter == 4:
        return '{}.09.13..{}.12.01'.format(year, year)
    
years = np.arange(2000,2021)

global extractor
global next_extractor

# function that instantiate an extractor (download a raster tile for extraction of patch) 
def prep_extractor(raster_url, grid_cell, cell_size, patch_size, margin):
    global next_extractor
    next_extractor = Extractor(raster_url, grid_cell, cell_size, patch_size, margin=margin, convert_uint8=False)

# read the csv file and select only the columns 'timeSerieID', 'lon', 'lat'
df_po = pd.read_csv('./occurrences/PO_anonymised_filtered_tsID.csv', usecols=['timeSerieID', 'lon', 'lat', 'x_EPSG3035', 'y_EPSG3035', 'year', 'year_ecodatacube_quarter'], delimiter=';')
df_pa = pd.read_csv('./occurrences/PA_public_anonymised_tsID.csv', usecols=['timeSerieID', 'lon', 'lat', 'x_EPSG3035', 'y_EPSG3035', 'year', 'year_ecodatacube_quarter'], delimiter=';')
df = pd.concat([df_po, df_pa])

# keep only one row for each unique timeSerieID
df = df.drop_duplicates(subset='timeSerieID', keep='first', ignore_index=True)

# Convert the x and y coordinates to integer values representing the grid cells
df['grid_x'] = (df['x_EPSG3035'] // cell_size).astype(int)
df['grid_y'] = (df['y_EPSG3035'] // cell_size).astype(int)

# Create a column with the cell indices as a tuple
df['grid_cell'] = list(zip(df['grid_y'], df['grid_x']))
grid_cells = df['grid_cell'].unique()

values = np.full((df.shape[0], 84), 65535, dtype=np.uint16)
header = []

time_stamp = 0
for year in years:
    for quarter in range(1,5):
        raster_path = '/data/zenith/share/GeoLifeCLEF/2023/data/landsat/lcv_{}_landsat.glad.ard_p50_30m_0..0cm_{}_eumap_epsg3035_v1.1.tif'.format(band, get_quarter_dates(year=year, quarter=quarter))
        print('Extraction of quarter {} {}!'.format(year, quarter))
        header.append('{}_{}'.format(year, quarter))
        try:
            # Prepare the first tile
            next_extractor = Extractor(raster_path, grid_cells[0], cell_size, patch_size, margin=margin, convert_uint8=False)

            # loop on grid_cells (= tiles) to extract patches on each tiles
            #for i, grid_cell in tqdm(enumerate(grid_cells), desc='Total extraction', total=grid_cells.shape[0]):
            for i, grid_cell in enumerate(grid_cells):
                # next_extractor become current extractor
                extractor = next_extractor
                next_extractor = None

                # create new threads to download in parallel next tiles and instantiate the next extractor
                if i < (len(grid_cells)-1):
                    thread_next = threading.Thread(target=prep_extractor, args=(raster_path, grid_cells[i+1], cell_size, patch_size, margin))
                    # start the new thread
                    thread_next.start()

                # get all occurrences that are within the tile
                df_cell = df.loc[df['grid_cell'] == grid_cell]

                # iterate through each row in the dataframe
                #for index, row in tqdm(df_cell.iterrows(), desc='Extraction of tile '+str(grid_cell), total=df_cell.shape[0]):
                for index, row in df_cell.iterrows():
                    if row['year'] > year or (row['year'] == year and row['year_ecodatacube_quarter'] >= quarter):
                        try:
                            timeSerieID = int(row['timeSerieID'])
                            lat = row['lat']
                            lon = row['lon']

                            # call the patch extractor to get the value
                            values[index, time_stamp] = extractor[(lat, lon)]
                            
                        except Exception as e:
                            # if an error occurs for a patch, save the id of the occurrence
                            with open(output+"errors.txt", "a") as f:
                                f.write('{};{};{};{}\n'.format(timeSerieID, year, quarter, e))
                # wait for the next extractor to be ready and get it
                if i < (len(grid_cells)-1):
                    thread_next.join()
        except Exception as e:
            print(e)
            print('error during extraction of quarter {} {}!'.format(year, quarter))
        time_stamp += 1

ids = df['timeSerieID'].to_numpy()
df = pd.DataFrame(data = values, index=ids, columns=header)
df = df.replace(65535, 'eos')

# Save time serie csv
df.to_csv('{}time_series_{}.csv'.format(output, band), index=True, sep=';', index_label='timeSerieID')

print('finished !')
# print total execution time
t = str(datetime.timedelta(seconds = (time.time() - start_time)))
print("--- "+t+" ---")
