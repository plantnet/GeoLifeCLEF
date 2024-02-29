"""
File: download.py
Author: maximiliense
Date: 2024-02-08
Description: This script downloads files from the GeoLifeCLEF 2024 dataset.
"""
import re
import os
import argparse
import requests

from tqdm import tqdm

variables = [
    'climate',
    'elevation',
    'humanfootprint',
    'landcover',
    'soilgrids',
    'satellitepatches',
    'satellitetimeseries'
]

rasters = {
    'climate': 'EnvironmentalRasters/Climate.zip',
    'elevation': 'EnvironmentalRasters/Elevation.zip',
    'humanfootprint': 'EnvironmentalRasters/HumanFootprint.zip',
    'landcover': 'EnvironmentalRasters/LandCover.zip',
    'soilgrids': 'EnvironmentalRasters/Soilgrids.zip'
}

presence_only = {
    'climate': {
        'not_cube': [
            'EnvironmentalRasters/Climate/GLC24-PO-train-bioclimatic-average.csv',
            'EnvironmentalRasters/Climate/GLC24-PO-train-bioclimatic-monthly.csv'],
        'cube': [
            'EnvironmentalRasters/Climate/GLC24-PO-train-bioclimatic-average.csv',
            'EnvironmentalRasters/Climate/Climatic_Monthly_2000-2019_cubes/GLC24-PO-train-bioclimatic-monthly.zip'
        ]
    },
    'elevation': {
        'not_cube': ['EnvironmentalRasters/Elevation/GLC24-PO-train-elevation.csv']
        },
    'humanfootprint': {
        'not_cube': ['EnvironmentalRasters/HumanFootprint/GLC24-PO-train-human-footprint.csv']
        },
    'landcover': {
        'not_cube': ['EnvironmentalRasters/LandCover/GLC24-PO-train-landcover.csv']},
    'soilgrids': {
        'not_cube': ['EnvironmentalRasters/Soilgrids/GLC24-PO-train-soilgrids.csv']
        },
    'satellitepatches': {
        'not_cube': [
            'SatellitePatches/PO-Train-SatellitePatches-NIR.zip',
            'SatellitePatches/PO-Train-SatellitePatches-RGB.zip'
        ]
    },
    'satellitetimeseries': {
        'not_cube': [
            'SatelliteTimeSeries/GLC24-PO-train-landsat-time-series-swir2.csv',
            'SatelliteTimeSeries/GLC24-PO-train-landsat-time-series-swir1.csv',
            'SatelliteTimeSeries/GLC24-PO-train-landsat-time-series-red.csv',
            'SatelliteTimeSeries/GLC24-PO-train-landsat-time-series-nir.csv',
            'SatelliteTimeSeries/GLC24-PO-train-landsat-time-series-green.csv',
            'SatelliteTimeSeries/GLC24-PO-train-landsat-time-series-blue.csv',
        ],
        'cube': [
            'SatelliteTimeSeries/cubes/GLC24-PO-train-landsat-time-series.zip'
        ]
    }
}

presence_absence = {
    'climate': {
        'not_cube': [
            'EnvironmentalRasters/Climate/GLC24-PA-test-bioclimatic-average.csv',
            'EnvironmentalRasters/Climate/GLC24-PA-train-bioclimatic-average.csv',
            'EnvironmentalRasters/Climate/GLC24-PA-test-bioclimatic-monthly.csv',
            'EnvironmentalRasters/Climate/GLC24-PA-train-bioclimatic-monthly.csv'],
        'cube': [
            'EnvironmentalRasters/Climate/GLC24-PA-test-bioclimatic-average.csv',
            'EnvironmentalRasters/Climate/GLC24-PA-train-bioclimatic-average.csv',
            'EnvironmentalRasters/Climate/Climatic_Monthly_2000-2019_cubes/GLC24-PA-test-bioclimatic-monthly.zip',
            'EnvironmentalRasters/Climate/Climatic_Monthly_2000-2019_cubes/GLC24-PA-train-bioclimatic-monthly.zip',
        ]
    },
    'elevation': {
        'not_cube': [
            'EnvironmentalRasters/Elevation/GLC24-PA-train-elevation.csv',
            'EnvironmentalRasters/Elevation/GLC24-PA-test-elevation.csv']
        },
    'humanfootprint': {
        'not_cube': [
            'EnvironmentalRasters/HumanFootprint/GLC24-PA-train-human-footprint.csv',
            'EnvironmentalRasters/HumanFootprint/GLC24-PA-test-human-footprint.csv'
            ]
        },
    'landcover': {
        'not_cube': [
            'EnvironmentalRasters/LandCover/GLC24-PA-train-landcover.csv',
            'EnvironmentalRasters/LandCover/GLC24-PA-test-landcover.csv'
            ]
    },
    'soilgrids': {
        'not_cube': [
            'EnvironmentalRasters/Soilgrids/GLC24-PA-train-soilgrids.csv',
            'EnvironmentalRasters/Soilgrids/GLC24-PA-test-soilgrids.csv'
        ]
    },
    'satellitepatches': {
        'not_cube': [
            'SatellitePatches/PA-Train-SatellitePatches-NIR.zip',
            'SatellitePatches/PA-Train-SatellitePatches-RGB.zip',
            'SatellitePatches/PA-Test-SatellitePatches-NIR.zip',
            'SatellitePatches/PA-Test-SatellitePatches-RGB.zip'
        ]
    },
    'satellitetimeseries': {
        'not_cube': [
            'SatelliteTimeSeries/GLC24-PA-train-landsat-time-series-swir2.csv',
            'SatelliteTimeSeries/GLC24-PA-train-landsat-time-series-swir1.csv',
            'SatelliteTimeSeries/GLC24-PA-train-landsat-time-series-red.csv',
            'SatelliteTimeSeries/GLC24-PA-train-landsat-time-series-nir.csv',
            'SatelliteTimeSeries/GLC24-PA-train-landsat-time-series-green.csv',
            'SatelliteTimeSeries/GLC24-PA-train-landsat-time-series-blue.csv',
            'SatelliteTimeSeries/GLC24-PA-test-landsat-time-series-swir2.csv',
            'SatelliteTimeSeries/GLC24-PA-test-landsat-time-series-swir1.csv',
            'SatelliteTimeSeries/GLC24-PA-test-landsat-time-series-red.csv',
            'SatelliteTimeSeries/GLC24-PA-test-landsat-time-series-nir.csv',
            'SatelliteTimeSeries/GLC24-PA-test-landsat-time-series-green.csv',
            'SatelliteTimeSeries/GLC24-PA-test-landsat-time-series-blue.csv'
        ],
        'cube': [
            'SatelliteTimeSeries/cubes/GLC24-PA-train-landsat-time-series.zip',
            'SatelliteTimeSeries/cubes/GLC24-PA-test-landsat-time-series.zip'
        ]
    }
}

metadata = {
    'po': [
        'PresenceOnlyOccurrences/GLC24-PO-metadata-train.csv'
        ],
    'pa': [
        'PresenceAbsenceSurveys/GLC24-PA-metadata-test.csv',
        'PresenceAbsenceSurveys/GLC24-PA-metadata-train.csv'
    ]
}
repository = 'https://lab.plantnet.org/seafile/d/bdb829337aa44a9489f6'
url_struct = f'{repository}/files/?p=/{{}}'

def find_url(file):
    """Given the folder and the file, return the url to download the file

    Args:
        category (str): folder containing the file
        file (str): file to download

    Raises:
        requests.exceptions.HTTPError: exception raised when the url is not found

    Returns:
        str: the url to direct download the file
    """
    response = requests.get(url_struct.format(file), timeout=60)
    url_key = "rawPath"
    pattern = f"{url_key}: '([^']+)'"
    # Find the first match
    match = re.search(pattern, response.text)

    if match:
        return match.group(1).replace('\\u002D', '-') + '?raw=1'

    raise requests.exceptions.HTTPError(f'Failed to find url for {file}')


def check_if_file_complete(path, response):
    """Given a path and a response from a request returns True or False
    if the file needs to be re-downloaded...

    Args:
        path (str): the location of the file on disk
        response (requests.models.Response): the response to the file request

    Returns:
        bool: True if the file needs to be downloaded, False otherwise.
    """
    # check if file exists
    if os.path.exists(path):
        # get file size from response
        file_size_server = int(response.headers.get('Content-Length', 0))
        file_size_downloaded = os.path.getsize(path)
        if file_size_server != file_size_downloaded:
            return True
        else:
            print(f'{path} already downloaded and complete.')
            return False
    else:
        return True


def download_file(url, filename):
    """Download a file given an url.

    Args:
        url (str): the url of the file to download
        filename (str): the location where to write the file
    """
    # Send a HTTP request to the URL of the zipfile
    response = requests.get(url, timeout=60, stream=True)
    # Make sure the request was successful
    assert response.status_code == 200, 'Failed to download file'

    # Download file only if needs to be
    if check_if_file_complete(filename, response):
        total_size_in_bytes = int(response.headers.get('content-length', 0))
        block_size = 1024  # 1 Kilobyte

        # Create dest directory if required
        if not os.path.exists(os.path.dirname(filename)):
            print(f'mkdir {os.path.dirname(filename)}')
            os.makedirs(os.path.dirname(filename))

        progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)

        with open(filename, 'wb') as file:
            for data in response.iter_content(block_size):
                progress_bar.update(len(data))
                file.write(data)

        progress_bar.close()

        if total_size_in_bytes not in (0, progress_bar.n):
            raise requests.exceptions.RequestException('Error.. size do not match...')


def process_download(file, data):
    """Download the file in category cat into the folder data. The function
    handles the exception that can rise from the sub-methods.

    Args:
        cat (str): the folder containing the file
        file (str): the file itself
        data (str): the folder to save the data in
    """
    try:
        u = find_url(file)
        print(f'Downloading {u} ({file})')
        download_file(u, f'{data}/{file}')
    except requests.exceptions.HTTPError:
        print(f'Failed to find url for {file}')
    except (AssertionError, requests.exceptions.RequestException):
        print(f'Failed to download file {file}\n\t{u}')

def process_option(struct, cube, variable, output):
    if variable in struct:
        if cube and 'cube' in struct[variable]:
            l = struct[variable]['cube']
        else:
            l = struct[variable]['not_cube']
        for f in l:
            process_download(f, output)
    
    
if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(
        description=('GeoLifeCLEF 2024 file downloader. '
                     'The script permits to download individual '
                     'files or all files at once.'))

    parser.add_argument(
        '--data', default='data', help='Destination for the downloads (default: data)')

    group = parser.add_argument_group('Data type')
    group.add_argument('--raster', action='store_true', help='Download full rasters')
    group.add_argument('--presence-only', action='store_true',
                       help='Download PO metadata and pre-extracted data if option set to true')
    group.add_argument('--presence-absence', action='store_true',
                       help='Download PA metadata and pre-extracted data if option set to true')
    group.add_argument('--pre-extracted', action='store_true',
                       help='Download pre-extracted data of PO and/or PA depending on the options set')

    group.add_argument('--cube', action='store_true',
                       help='Download cube instead of csv when available')

    group = parser.add_argument_group('Variables')
    for v in variables:
        group.add_argument(f'--{v}', action='store_true',
                           help=f'Download {v} raster or pre-extracted data')
    group.add_argument('--all-variables', action='store_true',
                       help='Select all variables')

    # Parse the arguments
    args = parser.parse_args()

    # Create output folder
    if not os.path.exists(args.data):
        print(f'mkdir {args.data}')
        os.mkdir(args.data)

    # download metadata
    if args.presence_only:
        process_download(metadata['po'][0], args.data)

    if args.presence_absence:
        for f in metadata['pa']:
            process_download(f, args.data)

    for v in variables:
        # process raster
        if args.all_variables or getattr(args, v):
            if args.raster:
                if v in rasters:
                    process_download(rasters[v], args.data)
            # process presence only
            if args.pre_extracted and args.presence_only:
                process_option(presence_only, args.cube, v, args.data)
            # process presence absence
            if args.pre_extracted and args.presence_absence:
                process_option(presence_absence, args.cube, v, args.data)
