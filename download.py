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

download_struct = {
    'EnvironmentalRasters': [
        ('Climate.zip', '26.7GB'),
        ('Elevation.zip', '13.1GB')
        ],
    'PresenceAbsenceSurveys': [
        ('GLC24_PA_metadata_test.csv', '5.4MB'),
        ('GLC24_PA_metadata_train.csv', '97.7MB')
    ],
    'PresenceOnlyOccurrences': [
        ('GLC24_PO_metadata_train.csv', '376.8MB'),
        #'po_train_patches_nir.zip',
        #'po_train_patches_rgb.zip'
    ],
    'SatellitePatches': [
        ('PA_Test_SatellitePatches_NIR.zip', '20.1MB'),
        ('PA_Test_SatellitePatches_RGB.zip', '19.6MB'),
        ('PA_Train_SatellitePatches_NIR.zip', '374.7MB'),
        ('PA_Train_SatellitePatches_RGB.zip', '353.4MB')
    ]
}

download_struct = {
    'EnvironmentalRasters': {
        'climate': ('Climate.zip', '26.7GB'),
        'elevation': ('Elevation.zip', '13.1GB')
    },
    'PresenceAbsenceSurveys': {
        'pa_metadata_test': ('GLC24_PA_metadata_test.csv', '5.4MB'),
        'pa_metadata_train': ('GLC24_PA_metadata_train.csv', '97.7MB')
    },
    'PresenceOnlyOccurrences': {
        'po_metadata_train': ('GLC24_PO_metadata_train.csv', '376.8MB')
    },
    'SatellitePatches': {
        'pa_test_satellite_patches_nir': ('PA_Test_SatellitePatches_NIR.zip', '20.1MB'),
        'pa_test_satellite_patches_rgb': ('PA_Test_SatellitePatches_RGB.zip', '19.6MB'),
        'pa_train_satellite_patches_nir': ('PA_Train_SatellitePatches_NIR.zip', '374.7MB'),
        'pa_train_satellite_patches_rgb': ('PA_Train_SatellitePatches_RGB.zip', '353.4MB')
    }
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
    if the file needs to be redownloaded...

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
        if file_size_server != file_size_downloaded: # TODO add epsilon
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

        progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
        if not os.path.exists(os.path.dirname(filename)):
            print(f'mkdir {os.path.dirname(filename)}')
            os.makedirs(os.path.dirname(filename))
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

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(
        description=('GeoLifeCLEF 2024 file downloader. '
                     'The script permits to download individual '
                     'files or all files at once.'))

    parser.add_argument(
        '--data', default='data', help='Destination for the downloads (default: data)')

    parser.add_argument('--all', action='store_true', help='Download all files')

    urls_dict = {}

    for k, v in download_struct.items():
        all_url_in_group = []
        group = parser.add_argument_group(k)
        for opt, f in v.items():
            urls_dict[opt] = f'{k}/{f[0]}'
            all_url_in_group.append(urls_dict[opt])
            group.add_argument(f'--{opt}', action='store_true', help=f[1])
        urls_dict[f'All{k}'] = all_url_in_group
        group.add_argument(f'--All{k}', action='store_true',
                           help='Download all files in this group')

    # Parse the arguments
    args = parser.parse_args()

    # Create output folder
    if not os.path.exists(args.data):
        print(f'mkdir {args.data}')
        os.mkdir(args.data)

    # select files to download (those where the corresponding args is true)
    for k, v in urls_dict.items():
        v = [v] if isinstance(v, str) else v
        if getattr(args, f'{k}') or args.all:
            for f in v:
                process_download(f, args.data)
