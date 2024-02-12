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
repository = 'https://lab.plantnet.org/seafile/d/bdb829337aa44a9489f6'
url_struct = f'{repository}/files/?p=/{{}}/{{}}'

def find_url(category, file):
    """Given the folder and the file, return the url to download the file

    Args:
        category (str): folder containing the file
        file (str): file to download

    Raises:
        requests.exceptions.HTTPError: exception raised when the url is not found

    Returns:
        str: the url to direct download the file
    """
    response = requests.get(url_struct.format(category, file), timeout=60)
    url_key = "rawPath"
    pattern = f"{url_key}: '([^']+)'"
    # Find the first match
    match = re.search(pattern, response.text)

    if match:
        return match.group(1).replace('\\u002D', '-') + '?raw=1'

    raise requests.exceptions.HTTPError(f'Failed to find url for {category}/{file}')


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

    total_size_in_bytes = int(response.headers.get('content-length', 0))
    block_size = 1024  # 1 Kilobyte

    progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)

    with open(filename, 'wb') as file:
        for data in response.iter_content(block_size):
            progress_bar.update(len(data))
            file.write(data)

    progress_bar.close()

    if total_size_in_bytes not in (0, progress_bar.n):
        raise requests.exceptions.RequestException('Error.. size do not match...')


def process_download(cat, file, data):
    """Download the file in category cat into the folder data. The function
    handles the exception that can rise from the sub-methods.

    Args:
        cat (str): the folder containing the file
        file (str): the file itself
        data (str): the folder to save the data in
    """
    try:
        u = find_url(cat, file)
        print(f'Downloading {u} ({file})')
        download_file(u, f'{data}/{file}')
    except requests.exceptions.HTTPError:
        print(f'Failed to find url for {cat}/{file}')
    except (AssertionError, requests.exceptions.RequestException):
        print(f'Failed to download file {cat}/{file}\n\t{u}')

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(
        description=('GeoLifeCLEF 2024 file downloader. '
                     'The script permits to download individual '
                     'files or all files at once.'))

    parser.add_argument(
        '--data', default='data', help='Destination for the downloads (default: data)')

    parser.add_argument('-all', action='store_true', help='Download all files')

    for k, tab in download_struct.items():
        group = parser.add_argument_group(k)
        for v, s in tab:
            group.add_argument(f'--{v}', action='store_true', help=s)

    # Parse the arguments
    args = parser.parse_args()

    # Create output folder
    if not os.path.exists(args.data):
        print(f'mkdir {args.data}')
        os.mkdir(args.data)

    # select files to download (those where the corresponding args is true)
    for k, tab in download_struct.items():
        for item, _ in tab:
            if getattr(args, f'{item}') or args.all:
                process_download(k, item, args.data)
