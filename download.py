"""
File: download.py
Author: maximiliense
Date: 2024-02-08
Description: This script downloads files from the GeoLifeCLEF 2024 dataset.
"""
# %% imports
import re
import os
import argparse
import requests

download_struct = {
    'EnvironmentalRasters': [
        'Climate.zip',
        'Elevation.zip'
        ],
    'PresenceAbsenceSurveys': [
        'GLC24_PA_metadata_test.csv',
        'GLC24_PA_metadata_train.csv'
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
    response = requests.get(url, timeout=60)

    # Make sure the request was successful
    assert response.status_code == 200, 'Failed to download file'

    # Write the content of the response to a zipfile
    with open(filename, 'wb') as f:
        f.write(response.content)

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description='GeoLifeCLEF 2024 file downloader')

    parser.add_argument(
        '--data', default='data', help='Destination for the downloads (default: data)')

    for k, tab in download_struct.items():
        group = parser.add_argument_group(k)
        for v in tab:
            group.add_argument(f'--{v}', action='store_true', help='Set the flag')

    # Parse the arguments
    args = parser.parse_args()

    # Create output folder
    if not os.path.exists(args.data):
        print(f'mkdir {args.data}')
        os.mkdir(args.data)

    # select files to download (those where the corresponding args is true)
    for k, tab in download_struct.items():
        for item in tab:
            if getattr(args, f'{item}'):
                try:
                    u = find_url(k, item)
                    print(f'Downloading {u} ({item})')
                    download_file(u, f'{args.data}/{item}')
                except requests.exceptions.HTTPError:
                    print(f'Failed to find url for {k}/{item}')
                except AssertionError:
                    print(f'Failed to download file {k}/{item}\n\t{u}')
