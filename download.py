# %% imports
import re
import argparse
import requests

# Create the parser
parser = argparse.ArgumentParser(description='GeoLifeCLEF 2024 file downloader')

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

for k, tab in download_struct.items():
    group = parser.add_argument_group(k)
    for v in tab:
        group.add_argument(f'--{v}', action='store_true', help='Set the flag')

# Parse the arguments
args = parser.parse_args()

# select files to download
repository = 'https://lab.plantnet.org/seafile/d/bdb829337aa44a9489f6'
url_struct = f'{repository}/files/?p=/{{}}/{{}}'

def find_url(category, file):
    response = requests.get(url_struct.format(category, file), timeout=60)
    url_key = "rawPath"
    pattern = f"{url_key}: '([^']+)'"
    # Find the first match
    match = re.search(pattern, response.text)

    if match:
        return match.group(1).replace('\\u002D', '-') + '?raw=1'
    else:
        print(f'Failed to find url for {category}/{file}')


def download_file(url, filename):
    # Send a HTTP request to the URL of the zipfile
    response = requests.get(url)

    # Make sure the request was successful
    assert response.status_code == 200, 'Failed to download file'

    # Write the content of the response to a zipfile
    with open(filename, 'wb') as f:
        f.write(response.content)


for k, tab in download_struct.items():
    for item in tab:
        if getattr(args, f'{item}'):
            url = find_url(k, item)
            print(url)
            download_file(url, item)
