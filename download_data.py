import argparse
import hashlib
import os
import urllib.request
from pathlib import Path

import tqdm


BASE_URL = "https://lilacdn.azureedge.net/geolifeclef-2020/"


def download_file(url, filename, resume=True):
    """
    Download a file from URL with the option to resume previously started download.
    """
    filename = Path(filename)

    with urllib.request.urlopen(url) as response:
        total_size = int(response.getheader("Content-Length"))

    headers = {}

    if resume and filename.exists():
        downloaded_size = filename.lstat().st_size

        if downloaded_size > 0:
            if downloaded_size >= total_size:
                print("CACHED", flush=True)
                return

            headers["Range"] = "bytes={}-".format(downloaded_size)
            open_mode = "ab"
            print("RESUMING... ", end="", flush=True)
    else:
        open_mode = "wb"
        downloaded_size = 0

    request = urllib.request.Request(url, headers=headers)

    print("")
    pbar = tqdm.tqdm(
        total=total_size,
        unit="B",
        unit_scale=True,
        unit_divisor=1024,
        mininterval=0.5,
        maxinterval=2,
        miniters=0,
    )

    with open(filename, open_mode) as fp:
        with urllib.request.urlopen(request) as response:
            if response.getcode() == 206:
                pbar.update(downloaded_size)

            for line in response:
                fp.write(line)
                pbar.update(len(line))

    pbar.close()


def check_file_md5sum(filename, md5sum):
    """
    Check the integrity of a file given an MD5 hash
    """
    hash_md5 = hashlib.md5()
    with open(filename, "rb") as fp:
        for chunk in iter(lambda: fp.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest() == md5sum


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Downloads the GeoLifeCLEF dataset",
    )
    parser.add_argument(
        "--force",
        dest="resume",
        action="store_false",
        help="forces to re-download everything (no download resume)",
    )
    parser.add_argument(
        "output_path",
        type=Path,
        help="output path to download and extract the data",
    )
    args = parser.parse_args()

    output_path = args.output_path
    os.makedirs(output_path, exist_ok=True)

    # Download files list with hashes
    filename = "hashes.md5"
    print("Downloading {}... ".format(filename), end="", flush=True)
    download_file(BASE_URL + filename, output_path / filename, resume=False)
    print("SUCCESS", flush=True)

    # Read files list
    filenames, hashes = [], []
    with open(output_path / filename, "r") as fp:
        for line in fp.readlines():
            md5sum, filename = line.split()
            filenames.append(filename)
            hashes.append(md5sum)

    # Download and check integrity of files
    n_files = len(filenames)
    for i, (filename, md5sum) in enumerate(zip(filenames, hashes)):
        print(
            "Downloading {} [{}/{}]... ".format(filename, i + 1, n_files),
            end="",
            flush=True,
        )
        download_file(BASE_URL + filename, output_path / filename, resume=args.resume)

        print("Checking integrity... ", end="", flush=True)
        valid = check_file_md5sum(output_path / filename, md5sum)

        if not valid:
            print("INTEGRITY ERROR", flush=True)
        else:
            print("SUCCESS", flush=True)
