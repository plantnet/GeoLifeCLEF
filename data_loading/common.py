from pathlib import Path

import numpy as np
from PIL import Image
import tifffile


def load_patch(patch_id, data_path, landcover_mapping=None, return_arrays=True):
    patch_id = str(patch_id)

    region_id = patch_id[0]
    if region_id == "1":
        region = "fr"
    elif region_id == "2":
        region = "us"
    else:
        raise ValueError("Incorrect 'patch_id' {}, can not extract region id from it".format(patch_id))

    subfolder1 = patch_id[-2:]
    subfolder2 = patch_id[-4:-2]

    filename = Path(data_path) / "patches" / region / subfolder1 / subfolder2 / patch_id

    rgb_filename = filename.with_name(filename.stem + "_rgb.jpg")
    rgb_patch = Image.open(rgb_filename)

    near_ir_filename = filename.with_name(filename.stem + "_near_ir.jpg")
    near_ir_patch = Image.open(near_ir_filename)

    landcover_filename = filename.with_name(filename.stem + "_landcover.tif")
    landcover_patch = tifffile.imread(landcover_filename)
    if landcover_mapping is not None:
        landcover_patch = landcover_mapping[landcover_patch]

    altitude_filename = filename.with_name(filename.stem + "_altitude.tif")
    altitude_patch = tifffile.imread(altitude_filename)

    if return_arrays:
        rgb_patch = np.asarray(rgb_patch)
        near_ir_patch = np.asarray(near_ir_patch)

    return (rgb_patch, near_ir_patch, altitude_patch, landcover_patch)
