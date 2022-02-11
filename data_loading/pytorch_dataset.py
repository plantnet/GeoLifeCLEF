from pathlib import Path

import numpy as np
import pandas as pd

from torch.utils.data import Dataset
import torch

from .environmental_raster import PatchExtractor
from .common import load_patch


class GeoLifeCLEF2021Dataset(Dataset):
    """Pytorch dataset handler for GeoLifeCLEF 2021 dataset.

    Parameters
    ----------
    root : string or pathlib.Path
        Root directory of dataset.
    subset : string, either "train", "val", "train+val" or "test"
        Use the given subset ("train+val" is the complete training data).
    use_rasters : boolean (optional)
        If True, extracts patches from rasters.
    patch_extractor : PatchExtractor object (optional)
        Patch extractor to use if rasters are used.
    transform : callable (optional)
        A function/transform that takes a list of arrays and returns a transformed version.
    target_transform : callable (optional)
        A function/transform that takes in the target and transforms it.
    """
    def __init__(self, root, subset, use_rasters=True, patch_extractor=None, transform=None, target_transform=None):
        self.root = Path(root)
        self.subset = subset
        self.transform = transform
        self.target_transform = target_transform

        possible_subsets = ["train", "val", "train+val", "test"]
        if subset not in possible_subsets:
            raise ValueError("Possible values for 'subset' are: {} (given {})".format(possible_subsets, subset))

        if subset == "test":
            subset_file_suffix = "test"
            self.training_data = False
        else:
            subset_file_suffix = "train"
            self.training_data = True

        df_fr = pd.read_csv(
            self.root / "observations" / "observations_fr_{}.csv".format(subset_file_suffix),
            sep=";",
            index_col="observation_id"
        )
        df_us = pd.read_csv(
            self.root / "observations" / "observations_us_{}.csv".format(subset_file_suffix),
            sep=";",
            index_col="observation_id"
        )
        df = pd.concat((df_fr, df_us))

        if self.training_data and subset != "train+val":
            ind = df.index[df["subset"] == subset]
            df = df.loc[ind]

        self.observation_ids = df.index
        self.coordinates = df[["latitude", "longitude"]].values

        if self.training_data:
            self.targets = df["species_id"].values
        else:
            self.targets = None

        # FIXME: add back landcover one hot encoding?
        # self.one_hot_size = 34
        # self.one_hot = np.eye(self.one_hot_size)

        if use_rasters:
            if patch_extractor is None:
                # 256 is mandatory as images have been extracted in 256 and will be stacked in the __getitem__ method
                patch_extractor = PatchExtractor(self.root / "rasters", size=256)
                patch_extractor.add_all_rasters()

            self.patch_extractor = patch_extractor
        else:
            self.patch_extractor = None

    def __len__(self):
        return len(self.observation_ids)

    def __getitem__(self, index):
        latitude = self.coordinates[index][0]
        longitude = self.coordinates[index][1]
        observation_id = self.observation_ids[index]

        patches = load_patch(observation_id, self.root / "patches")

        # FIXME: add back landcover one hot encoding?
        # lc = patches[3]
        # lc_one_hot = np.zeros((self.one_hot_size,lc.shape[0], lc.shape[1]))
        # row_index = np.arange(lc.shape[0]).reshape(lc.shape[0], 1)
        # col_index = np.tile(np.arange(lc.shape[1]), (lc.shape[0], 1))
        # lc_one_hot[lc, row_index, col_index] = 1

        # Extracting patch from rasters
        if self.patch_extractor is not None:
            environmental_patches = self.patch_extractor[(latitude, longitude)]
            patches = patches + tuple(environmental_patches)

        # Concatenate all patches into a single tensor
        patches = np.atleast_3d(*patches)
        patches = np.concatenate(patches, axis=-1, dtype=np.float32)

        # Transpose data to (CHANNELS, WIDTH, HEIGHT)
        patches = np.transpose(patches, (2, 0, 1))

        # Convert patches to Torch array
        patches = torch.from_numpy(patches)

        if self.transform:
            patches = self.transform(patches)

        if self.training_data:
            target = self.targets[index]

            if self.target_transform:
                target = self.target_transform(target)

            return patches, target
        else:
            return patches
