from torch.utils.data import Dataset
import torch

import numpy as np

from .environmental_raster_tools import PatchExtractor
from .common import load_patch


class DatasetGLC20(Dataset):
    def __init__(self, labels, dataset, ids, patches, rasters=None, patch_extractor=None, use_rasters=True):
        """
        :param labels: list of labels
        :param dataset: list of (latitude, longitude)
        :param ids: list of identifiers
        :param rasters: path to the rasters root
        :param patches: path to the patches root
        """
        self.labels = labels
        self.dataset = dataset
        self.ids = ids

        self.one_hot_size = 34
        self.one_hot = np.eye(self.one_hot_size)

        self.rasters = rasters
        self.patches = patches

        if patch_extractor is None and rasters is not None and use_rasters:
            # 256 is mandatory as images have been extracted in 256 and will be stacked in the __getitem__ method
            patch_extractor = PatchExtractor(rasters, size=256, verbose=True)
            patch_extractor.add_all_rasters()

        self.extractor = patch_extractor

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        latitude = self.dataset[idx][0]
        longitude = self.dataset[idx][1]
        id_ = str(self.ids[idx])

        patches = load_patch(id_, self.patches)

        # FIXME add back landcover one hot encoding?
        # lc = patches[3]
        # lc_one_hot = np.zeros((self.one_hot_size,lc.shape[0], lc.shape[1]))
        # row_idx = np.arange(lc.shape[0]).reshape(lc.shape[0], 1)
        # col_idx = np.tile(np.arange(lc.shape[1]), (lc.shape[0], 1))
        # lc_one_hot[lc, row_idx, col_idx] = 1

        # extracting patch from rasters
        if self.extractor is not None:
            environmental_patches = self.extractor[(latitude, longitude)]
            patches = patches + tuple(environmental_patches)

        # Concatenate all patches into a single tensor
        patches = np.concatenate(patches)

        # Transpose data to (CHANNELS, WIDTH, HEIGHT)
        patches = np.transpose(patches, (2, 0, 1))

        return torch.from_numpy(patches).float(), self.labels[idx]
