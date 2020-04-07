from torch.utils.data import Dataset
import os
import torch

import numpy as np

from environmental_raster_glc import PatchExtractor

patch_extractor = None


class DatasetGLC20(Dataset):
    def __init__(self, labels, dataset, ids, rasters, patches, use_rasters=True):
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
        global patch_extractor
        if patch_extractor is None and rasters is not None and use_rasters:
            # 256 is mandatory as images have been extracted in 256 and will be stacked in the __getitem__ method
            patch_extractor = PatchExtractor(rasters, size=256, verbose=True)
            patch_extractor.add_all()

        self.extractor = patch_extractor

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        latitude = self.dataset[idx][0]
        longitude = self.dataset[idx][1]
        id_ = str(self.ids[idx])

        # folders that contain patches
        folder_1 = id_[-2:]
        folder_2 = id_[-4:-2]

        # path to patches
        path = os.path.join(self.patches, folder_1, folder_2, id_)
        path_alti = path + '_alti.npy'
        path_rgb_ir_lc = path + '.npy'

        # extracting patch from rasters
        tensor = None
        if self.extractor is not None:
            tensor = self.extractor[(latitude, longitude)]

        # extracting altitude patch
        alti = np.load(path_alti)

        # extracting rgb infra-red and land cover (5, 256, 256)
        rgb_ir_lc = np.load(path_rgb_ir_lc).transpose((2, 0, 1))

        # transforming landcover in one hot encoding
        lc = rgb_ir_lc[4:5]
        lc_reshaped = lc.reshape((lc.shape[1] * lc.shape[2],))
        lc_one_hot = self.one_hot[lc_reshaped]
        lc_one_hot = lc_one_hot.reshape((self.one_hot_size, lc.shape[1], lc.shape[2]))

        # concatenating all patches
        if tensor is not None:
            tensor = np.concatenate((tensor,
                                     alti.reshape((1, alti.shape[0], alti.shape[1])),
                                     rgb_ir_lc[:4],
                                     lc_one_hot))
        else:
            tensor = np.concatenate((alti.reshape((1, alti.shape[0], alti.shape[1])),
                                     rgb_ir_lc[:4],
                                     lc_one_hot))

        return torch.from_numpy(tensor).float(), self.labels[idx]
