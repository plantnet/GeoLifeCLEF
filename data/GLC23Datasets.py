# Author: Benjamin Deneu <benjamin.deneu@inria.fr>
#         Theo Larcher <theo.larcher@inria.fr>
#
# License: GPLv3
#
# Python version: 3.10.6

from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset
import os

from data.GLC23PatchesProviders import MetaPatchProvider
from data.GLC23TimeSeriesProviders import MetaTimeSeriesProvider

class PatchesDataset(Dataset):
    def __init__(
        self,
        occurrences,
        providers,
        transform=None,
        target_transform=None,
        id_name="glcID",
        label_name="speciesId",
        item_columns=['lat', 'lon', 'patchID'],
    ):
        self.occurences = Path(occurrences)
        self.base_providers = providers
        self.transform = transform
        self.target_transform = target_transform
        self.provider = MetaPatchProvider(self.base_providers, self.transform)

        df = pd.read_csv(self.occurences, sep=";", header='infer', low_memory=False)

        self.observation_ids = df[id_name].values
        self.items = df[item_columns]
        self.targets = df[label_name].values

    def __len__(self):
        return len(self.observation_ids)

    def __getitem__(self, index):
        item = self.items.iloc[index].to_dict()

        patch = self.provider[item]

        target = self.targets[index]

        if self.target_transform:
            target = self.target_transform(target)

        return torch.from_numpy(patch).float(), target
    
    def plot_patch(self, index):
        item = self.items.iloc[index].to_dict()
        self.provider.plot_patch(item)


class PatchesDatasetOld(Dataset):
    def __init__(
        self,
        occurrences,
        providers,
        transform=None,
        target_transform=None,
        id_name="glcID",
        label_name="speciesId",
        item_columns=['lat', 'lon', 'patchID'],
    ):
        self.occurences = Path(occurrences)
        self.providers = providers
        self.transform = transform
        self.target_transform = target_transform

        df = pd.read_csv(self.occurences, sep=";", header='infer', low_memory=False)

        self.observation_ids = df[id_name].values
        self.items = df[item_columns]
        self.targets = df[label_name].values

    def __len__(self):
        return len(self.observation_ids)

    def __getitem__(self, index):
        item = self.items.iloc[index].to_dict()

        patches = []
        for provider in self.providers:
            patches.append(provider[item])

        # Concatenate all patches into a single tensor
        if len(patches) == 1:
            patches = patches[0]
        else:
            patches = np.concatenate(patches, axis=0)

        if self.transform:
            patches = self.transform(patches)

        target = self.targets[index]

        if self.target_transform:
            target = self.target_transform(target)

        return torch.from_numpy(patches).float(), target

class TimeSeriesDataset(Dataset):
    def __init__(
        self,
        occurrences,
        providers,
        transform=None,
        target_transform=None,
        id_name="glcID",
        label_name="speciesId",
        item_columns=['timeSerieID'],
    ):
        self.occurences = Path(occurrences)
        self.base_providers = providers
        self.transform = transform
        self.target_transform = target_transform
        self.provider = MetaTimeSeriesProvider(self.base_providers, self.transform)

        df = pd.read_csv(self.occurences, sep=";", header='infer', low_memory=False)

        self.observation_ids = df[id_name].values
        self.items = df[item_columns]
        self.targets = df[label_name].values
    
    def __len__(self):
        return len(self.observation_ids)

    def __getitem__(self, index):
        item = self.items.iloc[index].to_dict()

        patch = self.provider[item]

        target = self.targets[index]

        if self.target_transform:
            target = self.target_transform(target)

        return torch.from_numpy(patch).float(), target

    def plot_ts(self, index):
        item = self.items.iloc[index].to_dict()
        self.provider.plot_ts(item)