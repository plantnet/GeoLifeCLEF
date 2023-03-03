# Author: Benjamin Deneu <benjamin.deneu@inria.fr>
#         Theo Larcher <theo.larcher@inria.fr>
#
# License: GPLv3
#
# Python version: 3.10.6

import itertools
import logging
import math
import os
from abc import abstractmethod
from logging import warning

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split


class TimeSeriesProvider(object):
    def __init__(self, root_path, normalize=False) -> None:
        self.root_path = root_path
        self.normalize = normalize
        self.nb_layers = 0
        self.min_sequence = 0
        self.max_sequence = 0
        
    @abstractmethod
    def __getitem__(self, item):
        pass
    
    def __repr__(self):
        return self.__str__()
    
    @abstractmethod
    def __str__(self):
        result = f"{'-' * 50} \n"
        result += f'nb_layers: {self.nb_layers}\n'
        result += f'min_sequence: {self.min_sequence}\n'
        result += f'max_sequence: {self.max_sequence}\n'
        result += f'bands_names: {self.bands_names}\n'
        result += '-' * 50
        return result
    
    def __len__(self):
        return self.nb_layers
    
    def plot_ts(self, item):
        # patch = self[item]
        # if self.nb_layers==1:
        #     eos_start_ind = np.where(patch[0,0]==self.eos_replace_value)[0][0]
        #     plt.figure(figsize=(20, 20))
        #     plt.plot(range(eos_start_ind), patch[0,0,:eos_start_ind], '-.', c='blue', marker='+')
        #     plt.plot(range(eos_start_ind-1, patch.shape[2]), patch[0,0,eos_start_ind-1:], '', c='red', marker=['']+['+']*(patch.shape[2]-(eos_start_ind-1)))
        #     plt.title(item)
        #     #plt.xticks(patch[0], )
        # else:
        #     # calculate the number of rows and columns for the subplots grid
        #     rows = int(math.ceil(math.sqrt(self.nb_layers)))
        #     cols = int(math.ceil(self.nb_layers / rows))

        #     # create a figure with a grid of subplots
        #     fig, axs = plt.subplots(rows, cols, figsize=(20, 20))

        #     # flatten the subplots array to easily access the subplots
        #     axs = axs.flatten()

        #     # loop through the layers of patch data
        #     for i in range(self.nb_layers):
        #         # display the layer on the corresponding subplot
        #         axs[i].plot(patch[i], '--', c='blue', marker='+')
        #         axs[i].set_title('layer_{}: {}'.format(i, self.bands_names[i]))
        #         axs[i].axis('off')

        #     # remove empty subplots
        #     for i in range(self.nb_layers, rows*cols):
        #         fig.delaxes(axs[i])

        # # show the plot
        # plt.show()
        pass

class MetaTimeSeriesProvider(TimeSeriesProvider):
    def __init__(self, providers, transform=None):
        self.providers = providers
        self.layers_length = [provider.nb_layers for provider in self.providers]
        self.nb_layers = sum(self.layers_length)
        self.bands_names = [provider.bands_names for provider in self.providers]
        self.transform = transform
        self.eos_replace_value = [provider.eos_replace_value for provider in self.providers]
    
    def __getitem__(self, item):
        patch = np.concatenate([provider[item] for provider in self.providers], axis=1)
        if self.transform:
            patch = self.transform(patch)
        return patch
    
    def __str__(self):
        result = 'Providers:\n'
        for provider in self.providers:
            result += str(provider)
            result += '\n'
        return result

class CSVTimeSeriesProvider(TimeSeriesProvider):
    def __init__(self, 
                 ts_data_path,
                 normalize=False,
                 ts_id='timeSerieID',
                 features_col=[],
                 eos_replace_value=-1) -> None:
        super().__init__(normalize)
        self.ts_id = ts_id
        self.ts_data_path = ts_data_path
        self.ts_data = pd.read_csv(ts_data_path, sep=';', nrows=100).set_index(self.ts_id, drop=False)
        self.ts_data = self.ts_data.replace('eos', eos_replace_value).astype(np.int16)
        self.eos_replace_value = eos_replace_value
        self.max_sequence, self.min_sequence = self.get_min_max_sequence(eos_replace_value)
        self.nb_layers = 1
        if not features_col:
            self.features_col = list(self.ts_data.columns[1:])
        elif len(set(features_col).intersection(self.ts_data.columns)) == len(self.ts_data.columns):
            self.features_col = features_col
        else:
            raise KeyError('Some values in `features_col` do not match the `ts_data` column names.')
        self.bands_names = os.path.basename(os.path.splitext(ts_data_path)[0])
            
    def get_min_max_sequence(self, eos_replace_value):
        min_seq = len(self.ts_data.columns)
        max_row = self.ts_data.loc[(self.ts_data == eos_replace_value).sum(axis=1).idxmax()]
        max_seq = (max_row != eos_replace_value).sum()
        return min_seq, max_seq
        
    def __getitem__(self, item):
        # coordinates must be GPS-like, in the 4326 CRS
        tensor = np.array([self.ts_data.loc[item[self.ts_id], self.features_col]])
        tensor = np.expand_dims(tensor, axis=0)
        return tensor
    
    def __len__(self):
        return self.max_sequence

    
class MultipleCSVTimeSeriesProvider(TimeSeriesProvider):
    # Be careful not to place the label .csv with the data .csv and leaving select=None as the provider would then list all .csv files as data including the label file.
    def __init__(self, root_path,
                 select=[], # ['red', 'green', 'blue', 'ir', 'swir1', 'swir2']
                 normalize=False,
                 ts_id='timeSerieID',
                 features_col=[],
                 eos_replace_value=-1) -> None:
        super().__init__(root_path, normalize)
        self.root_path = root_path
        self.ts_id = ts_id
        self.select = [c.lower() for c in select]
        
        files = os.listdir(root_path)
        ts_paths = [f for f in files if f.endswith('.csv')]
        if select:
            select = [f'time_series_{r}.csv' for r in select]
            ts_paths = [r for r in files if r in select]
            if len(ts_paths)!=len(select):
                logging.warn(f'Could not find all files based on `select`. Loading only the ones which names match... (see the `ts_paths` attribute for the complete list)')
        self.ts_paths = ts_paths
        self.ts_providers = [CSVTimeSeriesProvider(root_path+path, normalize=normalize, ts_id=ts_id, features_col=features_col, eos_replace_value=eos_replace_value) for path in ts_paths]
        self.nb_layers = len(self.ts_providers)
        self.bands_names = [ts_.bands_names for ts_ in self.ts_providers]
        self.min_sequence = min(list(itertools.chain.from_iterable([[ts_.min_sequence] for ts_ in self.ts_providers])))
        self.max_sequence = max(list(itertools.chain.from_iterable([[ts_.max_sequence] for ts_ in self.ts_providers])))
        
    def __getitem__(self, item):
        return np.concatenate([ts_[item] for ts_ in self.ts_providers], axis=1)
