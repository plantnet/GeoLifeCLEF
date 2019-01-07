"""
"
"   Author: Maximilien Servajean - mservajean
"   Mail: servajean@lirmm.fr
"   Date: 04/01/2019
"
"   Description: The code to extract environmental tensors and environmental vectors given some environmental rasters.
"
"""
import numpy as np
import rasterio
import re
import os
import warnings
import pandas as pd
import matplotlib.pyplot as plt

MIN_ALLOWED_VALUE = -10000
EPS = 1

# metadata used to setup some rasters
raster_metadata = {
    'alti': {'nan': 0},
    'awc_top': {'nan': -1},
    'bs_top': {'nan': 30},
    'cec_top': {'nan': 0},
    'chbio_1': {'min_val': -10.7, 'max_val': 18.4, 'nan': -11.},
    'chbio_2': {'min_val': 7.8, 'max_val': 21.0, 'nan': 7.},
    'chbio_3': {'min_val': 41.1, 'max_val': 60., 'nan': 40.},
    'chbio_4': {'min_val': 302.7, 'max_val': 777.8, 'nan': 300.},
    'chbio_5': {'min_val': 6.1, 'max_val': 36.6, 'nan': 5.},
    'chbio_6': {'min_val': -28.3, 'max_val': 5.4, 'nan': -29.},
    'chbio_7': {'min_val': 16.7, 'max_val': 42., 'nan': 16.},
    'chbio_8': {'min_val': -14.2, 'max_val': 23.0, 'nan': -15.},
    'chbio_9': {'min_val': -17.7, 'max_val': 26.5, 'nan': -19.},
    'chbio_10': {'min_val': -2.8, 'max_val': 26.5, 'nan': -4.},
    'chbio_11': {'min_val': -17.7, 'max_val': 11.8, 'nan': -19.},
    'chbio_12': {'min_val': 318.3, 'max_val': 2543.3, 'nan': 317.},
    'chbio_13': {'min_val': 43., 'max_val': 285.5, 'nan': 42.},
    'chbio_14': {'min_val': 3.0, 'max_val': 135.6, 'nan': 2.},
    'chbio_15': {'min_val': 8.2, 'max_val': 26.5, 'nan': 7.},
    'chbio_16': {'min_val': 121.6, 'max_val': 855.6, 'nan': 120.},
    'chbio_17': {'min_val': 19.8, 'max_val': 421.3, 'nan': 19.},
    'chbio_18': {'min_val': 19.8, 'max_val': 851.7, 'nan': 19.},
    'chbio_19': {'min_val': 60.5, 'max_val': 520.4, 'nan': 60.},
    'clc': {'attrib_column': 'CLC_CODE', 'one_hot': True, 'nan': 0},
    'crusting': {'nan': -1.},
    'dgh': {'nan': 10},
    'dimp': {'nan': 50},
    'erodi': {'nan': -1},
    'etp': {'nan': 132},
    'oc_top': {'nan': 0},
    'pd_top': {'nan': 0},
    'proxi_eau_fast': {'attrib_column': None, 'nan': -1},
    'text': {'nan': -1},
}


class Raster(object):
    """
    Raster is dedicated to a single raster management...
    """
    def __init__(self, path, nan=None, normalized=False, transform=None, size=64, one_hot=False,
                 attrib_column='QUANTI', max_val=255, min_val=0):
        """
        Loads a tiff file describing an environmental raster into a numpy array and...

        :param path: the path of the raster (the directory)
        :param nan: the value to use when NaN number are present. If False, then default values will be used
        :param normalized: if True the raster will be normalized (minus the mean and divided by std)
        :param transform: if a function is given, it will be applied on each patch.
        :param size: the size of a patch (size x size)
        :param one_hot: if True, each patch will have a one hot encoding representation given the values in the raster.
        :param attrib_column: the name of the column that contains the correct value to use in the raster
        :param max_val: the maximum value within the raster (used to reconstruct correct values in the raster)
        :param min_val: the minimum value within the raster (used to reconstruct correct values in the raster)
        """
        self.path = path
        self.no_data = nan
        self.normalized = normalized
        self.transform = transform
        self.size = size
        self.one_hot = one_hot

        path = re.sub(r'/\/+/', '/', path)

        self.name = path.split('/')[-1] if path[-1] != '/' else path.split('/')[-2]

        # src.meta
        # to avoid the annoying corresponding warning, temporary warning disabling...
        warnings.filterwarnings("ignore")
        src = rasterio.open(path + '/' + self.name + '.tif', nodata=nan)
        warnings.filterwarnings("default")

        if src.meta['crs'] is None:

            with open(path + '/' + 'GeoMetaData.csv') as f:
                metadata = f.read()

            m_split = metadata.split('\n')[1].split(';')

            # loading file data
            self.x_min = float(m_split[1])
            self.y_min = float(m_split[2])
            self.x_resolution = float(m_split[5])
            self.y_resolution = float(m_split[6])
            self.n_rows = int(m_split[3])
            self.n_cols = int(m_split[4])
        else:
            self.x_min = src.bounds.left
            self.y_min = src.bounds.bottom
            self.x_resolution = src.res[0]
            self.y_resolution = src.res[1]
            self.n_rows = src.height
            self.n_cols = src.width

        # some tiff do not contain geo data (stored in the file GeoMetaData)
        # loading the raster
        self.raster = np.squeeze(src.read())
        src.close()

        # raster type is float
        self.raster = self.raster.astype(np.float)

        # value bellow min_value are considered incorrect and therefore no_data
        self.raster[self.raster < MIN_ALLOWED_VALUE] = nan
        self.raster[np.isnan(self.raster)] = nan
        # if the file exists, then it contains the correct values
        if os.path.isfile(path + '/attrib_' + self.name + '.csv'):
            if attrib_column is not None:

                df = pd.read_csv(path + '/attrib_' + self.name + '.csv', header='infer', sep=";")

                for line in df.iterrows():
                    if np.isnan(line[1][attrib_column]):
                        self.raster[self.raster == line[1]['storage_8bit']] = nan
                    else:
                        self.raster[self.raster == line[1]['storage_8bit']] = line[1][attrib_column]
        # if the file does not exist, then correct values must be reconstructed....
        else:
            self.raster = min_val + (max_val - min_val)*((self.raster/255) - 0.1) / 0.8

        if normalized:
            # normalizing the whole raster given available data (therefore avoiding no_data)...
            selected_cell = self.raster != nan
            self.raster[selected_cell] = (self.raster[selected_cell] - self.raster[selected_cell].mean()) \
                / self.raster[selected_cell].std()  # TODO all raster with nan or without nan

        if self.one_hot:
            # unique values for 1 hot encoding
            self.unique_values = np.unique(self.raster[self.raster != nan])

        # setting the shape of the raster
        self.shape = self.raster.shape

    def _get_patch(self, item, cancel_one_hot=False):
        """
        Avoid using this method directly

        :param item: the GPS position (latitude, longitude)
        :param cancel_one_hot: if True, one hot encoding will not be used
        :return: a patch
        """
        row_num = int(self.n_rows - (item[0] - self.y_min) / self.y_resolution)
        col_num = int((item[1] - self.x_min) / self.x_resolution)

        # environmental vector
        if self.size == 1:
            patch = self.raster[row_num, col_num]
            if self.one_hot and not cancel_one_hot:
                patch = np.array([(patch == i).astype(float) for i in self.unique_values])
            else:
                patch = patch[np.newaxis]
        # environmental tensor
        else:
            half_size = int(self.size/2)
            patch = self.raster[row_num-half_size:row_num+half_size, col_num - half_size:col_num+half_size]
            if self.one_hot and not cancel_one_hot:
                patch = np.array([(patch == i).astype(float) for i in self.unique_values])
            else:
                patch = patch[np.newaxis]

        return patch

    def __len__(self):
        """
        :return: the depth of the tensor/vector...
        """
        if self.one_hot:
            return int(self.unique_values.shape[0])
        else:
            return 1

    def __getitem__(self, item, cancel_one_hot=False):
        """
        The method to use to retrieve a patch.

        :param item: GPS position (latitude, longitude)
        :param cancel_one_hot: if true the one hot encoding representation will be disabled
        :return: the extracted patch with eventually some transformations
        """
        # item is a tuple of (latitude, longitude)
        patch = self._get_patch(item, cancel_one_hot)
        if self.transform:
            patch = self.transform(patch)

        return patch


class PatchExtractor(object):
    """
    PatchExtractor enables the extraction of an environmental tensor from multiple rasters given a GPS
    position.
    """
    def __init__(self, root_path, size=64, verbose=False):
        self.root_path = root_path
        self.size = size

        self.verbose = verbose

        self.rasters = []

    def add_all(self, normalized=False, transform=None):
        """
        Add all variables (rasters) available at root_path

        :param normalized: if True, each raster will be normalized
        :param transform: a function to apply on each patch
        """
        for key in sorted(raster_metadata.keys()):
            if 'ignore' not in raster_metadata[key]:
                self.append(key, normalized=normalized, transform=transform)

    def append(self, raster_name, **kwargs):
        """
        This method append a new raster given its name

        :param raster_name:
        :param kwargs: nan, normalized, transform
        """
        # you may want to add rasters one by one if specific configuration are required on a per raster
        # basis
        if self.verbose:
            print('Adding: '+raster_name)
        params = {**raster_metadata[raster_name]}
        for k in kwargs.keys():
            if kwargs[k] != 'default':
                params[k] = kwargs[k]

        self.rasters.append(Raster(self.root_path + '/' + raster_name, size=self.size, **params))

    def clean(self):
        """
        Remove all rasters from the extractor.
        """
        if self.verbose:
            print('Removing all rasters...')
        self.rasters = []

    def __getitem__(self, item, cancel_one_hot=False):
        """
        :param item: the GPS location (latitude, longitude)
        :return: return the environmental tensor or vector (size>1 or size=1)
        """
        return np.concatenate([r.__getitem__(item, cancel_one_hot) for r in self.rasters])

    def __len__(self):
        """
        :return: the number of variables (not the size of the tensor when some variables have a one hot encoding
                 representation)
        """
        return len(self.rasters)

    def plot(self, item, cancel_one_hot=True, return_fig=False, style='fivethirtyeight'):
        """
        Plot an environmental tensor (size > 1)...

        :param item: the GPS location (latitude, longitude)
        :param cancel_one_hot: if False, the variables that have to be display with a one hot encoding approach will
                               be displayed as such. If True, all variables will have only one dimension.
        :param return_fig: if True, the matplotlib fig will be returned, if False, it will be displayed
        :param style: style of the chart
        """
        if self.size > 1:
            with plt.style.context(style):
                metadata = [(r.name,
                             [item[1] - self.size // 2 * r.x_resolution,
                              item[1] + self.size // 2 * r.x_resolution,
                              item[0] - self.size // 2 * r.y_resolution,
                              item[0] + self.size // 2 * r.y_resolution]
                             ) for r in self.rasters for _ in range(1 if cancel_one_hot else len(r))]
                # metadata are the name of the variable and the bounding box in latitude-longitude coordinates

                # retrieve the patch... Eventually disabling the one hot encoding variables
                patch = self.__getitem__(item, cancel_one_hot)

                # computing number of rows and columns...
                nb_rows = (patch.shape[0] + 4) // 5
                nb_cols = 5

                fig = plt.figure(figsize=(nb_cols * 5, nb_rows * 3.8))
                for k, i in zip(metadata, range(patch.shape[0])):
                    plt.subplot(nb_rows, nb_cols, i + 1)
                    plt.title(k[0], fontsize=20)
                    plt.imshow(patch[i], extent=k[1], aspect='auto')
                    plt.colorbar()
                fig.tight_layout()
            if return_fig:
                return fig
            else:
                fig.show()
                plt.close(fig)
        else:
            raise ValueError('Plot works only for tensors: size must be > 1...')
