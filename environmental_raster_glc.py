"""
"
"   Author: Maximilien Servajean - mservajean
"   Mail: servajean@lirmm.fr
"   Date: 04/01/2019
"
"   Description: The code to extract environmental tensors and environmental vectors given some environmental rasters.
"
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import rasterio


# metadata used to setup some rasters
raster_metadata = {
    'bio_1': {'min_val': -116, 'max_val': 259, 'nan': -2147483647, 'new_nan': -117, 'mu': 101, 'sigma': 58},
    'bio_2': {'min_val': -53, 'max_val': 361, 'nan': -2147483647, 'new_nan': -54, 'mu': 131, 'sigma': 28},
    'bio_3': {'min_val': 19, 'max_val': 69, 'nan': -2147483647, 'new_nan': 18, 'mu': 36, 'sigma': 8},
    'bio_4': {'min_val': 1624, 'max_val': 13302, 'nan': -2147483647, 'new_nan': 1623, 'mu': 8267, 'sigma': 2152},
    'bio_5': {'min_val': -25, 'max_val': 457, 'nan': -2147483647, 'new_nan': -26, 'mu': 289, 'sigma': 48},
    'bio_6': {'min_val': -276, 'max_val': 183, 'nan': -2147483647, 'new_nan': -277, 'mu': -78, 'sigma': 83},
    'bio_7': {'min_val': 117, 'max_val': 515, 'nan': -2147483647, 'new_nan': 116, 'mu': 367, 'sigma': 72},
    'bio_8': {'min_val': -169, 'max_val': 332, 'nan': -2147483647, 'new_nan': -170, 'mu': 149, 'sigma': 82},
    'bio_9': {'min_val': -181, 'max_val': 331, 'nan': -2147483647, 'new_nan': -182, 'mu': 54, 'sigma': 114},
    'bio_10': {'min_val': -53, 'max_val': 361, 'nan': -2147483647, 'new_nan': -54, 'mu': 205, 'sigma': 47},
    'bio_11': {'min_val': -186, 'max_val': 220, 'nan': -2147483647, 'new_nan': -187, 'mu': -7, 'sigma': 80},
    'bio_12': {'min_val': -35, 'max_val': 3385, 'nan': -2147483647, 'new_nan': -36, 'mu': 746, 'sigma': 383},
    'bio_13': {'min_val': 7, 'max_val': 570, 'nan': -2147483647, 'new_nan': 6, 'mu': 98, 'sigma': 47},
    'bio_14': {'min_val': 0, 'max_val': 184, 'nan': -2147483647, 'new_nan': -1, 'mu': 34, 'sigma': 26},
    'bio_15': {'min_val': 5, 'max_val': 140, 'nan': -2147483647, 'new_nan': 4, 'mu': 38, 'sigma': 23},
    'bio_16': {'min_val': 19, 'max_val': 1546, 'nan': -2147483647, 'new_nan': 18, 'mu': 265, 'sigma': 132},
    'bio_17': {'min_val': 0, 'max_val': 612, 'nan': -2147483647, 'new_nan': -1, 'mu': 117, 'sigma': 84},
    'bio_18': {'min_val': 1, 'max_val': 777, 'nan': -2147483647, 'new_nan': 0, 'mu': 213, 'sigma': 107},
    'bio_19': {'min_val': 5, 'max_val': 1485, 'nan': -2147483647, 'new_nan': 4, 'mu': 163, 'sigma': 137},
    'bdticm': {'min_val': 0, 'max_val': 112467, 'nan': -2147483647, 'new_nan': -1, 'mu': 2579, 'sigma': 3058},
    'bldfie': {'min_val': 93, 'max_val': 1828, 'nan': -32768, 'new_nan': 92, 'mu': 1372, 'sigma': 137},
    'cecsol': {'min_val': 0, 'max_val': 385, 'nan': -32768, 'new_nan': -1, 'mu': 20, 'sigma': 8},
    'clyppt': {'min_val': 0, 'max_val': 81, 'nan': -32768, 'new_nan': -1, 'mu': 22, 'sigma': 8},
    'orcdrc': {'min_val': 0, 'max_val': 524, 'nan': -32768, 'new_nan': -1, 'mu': 24, 'sigma': 21},
    'phihox': {'min_val': 32, 'max_val': 98, 'nan': -32768, 'new_nan': 31, 'mu': 64, 'sigma': 11},
    'sltppt': {'min_val': 0, 'max_val': 86, 'nan': -32768, 'new_nan': -1, 'mu': 37, 'sigma': 11},
    'sndppt': {'min_val': 0, 'max_val': 99, 'nan': -32768, 'new_nan': -1, 'mu': 42, 'sigma': 14},
}


class Raster(object):
    """
    Raster is dedicated to a single raster management...
    """
    def __init__(self, path, country='FR', normalized=False, transform=None, size=256, nan=None, new_nan=None, mu=0,
                 sigma=1, **kw):
        """
        Loads a tiff file describing an environmental raster into a numpy array and...

        :type new_nan:
        :param path: the path of the raster (the directory)
        :param nan: the value to use when NaN number are present. If False, then default values will be used
        :param normalized: if True the raster will be normalized (minus the mean and divided by std)
        :param transform: if a function is given, it will be applied on each patch.
        :param size: the size of a patch (size x size)
        """
        path = Path(path)
        if not path.exists():
            raise ValueError("path should be the path to a raster, given non-existant path: {}".format(path))

        self.path = path
        self.name = path.name
        self.no_data = new_nan
        self.normalized = normalized
        self.transform = transform
        self.size = size

        # FIXME nodata argument useless here?
        # Loading the raster
        filename = path / (self.name + '_' + country + '.tif')
        with rasterio.open(filename, nodata=nan) as src:
            self.x_min = src.bounds.left
            self.y_min = src.bounds.bottom
            self.x_resolution = src.res[0]
            self.y_resolution = src.res[1]
            self.n_rows = src.height
            self.n_cols = src.width
            self.raster = np.squeeze(src.read())

        # FIXME where is this test? (min_value = x_min or min_val?)
        # value bellow min_value are considered incorrect and therefore no_data
        self.raster[self.raster == nan] = new_nan
        self.raster[np.isnan(self.raster)] = new_nan

        # if asked, normalize the whole raster where data is available
        if normalized:
            # FIXME should be new_nan?
            selected_cell = self.raster != nan
            self.raster[selected_cell] = (self.raster[selected_cell] - mu) / sigma

        # setting the shape of the raster
        self.shape = self.raster.shape

    def _get_patch(self, item):
        """
        Avoid using this method directly

        :param item: GPS position as tuple (latitude, longitude)
        :return: a patch
        """
        row_num = int(self.n_rows - (item[0] - self.y_min) / self.y_resolution)
        col_num = int((item[1] - self.x_min) / self.x_resolution)

        # environmental vector
        if self.size == 1:
            patch = self.raster[row_num, col_num].astype(np.float)
        else:
            half_size = int(self.size/2)
            patch = self.raster[
                row_num-half_size:row_num+half_size,
                col_num - half_size:col_num+half_size
            ].astype(np.float)

        patch = patch[np.newaxis]
        return patch

    def __len__(self):
        """
        :return: the depth of the tensor/vector...
        """
        return 1

    def __getitem__(self, item):
        """
        The method to use to retrieve a patch.

        :param item: GPS position as a tuple (latitude, longitude)
        :return: the extracted patch with eventually some transformations
        """
        patch = self._get_patch(item)

        if self.transform:
            patch = self.transform(patch)

        return patch


class PatchExtractor(object):
    """
    PatchExtractor handles the extraction of an environmental tensor from multiple rasters given a GPS
    position.
    """
    def __init__(self, root_path, size=256, verbose=False, resolution=1.):
        self.root_path = Path(root_path)
        if not self.root_path.exists():
            raise ValueError(
                "root_path should be the directory containing the rasters, given a non-existant path: {}".format(root_path)
            )

        self.size = size

        self.verbose = verbose
        self.resolution = resolution

        self.rasters_fr = []
        self.rasters_us = []

    def add_all(self, **kwargs):
        """
        Add all variables (rasters) available at root_path

        :param kwargs: updates the default arguments passed to Raster (nan, normalized, transform, etc.)
        """
        for key in raster_metadata.keys():
            if "ignore" not in raster_metadata[key]:
                self.append(key, **kwargs)

    def append(self, raster_name, **kwargs):
        """
        This method appends a new raster given its name.
        Can be useful to load only a subset of rasters or to pass configurations specific to each raster.

        :param raster_name:
        :param kwargs: updates the default arguments passed to Raster (nan, normalized, transform, etc.)
        """
        params = dict(raster_metadata[raster_name])
        params.update(kwargs)

        r_us = Raster(self.root_path / raster_name, 'USA', size=self.size, **params)
        r_fr = Raster(self.root_path / raster_name, 'FR', size=self.size, **params)

        self.rasters_us.append(r_us)
        self.rasters_fr.append(r_fr)

    def clean(self):
        """
        Remove all rasters from the extractor.
        """
        self.rasters_fr = []
        self.rasters_us = []

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        str_ = ''

        def raster_str(r):
            result = ''
            result += '-' * 50 + '\n'
            result += 'title: ' + r.name + '\n'
            result += '\t x_min: ' + str(r.x_min) + '\n'
            result += '\t y_min: ' + str(r.y_min) + '\n'
            result += '\t x_resolution: ' + str(r.x_resolution) + '\n'
            result += '\t y_resolution: ' + str(r.y_resolution) + '\n'
            result += '\t n_rows: ' + str(r.n_rows) + '\n'
            result += '\t n_cols: ' + str(r.n_cols) + '\n'
            return result

        for r in self.rasters_fr:
            str_ += raster_str(r)

        for r in self.rasters_us:
            str_ += raster_str(r)

        return str_

    def __getitem__(self, item):
        """
        :param item: the GPS location (latitude, longitude)
        :return: return the environmental tensor or vector (size>1 or size=1)
        """
        rasters = self._raster(item)

        if len(rasters) > 1:
            return np.concatenate([r[item] for r in rasters])
        else:
            return np.array([rasters[0][item]])

    def __len__(self):
        """
        :return: the number of variables (not the size of the tensor when some variables have a one hot encoding
                 representation)
        """
        return len(self.rasters_fr)

    def _raster(self, item):
        return self.rasters_fr if item[1] > -10. else self.rasters_us

    def plot(self, item, return_fig=False, n_cols=5, alpha=1., fig=None):
        """
        Plot an environmental tensor (size > 1)...

        :param alpha:
        :param n_cols:
        :param item: the GPS location (latitude, longitude)
        :param return_fig: if True, the matplotlib fig will be returned, if False, it will be displayed
        :param style: style of the chart
        """
        if self.size <= 1:
            raise ValueError("Plot works only for tensors: size must be > 1")

        rasters = self._raster(item)

        # metadata are the name of the variable and the bounding box in latitude-longitude coordinates
        metadata = [
            (r.name, [
                item[1] - (self.size // 2) * r.x_resolution,
                item[1] + (self.size // 2) * r.x_resolution,
                item[0] - (self.size // 2) * r.y_resolution,
                item[0] + (self.size // 2) * r.y_resolution
            ]) for r in rasters
        ]

        # retrieve the patch... Eventually disabling the one hot encoding variables
        patch = self[item]

        # computing number of rows and columns
        n_rows = (patch.shape[0] + (n_cols-1)) // n_cols

        if fig is None:
            fig = plt.figure(figsize=(n_cols * 6.4 * self.resolution, n_rows * 4.8 * self.resolution))

        axes = fig.subplots(n_rows, n_cols)
        axes = axes.ravel()

        for i, (ax, k) in enumerate(zip(axes, metadata)):
            p = np.squeeze(patch[i])
            im = ax.imshow(p, extent=k[1], aspect="equal")

            ax.set_title(k[0], fontsize=20)
            fig.colorbar(im, ax=ax)

        # for ax in axes:
        for ax in axes[len(metadata):]:
            ax.axis("off")

        fig.tight_layout()
        fig.patch.set_alpha(alpha)

        if return_fig:
            return fig
