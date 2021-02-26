"""
"
"   Author: Maximilien Servajean - mservajean
"   Mail: servajean@lirmm.fr
"   Date: 04/01/2019
"
"   Description: The code to extract environmental tensors and environmental vectors given some environmental rasters.
"
"""
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import rasterio


# metadata used to setup some rasters
raster_metadata = {
    'bio_1': {'nan': -117, 'mu': 101, 'sigma': 58},
    'bio_2': {'nan': -54, 'mu': 131, 'sigma': 28},
    'bio_3': {'nan': 18, 'mu': 36, 'sigma': 8},
    'bio_4': {'nan': 1623, 'mu': 8267, 'sigma': 2152},
    'bio_5': {'nan': -26, 'mu': 289, 'sigma': 48},
    'bio_6': {'nan': -277, 'mu': -78, 'sigma': 83},
    'bio_7': {'nan': 116, 'mu': 367, 'sigma': 72},
    'bio_8': {'nan': -170, 'mu': 149, 'sigma': 82},
    'bio_9': {'nan': -182, 'mu': 54, 'sigma': 114},
    'bio_10': {'nan': -54, 'mu': 205, 'sigma': 47},
    'bio_11': {'nan': -187, 'mu': -7, 'sigma': 80},
    'bio_12': {'nan': -36, 'mu': 746, 'sigma': 383},
    'bio_13': {'nan': 6, 'mu': 98, 'sigma': 47},
    'bio_14': {'nan': -1, 'mu': 34, 'sigma': 26},
    'bio_15': {'nan': 4, 'mu': 38, 'sigma': 23},
    'bio_16': {'nan': 18, 'mu': 265, 'sigma': 132},
    'bio_17': {'nan': -1, 'mu': 117, 'sigma': 84},
    'bio_18': {'nan': 0, 'mu': 213, 'sigma': 107},
    'bio_19': {'nan': 4, 'mu': 163, 'sigma': 137},
    'bdticm': {'nan': -1, 'mu': 2579, 'sigma': 3058},
    'bldfie': {'nan': 92, 'mu': 1372, 'sigma': 137},
    'cecsol': {'nan': -1, 'mu': 20, 'sigma': 8},
    'clyppt': {'nan': -1, 'mu': 22, 'sigma': 8},
    'orcdrc': {'nan': -1, 'mu': 24, 'sigma': 21},
    'phihox': {'nan': 31, 'mu': 64, 'sigma': 11},
    'sltppt': {'nan': -1, 'mu': 37, 'sigma': 11},
    'sndppt': {'nan': -1, 'mu': 42, 'sigma': 14},
}


class Raster(object):
    """
    Raster is dedicated to a single raster management...
    """
    def __init__(self, path, country='FR', normalized=False, transform=None, size=256, nan=None, mu=0,
                 sigma=1, out_of_bounds="error", **kw):
        """
        Loads a tiff file describing an environmental raster into a numpy array and...

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
        self.normalized = normalized
        self.transform = transform
        self.size = size
        self.out_of_bounds = out_of_bounds
        self.nan = nan

        # Loading the raster
        filename = path / (self.name + '_' + country + '.tif')
        with rasterio.open(filename) as dataset:
            self.dataset = dataset
            self.raster = np.squeeze(dataset.read())
            nodata = self.dataset.nodata

        # Changes nodata to user specified value
        if nan:
            self.raster[self.raster == nodata] = nan
            self.raster[np.isnan(self.raster)] = nan

        # if asked, normalize the whole raster where data is available
        if normalized:
            non_nan = self.raster != nan
            self.raster[non_nan] = (self.raster[non_nan] - mu) / sigma

        # setting the shape of the raster
        self.shape = self.raster.shape

    def _get_patch(self, item):
        """
        Avoid using this method directly

        :param item: GPS position as tuple (latitude, longitude)
        :return: a patch
        """
        row, col = self.dataset.index(item[1], item[0])

        # environmental vector
        if self.size == 1:
            patch = self.raster[row, col].astype(np.float)
        else:
            # FIXME can it happen that part of the patch is outside the raster? (and how about the mask of the dataset?)
            half_size = int(self.size / 2)
            # FIXME only way to trigger an exception? (slices don't)
            self.raster[row, col]
            patch = self.raster[
                row - half_size:row + half_size,
                col - half_size:col + half_size
            ].astype(np.float)

        patch = patch[np.newaxis]
        return patch

    def __len__(self):
        """
        :return: the depth of the tensor/vector
        """
        return self.dataset.count

    def __getitem__(self, item):
        """
        The method to use to retrieve a patch.

        :param item: GPS position as a tuple (latitude, longitude)
        :return: the extracted patch with eventually some transformations
        """
        try:
            patch = self._get_patch(item)

            if self.transform:
                patch = self.transform(patch)

            return patch
        except IndexError as e:
            if self.out_of_bounds == "error":
                raise e
            else:
                if self.out_of_bounds == "warn":
                    warnings.warn("GPS position ({}, {}) out of bounds".format(*item))

                if self.size == 1:
                    return np.array([self.nan], dtype=np.float)
                else:
                    patch = np.empty((1, self.size, self.size), dtype=np.float)
                    patch.fill(self.nan)
                    return patch

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'title: ' + self.name + '\n'


class PatchExtractor(object):
    """
    PatchExtractor handles the extraction of an environmental tensor from multiple rasters given a GPS
    position.
    """
    def __init__(self, root_path, size=256, resolution=1., out_of_bounds="error"):
        self.root_path = Path(root_path)
        if not self.root_path.exists():
            raise ValueError(
                "root_path should be the directory containing the rasters, given a non-existant path: {}".format(root_path)
            )

        self.size = size

        self.resolution = resolution
        self.out_of_bounds = out_of_bounds

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
        params["out_of_bounds"] = self.out_of_bounds
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
        return str(self)

    def __str__(self):
        result = ''

        for rasters in [self.rasters_fr, self.rasters_us]:
            for raster in rasters:
                result += '-' * 50 + '\n'
                result += str(raster)

        return result

    def __getitem__(self, item):
        """
        :param item: the GPS location (latitude, longitude)
        :return: return the environmental tensor or vector (size>1 or size=1)
        """
        rasters = self._raster(item)
        return np.concatenate([r[item] for r in rasters])

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
            (raster.name, [
                item[1] - (self.size // 2) * raster.dataset.res[0],
                item[1] + (self.size // 2) * raster.dataset.res[0],
                item[0] - (self.size // 2) * raster.dataset.res[1],
                item[0] + (self.size // 2) * raster.dataset.res[1]
            ]) for raster in rasters
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
