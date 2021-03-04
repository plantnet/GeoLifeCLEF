import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import rasterio


raster_names = [
    "bio_1", "bio_2", "bio_3", "bio_4", "bio_5", "bio_6", "bio_7", "bio_8", "bio_9",
    "bio_10", "bio_11", "bio_12", "bio_13", "bio_14", "bio_15", "bio_16", "bio_17",
    "bio_18", "bio_19", "bdticm", "bldfie", "cecsol", "clyppt", "orcdrc", "phihox",
    "sltppt", "sndppt"
]


class Raster(object):
    """
    Handles the loading and patch extraction of a single raster
    """

    def __init__(self, path, country, size=256, nan=None, out_of_bounds="error"):
        """
        Loads a GeoTIFF file describing an environmental raster

        :param path: the path of the raster (the directory)
        :param nan: the value to use when NaN number are present. If False, then default values will be used
        :param size: the size of a patch (size x size)
        """
        path = Path(path)
        if not path.exists():
            raise ValueError(
                "path should be the path to a raster, given non-existant path: {}".format(
                    path
                )
            )

        self.path = path
        self.name = path.name
        self.size = size
        self.out_of_bounds = out_of_bounds
        self.nan = nan

        # Loading the raster
        filename = path / "{}_{}.tif".format(self.name, country)
        with rasterio.open(filename) as dataset:
            self.dataset = dataset
            self.raster = dataset.read(1).astype(np.float32)
            mask = self.dataset.read_masks(1).astype(np.bool)

        # Changes nodata to user specified value
        if nan:
            self.raster[np.isnan(self.raster)] = nan
            self.raster[~mask] = nan

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
            patch = self.raster[row, col]
        else:
            # FIXME can it happen that part of the patch is outside the raster? (and how about the mask of the dataset?)
            half_size = int(self.size / 2)
            # FIXME only way to trigger an exception? (slices don't)
            self.raster[row, col]
            patch = self.raster[
                row - half_size:row + half_size,
                col - half_size:col + half_size
            ]

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
        :return: the extracted patch
        """
        try:
            return self._get_patch(item)
        except IndexError as e:
            if self.out_of_bounds == "error":
                raise e
            else:
                if self.out_of_bounds == "warn":
                    warnings.warn("GPS position ({}, {}) out of bounds".format(*item))

                if self.size == 1:
                    return np.array([self.nan], dtype=np.float32)
                else:
                    patch = np.empty((1, self.size, self.size), dtype=np.float32)
                    patch.fill(self.nan)
                    return patch

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "name: " + self.name + "\n"


class PatchExtractor(object):
    """
    PatchExtractor handles the extraction of an environmental tensor from multiple rasters given a GPS
    position.
    """

    def __init__(self, root_path, size=256):
        self.root_path = Path(root_path)
        if not self.root_path.exists():
            raise ValueError(
                "root_path should be the directory containing the rasters, given a non-existant path: {}".format(
                    root_path
                )
            )

        self.size = size

        self.rasters_fr = []
        self.rasters_us = []

    def add_all(self, **kwargs):
        """
        Add all variables (rasters) available at root_path

        :param kwargs: updates the default arguments passed to Raster (nan, out_of_bounds, etc.)
        """
        for raster_name in raster_names:
            self.append(raster_name, **kwargs)

    def append(self, raster_name, **kwargs):
        """
        This method appends a new raster given its name.
        Can be useful to load only a subset of rasters or to pass configurations specific to each raster.

        :param raster_name:
        :param kwargs: updates the default arguments passed to Raster (nan, out_of_bounds, etc.)
        """
        r_us = Raster(self.root_path / raster_name, "USA", size=self.size, **kwargs)
        r_fr = Raster(self.root_path / raster_name, "FR", size=self.size, **kwargs)

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
        result = ""

        for rasters in [self.rasters_fr, self.rasters_us]:
            for raster in rasters:
                result += "-" * 50 + "\n"
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
        if item[1] > -10.0:
            return self.rasters_fr
        else:
            return self.rasters_us

    def plot(self, item, return_fig=False, n_cols=5, alpha=1.0, fig=None, resolution=1.0):
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
            (
                raster.name,
                [
                    item[1] - (self.size // 2) * raster.dataset.res[0],
                    item[1] + (self.size // 2) * raster.dataset.res[0],
                    item[0] - (self.size // 2) * raster.dataset.res[1],
                    item[0] + (self.size // 2) * raster.dataset.res[1],
                ],
            )
            for raster in rasters
        ]

        # retrieve the patch... Eventually disabling the one hot encoding variables
        patch = self[item]

        # computing number of rows and columns
        n_rows = (patch.shape[0] + (n_cols - 1)) // n_cols

        if fig is None:
            fig = plt.figure(figsize=(n_cols * 6.4 * resolution, n_rows * 4.8 * resolution))

        axes = fig.subplots(n_rows, n_cols)
        axes = axes.ravel()

        for i, (ax, k) in enumerate(zip(axes, metadata)):
            p = np.squeeze(patch[i])
            im = ax.imshow(p, extent=k[1], aspect="equal", interpolation="none")

            ax.set_title(k[0], fontsize=20)
            fig.colorbar(im, ax=ax)

        # for ax in axes:
        for ax in axes[len(metadata):]:
            ax.axis("off")

        fig.tight_layout()
        fig.patch.set_alpha(alpha)

        if return_fig:
            return fig
