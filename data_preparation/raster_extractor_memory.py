import rasterio
import pyproj
import numpy as np
from rasterio.windows import Window
import math
import matplotlib.pyplot as plt
from random import random
from rasterio.windows import from_bounds
from PIL import Image


class Extractor(object):
    def __init__(self, raster_path, grid_cell, cell_size, patch_size, margin=1000, gamma=2.5, convert_uint8=True, fill_zero_if_error=False):
        # create a pyproj transformer object to convert lat, lon to EPSG:32738
        self.patch_size = patch_size
        self.fill_zero_if_error = fill_zero_if_error
        self.raster_path = raster_path
        self.grid_cell = grid_cell
        self.gamma = gamma

        # open the tif file with rasterio
        with rasterio.open(raster_path) as src:
            # read the metadata of the file
            meta = src.meta
            meta.update(count=src.count) # update the count of the meta to match the number of layers

            self.transformer = pyproj.Transformer.from_crs("epsg:4326", src.crs, always_xy=True)

            # get the NoData value from the raster
            self.nodata_value = src.nodatavals
            self.x_resolution = src.res[0]
            self.y_resolution = src.res[1]

            x_min, x_max = int((grid_cell[0] * cell_size) - margin), int(((grid_cell[0] + 1) * cell_size) + margin)
            y_min, y_max = int((grid_cell[1] * cell_size) - margin), int(((grid_cell[1] + 1) * cell_size) + margin)

            self.x_min = x_min
            self.y_min = y_min

            window  = window = from_bounds(y_min, x_min, y_max, x_max, src.transform)

            #print(window)

            # self.data = src.read(1, window=window, masked=False, out_dtype='uint16', boundless=True, fill_value=10000)
            self.data = src.read(1, window=window, masked=False, out_dtype='uint16')
            
            #self.data = np.where(self.data == self.nodata_value, 0, self.data)

            if convert_uint8:
                self.data = np.clip(self.data/10000.0, a_min=0, a_max=1.0)
                self.data = (self.data**(1/self.gamma))*256
                self.data = self.data.astype(np.uint8)

            #print(self.data.shape)

            #image = Image.fromarray(self.data)
            #image = image.convert("L")
            #image.save("./test_tile_"+str(grid_cell[0])+"_"+str(grid_cell[1])+".jpeg", "JPEG", quality=90)
        src.close()


    def __getitem__(self, item):
        """
        :param item: the GPS location (latitude, longitude)
        :return: return the environmental tensor or vector (size>1 or size=1)
        """
        
        # convert the lat, lon coordinates to EPSG:32738
        lon, lat = self.transformer.transform(item[1], item[0])

        # calculate the x, y coordinate of the point of interest
        x = int(self.data.shape[0] - (lat - self.x_min) / self.x_resolution)
        y = int((lon - self.y_min) / self.y_resolution)

        # read the data of the patch from all layers
        if self.patch_size == 1:
            patch_data = self.data[x, y]
        else:
            patch_data = self.data[x - (self.patch_size // 2): x + (self.patch_size // 2), y - (self.patch_size // 2): y + (self.patch_size // 2)]
        
        tensor = patch_data
        if self.patch_size != 1 and self.fill_zero_if_error and tensor.shape != (self.patch_size, self.patch_size):
            tensor = np.zeros((self.patch_size, self.patch_size))
        if self.patch_size != 1 and tensor.shape != (self.patch_size, self.patch_size):
            print('error', tensor.shape, self.raster_path, self.grid_cell)
            tensor = np.nan
        return tensor
    
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        result = '-' * 50 + '\n'
        result += 'grid_cell: ' + str(self.grid_cell) + '\n'
        result += 'x_min: ' + str(self.x_min) + '\n'
        result += 'y_min: ' + str(self.y_min) + '\n'
        result += 'x_resolution: ' + str(self.x_resolution) + '\n'
        result += 'y_resolution: ' + str(self.y_resolution) + '\n'
        result += '-' * 50

        return result

    def __len__(self):
        """
        :return: the number of variables (not the size of the tensor when some variables have a one hot encoding
                 representation)
        """
        return 1

    
    def plot_patch(self, item):
        patch = self[item]

        # calculate the number of rows and columns for the subplots grid
        rows = int(math.ceil(math.sqrt(self.nb_layers)))
        cols = int(math.ceil(self.nb_layers / rows))

        # create a figure with a grid of subplots
        fig, axs = plt.subplots(rows, cols, figsize=(20, 20))

        # flatten the subplots array to easily access the subplots
        axs = axs.flatten()

        # loop through the layers of patch data
        for i in range(self.nb_layers):
            # display the layer on the corresponding subplot
            axs[i].imshow(patch[i])
            axs[i].set_title(f'Layer {i+1}')
            axs[i].axis('off')

        # remove empty subplots
        for i in range(self.nb_layers, rows*cols):
            fig.delaxes(axs[i])

        # show the plot
        plt.show()


if __name__ == "__main__":
    p = Extractor('data/data_raw/current_chelsa.tif')
    p.plot_patch((-13.331141, 48.267285))


