import random

from data.GLC23PatchesProviders import MultipleRasterPatchProvider, RasterPatchProvider, JpegPatchProvider
from data.GLC23Datasets import PatchesDataset, PatchesDatasetMultiLabel

data_path = 'data/sample_data/' # root path of the data

# configure providers
p_rgb = JpegPatchProvider(data_path+'SatelliteImages/', dataset_stats='jpeg_patches_sample_stats.csv') # take all sentinel imagery layer (r,g,b,nir)
p_hfp_d = MultipleRasterPatchProvider(data_path+'EnvironmentalRasters/HumanFootprint/detailed/') # take all rasters from human footprint detailed
p_bioclim = MultipleRasterPatchProvider(data_path+'EnvironmentalRasters/Climate/BioClimatic_Average_1981-2010/', select=['bio1', 'bio2']) # take only bio1 and bio2 from bioclimatic rasters
p_hfp_s = RasterPatchProvider(data_path+'EnvironmentalRasters/HumanFootprint/summarized/HFP2009_WGS84.tif') # take the human footprint 2009 summurized raster

# create dataset
dataset = PatchesDataset(occurrences=data_path+'Presence_only_occurrences/Presences_only_train_sample.csv', providers=[p_hfp_d, p_bioclim, p_hfp_s, p_rgb])
dataset_multi = PatchesDatasetMultiLabel(occurrences=data_path+'Presence_only_occurrences/Presences_only_train_sample.csv', providers=[p_hfp_d, p_bioclim, p_hfp_s, p_rgb])

# print random tensors from dataset
# ids = [random.randint(0, len(dataset)-1) for i in range(1)]
ids = [0]
for id in ids:
    tensor, label = dataset[id]
    label_multi = dataset_multi[id][1]
    print('Tensor type: {}, tensor shape: {}, label: {}, \nlabel_multi: {}'.format(type(tensor), tensor.shape, label, label_multi))
    dataset.plot_patch(id)