import random

from data.GLC23PatchesProviders import MultipleRasterPatchProvider, RasterPatchProvider
from data.GLC23Datasets import PatchesDataset

data_path = './data/' # root path of the data

# configure providers
p_soil = MultipleRasterPatchProvider(data_path+'EnvironmentalRasters/Soilgrids/') # take all soilgris rasters
p_bioclim = MultipleRasterPatchProvider(data_path+'EnvironmentalRasters/Climate/BioClimatic_Average_1981-2010/', select=['bio1', 'bio7']) # take only bio1 and bio7 from bioclimatic rasters
p_hfp = RasterPatchProvider(data_path+'EnvironmentalRasters/HumanFootprint/summarized/HFP2009_WGS84.tif') # take the human footprint 2009 summurized raster

# create dataset
dataset = PatchesDataset(occurrences=data_path+'Presence_only_occurrences/PO_anonymised_filtered.csv', providers=(p_soil, p_bioclim, p_hfp))


# print random tensors from dataset
ids = [random.randint(0, len(dataset)) for i in range(10)]
for id in ids:
    tensor = dataset[id][0]
    label = dataset[id][1]
    print('Tensor type: {}, tensor shape: {}, label: {}'.format(type(tensor), tensor.shape, label))
    dataset.plot_patch(id)