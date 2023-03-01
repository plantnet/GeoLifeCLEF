import random

from data.GLC23PatchesProviders import MultipleRasterPatchProvider, RasterPatchProvider
from data.GLC23Datasets import PatchesDataset

data_path = './data/' # root path of the data

p_soil = MultipleRasterPatchProvider(data_path+'EnvironmentalRasters/Soilgrids/')
p_bioclim = MultipleRasterPatchProvider(data_path+'EnvironmentalRasters/Climate/BioClimatic_Average_1981-2010/', select=['bio1', 'bio7'])
p_hfp = RasterPatchProvider(data_path+'EnvironmentalRasters/HumanFootprint/summarized/HFP2009_WGS84.tif')

dataset = PatchesDataset(occurrences=data_path+'Presence_only_occurrences/PO_anonymised_filtered.csv', providers=(p_soil, p_bioclim, p_hfp))

ids = [random.randint(0, len(dataset)) for i in range(10)]
for id in ids:
    print(dataset[id])