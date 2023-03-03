import random

from data.GLC23TimeSeriesProviders import MultipleCSVTimeSeriesProvider, CSVTimeSeriesProvider
from data.GLC23Datasets import TimeSeriesDataset

data_path = '/home/tlarcher/Documents/Pl@ntNet/Seafile/LIRMM - Malpolon/time_series/' # root path of the data

# configure providers
ts_red = CSVTimeSeriesProvider(data_path+'time_series_red.csv')
ts_multi = MultipleCSVTimeSeriesProvider(data_path, select=['red', 'blue'])
ts_all = MultipleCSVTimeSeriesProvider(data_path)

# create dataset
dataset = TimeSeriesDataset(occurrences=data_path+'../src/PO_anonymised_filtered.csv', providers=(ts_red))


# print random tensors from dataset
ids = [random.randint(0, len(dataset)) for i in range(10)]
for id in ids:
    tensor = dataset[id][0]
    label = dataset[id][1]
    print('Tensor type: {}, tensor shape: {}, label: {}'.format(type(tensor), tensor.shape, label))
    dataset.plot_ts(id)