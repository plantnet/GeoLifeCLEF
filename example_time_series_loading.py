import random

from data.GLC23TimeSeriesProviders import MultipleCSVTimeSeriesProvider, CSVTimeSeriesProvider
from data.GLC23Datasets import TimeSeriesDataset

data_path = '/home/tlarcher/Documents/Pl@ntNet/git/GLC/data/sample_data/' # root path of the data
# configure providers
ts_red = CSVTimeSeriesProvider(data_path+'SatelliteTimeSeries/time_series_red.csv')
ts_multi = MultipleCSVTimeSeriesProvider(data_path+'SatelliteTimeSeries/', select=['red', 'blue'])
ts_all = MultipleCSVTimeSeriesProvider(data_path+'SatelliteTimeSeries/')

# create dataset
dataset = TimeSeriesDataset(occurrences=data_path+'Presence_only_occurrences/Presences_only_train_sample.csv',
                            providers=[ts_red, ts_multi, ts_all])

# print random tensors from dataset
ids = [random.randint(0, len(dataset)-1) for i in range(5)]
for id in ids:
    tensor = dataset[id][0]
    label = dataset[id][1]
    print('Tensor type: {}, tensor shape: {}, label: {}'.format(type(tensor), tensor.shape, label))
    dataset.plot_ts(id)