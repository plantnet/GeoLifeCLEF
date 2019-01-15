import torch
from torch.utils.data import Dataset

from environmental_raster_glc import PatchExtractor


class GeoLifeClefDataset(Dataset):
    def __init__(self, extractor, dataset, labels):
        self.extractor = extractor
        self.labels = labels
        self.dataset = dataset

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        tensor = self.extractor[self.dataset[idx]]
        return torch.from_numpy(tensor).float(), self.labels[idx]


if __name__ == '__main__':
    patch_extractor = PatchExtractor('/data/rasters_GLC19', size=64, verbose=True)

    patch_extractor.append('chbio_1')
    patch_extractor.append('text')
    # patch_extractor.add_all()

    # example of dataset
    dataset_list = [(43.61, 3.88), (42.61, 4.88), (46.15, -1.1), (49.54, -1.7)]
    labels_list = [0, 1, 0, 1]

    dataset_pytorch = GeoLifeClefDataset(patch_extractor, dataset_list, labels_list)

    print(len(dataset_pytorch), 'elements in the dataset')

    # dataset_pytorch can now be used in a data_loader
    data_loader = torch.utils.data.DataLoader(dataset_pytorch, shuffle=True, batch_size=2)

    for batch in data_loader:
        data, label = batch
        print('[batch, channels, width, height]:', data.size())
        print('[batch]:', label)
        print('*' * 5)
