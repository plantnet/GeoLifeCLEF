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
    patch_extractor = PatchExtractor('/home/data/rasters_GLC19', size=64, verbose=True)

    patch_extractor.add_all()

    dataset_list = [(43.61, 3.88), (42.61, 4.88)]
    labels_list = [0, 1]

    dataset_pytorch = GeoLifeClefDataset(patch_extractor, dataset_list, labels_list)
    print(dataset_pytorch[0])
    print(len(dataset_pytorch))

    data_loader = test_loader = torch.utils.data.DataLoader(dataset_pytorch, shuffle=True, batch_size=8)
    # dataset_pytorch can now be used in a data_loader
