from pathlib import Path

import numpy as np
import pytest

from data_loading.environmental_raster import PatchExtractor
from data_loading.pytorch_dataset import GeoLifeCLEF2022Dataset


DATA_PATH = Path("../data")


SUBSET_SIZE = {
    "train": 1587395,
    "val": 40080,
    "train+val": 1627475,
    "test": 36421,
}


@pytest.mark.parametrize("subset", SUBSET_SIZE.keys())
def test_pytorch_load_only_patches(subset):
    dataset = GeoLifeCLEF2022Dataset(DATA_PATH, subset, use_rasters=False)

    assert len(dataset) == SUBSET_SIZE[subset]

    result = dataset[0]

    if subset == "test":
        assert len(result) == 4
    else:
        assert len(result) == 2
        assert len(result[0]) == 4
        assert type(result[1]) == np.int64


@pytest.mark.parametrize("subset", SUBSET_SIZE.keys())
def test_pytorch_load_one_raster(subset):
    patch_extractor = PatchExtractor(DATA_PATH / "rasters", size=256)
    patch_extractor.append("bio_1")
    dataset = GeoLifeCLEF2022Dataset(
        DATA_PATH, subset, use_rasters=True, patch_extractor=patch_extractor
    )

    assert len(dataset) == SUBSET_SIZE[subset]

    result = dataset[0]

    if subset == "test":
        assert len(result) == 5
    else:
        assert len(result) == 2
        assert len(result[0]) == 5
        assert len(result[0][-1]) == 1
        assert type(result[1]) == np.int64


def test_pytorch_load_all():
    subset = "train"
    dataset = GeoLifeCLEF2022Dataset(DATA_PATH, subset, use_rasters=True)

    assert len(dataset) == SUBSET_SIZE[subset]

    result = dataset[0]

    if subset == "test":
        assert len(result) == 5
    else:
        assert len(result) == 2
        assert len(result[0]) == 5
        assert len(result[0][-1]) == 27
        assert type(result[1]) == np.int64
