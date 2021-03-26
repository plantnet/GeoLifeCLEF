from pathlib import Path

import numpy as np
import pytest

from data_loading.pytorch_dataset import GeoLifeCLEF2021Dataset


DATA_PATH = Path("../data")


SUBSET_SIZE = {
    "train": 1833272,
    "val": 45446,
    "train+val": 1878718,
    "test": 42405,
}


@pytest.mark.parametrize("subset", SUBSET_SIZE.keys())
def test_pytorch_load_only_patches(subset):
    dataset = GeoLifeCLEF2021Dataset(DATA_PATH, subset, use_rasters=False)

    assert len(dataset) == SUBSET_SIZE[subset]

    result = dataset[0]

    if subset == "test":
        assert len(result) == 6
    else:
        assert len(result) == 2
        assert len(result[0]) == 6
        assert type(result[1]) == np.int64


def test_pytorch_load_all():
    subset = "train"
    dataset = GeoLifeCLEF2021Dataset(DATA_PATH, subset, use_rasters=True)

    assert len(dataset) == SUBSET_SIZE[subset]

    result = dataset[0]

    if subset == "test":
        assert len(result) == 33
    else:
        assert len(result) == 2
        assert len(result[0]) == 33
        assert type(result[1]) == np.int64
