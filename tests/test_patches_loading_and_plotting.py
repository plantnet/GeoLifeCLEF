from pathlib import Path

import numpy as np
import pytest

from data_loading.common import load_patch
from plotting import visualize_observation_patch


DATA_PATH = Path("../data")


@pytest.mark.parametrize("observation_id", (10561900, 22068100))
def test_load_patch(observation_id):
    patches = load_patch(observation_id, DATA_PATH / "patches", return_arrays=True)

    assert len(patches) == 4

    rgb_patch, near_ir_patch, altitude_patch, landcover_patch = patches

    assert rgb_patch.shape == (256, 256, 3)
    assert rgb_patch.dtype == np.uint8

    assert near_ir_patch.shape == (256, 256)
    assert near_ir_patch.dtype == np.uint8

    assert altitude_patch.shape == (256, 256)
    assert altitude_patch.dtype == np.int16

    assert landcover_patch.shape == (256, 256)
    assert landcover_patch.dtype == np.uint8


@pytest.mark.parametrize("observation_id", (10561900, 22068100))
def test_patch_plotting(observation_id):
    patch = load_patch(observation_id, DATA_PATH / "patches")
    visualize_observation_patch(patch)
