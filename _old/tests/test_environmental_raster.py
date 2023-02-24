from pathlib import Path

import numpy as np
import pytest

from data_loading.environmental_raster import PatchExtractor


DATA_PATH = Path("../data")


def test_patch_extractor_single_variable():
    extractor = PatchExtractor(DATA_PATH / "rasters", size=8)
    extractor.append("bio_1")

    patch = extractor[43.61, 3.88]
    expected_patch = np.array(
        [
            [
                [
                    14.6,
                    14.595834,
                    14.608334,
                    14.670834,
                    14.6375,
                    14.679167,
                    14.658334,
                    14.625,
                ],
                [
                    14.7,
                    14.625,
                    14.658333,
                    14.620833,
                    14.620833,
                    14.658334,
                    14.658334,
                    14.6625,
                ],
                [
                    14.7625,
                    14.775,
                    14.754167,
                    14.729167,
                    14.704166,
                    14.741667,
                    14.791667,
                    14.795834,
                ],
                [
                    14.775,
                    14.8125,
                    14.825,
                    14.825,
                    14.845833,
                    14.854167,
                    14.858334,
                    14.795834,
                ],
                [
                    14.770833,
                    14.795834,
                    14.845833,
                    14.849999,
                    14.845833,
                    14.891666,
                    14.908334,
                    14.866667,
                ],
                [
                    14.820833,
                    14.875,
                    14.883333,
                    14.904166,
                    14.916667,
                    14.9375,
                    14.954166,
                    14.9375,
                ],
                [
                    14.904166,
                    14.904166,
                    14.9625,
                    14.983334,
                    14.970834,
                    14.975,
                    15.0125,
                    15.004167,
                ],
                [
                    14.991667,
                    14.954166,
                    14.9625,
                    15.0125,
                    15.041667,
                    15.05,
                    15.045834,
                    15.008333,
                ],
            ]
        ],
        dtype=np.float32,
    )

    assert patch.shape == (1, 8, 8)
    np.testing.assert_allclose(patch, expected_patch)


def test_patch_extractor_several_variables():
    extractor = PatchExtractor(DATA_PATH / "rasters", size=256)
    extractor.append("bio_1")
    extractor.append("bio_2")
    extractor.append("bio_3")

    patch = extractor[43.61, 3.88]

    assert patch.shape == (3, 256, 256)


def test_patch_extractor_out_of_rasters_bounds():
    size = 2048
    extractor = PatchExtractor(DATA_PATH / "rasters", size=size)
    extractor.append("bio_1")

    patch = extractor[48.480688, -4.702032]

    assert patch.shape == (1, size, size)


def test_vector_extractor_single_variable():
    extractor = PatchExtractor(DATA_PATH / "rasters", size=1)
    extractor.append("bio_1")

    patch = extractor[43.61, 3.88]

    assert patch.shape == (1,)


@pytest.mark.parametrize("size", (1, 256))
def test_patch_extractor_out_of_bounds(size):
    extractor = PatchExtractor(DATA_PATH / "rasters", size=size)
    extractor.append("bio_1", out_of_bounds="ignore")
    extractor.append("bio_2", out_of_bounds="ignore")
    extractor.append("bio_3", out_of_bounds="ignore")
    patch = extractor[0.0, 0.0]
    if size == 1:
        assert patch.shape == (3,)
    else:
        assert patch.shape == (3, size, size)

    extractor = PatchExtractor(DATA_PATH / "rasters", size=size)
    extractor.append("bio_1", out_of_bounds="error")
    with pytest.raises(IndexError):
        extractor[0, 0]

    extractor = PatchExtractor(DATA_PATH / "rasters", size=size)
    extractor.append("bio_1", out_of_bounds="warn")
    with pytest.warns(UserWarning):
        patch = extractor[0, 0]
        if size == 1:
            assert patch.shape == (1,)
        else:
            assert patch.shape == (1, size, size)


def test_patch_plotting():
    extractor = PatchExtractor(DATA_PATH / "rasters", size=256)
    extractor.append("bio_1")
    extractor.append("bio_2")
    extractor.append("bio_3")

    extractor.plot((43.61, 3.88))
