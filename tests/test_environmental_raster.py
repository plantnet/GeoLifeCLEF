from pathlib import Path

import pytest

from data_loading.environmental_raster import PatchExtractor


DATA_PATH = Path("../data")


def test_patch_extractor_single_variable():
    extractor = PatchExtractor(DATA_PATH / "rasters", size=256)
    extractor.append("bio_1")

    patch = extractor[43.61, 3.88]

    assert patch.shape == (1, 256, 256)


def test_patch_extractor_several_variables():
    extractor = PatchExtractor(DATA_PATH / "rasters", size=256)
    extractor.append("bio_1")
    extractor.append("bio_2")
    extractor.append("bio_3")

    patch = extractor[43.61, 3.88]

    assert patch.shape == (3, 256, 256)


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
