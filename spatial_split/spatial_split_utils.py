from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from geopy import distance

if TYPE_CHECKING:
    import numpy.typing as npt

WGS84_EQUATOR_CIRCUMFERENCE = 40075.017  # km
WGS84_MERIDIAN_CIRCUMFERENCE = 40007.863  # km


def coords_to_bin(
    x: npt.NDArray,
    y: npt.NDArray,
    x_bin_width: float,
    y_bin_width: float,
) -> tuple[npt.NDArray[np.int_], npt.NDArray[np.int_]]:
    """
    x: list of positive east-west coordinates of some sort
    y: list of positive north-south coordinates of some sort
    x_bin_width: bin width for x
    y_bin_width: bin width for y
    """

    assert np.all(x > 0)
    assert np.all(y > 0)
    assert x_bin_width > 0
    assert y_bin_width > 0

    # Compute bins
    x_bin_list = np.array(np.floor(x / x_bin_width), dtype=int)
    y_bin_list = np.array(np.floor(y / y_bin_width), dtype=int)

    return (x_bin_list, y_bin_list)


def lon_to_global_easting(
    lon: npt.NDArray,
    origin: tuple[float, float] = (0.0, 0.0),
) -> npt.NDArray:
    # converts lon to a global (nonnegative) easting in km
    global_eastings_km = np.zeros_like(lon)
    for i in range(len(lon)):
        lon_ofs = lon[i] - origin[1]
        if lon_ofs < -180:
            lon_ofs = 180 - (lon_ofs + 180)
        elif lon_ofs > 180:
            lon_ofs = -180 - (lon_ofs - 180)
        global_eastings_km[i] = distance.geodesic(origin, (origin[0], lon[i])).km
        if lon_ofs < 0:
            # avoid negative eastings
            global_eastings_km[i] = WGS84_EQUATOR_CIRCUMFERENCE - global_eastings_km[i]

    return global_eastings_km


def lat_to_global_northing(
    lat: npt.NDArray,
    origin: tuple[float, float] = (0.0, 0.0),
) -> npt.NDArray:
    # converts lon to a global (nonnegative) northing in km
    global_northings_km = np.zeros_like(lat)
    for i in range(len(lat)):
        lat_ofs = lat[i] - origin[0]
        if lat_ofs < -90:
            lat_ofs = 90 - (lat_ofs + 90)
        elif lat_ofs > 90:
            lat_ofs = -90 - (lat_ofs - 90)
        global_northings_km[i] = distance.geodesic(origin, (lat[i], origin[1])).km
        if lat_ofs < 0:
            # avoid negative northings
            global_northings_km[i] = (
                WGS84_MERIDIAN_CIRCUMFERENCE - global_northings_km[i]
            )

    return global_northings_km


def assign_block_ids(
    lon: npt.NDArray,
    lat: npt.NDArray,
    east_west_bin_km: float,
    north_south_bin_km: float,
    origin: tuple[float, float] = (0.0, 0.0),
) -> npt.NDArray:
    assert np.max(origin) <= 180
    assert np.min(origin) >= -180

    easting_list_km = lon_to_global_easting(lon, origin)
    northing_list_km = lat_to_global_northing(lat, origin)

    easting_bin_list, northing_bin_list = coords_to_bin(
        easting_list_km, northing_list_km, east_west_bin_km, north_south_bin_km
    )

    block_id_list = []
    for i in range(len(easting_bin_list)):
        block_id_list.append(
            "e"
            + str(easting_bin_list[i]).zfill(4)
            + "n"
            + str(northing_bin_list[i]).zfill(4)
        )

    return np.array(block_id_list)
