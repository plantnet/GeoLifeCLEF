from __future__ import annotations
from typing import Optional, TYPE_CHECKING
from collections.abc import Collection

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

if TYPE_CHECKING:
    import numpy.typing as npt
    import pandas as pd


def plot_map(
    *,
    region: Optional[str] = None,
    extent: Optional[npt.ArrayLike] = None,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Plots a map to show the observations on

    Parameters
    ----------
    region: string, either "fr" or "us"
        Region to show, France or US.
    extent: array-like of form [longitude min, longitude max, latitude min, latitude max]
        Explicit extent of the area to show, e.g., for zooming.
    ax: plt.Axes
        Provide an Axes to use instead of creating one.

    Returns
    -------
    plt.Axes:
        Returns the used Axes.
    """
    if region == "fr":
        extent = [-5.5, 10, 41, 52]
    elif region == "us":
        extent = [-126, -66, 24, 50]
    elif region is None and extent is None:
        raise ValueError("Either region or extent must be set")

    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    if ax is None:
        ax = plt.axes(projection=ccrs.PlateCarree())

    ax.set_extent(extent, crs=ccrs.PlateCarree())

    states_provinces = cfeature.NaturalEarthFeature(
        category="cultural",
        name="admin_0_countries",
        scale="10m",
        facecolor="none",
    )

    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(states_provinces, edgecolor="gray")

    ax.gridlines(
        draw_labels=True, dms=True, x_inline=False, y_inline=False, linestyle="--"
    )
    ax.set_aspect(1.25)

    return ax


def visualize_observation_patch(
    patch: tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray],
    *,
    observation_data: Optional[pd.Series] = None,
    landcover_labels: Optional[Collection] = None,
    return_fig: bool = False,
) -> Optional[plt.Figure]:
    """Plots patch data

    Parameters
    ----------
    patch : tuple of size 4 containing 2d array-like objects
        Patch data as returned by `load_patch`.
    observation_data : pandas Series
        Row of the dataframe containing the data of the observation.
    landcover_labels : list
        Labels corresponding to the landcover codes.
    return_fig : boolean
        If True, returns the created plt.Figure object

    Returns
    -------
    fig : plt.Figure
        If return_fig is True, the used plt.Figure object    Returns
    """
    rgb_patch, near_ir_patch, altitude_patch, landcover_patch = patch

    if landcover_labels is None:
        n_labels = np.max(landcover_patch) + 1
        landcover_labels = np.arange(n_labels)

    cmap = plt.cm.get_cmap("viridis", len(landcover_labels))

    legend_elements = []
    for landcover_label, color in zip(landcover_labels, cmap.colors):
        legend_elements.append(Patch(color=color, label=landcover_label))

    if observation_data is not None:
        import cartopy.crs as ccrs

        localisation = observation_data[["latitude", "longitude"]]
        species_id = getattr(observation_data, "species_id", None)
        species_name = getattr(observation_data, "GBIF_species_name", None)
        kingdom_name = getattr(observation_data, "GBIF_kingdom_name", None)

        fig = plt.figure(figsize=(15, 10))

        gs = fig.add_gridspec(1, 2, width_ratios=[1, 2])

        gs1 = gs[0].subgridspec(1, 1)
        ax = fig.add_subplot(gs1[0], projection=ccrs.PlateCarree())
        region = "fr" if localisation[1] > -6 else "us"
        plot_map(region=region, ax=ax)
        ax.scatter(
            localisation[1],
            localisation[0],
            marker="o",
            s=100,
            transform=ccrs.PlateCarree(),
        )
        ax.set_title("Observation localisation")

        s = "Observation id: {}".format(observation_data.name)
        s += "\nLocalisation: {:.3f}, {:.3f}".format(*localisation)
        if species_id:
            s += "\nSpecies id: {}".format(species_id)
        if species_name:
            s += "\nSpecies name: {}".format(species_name)
        if kingdom_name:
            s += "\nKingdom: {}".format(kingdom_name)
        pos = ax.get_position()
        fig.text(
            pos.x1 - pos.x0, pos.y0 - 0.2 * (pos.y1 - pos.y0), s, ha="center", va="top"
        )

        gs2 = gs[1].subgridspec(2, 2)
        axes = np.asarray([fig.add_subplot(gs) for gs in gs2])
    else:
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(13, 8))

    axes = axes.ravel()
    axes_iter = iter(axes)

    ax = next(axes_iter)
    ax.imshow(rgb_patch)
    ax.set_title("RGB image")

    ax = next(axes_iter)
    ax.imshow(near_ir_patch, cmap="gray")
    ax.set_title("Near-IR image")

    ax = next(axes_iter)
    vmin = round(altitude_patch.min(), -1)
    vmax = round(altitude_patch.max(), -1) + 10
    ax.imshow(altitude_patch)
    CS2 = ax.contour(
        altitude_patch,
        levels=np.arange(vmin, vmax, step=10),
        colors="w",
    )
    ax.clabel(CS2, inline=True, fontsize=10)
    ax.set_aspect("equal")
    ax.set_title("Altitude (in meters)")

    ax = next(axes_iter)
    ax.imshow(
        landcover_patch,
        interpolation="none",
        cmap=cmap,
        vmin=0,
        vmax=len(legend_elements),
    )
    ax.set_title("Land cover")
    visible_landcover_categories = np.unique(landcover_patch)
    legend = [legend_elements[i] for i in visible_landcover_categories]
    ax.legend(
        handles=legend, handlelength=0.75, bbox_to_anchor=(1, 0.5), loc="center left"
    )

    for ax in axes:
        ax.axis("off")

    if observation_data is None:
        fig.tight_layout()

    if return_fig:
        return fig
    else:
        return None
