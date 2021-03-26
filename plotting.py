import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def visualize_observation_patch(patch, landcover_labels=None, return_fig=False):
    """Plots patch data

    Parameters
    ----------
    patch : tuple of size 4 containing 2d array-like objects
        Patch data as returned by `load_patch`.
    landcover_labels : list of strings
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
        n_labels = np.max(landcover_patch)+1
        landcover_labels = np.arange(n_labels)

    cmap = plt.cm.get_cmap("viridis", len(landcover_labels))

    legend_elements = []
    for landcover_label, color in zip(landcover_labels, cmap.colors):
        legend_elements.append(
            Patch(color=color, label=landcover_label)
        )

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(9, 9))
    axes = axes.ravel()
    axes_iter = iter(axes)

    ax = next(axes_iter)
    ax.imshow(rgb_patch)
    ax.set_title("RGB satellite image")

    ax = next(axes_iter)
    ax.imshow(near_ir_patch, cmap="gray")
    ax.set_title("Near-IR satellite image")

    ax = next(axes_iter)
    im = ax.imshow(altitude_patch)
    ax.set_title("Altitude")
    fig.colorbar(im, ax=ax)

    ax = next(axes_iter)
    im = ax.imshow(landcover_patch, interpolation="none", cmap=cmap, vmin=0, vmax=len(legend_elements))
    ax.set_title("Land cover")
    visible_landcover_categories = np.unique(landcover_patch)
    legend = [legend_elements[i] for i in visible_landcover_categories]
    ax.legend(handles=legend, handlelength=.75, bbox_to_anchor=(1, 0.5), loc="center left")

    for ax in axes:
        ax.axis("off")

    fig.tight_layout()

    if return_fig:
        return fig
