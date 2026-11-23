"""Render paired Ghana terrain views from the completed analysis grid."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
import numpy as np
import rasterio


root = Path(__file__).resolve().parents[1]
cells = np.genfromtxt(
    root / "outputs/tables/ghana_intensity_comparison.csv",
    delimiter=",",
    skip_header=1,
    usecols=(1, 2, 4, 5, 9, 10),
)
x, y, observed, inclusion, x_m, y_m = cells.T
with rasterio.open(root / "data/processed/ghana_1km_elevation_m.tif") as raster:
    elevation = raster.read(1, masked=True).astype(float).filled(np.nan)
    rows, columns = rasterio.transform.rowcol(raster.transform, x, y)
    elevation = elevation[rows, columns]

if not np.isfinite(elevation).all():
    raise ValueError("Eligible cells must have finite elevation")

# Display bins use only eligible cells, just like the manuscript's intensity maps.
width_m = 8000
x_bin = np.floor(x_m / width_m).astype(int)
y_bin = np.floor(y_m / width_m).astype(int)
x_centres = np.arange(x_bin.min(), x_bin.max() + 1)
y_centres = np.arange(y_bin.min(), y_bin.max() + 1)
shape = (len(y_centres), len(x_centres))
bin_id = np.ravel_multi_index((y_bin - y_centres[0], x_bin - x_centres[0]), shape)
count = np.bincount(bin_id, minlength=np.prod(shape))


def average(values):
    total = np.bincount(bin_id, weights=values, minlength=count.size)
    return np.divide(total, count, out=np.full(count.shape, np.nan), where=count > 0).reshape(shape)


xx, yy = np.meshgrid((x_centres + 0.5) * width_m / 1000, (y_centres + 0.5) * width_m / 1000)
terrain = average(elevation)
fig = plt.figure(figsize=(12, 7.1), facecolor="white", layout="constrained")
panels = (
    ("Mapped cocoa occupancy", average(observed), "YlGn", Normalize(0, 100), "% of eligible cells"),
    (
        "Difference from environmental baseline",
        average(observed - inclusion),
        "BrBG",
        TwoSlopeNorm(vmin=-60, vcenter=0, vmax=60),
        "Observed minus mean simulated (percentage points)",
    ),
)
for index, (title, values, palette, scale, label) in enumerate(panels, start=1):
    axis = fig.add_subplot(1, 2, index, projection="3d", proj_type="ortho")
    colors = plt.get_cmap(palette)
    axis.plot_surface(
        xx, yy, terrain, facecolors=colors(scale(100 * values)),
        rstride=1, cstride=1, shade=False, linewidth=0, antialiased=False,
    )
    axis.set_title(title, fontsize=15, weight="bold", pad=3)
    axis.set_xlabel("Easting (km)", labelpad=8)
    axis.set_ylabel("Northing (km)", labelpad=9)
    axis.set_zlabel("Elevation (m)", labelpad=5)
    axis.set_xticks([150, 350, 550])
    axis.set_yticks([0, 250, 500, 750])
    axis.set_zlim(0, 800)
    axis.set_zticks([0, 400, 800])
    axis.set_box_aspect((np.ptp(xx), np.ptp(yy), 140))
    axis.view_init(elev=41, azim=-112)
    axis.tick_params(labelsize=8, pad=1)
    axis.grid(False)
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=scale, cmap=colors), ax=axis,
                       orientation="horizontal", shrink=0.82, pad=0.01,
                       extend="both" if index == 2 else "neither")
    bar.set_label(label)

fig.savefig(root / "outputs/figures/ghana_cocoa_terrain.png", dpi=180,
            bbox_inches="tight", pad_inches=0.2)
