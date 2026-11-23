"""Render Ghana occupancy and accessibility terrain views from the completed analysis grid."""

from contextlib import closing
from pathlib import Path
import sqlite3

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Patch
import numpy as np


root = Path(__file__).resolve().parents[1]
cells = np.genfromtxt(
    root / "outputs/tables/ghana_intensity_comparison.csv",
    delimiter=",",
    names=True,
    usecols=(
        "cell", "observed_presence", "spatial_block_20km_crossfit_simulated_inclusion",
        "x_m", "y_m",
    ),
)
cells.sort(order="cell")
observed = cells["observed_presence"]
inclusion = cells["spatial_block_20km_crossfit_simulated_inclusion"]
x_m, y_m = cells["x_m"], cells["y_m"]
grid_path = root / "data/processed/ghana_1km_candidate_points.gpkg"
with closing(sqlite3.connect(f"{grid_path.as_uri()}?mode=ro", uri=True)) as connection:
    grid = np.asarray(connection.execute("""
        SELECT cell, elevation_m, nearest_major_road_distance_m, nearest_town_distance_m
        FROM ghana_1km_candidate_points
        ORDER BY cell
    """).fetchall())

cell_ids, elevation, road_distance, town_distance = grid.T
if not np.array_equal(cell_ids, cells["cell"]) or np.unique(cell_ids).size != cells.size:
    raise ValueError("The intensity outputs and candidate grid must have the same unique cell IDs")
if not np.isfinite(grid).all():
    raise ValueError("Eligible cells must have finite terrain and distance values")

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
def render_pair(filename, panels):
    fig = plt.figure(figsize=(12, 7.1), facecolor="white", layout="constrained")
    for index, (title, values, palette, scale, label, extend) in enumerate(panels, start=1):
        axis = fig.add_subplot(1, 2, index, projection="3d", proj_type="ortho")
        colors = plt.get_cmap(palette)
        axis.plot_surface(
            xx, yy, terrain, facecolors=colors(scale(average(values))),
            rstride=1, cstride=1, shade=False, linewidth=0, antialiased=False,
        )
        axis.set_title(title, fontsize=15, weight="bold", pad=3)
        axis.set_xlabel("Easting (km)", labelpad=8)
        axis.set_ylabel("Northing (km)", labelpad=9)
        axis.set_zlabel("Elevation (m)", labelpad=5)
        axis.set_xticks([150, 350, 550])
        axis.set_yticks([0, 250, 500, 750])
        axis.set_zlim(0, max(800, np.nanmax(terrain)))
        axis.set_zticks([0, 400, 800])
        axis.set_box_aspect((np.ptp(xx), np.ptp(yy), 140))
        axis.view_init(elev=41, azim=-112)
        axis.tick_params(labelsize=8, pad=1)
        axis.grid(False)
        bar = fig.colorbar(
            plt.cm.ScalarMappable(norm=scale, cmap=colors), ax=axis,
            orientation="horizontal", shrink=0.82, pad=0.01, extend=extend,
        )
        bar.set_label(label)

    fig.savefig(root / "outputs/figures" / filename, dpi=180,
                bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)


render_pair("ghana_cocoa_terrain.png", (
    ("Mapped cocoa occupancy", 100 * observed, "YlGn",
     Normalize(0, 100), "% of eligible cells", "neither"),
    ("Difference from environmental baseline", 100 * (observed - inclusion), "BrBG",
     TwoSlopeNorm(vmin=-60, vcenter=0, vmax=60),
     "Observed minus mean simulated (percentage points)", "both"),
))

distance_scale = Normalize(0, 100)
render_pair("ghana_accessibility_terrain.png", (
    ("Distance to a major road", road_distance / 1000, "cividis",
     distance_scale, "Mean straight-line distance (km)", "max"),
    ("Distance to an urban centre", town_distance / 1000, "cividis",
     distance_scale, "Mean straight-line distance (km)", "max"),
))

quadrats = np.genfromtxt(
    root / "outputs/tables/ghana_quadrat_intensity.csv", delimiter=",", names=True,
)
quadrat_x = quadrats["x_m"] / 1000
quadrat_y = quadrats["y_m"] / 1000
observed_intensity = 100 * quadrats["observed_intensity"]
simulated_intensity = 100 * quadrats["simulated_median"]
residual = observed_intensity - simulated_intensity
intensity_scale = Normalize(0, 100)
residual_scale = Normalize(-100, 100)
intensity_palette = plt.get_cmap("YlGn")
residual_palette = plt.get_cmap("BrBG")
interval_labels = ("Below interval", "Within interval", "Above interval")
interval_palette = ("#a9521c", "#b6bac0", "#147a71")
interval_class = np.select(
    [quadrats["observed_intensity"] < quadrats["simulated_lower"],
     quadrats["observed_intensity"] > quadrats["simulated_upper"]],
    [0, 2], default=1,
)


def render_quadrats(filename, panels, limits, ticks, height_label, color_key):
    fig = plt.figure(figsize=(6.7 * len(panels), 6.5), layout="constrained")
    for index, (title, heights, colors) in enumerate(panels, start=1):
        axis = fig.add_subplot(1, len(panels), index, projection="3d", proj_type="ortho")
        axis.bar3d(
            quadrat_x - 11.5, quadrat_y - 11.5, np.minimum(heights, 0),
            23, 23, np.abs(heights), color=colors, shade=True,
            edgecolor="#ffffff", linewidth=0.15,
        )
        axis.set_title(title, fontsize=14, weight="bold")
        axis.set(xlabel="Easting (km)", ylabel="Northing (km)", zlabel=height_label,
                 xticks=[150, 350, 550], yticks=[0, 250, 500, 750],
                 zlim=limits, zticks=ticks)
        axis.set_box_aspect((np.ptp(quadrat_x), np.ptp(quadrat_y), 210))
        axis.view_init(elev=46, azim=-112)
        axis.tick_params(labelsize=8, pad=1)
        axis.grid(False)
    if color_key is not None:
        fig.legend(handles=color_key, loc="outside lower center", ncol=3, frameon=False)
    fig.savefig(root / "outputs/figures" / filename, dpi=180,
                bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)


render_quadrats(
    "ghana_quadrat_residual_3d.png",
    [("Observed minus simulated median", residual, residual_palette(residual_scale(residual)))],
    (-100, 100), [-100, -50, 0, 50, 100], "Difference (percentage points)",
    [Patch(color=residual_palette(value), label=label)
     for value, label in zip((0.15, 0.85), ("Below simulated median", "Above simulated median"))],
)
render_quadrats(
    "ghana_quadrat_intensity_3d.png",
    [("Observed", observed_intensity, intensity_palette(intensity_scale(observed_intensity))),
     ("Simulated median", simulated_intensity, intensity_palette(intensity_scale(simulated_intensity)))],
    (0, 100), [0, 25, 50, 75, 100], "Occupied eligible cells (%)", None,
)
render_quadrats(
    "ghana_quadrat_intervals_3d.png",
    [("Pointwise simulation intervals", observed_intensity,
      np.asarray(interval_palette)[interval_class])],
    (0, 100), [0, 25, 50, 75, 100], "Occupied eligible cells (%)",
    [Patch(color=color, label=label) for color, label in zip(interval_palette, interval_labels)],
)
