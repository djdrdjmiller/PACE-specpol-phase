"""Airborne context plots. These do not establish LTMP ground truth."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_matchup_context(
    sat,
    air,
    selected_matches,
    footprints,
    SELECTED_WINDOW_MIN,
    MAX_DISTANCE_KM,
    VARIANT,
    RATIO_CUTOFF,
    LI_CUTOFF,
    VALIDATION_DIR,
    tag,
):
    import cartopy.crs as ccrs

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(projection=ccrs.NorthPolarStereo(central_longitude=-45))
    ax.set_extent([-85, -10, 58, 87], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.scatter(
        sat.longitude,
        sat.latitude,
        s=1,
        color="0.75",
        transform=ccrs.PlateCarree(),
        label="Satellite samples",
    )
    ax.plot(
        air.longitude,
        air.latitude,
        color="#e66101",
        lw=1,
        transform=ccrs.PlateCarree(),
        label="G-III lidar flight",
    )
    if not selected_matches.empty:
        ax.scatter(
            selected_matches.air_longitude,
            selected_matches.air_latitude,
            s=4,
            color="#542788",
            transform=ccrs.PlateCarree(),
            label=f"Matches within {SELECTED_WINDOW_MIN} min",
        )
    ax.legend()
    ax.set_title("2024-06-10 | geographic overlap and accepted matches")
    plt.show()

    if not footprints.empty:
        palette = {
            "water_dominant": "#2166ac",
            "ice_dominant": "#b2182b",
            "unknown": "0.65",
        }
        fig, axes = plt.subplots(3, 1, figsize=(11, 11))
        for label, color in palette.items():
            q = footprints[footprints.lidar_top == label].sort_values("time")
            q = q[q.valid & np.isfinite(q.ratio)]
            axes[0].scatter(
                q.ratio, q.liquid_index, s=20, color=color, label=label, alpha=0.75
            )
            axes[1].scatter(q.time, q.liquid_index, s=15, color=color)
            axes[2].scatter(q.time, q.ratio, s=15, color=color)
        axes[0].axvline(RATIO_CUTOFF, color="k", ls="--")
        axes[0].axhline(LI_CUTOFF, color="k", ls="--")
        axes[0].set(xlabel="Normalized OCI 2260/1615 ratio", ylabel="Liquid index")
        axes[0].legend()
        axes[1].set(ylabel="Liquid index")
        axes[2].set(ylabel="Normalized ratio", xlabel="Satellite UTC time")
        axes[1].axhline(LI_CUTOFF, color="k", ls="--")
        axes[2].axhline(RATIO_CUTOFF, color="k", ls="--")
        fig.suptitle(
            f"{VARIANT} | {SELECTED_WINDOW_MIN} min / {MAX_DISTANCE_KM:g} km | lidar top evidence"
        )
        fig.tight_layout()
        fig.savefig(VALIDATION_DIR / f"{tag}_comparison.png", dpi=160)
        plt.show()


def plot_lidar_curtain(air, sat, hsrl_path, VALIDATION_DIR, tag):
    # Lidar curtain near the satellite overpass, including earlier profiles.
    import h5py
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch

    view = np.flatnonzero(
        (air.time >= sat.time.min() - pd.Timedelta(minutes=60))
        & (air.time <= sat.time.max() + pd.Timedelta(minutes=10))
    )
    if len(view):
        stride = max(1, int(np.ceil(len(view) / 2500)))
        rows = view[::stride]
        with h5py.File(hsrl_path) as f:
            z = np.asarray(f["DataProducts/CloudMask/Altitude"]).ravel() / 1000
            phase = f["DataProducts/CloudMask/Type"][rows, :]
        colors = ["#2166ac", "#fdae61", "#b2182b", "#762a83", "#cccccc"]
        cmap = ListedColormap(colors)
        cmap.set_bad("white")
        fig, ax = plt.subplots(figsize=(13, 5))
        ax.pcolormesh(
            air.time.iloc[rows],
            z,
            np.ma.masked_invalid(phase.T),
            cmap=cmap,
            norm=BoundaryNorm(np.arange(-0.5, 5), 5),
            shading="nearest",
        )
        ax.plot(
            air.time.iloc[rows],
            air.cloud_top_m.iloc[rows] / 1000,
            color="k",
            lw=0.6,
            label="Cloud top",
        )
        ax.set(
            ylim=(0, 10),
            ylabel="Altitude (km)",
            xlabel="G-III UTC time",
            title="Lidar phase context: earlier sampling, not a simultaneous satellite curtain",
        )
        ax.legend(
            handles=[
                Patch(color=c, label=l)
                for c, l in zip(
                    colors,
                    [
                        "Water-dominant",
                        "Mix: ambiguous",
                        "Ice-dominant",
                        "HOI",
                        "Aerosol",
                    ],
                )
            ],
            loc="upper right",
        )
        fig.tight_layout()
        fig.savefig(VALIDATION_DIR / f"{tag}_lidar_curtain.png", dpi=160)
        plt.show()
