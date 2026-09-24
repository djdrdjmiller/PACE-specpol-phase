"""Plots of individual matched sample centers."""

from pathlib import Path
import json
import numpy as np
import xarray as xr
from pace_specpol.reference import Reference
from pace_specpol.classification import sample_phase


def plot_cached_samples(
    cache_dir,
    field="liquid_index",
    reference=Reference(),
    ratio_threshold=1.27,
    li_threshold=0.3,
    point_size=0.3,
    vmin=None,
    vmax=None,
):
    """Plot every eligible cached sample at its center, without L3 averaging.

    Marker size is in points squared, not the physical footprint. Overlapping
    markers overplot in manifest order; this is not a unique-time mosaic.
    Reads local pair files only. Uses the same reference eligibility as the index.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    import cartopy.crs as ccrs

    if field not in ("liquid_index", "ratio", "phase"):
        raise ValueError("field must be liquid_index, ratio, or phase")
    source = Path(cache_dir)
    manifest = json.loads((source / "manifest.json").read_text())
    if not manifest.get("complete"):
        raise ValueError("Pair cache is incomplete")
    evaluate = reference.evaluator()
    fig, ax = plt.subplots(figsize=(15, 8), subplot_kw={"projection": ccrs.Robinson()})
    ax.set_global()
    if field == "phase":
        cmap = ListedColormap(["#4055df", "#298c42", "#d9363e"])
        style = {"cmap": cmap, "norm": BoundaryNorm([0.5, 1.5, 2.5, 3.5], 3)}
    else:
        defaults = (-0.5, 4.0) if field == "liquid_index" else (0.5, 3.5)
        style = {
            "cmap": "turbo",
            "vmin": defaults[0] if vmin is None else vmin,
            "vmax": defaults[1] if vmax is None else vmax,
        }
    artist = None
    for filename in manifest["files"]:
        with xr.open_dataset(source / filename, engine="h5netcdf") as ds:
            ds.load()
            baseline = evaluate(ds)
            ratio = np.full(ds.sizes["sample"], np.nan)
            np.divide(
                ds.raw_ratio.values,
                baseline,
                out=ratio,
                where=np.isfinite(baseline) & (baseline > 0),
            )
            li = ds.liquid_index.values
            good = np.isfinite(li) & np.isfinite(ratio)
            if not good.any():
                continue
            values = li if field == "liquid_index" else ratio
            if field == "phase":
                values = sample_phase(ratio, li, ratio_threshold, li_threshold)
            artist = ax.scatter(
                ds.longitude.values[good],
                ds.latitude.values[good],
                c=values[good],
                s=point_size,
                linewidths=0,
                transform=ccrs.PlateCarree(),
                rasterized=True,
                **style,
            )
    if artist is None:
        plt.close(fig)
        raise ValueError("No eligible samples")
    ax.coastlines(resolution="110m", linewidth=0.5, zorder=10)
    label = {
        "liquid_index": "Liquid index",
        "ratio": f"OCI ratio ({reference.mode})",
        "phase": "Sample phase candidate",
    }[field]
    bar = fig.colorbar(
        artist,
        ax=ax,
        orientation="horizontal",
        pad=0.04,
        extend="neither" if field == "phase" else "both",
        label=label,
    )
    if field == "phase":
        bar.set_ticks([1, 2, 3], labels=["Liquid", "Ice", "LTMP"])
    ax.set_title(
        label + " | cached sample centers; no geographic averaging\n"
        "Overlaps overplot in manifest order; marker size is not a footprint"
    )
    return fig, ax
