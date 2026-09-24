"""Interactive analysis dashboards; optional visualization dependencies."""

from pathlib import Path
import json
import numpy as np
from .aggregation import cell_ids
from pace_specpol.workflow import VariantReference


def dashboard(
    index,
    threshold=1.27,
    min_count=1,
    output_dir="mixed_phase_results",
    preview_width=900,
    ratio_clim=(0.5, 2.0),
    li_clim=(-0.5, 4.0),
    map_extent=None,
    layout="classic",
    map_framing="extent",
    graticules=True,
    show=True,
):
    """Widget UI; threshold release updates local sparse counts, with no NASA I/O.

    The preview uses nearest-cell sampling on Robinson or regional polar maps.
    layout="slide" exports 16:9; map_framing="granule" adjusts only regional display
    bounds. NetCDF export always preserves the full native grid.
    show=False updates the same figure without displaying widgets, for batch PNGs.
    """
    import ipywidgets as widgets
    from IPython.display import display, clear_output
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm, LogNorm, TwoSlopeNorm
    from matplotlib.patches import Rectangle, Patch
    import cartopy.crs as ccrs
    from html import escape

    if layout not in ("classic", "slide") or map_framing not in ("extent", "granule"):
        raise ValueError("layout: classic/slide; map_framing: extent/granule")
    slide = layout == "slide"
    index.threshold_index(threshold)
    if min_count < 1:
        raise ValueError("min_count must be positive")
    colors = ["#4055df", "#298c42", "#d9363e"]
    raw = index.metadata["reference"]["mode"] == "raw"
    corrected = index.metadata["reference"]["mode"] == "corrected"
    mode_label = (
        "RAW ratio; unnormalized — phase candidates only"
        if raw
        else "Atmospherically corrected ratio — provisional phase candidates"
        if corrected
        else "Normalized ratio — provisional phase candidates"
    )
    if "variant" in index.metadata["reference"]:
        mode_label += (
            f" | {index.metadata['reference']['variant']}"
            f" | {index.metadata['reference'].get('population', 'own')} population"
        )
    # Precompute projection lookup, so no new map reprojection is needed per slide.
    crs = (
        ccrs.Robinson()
        if map_extent is None
        else ccrs.NorthPolarStereo(central_longitude=-45)
    )
    xmin, xmax = crs.x_limits
    ymin, ymax = crs.y_limits
    if map_extent is not None:
        west, east, south, north = map_extent
        if not (-180 <= west < east <= 180 and 0 <= south < north < 90):
            raise ValueError(
                "Arctic map_extent needs ordered lon bounds and 0 <= south < north < 90"
            )
        boundary_lon = np.r_[
            np.linspace(west, east, 200),
            np.linspace(west, east, 200),
            np.full(200, west),
            np.full(200, east),
        ]
        boundary_lat = np.r_[
            np.full(200, south),
            np.full(200, north),
            np.linspace(south, north, 200),
            np.linspace(south, north, 200),
        ]
        boundary = crs.transform_points(ccrs.PlateCarree(), boundary_lon, boundary_lat)
        xmin, xmax = boundary[:, 0].min(), boundary[:, 0].max()
        ymin, ymax = boundary[:, 1].min(), boundary[:, 1].max()
    if map_framing == "granule":
        if map_extent is None:
            raise ValueError("granule framing requires a regional map_extent")
        _, occupied_counts = index.mean_fields()
        cells = np.flatnonzero(occupied_counts)
        ny, nx = index.config.shape
        lons = (cells % nx + 0.5) * index.config.resolution - 180
        lats = (cells // nx + 0.5) * index.config.resolution - 90
        inside = (lons >= west) & (lons <= east) & (lats >= south) & (lats <= north)
        if not inside.any():
            raise ValueError("No occupied cells within map_extent for granule framing")
        xy = crs.transform_points(ccrs.PlateCarree(), lons[inside], lats[inside])
        ylo, yhi = np.nanmin(xy[:, 1]), np.nanmax(xy[:, 1])
        xlo, xhi = np.nanmin(xy[:, 0]), np.nanmax(xy[:, 0])
        cy = (ylo + yhi) / 2
        cx = (xmin + xmax) / 2
        # Match the slide panel aspect without stretching the map projection.
        aspect = (0.423 * (40 / 3)) / (0.295 * 7.5)
        span_y = max((yhi - ylo) * 1.14, 100000.0)
        span_x = max(xmax - xmin, (xhi - xlo) * 1.14, span_y * aspect)
        span_y = span_x / aspect
        xmin, xmax = cx - span_x / 2, cx + span_x / 2
        ymin, ymax = cy - span_y / 2, cy + span_y / 2
    width = int(preview_width)
    height = max(1, round(width * (ymax - ymin) / (xmax - xmin)))
    xx, yy = np.meshgrid(
        np.linspace(xmin, xmax, width, endpoint=False) + (xmax - xmin) / (2 * width),
        np.linspace(ymin, ymax, height, endpoint=False) + (ymax - ymin) / (2 * height),
    )
    points = ccrs.PlateCarree().transform_points(crs, xx, yy)
    lat, lon = points[..., 1], points[..., 0]
    valid = np.isfinite(lat) & np.isfinite(lon) & (abs(lat) <= 90) & (abs(lon) <= 180)
    ids = cell_ids(lat[valid], lon[valid], index.config.resolution)
    unique, inverse = np.unique(ids, return_inverse=True)
    preview_low, preview_high = index.low[:, unique], index.high[:, unique]
    means, sample_total = index.mean_fields()
    sampling = index.sampling_summary()
    for bounds in (ratio_clim, li_clim):
        if len(bounds) != 2 or not np.isfinite(bounds).all() or bounds[0] >= bounds[1]:
            raise ValueError("Color limits must be two finite increasing numbers")
    if not li_clim[0] < index.config.li_threshold < li_clim[1]:
        raise ValueError("li_clim must bracket the liquid-index threshold")

    def ratio_norm(t):
        # Preserve the requested limits except when necessary to bracket the slider.
        pad = 0.05 * (ratio_clim[1] - ratio_clim[0])
        return TwoSlopeNorm(
            vmin=min(ratio_clim[0], t - pad),
            vcenter=t,
            vmax=max(ratio_clim[1], t + pad),
        )

    with plt.ioff():
        fig = plt.figure(figsize=((40 / 3, 7.5) if slide else (16, 12)))
        if slide:
            # Explicit panel/colorbar slots prevent colorbars from shrinking GeoAxes.
            slots = [
                (0.065, 0.585, 0.423, 0.295),
                (0.55, 0.585, 0.423, 0.295),
                (0.065, 0.13, 0.423, 0.295),
                (0.55, 0.13, 0.423, 0.295),
            ]
            axh = fig.add_axes(slots[0])
            axm = fig.add_axes(slots[1], projection=crs)
        else:
            gs = fig.add_gridspec(2, 2)
            axh = fig.add_subplot(gs[0, 0])
            axm = fig.add_subplot(gs[0, 1], projection=crs)

        def decorate_map(ax):
            ax.coastlines(resolution="110m", linewidth=0.45, zorder=4)
            if graticules:
                gl = ax.gridlines(
                    draw_labels=slide,
                    xlocs=np.arange(-180, 181, 20 if map_extent else 60),
                    ylocs=np.arange(60, 91, 10)
                    if map_extent
                    else np.arange(-60, 61, 30),
                    linewidth=0.45,
                    color="#adb7bf",
                    alpha=0.65,
                    zorder=0.5,
                    x_inline=False,
                    y_inline=False,
                )
                if slide:
                    gl.top_labels = False
                    gl.right_labels = False
                    gl.xlabel_style = {"size": 7, "color": "#7b858d"}
                    gl.ylabel_style = {"size": 7, "color": "#7b858d"}
                    gl.xpadding = -9
                    gl.ypadding = -9
                    gl.rotate_labels = False

        def map_bar(im, ax):
            if slide:
                # Read the actual projected map box (may be narrower than its slot).
                box = ax.get_position()
                # Reserve extension triangles inside the map width as well.
                cax = fig.add_axes(
                    [box.x0 + box.width * 0.025, 0.073, box.width * 0.95, 0.017]
                )
                return fig.colorbar(
                    im,
                    cax=cax,
                    orientation="horizontal",
                    extend="both",
                    extendfrac=0.025,
                )
            return fig.colorbar(
                im, ax=ax, orientation="horizontal", pad=0.07, extend="both"
            )

        xminh, xmaxh = index.edges[[0, -1]]
        yminh, ymaxh = index.li_edges[[0, -1]]
        lt = index.config.li_threshold
        p_liq = Rectangle(
            (xminh, yminh),
            threshold - xminh,
            ymaxh - yminh,
            color=colors[0],
            alpha=0.14,
        )
        p_ice = Rectangle(
            (threshold, yminh),
            xmaxh - threshold,
            max(0, min(lt, ymaxh) - yminh),
            color=colors[1],
            alpha=0.14,
        )
        p_mix = Rectangle(
            (threshold, max(lt, yminh)),
            xmaxh - threshold,
            max(0, ymaxh - max(lt, yminh)),
            color=colors[2],
            alpha=0.14,
        )
        for p in (p_liq, p_ice, p_mix):
            axh.add_patch(p)
        # Density is normalized by ALL indexed samples, including off-plot values.
        area = np.diff(index.edges)[:, None] * np.diff(index.li_edges)[None, :]
        density = index.histogram / (index.metadata["totals"]["indexed_samples"] * area)
        positive = density[density > 0]
        if positive.size:
            norm = LogNorm(
                vmin=float(positive.min()),
                vmax=float(max(positive.max(), positive.min() * 1.01)),
            )
            h = axh.pcolormesh(
                index.edges,
                index.li_edges,
                np.ma.masked_where(density.T <= 0, density.T),
                cmap="turbo",
                norm=norm,
                shading="flat",
            )
            if slide:
                cb = fig.colorbar(
                    h,
                    cax=fig.add_axes([0.065, 0.492, 0.423, 0.017]),
                    orientation="horizontal",
                )
                cb.set_label("")
                cb.ax.tick_params(labelsize=8, pad=1)
            else:
                fig.colorbar(
                    h,
                    ax=axh,
                    orientation="horizontal",
                    pad=0.16,
                    label="Joint probability density (log color scale)",
                )
        else:
            axh.text(
                0.5,
                0.5,
                "No samples within histogram display limits",
                transform=axh.transAxes,
                ha="center",
            )
        vline = axh.axvline(threshold, color="k", lw=1.2)
        (hline,) = axh.plot([threshold, xmaxh], [lt, lt], color="k", lw=1.2)
        axh.set(
            xlim=(xminh, xmaxh),
            ylim=(yminh, ymaxh),
            ylabel="Liquid Index",
            xlabel=(
                "OCI 2260 / 1615 TOA reflectance ratio"
                if raw
                else "OCI 2260 / 1615 corrected reflectance ratio"
                if corrected
                else "OCI reflectance ratio / simulated all-liquid ratio"
            ),
        )
        axh.set_title(
            "Joint histogram | log probability density"
            if slide
            else "Joint histogram of matched observations"
        )
        labels = [
            axh.text(
                0.05, 0.96, "Liquid", color=colors[0], transform=axh.transAxes, va="top"
            ),
            axh.text(
                0.7, 0.96, "LTMP", color=colors[2], transform=axh.transAxes, va="top"
            ),
            axh.text(
                0.82, 0.03, "Ice", color=colors[1], transform=axh.transAxes, va="bottom"
            ),
        ]
        cmap = ListedColormap(colors)
        cmap.set_bad((1, 1, 1, 0))
        map_image = axm.imshow(
            np.ma.masked_all((height, width)),
            origin="lower",
            extent=(xmin, xmax, ymin, ymax),
            cmap=cmap,
            norm=BoundaryNorm([0.5, 1.5, 2.5, 3.5], 3),
            interpolation="nearest",
            zorder=2,
        )
        axm.set_xlim(xmin, xmax)
        axm.set_ylim(ymin, ymax)
        decorate_map(axm)
        axm.legend(
            handles=[
                Patch(color=c, label=n)
                for c, n in zip(colors, ["Liquid", "Ice", "LTMP"])
            ],
            loc="lower center",
            bbox_to_anchor=(0.5, -0.24 if slide else -0.12),
            ncol=3,
            frameon=False,
            fontsize=9 if slide else None,
        )
        mean_images = []
        mean_bars = []
        threshold_annotations = []
        for j, (title, limits) in enumerate(
            (
                ("Mean liquid index", li_clim),
                (
                    "Mean raw TOA ratio"
                    if raw
                    else "Mean corrected ratio"
                    if corrected
                    else "Mean normalized ratio",
                    ratio_clim,
                ),
            )
        ):
            ax = (
                fig.add_axes(slots[2 + j], projection=crs)
                if slide
                else fig.add_subplot(gs[1, j], projection=crs)
            )
            continuous = plt.get_cmap("RdBu_r").copy()
            continuous.set_bad((1, 1, 1, 0))
            continuous.set_under("#053061")
            continuous.set_over("#67001f")
            im = ax.imshow(
                np.ma.masked_all((height, width)),
                origin="lower",
                extent=(xmin, xmax, ymin, ymax),
                cmap=continuous,
                norm=(
                    TwoSlopeNorm(vmin=limits[0], vcenter=lt, vmax=limits[1])
                    if j == 0
                    else ratio_norm(threshold)
                ),
                interpolation="nearest",
                zorder=2,
            )
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            decorate_map(ax)
            ax.set_title(
                title if slide else title + " | eligible matched samples",
                y=1.02 if slide else None,
            )
            bar = map_bar(im, ax)
            if slide:
                bar.ax.plot(
                    [0.5, 0.5],
                    [0.55, 1.75],
                    transform=bar.ax.transAxes,
                    color="#333333",
                    lw=0.8,
                    clip_on=False,
                )
                text = bar.ax.text(
                    0.5,
                    1.9,
                    "",
                    transform=bar.ax.transAxes,
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )
                threshold_annotations.append(text)
                bar.ax.tick_params(labelsize=8, pad=2)
            else:
                bar.set_label(
                    "Midpoint = phase threshold; unequal scales on either side"
                )
            mean_bars.append(bar)
            mean_images.append(im)
        fig.suptitle(
            mode_label + "\n" + str(index.metadata.get("data_scope", "")),
            fontsize=10 if slide else 12,
            y=0.98 if slide else 0.98,
        )
        if slide:
            for ax in [axh, axm] + [im.axes for im in mean_images]:
                ax.title.set_fontsize(10)
                ax.tick_params(labelsize=9)
                ax.xaxis.label.set_size(9)
                ax.yaxis.label.set_size(9)
            fig.text(
                0.5,
                0.018,
                "White: no eligible value or phase tie  |  Colors centered on thresholds; unequal scales on either side",
                ha="center",
                fontsize=8,
                color="#555555",
            )
        else:
            fig.subplots_adjust(top=0.9, bottom=0.08, wspace=0.23, hspace=0.4)
    slider = widgets.FloatSlider(
        value=threshold,
        min=float(index.edges[0]),
        max=float(index.edges[-1]),
        step=index.config.ratio_step,
        description="Ratio cutoff",
        continuous_update=False,
        readout_format=".2f",
        layout=widgets.Layout(width="550px"),
    )
    minimum = widgets.BoundedIntText(
        value=min_count, min=1, max=1000000, description="Min samples"
    )
    export = widgets.Button(
        description="Save current PNG + NetCDF",
        icon="save",
        layout=widgets.Layout(width="250px"),
    )
    output = widgets.Output()
    status = widgets.HTML()

    def update(*_):
        # Snap widget arithmetic to exact stored threshold grid.
        t = float(index.edges[np.argmin(abs(index.edges - slider.value))])
        phase = index.dominant(
            index.counts(t, preview_low, preview_high), minimum.value
        )
        raster = np.zeros((height, width), dtype="int8")
        raster[valid] = phase[inverse]
        map_image.set_data(np.ma.masked_where(raster == 0, raster))
        for j, im in enumerate(mean_images):
            field = np.full((height, width), np.nan)
            values = means[j, unique].copy()
            values[sample_total[unique] < minimum.value] = np.nan
            field[valid] = values[inverse]
            im.set_data(np.ma.masked_invalid(field))
        mean_images[1].set_norm(ratio_norm(t))
        for j, (im, bar) in enumerate(zip(mean_images, mean_bars)):
            bounds = im.norm
            if slide:
                bar.set_ticks(
                    [
                        bounds.vmin,
                        (bounds.vmin + bounds.vcenter) / 2,
                        bounds.vcenter,
                        (bounds.vcenter + bounds.vmax) / 2,
                        bounds.vmax,
                    ]
                )
                from matplotlib.ticker import FormatStrFormatter

                bar.ax.xaxis.set_major_formatter(FormatStrFormatter("%.3g"))
                threshold_annotations[j].set_text(
                    ("LI threshold = " if j == 0 else "Ratio threshold = ")
                    + f"{bounds.vcenter:g}"
                )
            else:
                bar.set_ticks(
                    np.r_[
                        np.linspace(bounds.vmin, bounds.vcenter, 4),
                        np.linspace(bounds.vcenter, bounds.vmax, 4)[1:],
                    ]
                )
        selected = means[1, sample_total >= minimum.value]
        ratio_stats = (
            f"Mapped-cell ratio means: median {np.median(selected):.3g}, "
            f"P99 {np.percentile(selected, 99):.3g}, max {selected.max():.3g}; "
            f"{np.mean(selected > ratio_clim[1]):.2%} above requested color maximum {ratio_clim[1]:g}. "
            if selected.size
            else "No cells meet the minimum count. "
        )
        vline.set_xdata([t, t])
        hline.set_xdata([t, xmaxh])
        p_liq.set_width(t - xminh)
        for p in (p_ice, p_mix):
            p.set_x(t)
            p.set_width(xmaxh - t)
        # Region labels follow their regions as the dividing line moves.
        labels[0].set_position(((t - xminh) / (xmaxh - xminh) / 2, 0.96))
        labels[1].set_position((((t + xmaxh) / 2 - xminh) / (xmaxh - xminh), 0.96))
        labels[2].set_position((((t + xmaxh) / 2 - xminh) / (xmaxh - xminh), 0.03))
        for label in labels:
            label.set_ha("center")
        axm.set_title(
            (
                f"Modal phase | ratio cutoff {t:.2f}; LI {lt:g} | {index.config.resolution:g}°"
                if slide
                else f"Most frequent sample phase | ratio cut {t:.2f}, LI cut {lt:g}\n"
                f"{index.config.resolution:g}° grid; white = missing, tied, or insufficient samples"
            ),
            fontsize=10 if slide else None,
            y=1.02 if slide else None,
        )
        totals = index.metadata["totals"]
        status.value = (
            f"<b>{escape(mode_label)}</b><br>{totals['indexed_samples']:,} eligible samples; "
            f"{totals['histogram_outside_samples']:,} outside histogram display; "
            f"{totals['reference_invalid_samples']:,} excluded by reference coverage. "
            "Threshold changes use the cached index only. Mean fields do not depend on the phase cutoff.<br>"
            f"Per occupied native cell: mean {sampling['mean_samples_per_occupied_cell']:.2f}, "
            f"median {sampling['median_samples']:g}, P90 {sampling['p90_samples']:g}, "
            f"max {sampling['maximum_samples']}; "
            f"{sampling['single_sample_cell_fraction']:.1%} contain one sample.<br>"
            + ratio_stats
            + "<br>Ratio colors recenter on the cutoff; limits expand only if needed to bracket it. "
            "Preview uses nearest-cell sampling and can miss small patches; export preserves every native grid cell."
        )
        if not show:
            return
        with output:
            clear_output(wait=True)
            if slide:
                # Inline tight-bbox rendering can crop GeoAxes decorations. Render
                # the complete slide and scale its preview to the notebook width.
                import io, base64
                from IPython.display import HTML

                buffer = io.BytesIO()
                with plt.rc_context({"savefig.bbox": None}):
                    fig.savefig(buffer, format="png", dpi=120, bbox_inches=None)
                payload = base64.b64encode(buffer.getvalue()).decode("ascii")
                display(
                    HTML(
                        '<img alt="Phase dashboard" style="display:block;width:100%;max-width:1600px;height:auto" '
                        'src="data:image/png;base64,' + payload + '">'
                    )
                )
            else:
                display(fig)

    def save(_):
        t = float(index.edges[np.argmin(abs(index.edges - slider.value))])
        folder = Path(output_dir)
        folder.mkdir(parents=True, exist_ok=True)
        prefix = (
            folder
            / f"phase_{index.metadata['fingerprint'][:8]}_ratio_{t:.2f}_LI_{lt:g}_n{minimum.value}_{layout}_{map_framing}"
        )
        index.export(str(prefix) + ".nc", t, minimum.value)
        with plt.rc_context({"savefig.bbox": None}):
            fig.savefig(
                str(prefix) + ".png",
                dpi=240 if slide else 200,
                bbox_inches=None if slide else "tight",
            )
        status.value += f"<br>Saved {escape(str(prefix))}.nc and .png"

    slider.observe(update, names="value")
    minimum.observe(update, names="value")
    export.on_click(save)
    ui = widgets.VBox([widgets.HBox([slider, minimum]), export, status, output])
    if show:
        display(ui)
    update()
    return {
        "widget": ui,
        "slider": slider,
        "minimum": minimum,
        "figure": fig,
        "update": update,
        "export_button": export,
        "index": index,
        "mean_images": mean_images,
        "mean_colorbars": mean_bars,
        "map_axes": [axm] + [im.axes for im in mean_images],
        "phase_image": map_image,
    }


def reference_dashboard(cache_dir, index_config=None, **dashboard_options):
    """Lazy-build one sparse index per variant/population/ratio space; reuse on switch."""
    import ipywidgets as widgets
    from IPython.display import clear_output
    from . import aggregation as mp

    root = Path(cache_dir)
    manifest = json.loads((root / "manifest.json").read_text())
    if not manifest.get("complete"):
        raise ValueError("Finish reference cache first")
    variants = list(manifest["settings"]["variants"])
    variant = widgets.Dropdown(
        options=variants,
        value="oci_2260" if "oci_2260" in variants else variants[0],
        description="CER / COT:",
    )
    population = widgets.Dropdown(
        options=[
            ("Each mode’s valid samples", "own"),
            ("Common valid samples", "common"),
        ],
        description="Population:",
    )
    mode = widgets.Dropdown(
        options=[
            ("Normalized", "normalized"),
            ("Corrected", "corrected"),
            ("Raw TOA", "raw"),
        ],
        description="Ratio:",
    )
    button = widgets.Button(description="Show selection", button_style="primary")
    status = widgets.HTML(
        "Select a mode, then Show selection. First use builds a local index; no satellite reads."
    )
    output = widgets.Output()
    current = None

    def show(_):
        nonlocal current
        button.disabled = True
        try:
            with output:
                clear_output(wait=True)
                if current is not None:
                    import matplotlib.pyplot as plt

                    plt.close(current["figure"])
                    current["widget"].close()
                idx = mp.build_index(
                    root,
                    reference=VariantReference(
                        variant.value, population.value, mode.value
                    ),
                    config=index_config or mp.IndexConfig(),
                )
                status.value = f"{variant.value} | {population.value} | {mode.value} | {idx.metadata['totals']['indexed_samples']:,} samples"
                current = dashboard(idx, **dashboard_options)
        except Exception as e:
            status.value = "Selection failed; see error below."
            with output:
                print(type(e).__name__ + ": " + str(e))
        finally:
            button.disabled = False

    button.on_click(show)
    return widgets.VBox(
        [widgets.HBox([variant, population, mode]), button, status, output]
    )
