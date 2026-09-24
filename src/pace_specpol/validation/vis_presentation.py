"""Fixed-size dashboard frames and exactly aligned map crops for slide decks."""

from datetime import datetime, timezone
from io import BytesIO
from pathlib import Path

from pace_specpol.liquid_reference import sha256
from pace_specpol.oci_companion import atomic_json

from .vis_dashboard import dashboard
from .vis_truecolor import rgb_on_map


def export_presentation(
    index,
    output_dir,
    *,
    thresholds=(1.10, 1.27, 1.45),
    context_threshold=1.27,
    rgb=None,
    dpi=240,
    overwrite=False,
    rgb_black=0.0,
    rgb_white=1.0,
    rgb_gamma=2.2,
    rgb_max_distance_km=4.0,
    **dashboard_options,
):
    """Export one fixed 16:9 canvas per threshold and a true-color context frame.

    All map crops use one integer pixel box from the modal-phase axes. The same
    Figure/GeoAxes is reused for every frame, including RGB, preventing drift.
    Data means remain fixed; the ratio color scale recenters as in the live UI.
    """
    import matplotlib.pyplot as plt
    from PIL import Image

    thresholds = tuple(float(t) for t in thresholds)
    if not thresholds or len(set(thresholds)) != len(thresholds):
        raise ValueError("Provide distinct thresholds")
    for threshold in (*thresholds, context_threshold):
        index.threshold_index(threshold)
    if dpi < 60 or int(dpi) != dpi:
        raise ValueError("dpi must be an integer >= 60")
    if any(
        k in dashboard_options for k in ("layout", "show", "threshold", "output_dir")
    ):
        raise ValueError(
            "Presentation export fixes layout, show, threshold and output_dir"
        )
    root = Path(output_dir)
    if root.exists() and any(root.iterdir()) and not overwrite:
        raise FileExistsError(
            "Export folder is not empty; change EXPORT_LABEL or set OVERWRITE=True"
        )
    root.mkdir(parents=True, exist_ok=True)
    ui = dashboard(
        index,
        threshold=thresholds[0],
        layout="slide",
        show=False,
        output_dir=root,
        **dashboard_options,
    )
    fig = ui["figure"]
    manifest = dict(
        complete=False,
        created_utc=datetime.now(timezone.utc).isoformat(),
        index_fingerprint=index.metadata["fingerprint"],
        reference=index.metadata["reference"],
        data_scope=index.metadata.get("data_scope"),
        thresholds=list(thresholds),
        liquid_index_threshold=index.config.li_threshold,
        resolution_deg=index.config.resolution,
        dpi=dpi,
        dashboard_options=dashboard_options,
        files=[],
    )
    atomic_json(root / "presentation_manifest.json", manifest)
    crop_box = None

    def save_frame(stem, crop_stem):
        nonlocal crop_box
        buffer = BytesIO()
        # Explicit rc override defeats Jupyter's usual tight-bbox cropping.
        with plt.rc_context({"savefig.bbox": None}):
            fig.savefig(
                buffer, format="png", dpi=dpi, bbox_inches=None, facecolor="white"
            )
        buffer.seek(0)
        with Image.open(buffer) as rendered:
            if crop_box is None:
                box = ui["map_axes"][0].get_position()
                width, height = rendered.size
                crop_box = (
                    round(box.x0 * width),
                    round((1 - box.y1) * height),
                    round(box.x1 * width),
                    round((1 - box.y0) * height),
                )
                manifest.update(
                    canvas_pixels=list(rendered.size),
                    map_crop_box_pixels=list(crop_box),
                    map_pixels=[crop_box[2] - crop_box[0], crop_box[3] - crop_box[1]],
                    projection=str(ui["map_axes"][0].projection),
                    projected_xlim=list(ui["map_axes"][0].get_xlim()),
                    projected_ylim=list(ui["map_axes"][0].get_ylim()),
                )
            if list(rendered.size) != manifest["canvas_pixels"]:
                raise RuntimeError("Figure canvas changed between frames")
            full = root / f"{stem}.png"
            full.write_bytes(buffer.getvalue())
            crop = root / f"{crop_stem}.png"
            rendered.crop(crop_box).save(crop)
        for path in (full, crop):
            manifest["files"].append(dict(name=path.name, sha256=sha256(path)))

    try:
        for threshold in thresholds:
            ui["slider"].value = threshold
            tag = f"{threshold:.2f}".replace(".", "p")
            save_frame(f"dashboard_ratio_{tag}", f"phase_map_ratio_{tag}")
        if rgb is not None:
            ui["slider"].value = context_threshold
            ax = ui["map_axes"][0]
            phase = ui["phase_image"]
            rgba = rgb_on_map(
                rgb,
                ax,
                phase.get_array().shape,
                black=rgb_black,
                white=rgb_white,
                gamma=rgb_gamma,
                max_distance_km=rgb_max_distance_km,
            )
            phase.set_visible(False)
            ax.imshow(
                rgba,
                extent=phase.get_extent(),
                origin="lower",
                interpolation="nearest",
                zorder=2,
            )
            ax.get_legend().set_visible(False)
            ax.set_title(
                "OCI true color | visible TOA reflectance | 5 km L1C",
                fontsize=10,
                y=1.02,
            )
            for text in fig.texts:
                if text.get_text().startswith("White: no eligible"):
                    text.set_text(
                        "True color: full visible scene; white = no RGB coverage  |  Phase diagnostics: eligible matched samples"
                    )
            save_frame("dashboard_truecolor_context", "truecolor_map")
            manifest["truecolor"] = dict(
                source=dict(rgb.attrs),
                context_threshold=context_threshold,
                black=rgb_black,
                white=rgb_white,
                gamma=rgb_gamma,
                max_distance_km=rgb_max_distance_km,
                stretch="clip((rho-black)/(white-black), 0, 1) ** (1/gamma)",
            )
        manifest["complete"] = True
        atomic_json(root / "presentation_manifest.json", manifest)
    finally:
        plt.close(fig)
        ui["widget"].close()
    return manifest
