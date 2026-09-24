"""Radiometry, cache integrity and slide registration checks without NASA access."""

import json
from dataclasses import asdict

import numpy as np
import pytest
import xarray as xr

from pace_specpol import matching as mp
from pace_specpol.liquid_reference import sha256
from pace_specpol.validation.vis_truecolor import (
    cache_rgb,
    presentation_inputs,
    read_rgb,
)
from pace_specpol.workflow import digest

from .test_matching import make_pair


def visible_file(tmp_path):
    _, path, tree = make_pair(tmp_path)
    wave = np.array([[645.0, 555.0, 470.0]] * 2)
    f0 = np.array([[100.0, 200.0, 300.0]] * 2)
    tree["sensor_views_bands"] = xr.DataTree(
        xr.Dataset(
            {
                "intensity_wavelength": ((mp.VIEW, mp.BAND), wave),
                "intensity_f0": ((mp.VIEW, mp.BAND), f0),
            }
        )
    )
    rho = np.broadcast_to([0.4, 0.3, 0.2], (2, 2, 2, 3)).copy()
    # A pixel whose channels are valid only in different views must be absent.
    rho[0, 0, 0, 0] = np.nan
    rho[0, 0, 1, 1] = np.nan
    # Another pixel is valid in view 1 only.
    rho[0, 1, 0, 0] = np.nan
    radiance = rho * f0 * np.cos(np.deg2rad(30)) / np.pi
    tree["observation_data"] = xr.DataTree(
        xr.Dataset(
            {
                "i": (
                    mp.DIMS + (mp.VIEW, mp.BAND),
                    radiance,
                    {"units": "W m^-2 sr^-1 um^-1"},
                ),
                "qc": (mp.DIMS + (mp.VIEW, mp.BAND), np.full_like(radiance, np.nan)),
                "number_of_observations": (mp.DIMS + (mp.VIEW,), np.ones((2, 2, 2))),
            }
        )
    )
    tree.to_netcdf(path, engine="h5netcdf")
    return path, tree


def test_visible_radiometry_and_view_selection(tmp_path):
    path, tree = visible_file(tmp_path)
    rgb = read_rgb(path)
    assert rgb.sizes["pixel"] == 3
    np.testing.assert_allclose(rgb.reflectance, [[0.4, 0.3, 0.2]] * 3, rtol=1e-6)
    np.testing.assert_array_equal(rgb.oci_view, [1, 0, 0])
    assert rgb.attrs["unknown_qc_pixels"] == 3
    with pytest.raises(ValueError, match="No visible RGB"):
        read_rgb(path, qc_policy="require_zero")
    with pytest.raises(ValueError, match="Wrong OCI granule"):
        read_rgb(path, expected_product="wrong.nc")
    tree["sensor_views_bands"]["intensity_wavelength"].values[:, 0] = 1615
    tree.to_netcdf(path, engine="h5netcdf")
    with pytest.raises(ValueError, match="distinct visible"):
        read_rgb(path)


def test_rgb_cache_reuse_and_changed_settings(tmp_path, monkeypatch):
    import pace_specpol.validation.vis_truecolor as vis

    path, tree = visible_file(tmp_path)
    desc = {"name": tree.attrs["product_name"], "url": "s3://test/exact.nc"}
    config = asdict(mp.MatchConfig())
    target = tmp_path / "rgb.nc"
    first = cache_rgb(desc, target, config, local_source=path)

    def unexpected(*args, **kwargs):
        raise AssertionError("A cached RGB must not read the source again")

    monkeypatch.setattr(vis, "read_rgb", unexpected)
    second = cache_rgb(desc, target, config, local_source=path)
    xr.testing.assert_equal(first, second)
    config["oci_qc_policy"] = "require_zero"
    with pytest.raises(ValueError, match="settings changed"):
        cache_rgb(desc, target, config, local_source=path)


def test_presentation_provenance_rejects_changed_source(tmp_path):
    stamp = "20240610T154205"
    name = stamp + ".nc"
    product = f"PACE_OCI.{stamp}.L1C.V3.5km.nc"
    roots = [tmp_path / k for k in ("pairs", "companions", "references")]
    for root in roots:
        root.mkdir()
    xr.Dataset(attrs={"report_json": json.dumps({"oci_product": product})}).to_netcdf(
        roots[0] / name
    )
    manifest = dict(
        complete=True,
        files=[name],
        pairs=[dict(stamp=stamp, oci={"name": product})],
        config={},
    )
    (roots[0] / "manifest.json").write_text(json.dumps(manifest))
    for i in (1, 2):
        settings = {"stage": i}
        fingerprint = digest(settings)
        xr.Dataset(
            attrs={
                "source_sha256": sha256(roots[i - 1] / name),
                "stage_fingerprint": fingerprint,
            }
        ).to_netcdf(roots[i] / name)
        (roots[i] / "manifest.json").write_text(
            json.dumps(
                dict(
                    complete=True,
                    files=[name],
                    settings=settings,
                    fingerprint=fingerprint,
                )
            )
        )
    desc, _ = presentation_inputs(*roots, stamp)
    assert desc["name"] == product
    # Manifest bookkeeping changes alone do not invalidate identical source data.
    manifest["note"] = "Relocated cache"
    (roots[0] / "manifest.json").write_text(json.dumps(manifest))
    presentation_inputs(*roots, stamp)
    xr.Dataset(attrs={"changed": "data"}).to_netcdf(roots[0] / name)
    with pytest.raises(ValueError, match="Source data changed"):
        presentation_inputs(*roots, stamp)


def test_aligned_exports_and_gap_mask(tmp_path, monkeypatch):
    import matplotlib

    matplotlib.use("Agg")
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    from cartopy.mpl.geoaxes import GeoAxes
    from PIL import Image

    from pace_specpol.validation.aggregation import IndexConfig, build_index
    from pace_specpol.validation.vis_presentation import export_presentation
    from pace_specpol.validation.vis_truecolor import rgb_on_map

    # Avoid downloading Natural Earth in unit tests; visual QA uses real coastlines.
    monkeypatch.setattr(GeoAxes, "coastlines", lambda *a, **kw: None)
    source = tmp_path / "source"
    source.mkdir()
    lat, lon = np.meshgrid(np.linspace(65, 80, 12), np.linspace(-65, -35, 12))
    xr.Dataset(
        {
            "latitude": ("sample", lat.ravel()),
            "longitude": ("sample", lon.ravel()),
            "liquid_index": ("sample", np.linspace(0, 4, lat.size)),
            "raw_ratio": ("sample", np.linspace(0.7, 2.0, lat.size)),
        }
    ).to_netcdf(source / "samples.nc", engine="h5netcdf")
    (source / "manifest.json").write_text(
        json.dumps(dict(complete=True, files=["samples.nc"], config={}))
    )
    index = build_index(source, config=IndexConfig(resolution=1))
    rgb = xr.Dataset(
        {
            "latitude": ("pixel", lat.ravel()),
            "longitude": ("pixel", lon.ravel()),
            "reflectance": (
                ("pixel", "channel"),
                np.tile([0.4, 0.3, 0.2], (lat.size, 1)),
            ),
        }
    )
    # Explicit distance cutoff keeps the far side of the map unfilled.
    fig = plt.figure()
    ax = fig.add_subplot(projection=ccrs.NorthPolarStereo(central_longitude=-45))
    ax.set_extent([-85, -10, 58, 87], crs=ccrs.PlateCarree())
    rgba = rgb_on_map(rgb, ax, (200, 400))
    assert np.any(rgba[..., 3] == 1) and np.any(rgba[..., 3] == 0)
    np.testing.assert_allclose(
        rgba[rgba[..., 3] == 1, :3][0],
        np.array([0.4, 0.3, 0.2]) ** (1 / 2.2),
        rtol=1e-6,
    )
    plt.close(fig)
    out = tmp_path / "exports"
    result = export_presentation(
        index,
        out,
        rgb=rgb,
        dpi=60,
        preview_width=300,
        map_extent=(-85, -10, 58, 87),
        map_framing="granule",
    )
    assert result["complete"]
    assert result["canvas_pixels"] == [800, 450]
    assert len(result["files"]) == 8
    crop_box = result["map_crop_box_pixels"]
    for item in result["files"]:
        assert sha256(out / item["name"]) == item["sha256"]
        with Image.open(out / item["name"]) as im:
            assert list(im.size) == (
                result["canvas_pixels"]
                if item["name"].startswith("dashboard_")
                else result["map_pixels"]
            )
    with (
        Image.open(out / "dashboard_truecolor_context.png") as full,
        Image.open(out / "truecolor_map.png") as crop,
    ):
        np.testing.assert_array_equal(np.asarray(full.crop(crop_box)), np.asarray(crop))
    a = np.asarray(Image.open(out / "dashboard_ratio_1p10.png"))
    b = np.asarray(Image.open(out / "dashboard_ratio_1p45.png"))
    # Histogram density and lower-left LI map are fixed; title/phase/ratio colors vary.
    np.testing.assert_array_equal(a[280:385, 60:360], b[280:385, 60:360])
    assert np.any(a != b)
    with pytest.raises(FileExistsError):
        export_presentation(index, out)
