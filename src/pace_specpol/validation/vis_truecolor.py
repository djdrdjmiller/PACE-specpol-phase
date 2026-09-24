"""Visible OCI L1C TOA reflectance for presentation context, not phase retrieval."""

import json
import re
from contextlib import ExitStack
from pathlib import Path

import numpy as np
import xarray as xr

from pace_specpol import matching as mp
from pace_specpol.liquid_reference import sha256
from pace_specpol.workflow import digest, write_dataset


def presentation_inputs(pair_dir, companion_dir, reference_dir, stamp):
    """Verify the completed per-file provenance chain without rewriting manifests."""
    roots = [Path(p) for p in (pair_dir, companion_dir, reference_dir)]
    manifests = [json.loads((p / "manifest.json").read_text()) for p in roots]
    name = f"{stamp}.nc"
    for root, manifest in zip(roots, manifests):
        if not manifest.get("complete") or manifest.get("files") != [name]:
            raise ValueError(
                f"Expected a completed single-granule cache for {stamp}: {root}"
            )
    for i in (1, 2):
        manifest = manifests[i]
        if digest(manifest["settings"]) != manifest["fingerprint"]:
            raise ValueError(f"Manifest fingerprint mismatch: {roots[i]}")
        with xr.open_dataset(roots[i] / name, engine="h5netcdf") as ds:
            if ds.attrs.get("stage_fingerprint") != manifest["fingerprint"]:
                raise ValueError(f"Stage fingerprint mismatch: {roots[i] / name}")
            if ds.attrs.get("source_sha256") != sha256(roots[i - 1] / name):
                raise ValueError(f"Source data changed: {roots[i - 1] / name}")
    pairs = manifests[0]["pairs"]
    if len(pairs) != 1 or pairs[0]["stamp"] != stamp:
        raise ValueError(
            "The raw cache catalog does not match the presentation granule"
        )
    descriptor = pairs[0]["oci"]
    with xr.open_dataset(roots[0] / name, engine="h5netcdf") as raw:
        report = json.loads(raw.attrs["report_json"])
        if report["oci_product"] != descriptor["name"]:
            raise ValueError("The raw paired file and OCI catalog name disagree")
    return descriptor, manifests[0]["config"]


def read_rgb(
    path,
    *,
    fs=None,
    expected_product=None,
    wavelengths=(645.0, 555.0, 469.0),
    max_solar_zenith=75.0,
    max_sensor_zenith=65.0,
    max_view_offset_s=200.0,
    qc_policy="ignore",
):
    """Read three visible bands, selecting a single valid view per L1C grid cell.

    Same radiance-to-TOA conversion as matching.py. Full scene, without the HARP2
    cloud/LI or liquid-reference masks. No gas, Rayleigh, or aerosol correction.
    """
    if qc_policy not in ("ignore", "require_zero"):
        raise ValueError("qc_policy must be ignore or require_zero")
    targets = np.asarray(wavelengths, dtype=float)
    if targets.shape != (3,) or not np.all(np.isfinite(targets)):
        raise ValueError("Specify three finite RGB wavelengths in nm")
    with ExitStack() as stack:
        tree = mp._open_tree(stack, path, fs)
        product = str(tree.attrs.get("product_name", ""))
        if not re.fullmatch(mp.PATTERNS["oci"], product):
            raise ValueError("Expected refined OCI L1C V3 5 km data")
        if expected_product is not None and product != expected_product:
            raise ValueError(
                f"Wrong OCI granule: {product}; expected {expected_product}"
            )
        sensor = tree["sensor_views_bands"]
        wave = sensor["intensity_wavelength"].transpose(mp.VIEW, mp.BAND).values
        selected = [int(np.nanargmin(abs(wave[0] - w))) for w in targets]
        actual = wave[:, selected]
        if (
            len(set(selected)) != 3
            or not np.isfinite(actual).all()
            or np.any(abs(actual - targets) > 10)
        ):
            raise ValueError(
                "Three distinct visible RGB bands within 10 nm are required"
            )
        f0 = (
            sensor["intensity_f0"]
            .isel({mp.BAND: selected})
            .transpose(mp.VIEW, mp.BAND)
            .values
        )
        if not np.all(np.isfinite(f0) & (f0 > 0)):
            raise ValueError("Invalid visible-band solar irradiance")
        obs, geo = tree["observation_data"], tree["geolocation_data"]
        if "W" not in str(obs["i"].attrs.get("units")) or "sr" not in str(
            obs["i"].attrs.get("units")
        ):
            raise ValueError("Expected spectral radiance in observation_data/i")
        radiance = (
            obs["i"]
            .isel({mp.BAND: selected})
            .transpose(*mp.DIMS, mp.VIEW, mp.BAND)
            .values
        )
        qc = (
            obs["qc"]
            .isel({mp.BAND: selected})
            .transpose(*mp.DIMS, mp.VIEW, mp.BAND)
            .values
        )
        sza = geo["solar_zenith_angle"].transpose(*mp.DIMS, mp.VIEW).values
        vza = geo["sensor_zenith_angle"].transpose(*mp.DIMS, mp.VIEW).values
        offsets = (
            tree["bin_attributes"]["view_time_offsets"]
            .transpose(*mp.DIMS, mp.VIEW)
            .values
        )
        nobs = obs["number_of_observations"].transpose(*mp.DIMS, mp.VIEW).values
        distance = float(tree.attrs["sun_earth_distance"])
        if not 0.95 < distance < 1.05:
            raise ValueError("Unexpected Earth-Sun distance")
        with np.errstate(divide="ignore", invalid="ignore"):
            rho = (
                np.pi
                * radiance
                * distance**2
                / (f0 * np.cos(np.deg2rad(sza))[..., None])
            )
        valid = (
            np.all(np.isfinite(rho) & (rho >= 0) & (rho <= 3), axis=-1)
            & np.isfinite(sza)
            & (sza >= 0)
            & (sza <= max_solar_zenith)
            & np.isfinite(vza)
            & (vza >= 0)
            & (vza <= max_sensor_zenith)
            & np.isfinite(offsets)
            & (abs(offsets) <= max_view_offset_s)
            & np.isfinite(nobs)
            & (nobs > 0)
        )
        if qc_policy == "require_zero":
            valid &= np.all(np.isfinite(qc) & (qc == 0), axis=-1)
        view = np.argmin(np.where(valid, abs(offsets), np.inf), axis=-1)
        lat = geo["latitude"].transpose(*mp.DIMS).values
        lon = geo["longitude"].transpose(*mp.DIMS).values
        good = (
            valid.any(axis=-1)
            & np.isfinite(lat)
            & np.isfinite(lon)
            & (abs(lat) <= 90)
            & (abs(lon) <= 180)
        )
        iy, ix = np.nonzero(good)
        iv = view[iy, ix]
        if not iy.size:
            raise ValueError("No visible RGB pixels passed screening")
        result = xr.Dataset(
            {
                "latitude": ("pixel", lat[good]),
                "longitude": ("pixel", lon[good]),
                "reflectance": (
                    ("pixel", "channel"),
                    rho[iy, ix, iv].astype("float32"),
                ),
                "oci_view": ("pixel", iv.astype("int8")),
                "view_offset_s": ("pixel", offsets[iy, ix, iv]),
            },
            coords={"channel": ["red", "green", "blue"]},
        )
        result.attrs.update(
            product_name=product,
            time_coverage_start=str(tree.attrs.get("time_coverage_start", "")),
            time_coverage_end=str(tree.attrs.get("time_coverage_end", "")),
            actual_wavelengths_nm=json.dumps(actual.tolist()),
            radiometry="rho = pi * I * sun_earth_distance^2 / (F0 * cos(solar_zenith))",
            atmosphere="TOA; no atmospheric correction",
            view_selection="Smallest absolute nadir time offset with all RGB bands valid",
            population="Full visible L1C scene; no HARP2 or reference-validity mask",
            qc_policy=qc_policy,
            unknown_qc_pixels=int((~np.isfinite(qc[iy, ix, iv]).all(axis=-1)).sum()),
        )
        return result


def cache_rgb(descriptor, cache_path, match_config, *, local_source=None):
    """Cache only the visible reflectances; Earthdata/S3 is used on first AWS read."""
    options = dict(
        max_solar_zenith=match_config["max_solar_zenith"],
        max_sensor_zenith=match_config["max_sensor_zenith"],
        max_view_offset_s=match_config["max_view_offset_s"],
        qc_policy=match_config["oci_qc_policy"],
    )
    source = Path(local_source).expanduser() if local_source is not None else None
    settings = dict(
        schema=1,
        descriptor=descriptor,
        options=options,
        wavelengths_nm=[645.0, 555.0, 469.0],
        local_source_sha256=sha256(source) if source else None,
    )
    fingerprint = digest(settings)
    cache_path = Path(cache_path)
    if cache_path.exists():
        with xr.open_dataset(cache_path, engine="h5netcdf") as ds:
            if ds.attrs.get("rgb_fingerprint") != fingerprint:
                raise ValueError(
                    "RGB source/settings changed; choose a new RGB_CACHE path"
                )
            return ds.load()
    if source is not None:
        rgb = read_rgb(source, expected_product=descriptor["name"], **options)
    else:
        import earthaccess

        earthaccess.login()
        found = earthaccess.search_data(
            short_name="PACE_OCI_L1C_SCI",
            version="3",
            granule_name=descriptor["name"],
            count=-1,
        )
        exact = [
            g for g in found if mp._descriptor(g, "oci")["url"] == descriptor["url"]
        ]
        if len(exact) != 1:
            raise ValueError(
                "Could not identify exactly the cached OCI L1C source in CMR"
            )
        fs = earthaccess.get_s3_filesystem(results=exact)
        rgb = read_rgb(
            descriptor["url"], fs=fs, expected_product=descriptor["name"], **options
        )
    rgb.attrs.update(rgb_fingerprint=fingerprint, rgb_settings=json.dumps(settings))
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    write_dataset(rgb, cache_path)
    return rgb


def rgb_on_map(rgb, ax, shape, *, black=0.0, white=1.0, gamma=2.2, max_distance_km=4.0):
    """Nearest 5 km L1C centers on the dashboard raster, with a finite gap cutoff."""
    import cartopy.crs as ccrs
    from scipy.spatial import cKDTree

    if not np.isfinite([black, white, gamma, max_distance_km]).all() or not (
        0 <= black < white and gamma > 0 and 0 < max_distance_km <= 20
    ):
        raise ValueError(
            "Require 0 <= black < white, gamma > 0 and 0 < distance <= 20 km"
        )
    height, width = shape
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xx, yy = np.meshgrid(
        np.linspace(xmin, xmax, width, endpoint=False) + (xmax - xmin) / (2 * width),
        np.linspace(ymin, ymax, height, endpoint=False) + (ymax - ymin) / (2 * height),
    )
    points = ccrs.PlateCarree().transform_points(ax.projection, xx, yy)
    lon, lat = points[..., 0], points[..., 1]
    valid = np.isfinite(lat) & np.isfinite(lon) & (abs(lat) <= 90)

    def xyz(lat, lon):
        lat, lon = np.deg2rad(lat), np.deg2rad(lon)
        return np.column_stack(
            [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]
        )

    tree = cKDTree(xyz(rgb.latitude.values, rgb.longitude.values))
    distance, ids = tree.query(
        xyz(lat[valid], lon[valid]),
        distance_upper_bound=2 * np.sin(max_distance_km / (2 * 6371.0)),
    )
    hit = np.isfinite(distance)
    rgba = np.zeros((height, width, 4), dtype="float32")
    values = np.zeros((valid.sum(), 4), dtype="float32")
    values[hit, :3] = np.clip(
        (rgb.reflectance.values[ids[hit]] - black) / (white - black), 0, 1
    ) ** (1 / gamma)
    values[hit, 3] = 1
    rgba[valid] = values
    if not hit.any():
        raise ValueError("No RGB pixels fall within the dashboard framing")
    return rgba
