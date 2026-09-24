"""Discover, align, and cache matched HARP2/OCI observations.

HARP2 GPC V4.0 and OCI L1C V3 are paired by timestamp, then checked for identical
geolocation and nadir timing. This stage preserves TOA reflectances; liquid
reference reconstruction and atmospheric correction are separate stages.
"""

from contextlib import ExitStack
from dataclasses import dataclass, asdict
from datetime import date, timedelta
from pathlib import Path
import hashlib
import json
import re
import warnings

import numpy as np
import xarray as xr
import dask

CODE_VERSION = "1.1"
LI = "cloud_liquid_index"
DIMS = ("bins_along_track", "bins_across_track")
VIEW = "number_of_views"
BAND = "intensity_bands_per_view"
COLLECTIONS = {
    "harp": ("PACE_HARP2_L2_CLOUD_GPC", "4.0"),
    "oci": ("PACE_OCI_L1C_SCI", "3"),
}
PATTERNS = {
    "harp": r"PACE_HARP2\.(\d{8}T\d{6})\.L2\.CLOUD_GPC\.V4_0\.nc",
    "oci": r"PACE_OCI\.(\d{8}T\d{6})\.L1C\.V3\.5km\.nc",
}


def _digest(obj):
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, default=str).encode()
    ).hexdigest()


def _write_json(path, obj):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, default=str))
    tmp.replace(path)


@dataclass(frozen=True)
class MatchConfig:
    start_date: str = "2024-09-03"
    end_date: str = "2024-09-05"  # inclusive UTC
    cache_dir: str = "mixed_phase_cache"
    workers: int = 2
    batch_size: int = 4
    geo_tolerance_deg: float = 1e-4
    nadir_time_tolerance_s: float = 1.0
    max_view_offset_s: float = 200.0
    max_solar_zenith: float = 75.0
    max_sensor_zenith: float = 65.0
    min_rho1615: float = 0.01
    max_reflectance: float = 3.0
    cloud_reject_bits: int = (
        1  # reject cloud-presence discrepancy, not liquid retrieval failures
    )
    min_cloud_optical_thickness: float | None = None
    oci_qc_policy: str = (
        "ignore"  # tested OCI qc is entirely fill; 'require_zero' is available
    )
    allow_unpaired: bool = False

    def validate(self):
        if date.fromisoformat(self.end_date) < date.fromisoformat(self.start_date):
            raise ValueError("end_date must not precede start_date")
        if self.oci_qc_policy not in ("ignore", "require_zero"):
            raise ValueError("oci_qc_policy must be ignore or require_zero")
        if self.workers < 1 or self.batch_size < 1:
            raise ValueError("workers and batch_size must be positive")
        if not 0 <= self.cloud_reject_bits <= 255:
            raise ValueError("cloud_reject_bits must be 0..255")
        if not 0 < self.min_rho1615 < self.max_reflectance:
            raise ValueError("invalid reflectance bounds")
        if not 0 < self.max_solar_zenith < 90 or not 0 < self.max_sensor_zenith < 90:
            raise ValueError("zenith limits must be between 0 and 90")
        if (
            self.geo_tolerance_deg <= 0
            or self.nadir_time_tolerance_s <= 0
            or self.max_view_offset_s <= 0
        ):
            raise ValueError("alignment/time tolerances must be positive")
        if (
            self.min_cloud_optical_thickness is not None
            and self.min_cloud_optical_thickness < 0
        ):
            raise ValueError("minimum optical thickness must be nonnegative")

    @property
    def bounds(self):
        self.validate()
        return (
            np.datetime64(self.start_date, "ns"),
            np.datetime64(date.fromisoformat(self.end_date) + timedelta(days=1), "ns"),
        )

    @property
    def fingerprint(self):
        settings = asdict(self)
        for key in ("cache_dir", "workers", "batch_size"):
            settings.pop(key)
        return _digest({"code": CODE_VERSION, "settings": settings})


def _descriptor(g, kind):
    u = g["umm"]
    short, version = COLLECTIONS[kind]
    c = u["CollectionReference"]
    if c["ShortName"] != short or str(c["Version"]) != version:
        raise ValueError(f"Unexpected {kind} collection/version")
    urls = [
        x["URL"]
        for x in u.get("RelatedUrls", [])
        if x.get("Type") == "GET DATA VIA DIRECT ACCESS"
        and x["URL"].startswith("s3://")
        and x["URL"].endswith(".nc")
    ]
    if len(urls) != 1:
        raise ValueError(f"Expected one clean S3 link, found {len(urls)}")
    name = urls[0].rsplit("/", 1)[-1]
    match = re.fullmatch(PATTERNS[kind], name)
    if match is None:
        raise ValueError(f"Unexpected non-NRT filename: {name}")
    return {
        "name": name,
        "stamp": match[1],
        "url": urls[0],
        "concept_id": g.get("meta", {}).get("concept-id"),
        "revision_id": g.get("meta", {}).get("revision-id"),
    }


def discover_pairs(cfg):
    """No row-order pairing and no fallback to a different product version."""
    import earthaccess

    cfg.validate()
    lower, upper = cfg.bounds
    temporal = (str(lower) + "Z", str(upper - np.timedelta64(1, "ns")) + "Z")
    catalogs, raw = {}, {}
    for kind, (short, version) in COLLECTIONS.items():
        raw[kind] = earthaccess.search_data(
            short_name=short, version=version, temporal=temporal, count=-1
        )
        catalog = {}
        for g in raw[kind]:
            entry = _descriptor(g, kind)
            if entry["stamp"] in catalog and catalog[entry["stamp"]] != entry:
                raise ValueError(f"Ambiguous {kind} granules at {entry['stamp']}")
            catalog[entry["stamp"]] = entry
        catalogs[kind] = catalog
    if not catalogs["harp"]:
        raise RuntimeError("No HARP2 V4.0 non-NRT granules found")
    missing = sorted(set(catalogs["harp"]) - set(catalogs["oci"]))
    if missing and not cfg.allow_unpaired:
        raise RuntimeError(
            f"{len(missing)} HARP2 granules lack exact OCI L1C partners. "
            "Review coverage before setting allow_unpaired=True."
        )
    pairs = [
        {"stamp": t, "harp": catalogs["harp"][t], "oci": catalogs["oci"][t]}
        for t in sorted(set(catalogs["harp"]) & set(catalogs["oci"]))
    ]
    if not pairs:
        raise RuntimeError("No matched timestamps")
    report = {
        "harp_granules": len(catalogs["harp"]),
        "oci_granules": len(catalogs["oci"]),
        "pairs": len(pairs),
        "unpaired_harp_stamps": missing,
        "unused_oci_stamps": sorted(set(catalogs["oci"]) - set(catalogs["harp"])),
    }
    print(report)
    return pairs, raw, report


def _open_tree(stack, path, fs):
    if fs is not None:
        path = stack.enter_context(
            fs.open(
                path,
                "rb",
                block_size=8 * 2**20,
                cache_type="blockcache",
                cache_options={"maxblocks": 16},
            )
        )
    return stack.enter_context(
        xr.open_datatree(path, engine="h5netcdf", mask_and_scale=True)
    )


def _harp_times(tree):
    t = tree["bin_attributes"]["nadir_view_time"]
    if np.issubdtype(t.dtype, np.datetime64):
        return t.values.astype("datetime64[ns]")
    if t.attrs.get("units") != "seconds from UTC midnight":
        raise ValueError("Unrecognized HARP2 time units")
    start = np.datetime64(str(tree.attrs["time_coverage_start"]).rstrip("Z"), "ns")
    midnight = start.astype("datetime64[D]").astype("datetime64[ns]")
    s0 = (start - midnight) / np.timedelta64(1, "s")
    seconds = t.values.astype(float)
    seconds += 86400 * np.round((s0 - seconds) / 86400)
    ok = np.isfinite(seconds)
    times = midnight + (np.where(ok, seconds, 0) * 1e9).astype("timedelta64[ns]")
    times[~ok] = np.datetime64("NaT", "ns")
    return times


def _check_grid(h, o, cfg):
    hg, og = h["geolocation_data"], o["geolocation_data"]
    valid = None
    for name in ("latitude", "longitude"):
        a, b = hg[name], og[name]
        if a.dims != DIMS or b.dims != DIMS or a.shape != b.shape:
            raise ValueError(
                "Different HARP2/OCI grids; explicit regridding is required"
            )
        av, bv = a.values, b.values
        finite = np.isfinite(av) & np.isfinite(bv)
        if not finite.any():
            raise ValueError("No common finite geolocation")
        diff = abs(av - bv) if name == "latitude" else abs((av - bv + 180) % 360 - 180)
        if np.max(diff[finite]) > cfg.geo_tolerance_deg:
            raise ValueError(
                f"{name} grids differ beyond tolerance; refusing positional pairing"
            )
        bound = 90 if name == "latitude" else 180
        valid = (
            finite & (abs(av) <= bound)
            if valid is None
            else valid & finite & (abs(av) <= bound)
        )
    ht = _harp_times(h)
    ot = o["bin_attributes"]["nadir_view_time"]
    if ot.dims != (DIMS[0],) or not np.issubdtype(ot.dtype, np.datetime64):
        raise ValueError("Expected CF-decoded OCI nadir times on along-track rows")
    ot = ot.values.astype("datetime64[ns]")
    if ht.shape != ot.shape:
        raise ValueError("Different row timing shapes")
    signed_delta = (ht - ot) / np.timedelta64(1, "s")
    good = np.isfinite(signed_delta)
    correction_s = 0
    if good.any() and np.max(abs(signed_delta[good])) > cfg.nadir_time_tolerance_s:
        # Some OCI V3 files after midnight retain the preceding day's CF epoch
        # while their seconds reset to zero. Verified for 20240905T000311 and
        # 20240905T000811. Repair only a uniform +/- one-day error supported by
        # BOTH the HARP2 times and OCI's own global coverage interval.
        shift_days = int(np.rint(np.median(signed_delta[good]) / 86400))
        if abs(shift_days) == 1:
            shift_s = shift_days * 86400
            candidate = ot + np.timedelta64(shift_s, "s")
            start = np.datetime64(str(o.attrs["time_coverage_start"]).rstrip("Z"), "ns")
            end = np.datetime64(str(o.attrs["time_coverage_end"]).rstrip("Z"), "ns")
            tolerance = np.timedelta64(round(cfg.nadir_time_tolerance_s * 1e9), "ns")
            uniform = np.all(
                abs(signed_delta[good] - shift_s) <= cfg.nadir_time_tolerance_s
            )
            inside = np.all(
                (candidate[good] >= start - tolerance)
                & (candidate[good] <= end + tolerance)
            )
            if uniform and inside:
                ot = candidate
                correction_s = shift_s
                warnings.warn(
                    "Corrected a uniform one-day OCI CF-time epoch error; "
                    "the offset and source units are recorded in each pair report.",
                    RuntimeWarning,
                )
    delta = abs((ot - ht) / np.timedelta64(1, "s"))
    good = np.isfinite(delta)
    if not good.any() or np.max(delta[good]) > cfg.nadir_time_tolerance_s:
        largest = float(np.max(delta[good])) if good.any() else None
        raise ValueError(
            f"HARP2/OCI nadir-time mismatch in {o.attrs.get('product_name')}; "
            f"maximum difference={largest} seconds. No safe day-offset correction found."
        )
    # In-memory audit metadata only; source files are never modified.
    o.attrs["_reader_nadir_time_correction_seconds"] = correction_s
    return valid & good[:, None], ht, ot


def _two_bands(tree):
    sensor = tree["sensor_views_bands"]
    wave = sensor["intensity_wavelength"].transpose(VIEW, BAND).values
    selected = [int(np.nanargmin(abs(wave[0] - x))) for x in (1615, 2260)]
    for i, target in zip(selected, (1615, 2260)):
        if not np.all(np.isfinite(wave[:, i])) or np.any(abs(wave[:, i] - target) > 15):
            raise ValueError(
                "Expected OCI ~1615/~2260 nm channels not found for every view"
            )
    if selected[0] == selected[1]:
        raise ValueError("SWIR channel selection is ambiguous")
    f0 = (
        sensor["intensity_f0"]
        .isel({BAND: selected})
        .transpose(VIEW, BAND)
        .values.astype("float64")
    )
    if not np.all(np.isfinite(f0) & (f0 > 0)):
        raise ValueError("Invalid solar irradiances")
    i = tree["observation_data"]["i"]
    if "W" not in str(i.attrs.get("units")) or "sr" not in str(i.attrs.get("units")):
        raise ValueError("Expected spectral radiance; do not treat i as reflectance")
    return selected, wave[:, selected], f0


def match_pair(harp_path, oci_path, cfg, fs_harp=None, fs_oci=None):
    """Return one sample per HARP2 cell with both SWIR bands from the SAME OCI view."""
    cfg.validate()
    with ExitStack() as stack:
        h, o = (
            _open_tree(stack, harp_path, fs_harp),
            _open_tree(stack, oci_path, fs_oci),
        )
        stamps = []
        for tree, kind in ((h, "harp"), (o, "oci")):
            match = re.fullmatch(
                PATTERNS[kind], str(tree.attrs.get("product_name", ""))
            )
            if match is None:
                raise ValueError(f"Unexpected {kind} product_name")
            stamps.append(match[1])
        if stamps[0] != stamps[1]:
            raise ValueError("Different granule timestamps")
        grid_ok, ht, ot = _check_grid(h, o, cfg)
        hp = h["geophysical_data"]
        li = hp[LI]
        if li.dims != DIMS or li.shape != grid_ok.shape:
            raise ValueError("LI is not on the verified grid")
        liv = li.values
        good = grid_ok & np.isfinite(liv) & (liv >= -10) & (liv <= 10)
        quality = hp["cloud_quality"].values
        if quality.shape != li.shape:
            raise ValueError("Cloud QA shape mismatch")
        qfinite = np.isfinite(quality)
        qi = np.where(qfinite, quality, 0).astype("uint16")
        good &= qfinite & ((qi & cfg.cloud_reject_bits) == 0)
        cot = hp["cloud_optical_thickness"].values
        if cfg.min_cloud_optical_thickness is not None:
            good &= np.isfinite(cot) & (cot >= cfg.min_cloud_optical_thickness)
        lower, upper = cfg.bounds
        good &= ((ht >= lower) & (ht < upper))[:, None]
        selected, wave, f0 = _two_bands(o)
        obs, geo = o["observation_data"], o["geolocation_data"]
        # Selection is explicit, but contiguous spectral storage can still require
        # substantial underlying S3 I/O. Read both bands in one selection.
        intensity = (
            obs["i"]
            .isel({BAND: selected})
            .transpose(*DIMS, VIEW, BAND)
            .values.astype("float64")
        )
        qc = obs["qc"].isel({BAND: selected}).transpose(*DIMS, VIEW, BAND).values

        def geometry(name):
            return geo[name].transpose(*DIMS, VIEW).values

        sza = geometry("solar_zenith_angle")
        vza = geometry("sensor_zenith_angle")
        raa = abs(
            (geometry("sensor_azimuth_angle") - geometry("solar_azimuth_angle") + 180)
            % 360
            - 180
        )
        offset = o["bin_attributes"]["view_time_offsets"].transpose(*DIMS, VIEW).values
        nobs = obs["number_of_observations"].transpose(*DIMS, VIEW).values
        distance = float(o.attrs["sun_earth_distance"])
        if not 0.95 < distance < 1.05:
            raise ValueError("Unexpected Earth-Sun distance")
        with np.errstate(divide="ignore", invalid="ignore"):
            rho = (
                np.pi
                * intensity
                * distance**2
                / (f0[None, None, :, :] * np.cos(np.deg2rad(sza))[..., None])
            )
        vok = (
            np.all(
                np.isfinite(intensity) & (intensity > 0) & (intensity <= 999), axis=-1
            )
            & np.all(
                np.isfinite(rho) & (rho > 0) & (rho <= cfg.max_reflectance), axis=-1
            )
            & (rho[..., 0] >= cfg.min_rho1615)
            & np.isfinite(sza)
            & (sza >= 0)
            & (sza <= cfg.max_solar_zenith)
            & np.isfinite(vza)
            & (vza >= 0)
            & (vza <= cfg.max_sensor_zenith)
            & np.isfinite(raa)
            & np.isfinite(offset)
            & (abs(offset) <= cfg.max_view_offset_s)
            & np.isfinite(nobs)
            & (nobs > 0)
        )
        qc_known = np.all(np.isfinite(qc), axis=-1)
        qc_good = qc_known & np.all(qc == 0, axis=-1)
        if cfg.oci_qc_policy == "require_zero":
            vok &= qc_good
        off_ns = (np.where(np.isfinite(offset), offset, 0) * 1e9).astype(
            "timedelta64[ns]"
        )
        otime = ot[:, None, None] - off_ns
        vok &= (otime >= lower) & (otime < upper)
        # Prefer smallest temporal separation if both tilt views are valid.
        view = np.argmin(np.where(vok, abs(offset), np.inf), axis=-1)
        good &= vok.any(axis=-1)
        iy, ix = np.nonzero(good)
        iv = view[iy, ix]

        def pick(array):
            return array[iy, ix, iv]

        p = pick(rho)
        fields = {
            "latitude": h["geolocation_data"]["latitude"].values[good],
            "longitude": h["geolocation_data"]["longitude"].values[good],
            "liquid_index": liv[good],
            "rho1615": p[:, 0],
            "rho2260": p[:, 1],
            "raw_ratio": p[:, 1] / p[:, 0],
            "sza": pick(sza),
            "vza": pick(vza),
            "raa": pick(raa),
            "cot": cot[good],
            "re": hp["cloud_bow_droplet_effective_radius"].values[good],
            "view_offset_s": pick(offset),
            "oci_view": iv.astype("int8"),
            "oci_qc_known": pick(qc_known).astype("int8"),
            "oci_qc_zero": pick(qc_good).astype("int8"),
            "harp_row": iy.astype("int32"),
            "harp_column": ix.astype("int32"),
        }
        ds = xr.Dataset(
            {
                k: (
                    "sample",
                    np.asarray(v, dtype="float32")
                    if np.asarray(v).dtype.kind == "f"
                    else v,
                )
                for k, v in fields.items()
            }
        )
        ds["harp_nadir_time"] = (
            "sample",
            np.broadcast_to(ht[:, None], good.shape)[good],
        )
        ds["oci_observation_time"] = ("sample", pick(otime))
        report = {
            "harp_product": h.attrs["product_name"],
            "oci_product": o.attrs["product_name"],
            "matched_samples": int(good.sum()),
            "finite_li_samples": int(np.isfinite(liv).sum()),
            "oci_qc_finite_band_values": int(np.isfinite(qc).sum()),
            "selected_wavelengths_nm": wave.tolist(),
            "selected_f0": f0.tolist(),
            "oci_qc_policy": cfg.oci_qc_policy,
            "oci_nadir_time_correction_seconds": o.attrs.get(
                "_reader_nadir_time_correction_seconds", 0
            ),
            "oci_nadir_time_source_units": str(
                o["bin_attributes"]["nadir_view_time"].encoding.get("units", "")
            ),
            "time_reader_revision": "2026-09-21-guarded-day-epoch-repair",
            "matched_samples_unknown_oci_qc": int((~pick(qc_known)).sum()),
            "alignment": "same timestamp; verified geolocation and nadir times; surface-grid pairing",
        }
        ds.attrs.update(
            report_json=json.dumps(report),
            match_config=json.dumps(asdict(cfg)),
            extraction_fingerprint=cfg.fingerprint,
            code_version=CODE_VERSION,
            radiometry="TOA reflectance = pi * I * earth_sun_distance^2 / (F0 * cos(sza))",
            ratio_definition="rho2260 / rho1615; no liquid normalization yet",
            cloud_parallax_correction="not applied",
            view_time_offset_convention="nadir_minus_observation",
        )
        return ds, report


def cache_pairs(pairs, raw, cfg, discovery_report=None):
    """Resumable per-pair NetCDF extraction; fail on any unread/changed pair."""
    import earthaccess
    import time

    cfg.validate()
    out = Path(cfg.cache_dir)
    out.mkdir(parents=True, exist_ok=True)
    manifest = {
        "code_version": CODE_VERSION,
        "config": asdict(cfg),
        "fingerprint": cfg.fingerprint,
        "pairs": pairs,
        "discovery": discovery_report,
        "complete": False,
        "scope": f"{cfg.start_date}–{cfg.end_date} UTC | {len(pairs)} matched granule pairs",
    }
    manifest_path = out / "manifest.json"
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text())
        if old["fingerprint"] != cfg.fingerprint or old["pairs"] != pairs:
            raise ValueError(
                "Cache configuration/catalog changed; choose a new cache_dir"
            )
    _write_json(manifest_path, manifest)
    records = []
    paths = []
    last_auth, hfs, ofs = -np.inf, None, None

    def extract(pair):
        target = out / (pair["stamp"] + ".nc")
        if target.exists():
            with xr.open_dataset(target, engine="h5netcdf") as saved:
                if saved.attrs.get("extraction_fingerprint") != cfg.fingerprint:
                    raise ValueError("Stale pair cache")
                return str(target), json.loads(saved.attrs["report_json"])
        ds, report = match_pair(pair["harp"]["url"], pair["oci"]["url"], cfg, hfs, ofs)
        tmp = target.with_suffix(".nc.tmp")
        ds.to_netcdf(
            tmp,
            engine="h5netcdf",
            encoding={k: {"zlib": True, "complevel": 3} for k in ds.data_vars},
        )
        tmp.replace(target)
        return str(target), report

    for i in range(0, len(pairs), cfg.batch_size):
        if time.monotonic() - last_auth > 2400:
            hfs = earthaccess.get_s3_filesystem(results=raw["harp"])
            ofs = earthaccess.get_s3_filesystem(results=raw["oci"])
            last_auth = time.monotonic()
        tasks = [dask.delayed(extract)(p) for p in pairs[i : i + cfg.batch_size]]
        result = dask.compute(*tasks, scheduler="threads", num_workers=cfg.workers)
        for path, report in result:
            paths.append(path)
            records.append(report)
        print(
            f"{len(paths)}/{len(pairs)} pairs cached; {sum(r['matched_samples'] for r in records):,} samples",
            flush=True,
        )
    manifest.update(complete=True, records=records, files=[Path(p).name for p in paths])
    _write_json(manifest_path, manifest)
    if not sum(r["matched_samples"] for r in records):
        raise RuntimeError(
            "No samples passed matching/screening; inspect QC and geometry settings"
        )
    return paths
