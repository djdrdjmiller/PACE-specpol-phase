"""OCI CLD companions and original OB.DAAC meteorology for cached HARP2/OCI pairs.

CLD is native ~1 km: nearest native pixel to each ~5 km L1C sample, with time/
distance checks. This is NOT NASA footprint aggregation or an L2 reretrieval.
"""

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
import json
import re
import shutil
import time
from urllib.parse import quote
import numpy as np
import xarray as xr
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator
import requests
from pace_specpol.liquid_reference import sha256


@dataclass(frozen=True)
class CompanionConfig:
    cloud_version: str = "3.1"
    max_distance_km: float = 3.0
    max_time_difference_s: float = 10.0
    ancillary_cache: str = "oci_ancillary_downloads"
    max_download_cache_mb: int = 512
    minimum_free_disk_mb: int = 512
    transport: str = "https"  # https (public cloud redirect) | earthaccess (AWS S3)


def atomic_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, default=str))
    tmp.replace(path)


def discover_cloud(start_date, end_date, version="3.1"):
    """Exact refined collection, paginated CMR; do not infer filenames from L1C."""
    end = (datetime.fromisoformat(end_date) + timedelta(days=1)).strftime("%Y-%m-%d")
    params = {
        "short_name": "PACE_OCI_L2_CLOUD",
        "version": version,
        "temporal": f"{start_date}T00:00:00Z,{end}T00:00:00Z",
        "page_size": 2000,
    }
    result = []
    page = 1
    while True:
        params["page_num"] = page
        r = requests.get(
            "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
            params=params,
            timeout=90,
        )
        r.raise_for_status()
        items = r.json()["items"]
        for item in items:
            u = item["umm"]
            c = u["CollectionReference"]
            if c["ShortName"] != "PACE_OCI_L2_CLOUD" or c["Version"] != version:
                raise ValueError("Wrong collection")
            urls = [
                v["URL"]
                for v in u["RelatedUrls"]
                if v["URL"].startswith("s3://") and v["URL"].endswith(".nc")
            ]
            if len(urls) != 1:
                raise ValueError("Expected one S3 granule URL")
            ext = u["TemporalExtent"]["RangeDateTime"]
            result.append(
                {
                    "url": urls[0],
                    "start": ext["BeginningDateTime"],
                    "end": ext["EndingDateTime"],
                    "concept_id": item["meta"]["concept-id"],
                }
            )
        if len(items) < 2000:
            break
        page += 1
    return sorted(result, key=lambda x: x["start"])


@contextmanager
def cloud_file(url, transport="https"):
    """Bounded remote reads, no full L1B/CLD download cache."""
    import h5netcdf

    if transport == "earthaccess":
        import earthaccess

        f = earthaccess.open(
            [url],
            provider="OB_CLOUD",
            show_progress=False,
            open_kwargs={
                "block_size": 2**20,
                "cache_type": "blockcache",
                "cache_options": {"maxblocks": 16},
            },
        )[0]
    elif transport == "https":
        import fsspec

        prefix = "s3://ob-cumulus-prod-public/"
        if not url.startswith(prefix):
            raise ValueError("Unexpected cloud bucket")
        public = url.replace(
            prefix, "https://obdaac-tea.earthdatacloud.nasa.gov/ob-cumulus-prod-public/"
        )
        with requests.get(
            public, headers={"Range": "bytes=0-15"}, timeout=90, stream=True
        ) as r:
            r.raise_for_status()
            resolved = r.url
        f = fsspec.open(
            resolved,
            block_size=2**20,
            cache_type="blockcache",
            cache_options={"maxblocks": 16},
        ).open()
    else:
        raise ValueError("transport must be https or earthaccess")
    try:
        with h5netcdf.File(f, "r") as nc:
            yield nc
    finally:
        f.close()


def read_valid(var):
    a = np.asarray(var[:], dtype=float)
    attrs = var.attrs
    if "_FillValue" in attrs:
        a[a == attrs["_FillValue"]] = np.nan
    if "valid_min" in attrs:
        a[a < attrs["valid_min"]] = np.nan
    if "valid_max" in attrs:
        a[a > attrs["valid_max"]] = np.nan
    return a * attrs.get("scale_factor", 1) + attrs.get("add_offset", 0)


def xyz(lat, lon):
    lat, lon = np.deg2rad(lat), np.deg2rad(lon)
    return np.column_stack(
        (np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat))
    )


class DownloadCache:
    """Disk-bounded cache. Evicts ONLY files listed in its own managed-files.json."""

    def __init__(self, config, session=None):
        self.config = config
        self.root = Path(config.ancillary_cache).expanduser()
        self.root.mkdir(parents=True, exist_ok=True)
        self.session = session or requests.Session()
        self.record = self.root / "managed-files.json"
        self.owned = json.loads(self.record.read_text()) if self.record.exists() else {}

    def get(self, name):
        name = Path(name).name
        if not re.fullmatch(r"GMAO_[A-Za-z0-9_.-]+\.(PROFILE|MET)\.nc", name):
            raise ValueError(f"Unsupported original ancillary filename: {name}")
        target = self.root / name
        if target.exists():
            with xr.open_dataset(target) as d:
                if d.attrs.get("product_name") != name:
                    raise ValueError("Ancillary product mismatch")
            if name in self.owned:
                self.owned[name]["used"] = time.time()
                atomic_json(self.record, self.owned)
            return target
        url = "https://oceandata.sci.gsfc.nasa.gov/ob/getfile/" + quote(name)
        # requests can use ~/.netrc; pass earthaccess.get_requests_https_session()
        # on Earthdata-authenticated systems when needed. No signed URLs persisted.
        with self.session.get(url, stream=True, timeout=(30, 120)) as r:
            r.raise_for_status()
            if "text/html" in r.headers.get("Content-Type", ""):
                raise RuntimeError(
                    "Ancillary download returned HTML; Earthdata login required"
                )
            size = int(r.headers.get("Content-Length", 0))
            limit = self.config.max_download_cache_mb * 2**20
            if size > limit:
                raise OSError("One ancillary exceeds configured download-cache limit")

            def used():
                return sum(
                    (self.root / n).stat().st_size
                    for n in self.owned
                    if (self.root / n).exists()
                )

            for n in sorted(self.owned, key=lambda n: self.owned[n]["used"]):
                if used() + size <= limit:
                    break
                (self.root / n).unlink(missing_ok=True)
                del self.owned[n]
            if (
                shutil.disk_usage(self.root).free
                < size + self.config.minimum_free_disk_mb * 2**20
            ):
                raise OSError("Insufficient disk headroom for ancillary download")
            tmp = target.with_suffix(".nc.part")
            written = 0
            try:
                with tmp.open("wb") as f:
                    for block in r.iter_content(2**20):
                        written += len(block)
                        if used() + written > limit:
                            raise OSError("Ancillary cache size limit reached")
                        f.write(block)
                with xr.open_dataset(tmp) as d:
                    if d.attrs.get("product_name") != name:
                        raise ValueError("Wrong ancillary product")
                tmp.replace(target)
            except Exception:
                tmp.unlink(missing_ok=True)
                raise
        self.owned[name] = {"used": time.time(), "sha256": sha256(target)}
        atomic_json(self.record, self.owned)
        return target


def sample_field(path, var, lat, lon):
    """Bilinear spatial interpolation, periodic longitude; no missing-data fill."""
    with xr.open_dataset(path, engine="h5netcdf") as ds:
        a = ds[var].sortby("lat").sortby("lon")
        if "lev" in a.dims:
            a = a.transpose("lat", "lon", "lev")
        else:
            a = a.transpose("lat", "lon")
        lats = a.lat.values
        lons = a.lon.values
        v = a.values
        if var == "QV" and a.attrs.get("units") not in ("kg kg-1", "kg/kg"):
            raise ValueError("Unexpected QV units")
        if var == "PS" and a.attrs.get("units") != "Pa":
            raise ValueError("Unexpected PS units")
        p = ds.lev.values.copy() if "lev" in a.dims else None
        tm = (
            np.datetime64(
                ds.attrs["time_coverage_start"].replace("Z", ""), "ns"
            ).astype("int64")
            / 1e9
        )
    v = np.concatenate((v[:, -1:], v, v[:, :1]), axis=1)
    lons = np.r_[lons[-1] - 360, lons, lons[0] + 360]
    result = RegularGridInterpolator(
        (lats, lons), v, bounds_error=False, fill_value=np.nan
    )(np.column_stack((lat, (lon + 180) % 360 - 180)))
    return result, p, tm


def sample_times(names, var, lat, lon, seconds, cache):
    arrays = []
    times = []
    pressure = None
    provenance = []
    for name in sorted(set(names)):
        path = cache.get(name)
        value, p, t = sample_field(path, var, lat, lon)
        if pressure is not None and not np.array_equal(pressure, p):
            raise ValueError("Profile pressure grids changed")
        pressure = p
        arrays.append(value)
        times.append(t)
        provenance.append({"file": name, "sha256": sha256(path)})
    order = np.argsort(times)
    times = np.array(times)[order]
    arrays = [arrays[i] for i in order]
    out = np.full_like(arrays[0], np.nan)
    if len(times) == 1:
        good = np.abs(seconds - times[0]) <= 1
        out[good] = arrays[0][good]
    else:
        for i in range(len(times) - 1):
            good = (seconds >= times[i]) & (seconds <= times[i + 1])
            w = (seconds[good] - times[i]) / (times[i + 1] - times[i])
            if out.ndim > 1:
                w = w[:, None]
            out[good] = (1 - w) * arrays[i][good] + w * arrays[i + 1][good]
    return out, pressure, provenance


def above_cloud_water(pressure_hpa, specific_humidity, ctp, surface_pressure_hpa):
    """CHIMAERA 101-level log-pressure humidity interpolation and PW integration.

    Reject profiles whose cloud-top integration involves surface-adjusted levels;
    avoids guessing the operational surface-level repair logic. This loss is
    explicit (NaN PW). Above those levels the source's humidity calculation is
    reproduced: q -> g/kg mixing ratio; 0.003 g/kg above top profile; integral /980.616.
    """
    p = np.asarray(pressure_hpa, float)
    order = np.argsort(p)
    p = p[order]
    q = np.asarray(specific_humidity, float)[:, order]
    l = np.arange(101, 0, -1, dtype=float)
    pp = (
        -1.550789414500298e-4 * l * l - 5.593654380586063e-2 * l + 7.451622227151780
    ) ** 3.5
    result = np.full(len(ctp), np.nan)
    for i, pc in enumerate(ctp):
        ps = surface_pressure_hpa[i]
        if not np.isfinite(pc + ps) or pc <= pp[0] or pc >= ps:
            continue
        # First 101-level endpoint above CTP needs valid unaffected bracketing data.
        end = min(np.searchsorted(pp, pc), 100)
        upper = p[min(np.searchsorted(p, pp[end]), len(p) - 1)]
        if upper >= ps or pp[end] > p[-1]:
            continue
        need = p <= upper
        qi = q[i, need]
        if not np.all(np.isfinite(qi) & (qi >= 0) & (qi < 1)):
            continue
        mixing = 1000 * qi / (1 - qi)
        ww = np.interp(np.log(pp[: end + 1]), np.log(p[need]), mixing, left=0.003)
        grid = np.r_[pp[:end], pc]
        last = np.interp(pc, pp[end - 1 : end + 1], ww[end - 1 : end + 1])
        values = np.r_[ww[:end], last]
        result[i] = np.sum(np.diff(grid) * (values[:-1] + values[1:]) * 0.5) / 980.616
    return result


class OCICompanion:
    def __init__(self, catalog, config=CompanionConfig(), session=None):
        self.catalog = catalog
        self.cfg = config
        self.cache = DownloadCache(config, session)
        if config.max_distance_km <= 0 or config.max_time_difference_s <= 0:
            raise ValueError("Positive match tolerances required")

    def __call__(self, ds):
        # Old extraction code used nadir + offset; the NASA L1C writer stores
        # offset = nadir - mean(native scan time). Repair without changing radiometry.
        offset = ds.view_offset_s.values.astype(float)
        if ds.attrs.get("view_time_offset_convention") == "nadir_minus_observation":
            seconds = (
                ds.oci_observation_time.values.astype("datetime64[ns]").astype("int64")
                / 1e9
            )
        elif "oci_observation_time" in ds:
            seconds = (
                ds.oci_observation_time.values.astype("datetime64[ns]").astype("int64")
                / 1e9
                - 2 * offset
            )
        else:
            seconds = (
                ds.harp_nadir_time.values.astype("datetime64[ns]").astype("int64") / 1e9
                - offset
            )

        def epoch(s):
            return np.datetime64(s.replace("Z", ""), "ns").astype("int64") / 1e9

        margin = self.cfg.max_time_difference_s
        candidates = [
            g
            for g in self.catalog
            if epoch(g["start"]) <= seconds.max() + margin
            and epoch(g["end"]) >= seconds.min() - margin
        ]
        if not candidates:
            raise ValueError(
                "No cloud granules overlap corrected OCI observation times"
            )
        n = ds.sizes["sample"]
        best = np.full(n, np.inf)
        out = ds.copy()
        owner = np.full(n, -1, int)
        keys = [
            "cer_21",
            "cot_21",
            "cld_phase_21",
            "cer_22",
            "cot_22",
            "cld_phase_22",
            "ctp",
            "cld_non_abs_band",
        ]
        values = {k: np.full(n, np.nan) for k in keys}
        tlat = np.full(n, np.nan)
        tlon = tlat.copy()
        native_time = tlat.copy()
        delta = tlat.copy()
        rows = np.full(n, -1, int)
        cols = rows.copy()
        sources = []
        targets = xyz(ds.latitude.values, ds.longitude.values)
        for gi, entry in enumerate(candidates):
            with cloud_file(entry["url"], self.cfg.transport) as f:
                pars = dict(
                    f.groups["processing_control"].groups["input_parameters"].attrs
                )
                if str(f.attrs["processing_version"]) != self.cfg.cloud_version:
                    raise ValueError("CLD processing version mismatch")
                if not re.fullmatch(
                    r"PACE_OCI\.\d{8}T\d{6}\.L1B\.V3\.nc", Path(str(pars["ifile"])).name
                ):
                    raise ValueError("Unexpected source L1B version")
                nav = f.groups["navigation_data"]
                geo = f.groups["geophysical_data"]
                lat = read_valid(nav.variables["latitude"])
                lon = read_valid(nav.variables["longitude"])
                good = (
                    np.isfinite(lat + lon) & (np.abs(lat) <= 90) & (np.abs(lon) <= 180)
                )
                flat = np.flatnonzero(good)
                if not len(flat):
                    sources.append(None)
                    continue
                tree = cKDTree(xyz(lat.ravel()[flat], lon.ravel()[flat]))
                dist, near = tree.query(targets, k=min(8, len(flat)))
                if dist.ndim == 1:
                    dist = dist[:, None]
                    near = near[:, None]
                source = flat[near]
                row, col = np.unravel_index(source, lat.shape)
                times = read_valid(f.groups["scan_line_attributes"].variables["time"])[
                    row
                ]
                km = 2 * 6371.0088 * np.arcsin(np.clip(dist / 2, 0, 1))
                allowed = (abs(times - seconds[:, None]) <= margin) & (
                    km <= self.cfg.max_distance_km
                )
                score = np.where(allowed, km, np.inf)
                j = np.argmin(score, axis=1)
                ix = np.arange(n)
                distance = score[ix, j]
                take = distance < best
                selected = source[ix, j]
                for k in keys:
                    a = read_valid(geo.variables[k]).ravel()[selected]
                    values[k][take] = a[take]
                best[take] = distance[take]
                owner[take] = gi
                rows[take] = row[ix, j][take]
                cols[take] = col[ix, j][take]
                tlat[take] = lat.ravel()[selected][take]
                tlon[take] = lon.ravel()[selected][take]
                native_time[take] = times[ix, j][take]
                delta[take] = native_time[take] - seconds[take]
                sources.append(
                    {
                        "cloud_file": str(f.attrs["product_name"]),
                        "cloud_concept_id": entry["concept_id"],
                        "processing_software": str(
                            f.groups["processing_control"].attrs.get(
                                "software_version", ""
                            )
                        ),
                        "input_parameters": {
                            k: str(v)
                            for k, v in pars.items()
                            if k.startswith(("anc_profile", "met", "cld_"))
                        },
                    }
                )
        for k, a in values.items():
            key = (
                k.replace("cld_phase", "oci_phase")
                if k.startswith("cld_phase")
                else ("oci_" + k if k.startswith(("cer", "cot")) else k)
            )
            out[key] = ("sample", a.astype("float32"))
        matched = owner >= 0
        # Enforce the requested date range using the repaired observation timestamps.
        config = json.loads(ds.attrs.get("match_config", "{}"))
        if "start_date" in config:
            lower = epoch(config["start_date"])
            upper = epoch(config["end_date"]) + 86400
            matched &= (seconds >= lower) & (seconds < upper)
        out["oci_ocean"] = ("sample", (out.cld_non_abs_band.values == 2).astype("int8"))
        out["oci_match_valid"] = ("sample", matched.astype("int8"))
        out["oci_match_distance_km"] = (
            "sample",
            np.where(np.isfinite(best), best, np.nan).astype("float32"),
        )
        out["oci_match_time_difference_s"] = ("sample", delta.astype("float32"))
        out["oci_l2_row"] = ("sample", rows.astype("int32"))
        out["oci_l2_column"] = ("sample", cols.astype("int32"))
        out["oci_l2_source_index"] = ("sample", owner.astype("int16"))
        out["oci_observation_time"] = (
            "sample",
            (seconds * 1e9).astype("datetime64[ns]"),
        )
        out.attrs["view_time_offset_convention"] = "nadir_minus_observation"
        pw = np.full(n, np.nan)
        for gi, provenance in enumerate(sources):
            take = (owner == gi) & matched
            if not take.any():
                continue
            pars = provenance["input_parameters"]
            profile = [
                Path(str(pars[k])).name
                for k in ("anc_profile1", "anc_profile2", "anc_profile3")
                if str(pars.get(k, ""))
            ]
            met = [
                Path(str(pars[k])).name
                for k in ("met1", "met2", "met3")
                if str(pars.get(k, ""))
            ]
            if not profile or not met:
                raise ValueError(
                    "CLD lacks original profile/MET provenance; no substitute selected"
                )
            q, p, qprov = sample_times(
                profile, "QV", tlat[take], tlon[take], native_time[take], self.cache
            )
            ps, _, pprov = sample_times(
                met, "PS", tlat[take], tlon[take], native_time[take], self.cache
            )
            pw[take] = above_cloud_water(p, q, out.ctp.values[take], ps / 100)
            provenance["ancillary_files"] = qprov + pprov
        out["above_cloud_water"] = ("sample", pw.astype("float32"))
        out.above_cloud_water.attrs["units"] = "g cm-2"
        out.ctp.attrs["units"] = "hPa"
        out.attrs["companion_provenance"] = json.dumps(
            {
                "sources": sources,
                "spatial_transfer": "nearest time-compatible native OCI L2 pixel, not footprint mean; no parallax correction",
                "profile_handling": "101-level log-pressure interpolation; near-surface affected profiles rejected",
                "time_handling": "nadir minus view offset; legacy cache timestamps repaired; adjacent granules searched",
            }
        )
        return out
