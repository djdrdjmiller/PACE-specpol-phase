"""ARCSIX exploratory matchups. No automatic lidar-to-LTMP truth conversion."""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import h5py
import xarray as xr
from scipy.spatial import cKDTree

COLLECTIONS = {
    "lidar": "ARCSIX_AircraftRemoteSensing_LaRC-G3_HALO_Data",
    "probes": "ARCSIX_Cloud_AircraftInSitu_P3B_Data",
    "navigation": "ARCSIX_MetNav_AircraftInSitu_P3B_Data",
    "marli": "ARCSIX_AircraftRemoteSensing_P3B_MARLi_Data",
}


def discover(day="2024-06-10"):
    import earthaccess

    return {
        k: earthaccess.search_data(
            short_name=v, temporal=(day + "T00:00:00Z", day + "T23:59:59Z"), count=-1
        )
        for k, v in COLLECTIONS.items()
    }


def name(g):
    return g["umm"]["GranuleUR"]


def download_exact(catalog, filename, folder):
    import earthaccess

    hits = [
        g for g in catalog if name(g) == filename or name(g).endswith("/" + filename)
    ]
    if len(hits) != 1:
        raise ValueError(
            f"Expected one catalog entry for {filename}; found {len(hits)}"
        )
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / filename
    if target.exists():
        return target
    paths = earthaccess.download(hits, local_path=str(folder))
    if len(paths) != 1 or not Path(paths[0]).is_file():
        raise RuntimeError("Download did not complete")
    return Path(paths[0])


def read_hsrl(
    path, day="2024-06-10", top_depth_m=90, min_bins=3, purity=0.75, chunk=512
):
    """Summarize the top-depth window, not an entire-column phase label.

    Code 1 is ambiguous ice/water OR ice/HOI; code 3 is HOI. Neither is LTMP truth.
    GPS navigation is used as a nadir approximation; no pointing/parallax correction.
    """
    with h5py.File(path) as f:
        print("\n".join(f["000_Readme"].asstr()[...].ravel()))
        phase = f["DataProducts/CloudMask/Type"]
        description = str(phase.attrs.get("Long_Name", ""))
        if "0 = Water Dominant" not in description or "3 = HOI" not in description:
            raise ValueError(
                "Unrecognized phase-code metadata; inspect before using reader"
            )
        z = np.asarray(f["DataProducts/CloudMask/Altitude"]).ravel()
        retrieval_top = np.asarray(f["DataProducts/Cloud/Cloud Top Height"]).ravel()
        aircraft_alt = np.asarray(f["Nav/gps_alt"]).ravel()
        surface_alt = np.asarray(f["DataProducts/Surface/altSurface"]).ravel()
        top = np.full(len(retrieval_top), np.nan)
        t = np.asarray(f["Nav/gps_time"]).ravel()
        result = pd.DataFrame(
            {
                "time": pd.Timestamp(day) + pd.to_timedelta(t, unit="h"),
                "latitude": np.asarray(f["Nav/gps_lat"]).ravel(),
                "longitude": np.asarray(f["Nav/gps_lon"]).ravel(),
                "retrieval_cloud_top_m": retrieval_top,
                "profile_id": np.arange(len(top)),
            }
        )
        counts = np.zeros((len(top), 4), int)
        below_ice = np.zeros(len(top), int)
        for lo in range(0, len(top), chunk):
            hi = min(lo + chunk, len(top))
            a = phase[lo:hi, :]
            # Property-retrieval CTH may target a lower water layer underneath ice.
            # Use the highest classified cloud bin visible below the aircraft instead.
            floor = np.where(np.isfinite(surface_alt[lo:hi]), surface_alt[lo:hi], 0)
            observed = (
                np.isin(a, [0, 1, 2, 3])
                & (z[None, :] <= aircraft_alt[lo:hi, None])
                & (z[None, :] >= floor[:, None])
            )
            candidate = np.max(np.where(observed, z[None, :], -np.inf), axis=1)
            top[lo:hi] = np.where(np.isfinite(candidate), candidate, np.nan)
            window = (
                observed
                & (z[None, :] <= top[lo:hi, None])
                & (z[None, :] >= top[lo:hi, None] - top_depth_m)
            )
            for code in range(4):
                counts[lo:hi, code] = np.sum(window & (a == code), axis=1)
            below = (
                observed
                & (z[None, :] < top[lo:hi, None] - top_depth_m)
                & (z[None, :] >= top[lo:hi, None] - 500)
            )
            below_ice[lo:hi] = np.sum(below & (a == 2), axis=1)
    result["cloud_top_m"] = top
    n = counts.sum(1)
    for k in range(4):
        result[f"top_fraction_{k}"] = np.divide(
            counts[:, k], n, out=np.full(len(n), np.nan), where=n > 0
        )
    result["top_cloud_bins"] = n
    result["lidar_top"] = "unknown"
    for k, label in [(0, "water_dominant"), (2, "ice_dominant")]:
        good = (
            (n >= min_bins)
            & (result[f"top_fraction_{k}"] >= purity)
            & (counts[:, 1] == 0)
            & (counts[:, 3] == 0)
        )
        result.loc[good, "lidar_top"] = label
    result.loc[(counts[:, 1] + counts[:, 3]) > 0, "lidar_top"] = "ambiguous_or_HOI"
    result["ice_bins_90_to_500m_below_top"] = (
        below_ice  # evidence only; may be another layer
    )
    result.attrs["phase_description"] = description
    return result


def read_satellite(folder, variant="oci_2260", population="own"):
    folder = Path(folder)
    m = json.loads((folder / "manifest.json").read_text())
    if not m["complete"]:
        raise ValueError("Reference cache incomplete")
    frames = []
    for filename in m["files"]:
        with xr.open_dataset(folder / filename) as d:
            if d.attrs.get("view_time_offset_convention") != "nadir_minus_observation":
                raise ValueError(
                    "Use a newly extracted Arctic cache with corrected offset convention"
                )
            keys = [
                "latitude",
                "longitude",
                "liquid_index",
                "raw_ratio",
                "sza",
                "vza",
                "raa",
                "rho1615",
                "rho2260",
            ]
            a = {k: d[k].values for k in keys}
            a.update(
                time=d.oci_observation_time.values,
                ratio=d[f"{variant}_normalized_ratio"].values,
                harp_time=d.harp_nadir_time.values,
                valid=(d[f"{variant}_reference_status"].values == 0),
            )
            if population == "common":
                a["valid"] &= d.common_valid.values == 1
            elif population != "own":
                raise ValueError("population must be own or common")
            a["sat_id"] = [
                f"{Path(filename).stem}:{r}:{c}:{v}"
                for r, c, v in zip(
                    d.harp_row.values, d.harp_column.values, d.oci_view.values
                )
            ]
            frames.append(pd.DataFrame(a))
    return pd.concat(frames, ignore_index=True)


def xyz(lat, lon):
    p, l = np.deg2rad(lat), np.deg2rad(lon)
    return np.column_stack([np.cos(p) * np.cos(l), np.cos(p) * np.sin(l), np.sin(p)])


def match_profiles(satellite, airborne, max_minutes=10, max_km=3):
    """Each profile gets one closest time-eligible sample; geometry only, before phase/ratio QC.
    Radius uses centers, NOT a true footprint polygon. No advection/parallax correction.
    """
    if max_minutes <= 0 or max_km <= 0:
        raise ValueError("Positive tolerances required")
    s = satellite.reset_index(drop=True)
    a = airborne.reset_index(drop=True)
    good = np.isfinite(s.latitude) & np.isfinite(s.longitude) & s.time.notna()
    ids = np.flatnonzero(good)
    tree = cKDTree(xyz(s.latitude[good], s.longitude[good]))
    st = pd.to_datetime(s.time).to_numpy(dtype="datetime64[ns]")
    at = pd.to_datetime(a.time).to_numpy(dtype="datetime64[ns]")
    records = []
    radius = 2 * np.sin(max_km / (2 * 6371.0088))
    for i in range(len(a)):
        if pd.isna(at[i]) or not np.isfinite(a.latitude.iloc[i] + a.longitude.iloc[i]):
            continue
        point = xyz([a.latitude.iloc[i]], [a.longitude.iloc[i]])[0]
        candidates = ids[tree.query_ball_point(point, radius)]
        dt = (at[i] - st[candidates]) / np.timedelta64(1, "s")
        eligible = abs(dt) <= max_minutes * 60
        candidates = candidates[eligible]
        dt = dt[eligible]
        if not len(candidates):
            continue
        chord = np.linalg.norm(
            xyz(s.latitude.iloc[candidates], s.longitude.iloc[candidates]) - point,
            axis=1,
        )
        best = np.lexsort((abs(dt), chord))[0]
        j = candidates[best]
        records.append(
            (
                i,
                j,
                float(2 * 6371.0088 * np.arcsin(np.clip(chord[best] / 2, 0, 1))),
                float(dt[best]),
            )
        )
    out = pd.DataFrame(
        records,
        columns=["air_index", "sat_index", "distance_km", "air_minus_sat_seconds"],
    )
    if len(out):
        out = out.join(a.add_prefix("air_"), on="air_index").join(
            s.add_prefix("sat_"), on="sat_index"
        )
    return out


def aggregate_matches(matches, min_profiles=3, purity=0.8):
    """One record per satellite sample, preserving all lidar label fractions."""
    rows = []
    if matches.empty:
        return pd.DataFrame()
    for sid, g in matches.groupby("sat_sat_id"):
        labels = g.air_lidar_top
        row = {
            k[4:]: g[k].iloc[0] for k in g if k.startswith("sat_") and k != "sat_index"
        }
        row.update(
            n_profiles=len(g),
            dt_median_min=g.air_minus_sat_seconds.median() / 60,
            distance_max_km=g.distance_km.max(),
            lidar_top="unknown",
        )
        for label in ["water_dominant", "ice_dominant", "ambiguous_or_HOI", "unknown"]:
            row["fraction_" + label] = float((labels == label).mean())
        for label in ["water_dominant", "ice_dominant"]:
            if len(g) >= min_profiles and row["fraction_" + label] >= purity:
                row["lidar_top"] = label
        rows.append(row)
    return pd.DataFrame(rows)


def binary_scores(y, p):
    y = np.asarray(y, bool)
    p = np.asarray(p, bool)
    tp = int(np.sum(y & p))
    fn = int(np.sum(y & ~p))
    fp = int(np.sum(~y & p))
    tn = int(np.sum(~y & ~p))
    div = lambda a, b: a / b if b else np.nan
    return dict(
        tp=tp,
        fn=fn,
        fp=fp,
        tn=tn,
        recall=div(tp, tp + fn),
        precision=div(tp, tp + fp),
        false_positive_rate=div(fp, fp + tn),
        balanced_accuracy=(div(tp, tp + fn) + div(tn, tn + fp)) / 2,
    )


def liquid_top_sweep(footprints, thresholds):
    if footprints.empty:
        return pd.DataFrame()
    q = footprints[
        footprints.lidar_top.isin(["water_dominant", "ice_dominant"])
        & np.isfinite(footprints.liquid_index)
    ]
    return pd.DataFrame(
        [
            dict(
                threshold=t,
                n=len(q),
                **binary_scores(q.lidar_top == "water_dominant", q.liquid_index >= t),
            )
            for t in thresholds
        ]
    )


def reviewed_ratio_sweep(footprints, labels, thresholds, li_threshold=0.3):
    """LTMP vs liquid conditional on independent liquid-top evidence and high LI.
    labels must be independently reviewed, not generated from the ratio.
    """
    required = {"sat_id", "reference_phase", "evidence", "split", "segment_id"}
    if not required.issubset(labels):
        raise ValueError(f"Required columns: {required}")
    if labels.sat_id.duplicated().any():
        raise ValueError("Duplicate reference sat_id")
    known = labels.reference_phase.isin(["liquid", "LTMP"])
    if (labels.loc[known, "evidence"].fillna("").str.strip() == "").any():
        raise ValueError("Document independent evidence for every label")
    q = footprints.merge(labels, on="sat_id", validate="one_to_one")
    q = q[
        (q.lidar_top == "water_dominant")
        & q.valid
        & (q.liquid_index >= li_threshold)
        & np.isfinite(q.ratio)
        & q.reference_phase.isin(["liquid", "LTMP"])
    ]
    if not len(q):
        return pd.DataFrame(), q
    if not q["split"].isin(["train", "test"]).all():
        raise ValueError("Use train/test splits")
    if q.segment_id.isna().any():
        raise ValueError("Assign a flight/cloud segment")
    if (q.groupby("segment_id")["split"].nunique() > 1).any():
        raise ValueError("A segment cannot be in both train and test")
    rows = []
    train = q[q["split"] == "train"]
    if train.reference_phase.nunique() != 2:
        raise ValueError("Training needs independently supported liquid and LTMP cases")
    for t in thresholds:
        rows.append(
            dict(
                threshold=t,
                **binary_scores(train.reference_phase == "LTMP", train.ratio >= t),
            )
        )
    return pd.DataFrame(rows), q


def read_icartt(path):
    """ICARTT FFI 1001 only; honors variable scale factors and missing markers."""
    lines = Path(path).read_text(errors="replace").splitlines()
    first = lines[0].replace(",", " ").split()
    n, ffi = map(int, first[:2])
    if ffi != 1001:
        raise ValueError("Only FFI 1001 supported; inspect this file separately")
    nv = int(lines[9].strip())
    sc = np.fromstring(lines[10], sep=",")
    missing = np.fromstring(lines[11], sep=",")
    if len(sc) != nv or len(missing) != nv:
        raise ValueError("Unexpected ICARTT scale/missing layout")
    names = [x.strip() for x in lines[n - 1].split(",")]
    a = np.loadtxt(path, delimiter=",", skiprows=n, ndmin=2)
    if a.shape[1] != nv + 1 or len(names) != nv + 1:
        raise ValueError("Unexpected ICARTT column count")
    limit_flags = []
    for line in lines[:n]:
        if line.startswith(("ULOD_FLAG:", "LLOD_FLAG:")):
            try:
                limit_flags.append(float(line.split(":", 1)[1].strip()))
            except ValueError:
                pass
    for j in range(nv):
        bad = (a[:, j + 1] == missing[j]) | np.isin(a[:, j + 1], limit_flags)
        a[:, j + 1] *= sc[j]
        a[bad, j + 1] = np.nan
    d = pd.DataFrame(a, columns=names)
    d.attrs.update(
        header="\n".join(lines[:n]), variable_descriptions=lines[12 : 12 + nv]
    )
    return d


def block_bootstrap_scores(reviewed, threshold, split="test", repeats=500, seed=42):
    """Descriptive uncertainty from whole review segments; no random-profile resampling."""
    q = reviewed[reviewed["split"] == split]
    groups = [g for _, g in q.groupby("segment_id")]
    if len(groups) < 3:
        return (
            pd.DataFrame()
        )  # insufficient independent blocks even for a pilot interval
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(repeats):
        b = pd.concat(
            [groups[i] for i in rng.integers(0, len(groups), len(groups))],
            ignore_index=True,
        )
        if b.reference_phase.nunique() != 2:
            continue
        rows.append(binary_scores(b.reference_phase == "LTMP", b.ratio >= threshold))
    return pd.DataFrame(rows)
