"""Optional sample aggregation for analysis; not an algorithm input grid."""

from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
import json
import numpy as np
import xarray as xr
from scipy import sparse
from pace_specpol.matching import CODE_VERSION, _digest, _write_json
from pace_specpol.reference import Reference


@dataclass(frozen=True)
class IndexConfig:
    resolution: float = 0.1
    li_threshold: float = 0.3
    ratio_min: float = 0.5
    ratio_max: float = 3.5
    ratio_step: float = 0.01
    histogram_li_min: float = 0.0
    histogram_li_max: float = 10.0
    histogram_li_bins: int = 100

    def validate(self):
        if not np.isfinite(self.resolution) or not 0 < self.resolution <= 180:
            raise ValueError("Invalid grid resolution")
        if any(
            not np.isclose(span / self.resolution, round(span / self.resolution))
            for span in (180, 360)
        ):
            raise ValueError("Grid resolution must divide 180 and 360")
        if not np.isfinite(self.li_threshold):
            raise ValueError("LI threshold must be finite")
        if not (
            np.isfinite(self.ratio_min)
            and np.isfinite(self.ratio_max)
            and np.isfinite(self.ratio_step)
            and 0 < self.ratio_step < self.ratio_max - self.ratio_min
        ):
            raise ValueError("Invalid ratio slider range")
        if not np.isclose(
            (self.ratio_max - self.ratio_min) / self.ratio_step,
            round((self.ratio_max - self.ratio_min) / self.ratio_step),
        ):
            raise ValueError("ratio_step must divide the slider range")
        if self.histogram_li_bins < 1 or self.histogram_li_max <= self.histogram_li_min:
            raise ValueError("Invalid histogram settings")

    @property
    def edges(self):
        self.validate()
        n = round((self.ratio_max - self.ratio_min) / self.ratio_step)
        return np.round(np.linspace(self.ratio_min, self.ratio_max, n + 1), 12)

    @property
    def shape(self):
        return round(180 / self.resolution), round(360 / self.resolution)


def cell_ids(lat, lon, resolution):
    ny, nx = round(180 / resolution), round(360 / resolution)
    iy = np.minimum(
        np.floor((np.asarray(lat) + 90) / resolution).astype("int64"), ny - 1
    )
    ix = np.minimum(
        np.floor(((np.asarray(lon) + 180) % 360) / resolution).astype("int64"), nx - 1
    )
    return iy * nx + ix


def build_index(cache_dir, reference=Reference(), config=IndexConfig(), index_dir=None):
    """Stream cached pairs once; sparse per-cell ratio bins retain exact slider cuts.

    Samples outside slider/histogram limits are retained in the classification.
    Histogram plot bounds do not define a science mask.
    """
    config.validate()
    source = Path(cache_dir)
    manifest = json.loads((source / "manifest.json").read_text())
    if not manifest.get("complete"):
        raise ValueError(
            "Pair extraction incomplete; finish or explicitly create a separate pilot cache"
        )
    paths = [source / f for f in manifest["files"]]
    meta = {
        "code_version": CODE_VERSION,
        "index_schema": 2,
        "reference": reference.metadata(),
        "index_config": asdict(config),
        "source_manifest_sha256": _digest(manifest),
        "match_config": manifest["config"],
        "data_scope": manifest.get("scope", "all paired granules in cache manifest"),
    }
    fingerprint = _digest(meta)
    target = Path(index_dir) if index_dir else source / ("index_" + fingerprint[:12])
    if (target / "metadata.json").exists():
        saved = json.loads((target / "metadata.json").read_text())
        if saved["fingerprint"] != fingerprint:
            raise ValueError(
                "Index directory belongs to another reference/configuration"
            )
        return PhaseIndex.load(target)
    evaluate = reference.evaluator()
    edges = config.edges
    yedges = np.linspace(
        config.histogram_li_min, config.histogram_li_max, config.histogram_li_bins + 1
    )
    shape = (len(edges) + 1, int(np.prod(config.shape)))
    low = sparse.csr_matrix(shape, dtype="int64")
    high = sparse.csr_matrix(shape, dtype="int64")
    histogram = np.zeros((len(edges) - 1, len(yedges) - 1), dtype="int64")
    sums = np.zeros((2, shape[1]), dtype="float64")
    totals = {
        "cached_samples": 0,
        "indexed_samples": 0,
        "reference_invalid_samples": 0,
        "histogram_outside_samples": 0,
    }
    for i, path in enumerate(paths):
        with xr.open_dataset(path, engine="h5netcdf") as ds:
            ds.load()
            n = ds.sizes["sample"]
            totals["cached_samples"] += n
            if n == 0:
                continue
            baseline = evaluate(ds)
            good = np.isfinite(baseline) & (baseline > 0)
            totals["reference_invalid_samples"] += int((~good).sum())
            ratio = np.full(n, np.nan, dtype="float64")
            np.divide(ds.raw_ratio.values, baseline, out=ratio, where=good)
            li = ds.liquid_index.values
            good &= np.isfinite(ratio) & np.isfinite(li)
            ratio, li = ratio[good], li[good]
            cells = cell_ids(
                ds.latitude.values[good], ds.longitude.values[good], config.resolution
            )
            np.add.at(sums[0], cells, li)
            np.add.at(sums[1], cells, ratio)
            # side=right ensures exactly-at-threshold samples go to ice/LTMP.
            bins = np.searchsorted(edges, ratio, side="right")
            is_high = li >= config.li_threshold
            for select, category in ((~is_high, "low"), (is_high, "high")):
                part = sparse.coo_matrix(
                    (
                        np.ones(int(select.sum()), dtype="int64"),
                        (bins[select], cells[select]),
                    ),
                    shape=shape,
                ).tocsr()
                if category == "low":
                    low = low + part
                else:
                    high = high + part
            h = np.histogram2d(ratio, li, bins=(edges, yedges))[0].astype("int64")
            histogram += h
            totals["indexed_samples"] += len(ratio)
            totals["histogram_outside_samples"] += int(len(ratio) - h.sum())
        if (i + 1) % 20 == 0 or i + 1 == len(paths):
            print(f"Indexed {i + 1}/{len(paths)} cached pairs", flush=True)
    if not totals["indexed_samples"]:
        raise RuntimeError(
            "No finite normalized samples: check LUT coverage/reference and matching settings"
        )
    meta.update(
        fingerprint=fingerprint,
        totals=totals,
        meaning="sample classifications aggregated by count; not classification of L3 means",
    )
    target.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(target / "low_li.npz", low)
    sparse.save_npz(target / "high_li.npz", high)
    np.savez_compressed(
        target / "histogram.npz", counts=histogram, ratio_edges=edges, li_edges=yedges
    )
    np.savez_compressed(target / "sums.npz", liquid_index=sums[0], ratio=sums[1])
    _write_json(target / "metadata.json", meta)  # completion marker written last
    return PhaseIndex(low, high, histogram, edges, yedges, meta, target, sums)


class PhaseIndex:
    def __init__(
        self, low, high, histogram, edges, li_edges, metadata, directory=None, sums=None
    ):
        self.low, self.high = low, high
        self.sums = sums
        self.histogram, self.edges, self.li_edges = histogram, edges, li_edges
        self.metadata, self.directory = metadata, directory
        self.config = IndexConfig(**metadata["index_config"])

    @classmethod
    def load(cls, directory):
        p = Path(directory)
        sums = None
        if (p / "sums.npz").exists():
            with np.load(p / "sums.npz") as z:
                sums = np.stack((z["liquid_index"], z["ratio"]))
        with np.load(p / "histogram.npz") as h:
            return cls(
                sparse.load_npz(p / "low_li.npz"),
                sparse.load_npz(p / "high_li.npz"),
                h["counts"],
                h["ratio_edges"],
                h["li_edges"],
                json.loads((p / "metadata.json").read_text()),
                p,
                sums,
            )

    def threshold_index(self, threshold):
        i = int(np.argmin(abs(self.edges - threshold)))
        if not np.isclose(self.edges[i], threshold, rtol=0, atol=1e-9):
            raise ValueError(
                "Threshold must be on the configured slider grid; rebuild index for other cuts"
            )
        return i

    def counts(self, threshold, low=None, high=None):
        low = self.low if low is None else low
        high = self.high if high is None else high
        stop = self.threshold_index(threshold) + 1

        def total(m):
            return np.asarray(m.sum(axis=0)).ravel().astype("int64")

        left_low, left_high = total(low[:stop]), total(high[:stop])
        return np.stack(
            (left_low + left_high, total(low) - left_low, total(high) - left_high)
        )

    @staticmethod
    def dominant(counts, min_count=1):
        if min_count < 1:
            raise ValueError("min_count must be >= 1")
        total = counts.sum(axis=0)
        maximum = counts.max(axis=0)
        ties = (counts == maximum).sum(axis=0) > 1
        phase = (np.argmax(counts, axis=0) + 1).astype("int8")
        phase[(total < min_count) | ties] = 0
        return phase

    def mean_fields(self, min_count=1):
        """Arithmetic sample means, including values outside histogram bounds."""
        if self.sums is None:
            raise ValueError(
                "Re-run build_index from the existing pair cache to add mean fields; no NASA downloads needed. Use a new index_dir if explicitly set."
            )
        total = (
            np.asarray(self.low.sum(axis=0)).ravel()
            + np.asarray(self.high.sum(axis=0)).ravel()
        )
        means = np.full(self.sums.shape, np.nan, dtype="float64")
        np.divide(
            self.sums, total[None, :], out=means, where=total[None, :] >= min_count
        )
        return means, total

    def sampling_summary(self):
        _, total = self.mean_fields()
        n = total[total > 0]
        return {
            "occupied_cells": int(n.size),
            "mean_samples_per_occupied_cell": float(n.mean()),
            "median_samples": float(np.median(n)),
            "p90_samples": float(np.percentile(n, 90)),
            "maximum_samples": int(n.max()),
            "single_sample_cell_fraction": float(np.mean(n == 1)),
        }

    def dataset(self, threshold=1.27, min_count=1):
        counts = self.counts(threshold)
        total = counts.sum(axis=0)
        fraction = np.full(counts.shape, np.nan, dtype="float32")
        np.divide(counts, total, out=fraction, where=total[None, :] > 0)
        phase = self.dominant(counts, min_count)
        ny, nx = self.config.shape
        res = self.config.resolution
        ds = xr.Dataset(
            {
                "dominant_phase": (("lat", "lon"), phase.reshape(ny, nx)),
                "sample_count": (("lat", "lon"), total.reshape(ny, nx)),
                "phase_count": (("phase", "lat", "lon"), counts.reshape(3, ny, nx)),
                "phase_fraction": (
                    ("phase", "lat", "lon"),
                    fraction.reshape(3, ny, nx),
                ),
            },
            coords={
                "lat": -90 + (np.arange(ny) + 0.5) * res,
                "lon": -180 + (np.arange(nx) + 0.5) * res,
                "phase": ["liquid", "ice", "LTMP"],
            },
        )
        if self.sums is not None:
            means, _ = self.mean_fields(min_count)
            for name, values, label in zip(
                ("mean_liquid_index", "mean_ratio"),
                means,
                (
                    "Arithmetic mean HARP2 liquid index",
                    "Arithmetic mean per-sample OCI ratio in configured reference mode",
                ),
            ):
                ds[name] = (("lat", "lon"), values.reshape(ny, nx).astype("float32"))
                ds[name].attrs.update(units="1", long_name=label)
        ds.dominant_phase.attrs.update(
            flag_values=np.array([0, 1, 2, 3], dtype="int8"),
            flag_meanings="missing_insufficient_or_tied liquid ice LTMP",
        )
        ds.lat.attrs.update(units="degrees_north", standard_name="latitude")
        ds.lon.attrs.update(units="degrees_east", standard_name="longitude")
        ds.phase_fraction.attrs.update(
            units="1",
            long_name="Fraction of eligible matched samples, not cloud volume or water fraction",
        )
        ds.attrs.update(
            title="Exploratory OCI/HARP2 cloud phase candidates",
            ratio_threshold=float(threshold),
            liquid_index_threshold=self.config.li_threshold,
            minimum_sample_count=int(min_count),
            reference=json.dumps(self.metadata["reference"]),
            provenance=json.dumps(self.metadata),
            algorithm_version=CODE_VERSION,
            tie_policy="blank; never resolve ties by class ordering",
            date_created=datetime.now(timezone.utc).isoformat(),
        )
        return ds

    def export(self, path, threshold=1.27, min_count=1):
        ds = self.dataset(threshold, min_count)
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        ds.to_netcdf(
            path,
            engine="h5netcdf",
            encoding={k: {"zlib": True, "complevel": 3} for k in ds.data_vars},
        )
        return path
