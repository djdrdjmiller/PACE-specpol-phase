"""Resumable companion/correction stages and variant selector for existing dashboard."""

from dataclasses import asdict, replace
from pathlib import Path
import json
import shutil
import numpy as np
import xarray as xr
from pace_specpol.liquid_reference import (
    LiquidLUT,
    TransmissionLUT,
    correct_samples,
    sha256,
    VERSION,
)
from pace_specpol.oci_companion import atomic_json


def digest(obj):
    import hashlib

    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, default=str).encode()
    ).hexdigest()


def write_dataset(ds, path):
    path = Path(path)
    tmp = path.with_suffix(".nc.tmp")
    if shutil.disk_usage(path.parent).free < ds.nbytes * 2 + 128 * 2**20:
        raise OSError("Insufficient disk headroom for output")
    try:
        ds.to_netcdf(
            tmp,
            engine="h5netcdf",
            encoding={k: {"zlib": True, "complevel": 3} for k in ds.data_vars},
        )
        tmp.replace(path)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise


def _prepare(source, target, settings, max_pairs):
    source = Path(source).expanduser()
    target = Path(target).expanduser()
    if source.resolve() == target.resolve():
        raise ValueError("Use a separate output directory")
    manifest = json.loads((source / "manifest.json").read_text())
    if not manifest.get("complete"):
        raise ValueError("Source cache is incomplete")
    if max_pairs is not None and max_pairs < 1:
        raise ValueError("max_pairs must be positive or None")
    files = manifest["files"] if max_pairs is None else manifest["files"][:max_pairs]
    if not files:
        raise ValueError("No source files")
    for name in files:
        if Path(name).name != name:
            raise ValueError("Manifest files must be basenames")
    settings = dict(
        settings,
        source_manifest_sha256=sha256(source / "manifest.json"),
        files=files,
        version=VERSION,
    )
    fingerprint = digest(settings)
    target.mkdir(parents=True, exist_ok=True)
    path = target / "manifest.json"
    if path.exists() and json.loads(path.read_text()).get("fingerprint") != fingerprint:
        raise ValueError("Output settings changed; choose a new output directory")
    out = dict(
        complete=False,
        config=manifest["config"],
        files=[],
        records=[],
        fingerprint=fingerprint,
        settings=settings,
        scope=manifest.get("scope", "")
        + (f" | PILOT {len(files)} pairs" if max_pairs else ""),
    )
    atomic_json(path, out)
    return source, target, files, out


def cache_companions(source_dir, output_dir, provider, max_pairs=None):
    source, target, files, manifest = _prepare(
        source_dir,
        output_dir,
        dict(
            stage="companions", provider=asdict(provider.cfg), catalog=provider.catalog
        ),
        max_pairs,
    )
    for i, name in enumerate(files):
        path = target / name
        source_hash = sha256(source / name)
        if path.exists():
            with xr.open_dataset(path, engine="h5netcdf") as ds:
                if (
                    ds.attrs.get("stage_fingerprint") != manifest["fingerprint"]
                    or ds.attrs.get("source_sha256") != source_hash
                ):
                    raise ValueError(f"Stale companion: {path}")
                record = json.loads(ds.attrs["stage_report"])
        else:
            with xr.open_dataset(source / name, engine="h5netcdf") as f:
                ds = f.load()
            if ds.sizes["sample"] == 0:
                # Process no remote data; explicit empty schema for correction stage.
                for k in [
                    "oci_cer_21",
                    "oci_cot_21",
                    "oci_phase_21",
                    "oci_cer_22",
                    "oci_cot_22",
                    "oci_phase_22",
                    "ctp",
                    "above_cloud_water",
                    "oci_ocean",
                    "oci_match_valid",
                ]:
                    ds[k] = ("sample", np.array([], dtype="float32"))
            else:
                ds = provider(ds)
            record = {
                "file": name,
                "samples": ds.sizes["sample"],
                "matched": int(ds.oci_match_valid.sum()),
                "water_valid": int(np.isfinite(ds.above_cloud_water).sum()),
            }
            ds.attrs.update(
                stage_fingerprint=manifest["fingerprint"],
                source_sha256=source_hash,
                stage_report=json.dumps(record),
            )
            write_dataset(ds, path)
        manifest["files"].append(name)
        manifest["records"].append(record)
        atomic_json(target / "manifest.json", manifest)
        print(
            f"Companion {i + 1}/{len(files)} | matched {record['matched']:,} | water {record['water_valid']:,}",
            flush=True,
        )
    manifest["complete"] = True
    atomic_json(target / "manifest.json", manifest)
    return target


def default_variants(config):
    return {
        f"{source}_{band}": replace(config, cer_source=source, oci_reference_band=band)
        for source in ("oci", "harp2")
        for band in (2260, 2130)
    }


def cache_references(companion_dir, output_dir, variants, max_pairs=None):
    if not variants:
        raise ValueError("At least one reference mode required")
    for name in variants:
        if not name.replace("_", "").isalnum():
            raise ValueError("Variant names must be alphanumeric/underscore")
    metadata = {k: v.metadata() for k, v in variants.items()}
    # Enforce same physical LUT set so a single in-memory band can be reused.
    first = next(iter(variants.values()))
    for c in variants.values():
        if (c.ms_path, c.phase_path, c.transmittance_path) != (
            first.ms_path,
            first.phase_path,
            first.transmittance_path,
        ):
            raise ValueError("Use one shared physical LUT set per comparison cache")
    source, target, files, manifest = _prepare(
        companion_dir,
        output_dir,
        dict(stage="references", variants=metadata),
        max_pairs,
    )
    lut = LiquidLUT(first)
    trans = TransmissionLUT(first.transmittance_path)
    for i, name in enumerate(files):
        path = target / name
        source_hash = sha256(source / name)
        if path.exists():
            with xr.open_dataset(path, engine="h5netcdf") as ds:
                if (
                    ds.attrs.get("stage_fingerprint") != manifest["fingerprint"]
                    or ds.attrs.get("source_sha256") != source_hash
                ):
                    raise ValueError(f"Stale reference: {path}")
                report = json.loads(ds.attrs["stage_report"])
        else:
            with xr.open_dataset(source / name, engine="h5netcdf") as f:
                ds = f.load()
            out = ds.copy()
            common = np.ones(ds.sizes["sample"], bool)
            report = {"file": name, "samples": ds.sizes["sample"], "variants": {}}
            for variant, cfg in variants.items():
                corrected = correct_samples(ds, cfg, lut, trans)
                for k in corrected.data_vars:
                    if k not in ds:
                        out[f"{variant}_{k}"] = corrected[k]
                valid = corrected.reference_status.values == 0
                common &= valid
                status = corrected.reference_status.values
                report["variants"][variant] = {
                    "valid": int(valid.sum()),
                    "excluded": int((~valid).sum()),
                    "ice_cer_as_liquid_valid": int(
                        (valid & (corrected.ice_cer_as_liquid.values == 1)).sum()
                    ),
                    "reason_counts": {
                        label: int(((status & bit) != 0).sum())
                        for bit, label in [
                            (1, "cer"),
                            (2, "cot"),
                            (4, "transmission"),
                            (8, "reflectance"),
                            (16, "surface"),
                            (32, "phase"),
                            (64, "match"),
                        ]
                    },
                }
            out["common_valid"] = ("sample", common.astype("int8"))
            report["common_valid"] = int(common.sum())
            out.attrs.update(
                stage_fingerprint=manifest["fingerprint"],
                source_sha256=source_hash,
                stage_report=json.dumps(report),
                reference_variants=json.dumps(metadata),
            )
            write_dataset(out, path)
        manifest["files"].append(name)
        manifest["records"].append(report)
        atomic_json(target / "manifest.json", manifest)
        print(
            f"Reference {i + 1}/{len(files)} | "
            + ", ".join(f"{k}: {v['valid']:,}" for k, v in report["variants"].items()),
            flush=True,
        )
    manifest["complete"] = True
    atomic_json(target / "manifest.json", manifest)
    return target


class VariantReference:
    def __init__(self, variant="oci_2260", population="own", mode="normalized"):
        self.variant = variant
        self.population = population
        self.mode = mode
        if population not in ("own", "common"):
            raise ValueError("population: own or common")
        if mode not in ("normalized", "corrected", "raw"):
            raise ValueError("mode: normalized, corrected or raw")

    def metadata(self):
        return dict(
            mode=self.mode,
            variant=self.variant,
            population=self.population,
            version=VERSION,
            description="265 K liquid reference, original OCI ancillary files",
        )

    def evaluator(self):
        def evaluate(ds):
            k = self.variant
            if self.mode == "normalized":
                v = ds[f"{k}_effective_reference"].values.copy()
            elif self.mode == "corrected":
                v = (ds[f"{k}_transmission2260"] / ds[f"{k}_transmission1615"]).values
            else:
                v = np.ones(ds.sizes["sample"])
            # Both own and common views use the selected reference population even
            # for raw/corrected diagnostics, so plots compare identical samples.
            valid = ds[f"{k}_reference_status"].values == 0
            if self.population == "common":
                valid &= ds.common_valid.values == 1
            return np.where(valid, v, np.nan)

        return evaluate
