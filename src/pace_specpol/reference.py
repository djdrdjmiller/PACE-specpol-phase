"""Raw, scalar, and externally supplied ratio-reference adapters."""

from dataclasses import dataclass, asdict
from pathlib import Path
import hashlib
import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator


@dataclass(frozen=True)
class Reference:
    mode: str = "raw"  # raw | scalar | lut
    scalar: float | None = None
    lut_path: str | None = None
    description: str = "No all-liquid simulation supplied; raw-ratio exploration only"

    def metadata(self):
        result = asdict(self)
        if self.mode not in ("raw", "scalar", "lut"):
            raise ValueError("Reference mode must be raw, scalar or lut")
        if self.mode == "scalar" and (
            self.scalar is None or not np.isfinite(self.scalar) or self.scalar <= 0
        ):
            raise ValueError("Supply a finite positive all-liquid simulated ratio")
        if self.mode != "raw" and (
            not self.description.strip() or self.description.startswith("No all-liquid")
        ):
            raise ValueError(
                "Describe the actual simulation, bandpasses, and assumptions"
            )
        if self.mode == "lut":
            if self.lut_path is None:
                raise ValueError("Supply lut_path")
            result["lut_sha256"] = hashlib.sha256(
                Path(self.lut_path).read_bytes()
            ).hexdigest()
        return result

    def evaluator(self):
        self.metadata()
        if self.mode == "raw":
            return lambda ds: np.ones(ds.sizes["sample"], dtype=float)
        if self.mode == "scalar":
            return lambda ds: np.full(ds.sizes["sample"], self.scalar)
        with xr.open_dataset(self.lut_path, engine="h5netcdf") as f:
            a = f["liquid_reflectance_ratio"].load()
            if (
                float(f.attrs.get("numerator_nominal_nm", 0)) != 2260
                or float(f.attrs.get("denominator_nominal_nm", 0)) != 1615
            ):
                raise ValueError(
                    "LUT must explicitly describe OCI 2260/1615 reflectance ratio"
                )
        if not a.dims or not set(a.dims) <= {"sza", "vza", "raa", "cot", "re"}:
            raise ValueError(
                "Supported LUT dimensions: sza, vza, raa, cot, re (any nonempty subset)"
            )
        for dim in a.dims:
            if a.sizes[dim] < 2 or not np.all(np.diff(a[dim].values) > 0):
                raise ValueError(
                    "LUT axes must have >=2 strictly increasing coordinates"
                )
        if not np.all(np.isfinite(a.values) & (a.values > 0)):
            raise ValueError("LUT ratios must be finite and positive")
        interpolation = RegularGridInterpolator(
            tuple(a[d].values for d in a.dims),
            a.values,
            bounds_error=False,
            fill_value=np.nan,
        )

        def evaluate(ds):
            points = np.column_stack([ds[d].values for d in a.dims])
            return interpolation(points)

        return evaluate
