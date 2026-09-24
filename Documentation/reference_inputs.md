# Optional external ratio references

This document describes the generic `Reference` adapter used by the paired-observation
notebook. The 265 K workflow reconstructs liquid reflectances and corrects atmospheric
transmission separately; see [algorithm details](algorithm_details.md). These are
distinct reference interfaces.

The normalized observable is:

R = (rho2260 / rho1615)_observed / (rho2260 / rho1615)_all-liquid,simulated

No simulated baseline is supplied or inferred from the observed distribution.
The default is explicitly labeled raw-ratio exploration. The supplied ACTIVATE
cutoff of 1.27 is not validated for raw OCI ratios or even automatically for
normalized OCI ratios: OCI's shorter band differs from RSP's 1590 nm channel.

## Scalar reference

In the notebook, set:

```python
from pace_specpol.reference import Reference

reference = Reference(
    mode="scalar",
    scalar=YOUR_SIMULATED_RATIO,
    description="Describe RT model, OCI bandpasses, atmosphere, geometry, cloud assumptions",
)
```

The scalar must be finite and positive. It represents a simplifying assumption
that one baseline is appropriate for all selected observations.

## Lookup table

Supply a NetCDF file containing `liquid_reflectance_ratio`. Supported axes are
any nonempty subset of:

| Axis | Meaning | Units |
|---|---|---|
| sza | OCI solar zenith | degrees |
| vza | OCI sensor zenith | degrees |
| raa | Absolute folded sensor-minus-solar azimuth | degrees, 0–180 |
| cot | GPC ancillary cloud optical thickness | dimensionless |
| re | GPC cloud-bow effective radius | micrometers |

Axes must be strictly increasing with at least two points each. The evaluator
uses multilinear interpolation and excludes samples outside the LUT domain or
with missing inputs. It never extrapolates. Ratios must be finite and positive.
Global attributes must include `numerator_nominal_nm=2260` and
`denominator_nominal_nm=1615`.

For example, this creates an intentionally incomplete template, NOT a simulation:

```python
import numpy as np
import xarray as xr

coords = {
    "sza": np.arange(0, 76, 5),
    "vza": np.arange(0, 66, 5),
    "raa": np.arange(0, 181, 10),
}
shape = tuple(len(values) for values in coords.values())
lut = xr.Dataset(
    {"liquid_reflectance_ratio": (tuple(coords), np.full(shape, np.nan))},
    coords=coords,
    attrs={
        "numerator_nominal_nm": 2260,
        "denominator_nominal_nm": 1615,
        "description": "REPLACE with simulation source and fixed cloud/atmosphere assumptions",
    },
)
# Populate this variable from your radiative-transfer calculations before use.
# All-NaN placeholders are rejected by the evaluator.
lut.to_netcdf("all_liquid_reference_TEMPLATE.nc", engine="h5netcdf")
```

After filling the table:

```python
from pace_specpol.reference import Reference
from pace_specpol.validation.aggregation import build_index
from pace_specpol.validation.vis_dashboard import dashboard

reference = Reference(
    mode="lut",
    lut_path="my_filled_all_liquid_reference.nc",
    description="Describe RT model, OCI response functions, cloud and atmosphere assumptions",
)
index = build_index(cfg.cache_dir, reference, index_cfg)
viewer = dashboard(index, threshold=1.27)
```

The lookup-table SHA-256 is recorded in the index provenance. A new reference
creates a new index from already-cached samples; no new satellite access occurs.
A LUT using `re` will omit missing cloud-bow retrievals, potentially biasing the
population against ice or otherwise difficult cases. No effective radius is
invented for those observations. A geometry-only LUT must document its fixed
microphysical assumptions instead.

The real OCI test file's selected centers were 1618.034 and 2258.429 nm. Use OCI
spectral-response-convolved simulations consistent with TOA reflectance and
atmospheric absorption; wavelength centers alone do not define band-integrated
reflectance. This generic adapter does not perform radiative transfer or gas
correction; its supplied reference must use conventions consistent with the
observed TOA ratio. The separate 265 K workflow does implement the documented
LUT reconstruction and transmission correction.
