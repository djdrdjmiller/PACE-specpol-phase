"""Small synthetic inputs; no mission data bundled."""

import numpy as np
import xarray as xr


def samples():
    vals = {
        "sza": [30] * 3,
        "vza": [20] * 3,
        "raa": [50] * 3,
        "raw_ratio": [1, 1.2, 1.4],
        "rho1615": [0.2] * 3,
        "rho2260": [0.2, 0.24, 0.28],
        "ctp": [700] * 3,
        "above_cloud_water": [1] * 3,
        "oci_cot_22": [10] * 3,
        "oci_cer_22": [10, 35, 10],
        "oci_phase_22": [2, 3, 3],
        "oci_cot_21": [10] * 3,
        "oci_cer_21": [10] * 3,
        "oci_phase_21": [2, 3, 3],
        "oci_ocean": [1, 1, 0],
        "oci_match_valid": [1] * 3,
        "re": [10, np.nan, 10],
        "latitude": [0] * 3,
        "longitude": [0] * 3,
        "liquid_index": [1] * 3,
    }
    return xr.Dataset({k: ("sample", v) for k, v in vals.items()})
