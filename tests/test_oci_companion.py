"""Offline numerical, selection, and resume tests for the liquid-reference extension."""

import numpy as np
import pytest
import xarray as xr
from pace_specpol.oci_companion import above_cloud_water, sample_field


def test_water_integral_and_invalid_profile():
    p = np.array([0.1, 1, 10, 100, 500, 700, 800, 900, 1000.0])
    q = np.full((2, len(p)), 0.01)
    result = above_cloud_water(p, q, [700, 990], [1010, 995])
    # Constant mixing ratio apart from tiny upper-atmosphere segment.
    assert result[0] == pytest.approx(
        (1000 * 0.01 / 0.99) * (700 - 0.1) / 980.616, rel=2e-4
    )
    assert np.isnan(result[1])  # near-surface convention not guessed
    q[0, 4] = np.nan
    assert np.isnan(above_cloud_water(p, q, [700, 700], [1010, 1010])[0])


def test_periodic_geos_sampling(tmp_path):
    p = tmp_path / "geo.nc"
    a = xr.Dataset(
        {
            "PS": (
                ("lat", "lon"),
                [[100, 200, 300, 400], [200, 300, 400, 500]],
                {"units": "Pa"},
            )
        },
        coords={"lat": [-1, 1], "lon": [-180, -90, 0, 90]},
        attrs={"time_coverage_start": "2024-09-03T00:00:00Z"},
    )
    a.to_netcdf(p, engine="h5netcdf")
    v, _, _ = sample_field(p, "PS", np.array([0, 0]), np.array([179, -181]))
    np.testing.assert_allclose(v[0], v[1])
