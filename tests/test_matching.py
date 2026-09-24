from dataclasses import replace
import numpy as np
import xarray as xr
import pytest
from pace_specpol import matching as m


def make_pair(tmp_path):
    dims = m.DIMS
    vdim = dims + (m.VIEW,)
    bdim = vdim + (m.BAND,)
    shape = (2, 2)
    geo = {"latitude": (dims, np.zeros(shape)), "longitude": (dims, np.zeros(shape))}
    attrs = {
        "time_coverage_start": "2024-09-03T00:12:41Z",
        "time_coverage_end": "2024-09-03T00:17:41Z",
    }
    h = xr.DataTree.from_dict(
        {
            "/": xr.Dataset(
                attrs=dict(
                    attrs,
                    product_name="PACE_HARP2.20240903T001241.L2.CLOUD_GPC.V4_0.nc",
                )
            ),
            "/geolocation_data": xr.Dataset(geo),
            "/bin_attributes": xr.Dataset(
                {
                    "nadir_view_time": (
                        "bins_along_track",
                        [761.0, 762.0],
                        {"units": "seconds from UTC midnight"},
                    )
                }
            ),
            "/geophysical_data": xr.Dataset(
                {
                    m.LI: (dims, np.full(shape, 1.0)),
                    "cloud_quality": (dims, np.zeros(shape)),
                    "cloud_optical_thickness": (dims, np.full(shape, 10.0)),
                    "cloud_bow_droplet_effective_radius": (dims, np.full(shape, 10.0)),
                }
            ),
        }
    )
    g = dict(geo)
    for name, value in [
        ("solar_zenith_angle", 30),
        ("sensor_zenith_angle", 20),
        ("solar_azimuth_angle", 40),
        ("sensor_azimuth_angle", 80),
    ]:
        g[name] = (vdim, np.full((2, 2, 2), float(value)))
    radiance = np.ones((2, 2, 2, 2)) * 2
    # At first cell each band exists only in a DIFFERENT view: must be rejected.
    radiance[0, 0, 0, 1] = np.nan
    radiance[0, 0, 1, 0] = np.nan
    offset = np.full((2, 2, 2), 20.0)
    offset[..., 0] = -10.0
    o = xr.DataTree.from_dict(
        {
            "/": xr.Dataset(
                attrs=dict(
                    attrs,
                    product_name="PACE_OCI.20240903T001241.L1C.V3.5km.nc",
                    sun_earth_distance=1.0,
                )
            ),
            "/geolocation_data": xr.Dataset(g),
            "/sensor_views_bands": xr.Dataset(
                {
                    "intensity_wavelength": ((m.VIEW, m.BAND), [[1618.0, 2258.0]] * 2),
                    "intensity_f0": ((m.VIEW, m.BAND), [[235.0, 74.0]] * 2),
                }
            ),
            "/bin_attributes": xr.Dataset(
                {
                    "nadir_view_time": (
                        "bins_along_track",
                        np.array(
                            ["2024-09-03T00:12:41", "2024-09-03T00:12:42"],
                            dtype="datetime64[ns]",
                        ),
                    ),
                    "view_time_offsets": (vdim, offset),
                }
            ),
            "/observation_data": xr.Dataset(
                {
                    "i": (bdim, radiance, {"units": "W m^-2 sr^-1 um^-1"}),
                    "qc": (bdim, np.full((2, 2, 2, 2), np.nan)),
                    "number_of_observations": (vdim, np.ones((2, 2, 2))),
                }
            ),
        }
    )
    hp, op = tmp_path / "h.nc", tmp_path / "o.nc"
    h.to_netcdf(hp, engine="h5netcdf")
    o.to_netcdf(op, engine="h5netcdf")
    return hp, op, o


def test_pair_alignment_and_same_view(tmp_path):
    hp, op, o = make_pair(tmp_path)
    cfg = m.MatchConfig()
    samples, report = m.match_pair(hp, op, cfg)
    assert len(samples.sample) == 3  # the cross-view-only pair was removed
    assert np.all(samples.oci_view == 0)  # smaller absolute offset, no double counting
    np.testing.assert_allclose(samples.raw_ratio, 235 / 74, rtol=1e-7)
    strict, _ = m.match_pair(hp, op, replace(cfg, oci_qc_policy="require_zero"))
    assert strict.sizes["sample"] == 0
    later, _ = m.match_pair(hp, op, replace(cfg, start_date="2024-09-04"))
    assert later.sizes["sample"] == 0
    o["geolocation_data"]["latitude"].values[0, 0] = 1
    o.to_netcdf(op, engine="h5netcdf")
    with pytest.raises(ValueError, match="grids differ"):
        m.match_pair(hp, op, cfg)


def test_guarded_midnight_epoch_repair(tmp_path):
    hp, op, o = make_pair(tmp_path)
    cfg = m.MatchConfig()
    # OCI global coverage still has the correct date, but CF-decoded rows are
    # uniformly one day early, as in the observed September 5 source granules.
    o["bin_attributes"]["nadir_view_time"] = o["bin_attributes"][
        "nadir_view_time"
    ] - np.timedelta64(1, "D")
    o.to_netcdf(op, engine="h5netcdf")
    with pytest.warns(RuntimeWarning, match="one-day OCI"):
        samples, report = m.match_pair(hp, op, cfg)
    assert samples.sizes["sample"] == 3
    assert report["oci_nadir_time_correction_seconds"] == 86400
    # NASA writes offset = nadir - observation, so a -10 s offset means later.
    expected = samples.harp_nadir_time.values + np.timedelta64(10, "s")
    np.testing.assert_array_equal(samples.oci_observation_time.values, expected)
    # A sub-day mismatch must still fail; never just loosen the tolerance.
    o["bin_attributes"]["nadir_view_time"] = o["bin_attributes"][
        "nadir_view_time"
    ] + np.timedelta64(2, "s")
    o.to_netcdf(op, engine="h5netcdf")
    with pytest.raises(ValueError, match="nadir-time mismatch"):
        m.match_pair(hp, op, cfg)


def test_day_correction_requires_independent_coverage(tmp_path):
    hp, op, o = make_pair(tmp_path)
    o["bin_attributes"]["nadir_view_time"] = o["bin_attributes"][
        "nadir_view_time"
    ] - np.timedelta64(1, "D")
    o.attrs["time_coverage_start"] = "2024-09-02T00:12:41Z"
    o.attrs["time_coverage_end"] = "2024-09-02T00:17:41Z"
    o.to_netcdf(op, engine="h5netcdf")
    with pytest.raises(ValueError, match="nadir-time mismatch"):
        m.match_pair(hp, op, m.MatchConfig())
