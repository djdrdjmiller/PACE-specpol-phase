import numpy as np
import xarray as xr
import pytest
from pace_specpol.reference import Reference


def test_lut_no_extrapolation(tmp_path):
    lut = xr.Dataset(
        {"liquid_reflectance_ratio": ("sza", [1.0, 3.0])},
        coords={"sza": [0.0, 90.0]},
        attrs={"numerator_nominal_nm": 2260, "denominator_nominal_nm": 1615},
    )
    path = tmp_path / "lut.nc"
    lut.to_netcdf(path, engine="h5netcdf")
    evaluate = Reference(
        mode="lut", lut_path=str(path), description="Synthetic test"
    ).evaluator()
    result = evaluate(xr.Dataset({"sza": ("sample", [0.0, 45.0, 90.0, 100.0, np.nan])}))
    np.testing.assert_allclose(result[:3], [1, 2, 3])
    assert np.isnan(result[3:]).all()
    with pytest.raises(ValueError):
        Reference(mode="scalar", scalar=0).evaluator()
