"""Offline numerical, selection, and resume tests for the liquid-reference extension."""

import json
from dataclasses import replace
import numpy as np
import pytest
import xarray as xr
from pace_specpol.workflow import cache_references, VariantReference


from .sample_data import samples


def test_cache_resume_and_common_population(config, tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    samples().to_netcdf(source / "x.nc", engine="h5netcdf")
    (source / "manifest.json").write_text(
        json.dumps({"complete": True, "config": {}, "files": ["x.nc"]})
    )
    variants = {"oci_2260": config, "harp2_2260": replace(config, cer_source="harp2")}
    out = cache_references(source, tmp_path / "out", variants)
    before = (out / "x.nc").stat().st_mtime_ns
    cache_references(source, tmp_path / "out", variants)
    assert (out / "x.nc").stat().st_mtime_ns == before
    with xr.open_dataset(out / "x.nc") as ds:
        np.testing.assert_array_equal(ds.common_valid, [1, 0, 0])
        v = VariantReference("oci_2260", "common").evaluator()(ds)
        assert np.isfinite(v).tolist() == [True, False, False]
    with pytest.raises(ValueError):
        cache_references(
            source, out, {"oci_2260": replace(config, oci_phase_policy="liquid")}
        )


def test_companion_atomic_resume(tmp_path):
    from pace_specpol.oci_companion import CompanionConfig
    from pace_specpol.workflow import cache_companions

    class Provider:
        cfg = CompanionConfig()
        catalog = []
        calls = 0

        def __call__(self, ds):
            self.calls += 1
            return ds

    source = tmp_path / "raw"
    source.mkdir()
    samples().to_netcdf(source / "x.nc", engine="h5netcdf")
    (source / "manifest.json").write_text(
        json.dumps({"complete": True, "config": {}, "files": ["x.nc"]})
    )
    provider = Provider()
    out = cache_companions(source, tmp_path / "companions", provider)
    cache_companions(source, out, provider)
    assert provider.calls == 1
    # A changed source file must never silently reuse the derived companion.
    ds = samples()
    ds.attrs["changed"] = "yes"
    ds.to_netcdf(source / "x.nc", engine="h5netcdf")
    with pytest.raises(ValueError, match="Stale companion"):
        cache_companions(source, out, provider)
