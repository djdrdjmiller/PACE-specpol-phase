import json
import numpy as np
import xarray as xr
import pytest
from pace_specpol.validation import aggregation as m
from pace_specpol.classification import sample_phase


def test_sparse_counts_and_reference(tmp_path):
    ratio = np.array([0.1, 1.26, 1.27, 1.27, 2, 4, 1.5, 1.5, 1.5])
    li = np.array([-0.2, 2, 0.29, 0.3, 2, 0, 0, 0, 1.0])
    ds = xr.Dataset(
        {
            "raw_ratio": ("sample", ratio),
            "liquid_index": ("sample", li),
            "latitude": ("sample", np.zeros(9)),
            "longitude": ("sample", [0] * 6 + [20] * 3),
            "sza": ("sample", np.full(9, 45.0)),
        }
    )
    ds.to_netcdf(tmp_path / "samples.nc", engine="h5netcdf")
    (tmp_path / "manifest.json").write_text(
        json.dumps({"complete": True, "files": ["samples.nc"], "config": {}})
    )
    index = m.build_index(tmp_path, config=m.IndexConfig(resolution=1))
    for cut in [0.5, 1.26, 1.27, 1.28, 2, 3.5]:
        classes = sample_phase(ratio, li, cut)
        cells = m.cell_ids(ds.latitude, ds.longitude, 1)
        expected = np.stack(
            [np.bincount(cells[classes == k], minlength=64800) for k in (1, 2, 3)]
        )
        np.testing.assert_array_equal(index.counts(cut), expected)
    assert index.metadata["totals"]["histogram_outside_samples"] == 2
    with pytest.raises(ValueError):
        index.counts(1.275)
    norm = m.build_index(
        tmp_path,
        reference=m.Reference(mode="scalar", scalar=2, description="Synthetic test"),
        config=m.IndexConfig(resolution=1),
    )
    np.testing.assert_array_equal(
        norm.counts(1.27).sum(axis=1),
        np.bincount(sample_phase(ratio / 2, li), minlength=4)[1:],
    )
    index.export(tmp_path / "result.nc")
    with xr.open_dataset(tmp_path / "result.nc") as reread:
        assert reread.sample_count.sum() == 9
        np.testing.assert_allclose(
            reread.phase_fraction.sum("phase").values[reread.sample_count.values > 0], 1
        )
    loaded = m.PhaseIndex.load(index.directory)
    np.testing.assert_array_equal(loaded.counts(1.27), index.counts(1.27))


def test_date_line_and_poles():
    ids = m.cell_ids([0, 0, 90, -90], [180, -180, 180, -180], 1)
    assert ids[0] == ids[1]
    assert ids[2] == 179 * 360
    assert ids[3] == 0


def test_continuous_means_and_reload(tmp_path):
    # Unequal ratios and negative LI expose clipping, sum/count and bin-center errors.
    ds = xr.Dataset(
        {
            "raw_ratio": ("sample", [0.1, 4.0, 2.0]),
            "liquid_index": ("sample", [-1.0, 5.0, 8.0]),
            "latitude": ("sample", [0.0, 0.0, 0.0]),
            "longitude": ("sample", [0.0, 0.0, 20.0]),
        }
    )
    ds.to_netcdf(tmp_path / "samples.nc", engine="h5netcdf")
    (tmp_path / "manifest.json").write_text(
        json.dumps({"complete": True, "files": ["samples.nc"], "config": {}})
    )
    index = m.build_index(
        tmp_path,
        reference=m.Reference(mode="scalar", scalar=2, description="Test"),
        config=m.IndexConfig(resolution=1),
    )
    loaded = m.PhaseIndex.load(index.directory)
    cell = m.cell_ids([0.0], [0.0], 1)[0]
    means, total = loaded.mean_fields()
    np.testing.assert_allclose(means[:, cell], [2.0, 1.025])
    assert total[cell] == 2
    assert loaded.sampling_summary()["single_sample_cell_fraction"] == 0.5
    result = loaded.dataset(min_count=2)
    assert np.isfinite(result.mean_ratio).sum() == 1
    assert float(result.mean_liquid_index.sel(lat=0.5, lon=0.5)) == 2.0
