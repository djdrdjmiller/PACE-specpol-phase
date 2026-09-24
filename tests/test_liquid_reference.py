"""Offline numerical, selection, and resume tests for the liquid-reference extension."""

from dataclasses import replace
import numpy as np
from pace_specpol.liquid_reference import (
    LiquidLUT,
    TransmissionLUT,
    correct_samples,
    scattering_angle,
)


from .sample_data import samples


def test_ss_normalization_truncation_and_node_interpolation(config):
    lut = LiquidLUT(config)
    # overhead sun and nadir: airmass 2; exp shortcut >10 ->0.
    full = 0.1 + 0.8 * (4 / 2) / (1 - 0.2 * 0.8) / 8
    result = lut.reflectance(
        [0, 0, 0], [0, 0, 0], [0, 0, 0], [10, 0, 5], [10, 10, 10], 2260
    )
    np.testing.assert_allclose(result, [full, 0.1, 0.1 + 0.5 * (full - 0.1)], rtol=1e-6)
    assert np.isnan(lut.reflectance([0], [0], [0], [10], [35], 2260))[0]


def test_geometry_convention():
    # Sun and sensor in the same outward direction -> backscatter.
    np.testing.assert_allclose(scattering_angle(30, 30, 180), 180, atol=1e-6)
    np.testing.assert_allclose(scattering_angle(30, 30, 0), 120, atol=1e-6)


def test_transmission_effective_path_and_nearest_pressure(config):
    t = TransmissionLUT(config.transmittance_path)
    result = t.evaluate([0, 0], [0, 0], [749, 750], [1, 1], 2260)
    np.testing.assert_allclose(
        result, [0.8 + 0.005 - 0.005 + 0.007, 0.8 + 0.005 - 0.005 + 0.008], rtol=1e-6
    )
    assert np.isnan(t.evaluate(0, 0, 1100, 1, 2260))
    assert np.isnan(t.evaluate(0, 0, 700, 9, 2260))
    # Do not invent data in unpopulated LUT corners.
    t.d["counts"][:] = 0
    assert np.isnan(t.evaluate(0, 0, 700, 1, 2260))


def test_ratio_orientation_and_status(config):
    ds = samples()
    r = correct_samples(ds, config)
    np.testing.assert_allclose(
        r.normalized_ratio[0],
        ds.raw_ratio[0]
        * r.transmission1615[0]
        / r.transmission2260[0]
        / r.liquid_ratio[0],
    )
    assert r.reference_status.values.tolist() == [0, 9, 16]
    r = correct_samples(ds, replace(config, cer_source="harp2"))
    assert np.isnan(r.normalized_ratio[1])
    r = correct_samples(
        ds, replace(config, oci_phase_policy="liquid", ocean_only=False)
    )
    assert r.reference_status.values[2] & 32
