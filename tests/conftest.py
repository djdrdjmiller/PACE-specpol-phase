"""Offline numerical, selection, and resume tests for the liquid-reference extension."""

import numpy as np
import pytest
from pyhdf.SD import SD, SDC
from pace_specpol.liquid_reference import LiquidConfig


def hdf(path, data):
    f = SD(str(path), SDC.WRITE | SDC.CREATE)
    try:
        for name, value in data.items():
            v = np.asarray(value, dtype="float32")
            s = f.create(name, SDC.FLOAT32, v.shape)
            s[:] = v
            s.endaccess()
    finally:
        f.end()


@pytest.fixture
def config(tmp_path):
    wave = [0.645, 0.865, 1.25, 1.616, 2.13, 2.26]
    shape = (2, 2, 2, 2, 6, 2)
    ms = tmp_path / "ocean_msr_test.hdf"
    pf = tmp_path / "phase.hdf"
    tr = tmp_path / "trans.hdf"
    data = {
        "Wavelengths": wave,
        "OpticalThickness": [0, 10],
        "ParticleRadius": [5, 15],
        "ReflectanceSolarZenith": [0.3, 1],
        "ReflectanceSensorZenith": [0.3, 1],
        "ReflectanceRelativeAzimuth": [0, 180],
        "ExtinctionCoefficient": np.ones((2, 6)),
        "SingleScatterAlbedo": np.full((2, 6), 0.8),
        "TruncationFactor": np.full((2, 6), 0.2),
        "PhaseFuncNormConstant": np.full((2, 6), 2),
        "MultiScatBDReflectance": np.full(shape, 0.1),
    }
    hdf(ms, data)
    hdf(
        pf,
        {
            "WaveLengths": wave,
            "ParticleRadiusWater": [5, 15],
            "ScatAnglesWater": [0, 180],
            "WaterPhaseFuncVals": np.full((2, 6, 2), 4),
        },
    )
    p = np.arange(100, 1001, 100)
    u = np.arange(1, 21) / 20
    w = np.r_[np.arange(11) * 0.02, 0.2 + np.arange(1, 43) * 0.2]
    t = np.empty((8, 20, 53, 10))
    for b in range(8):
        t[b] = (
            0.8
            + 0.01 * u[:, None, None]
            - 0.005 * w[None, :, None]
            + 0.00001 * p[None, None, :]
        )
    hdf(
        tr,
        {
            "Cloud_Top_Pressure": p,
            "View_Angle_Cosine": u,
            "Precipitable_Water": w,
            "Wavelengths": [0.64, 0.86, 0.94, 1.2, 1.38, 1.6, 2.1, 2.23],
            "Transmittance": t,
            "counts": np.ones(t.shape),
        },
    )
    return LiquidConfig(str(ms), str(pf), str(tr), ocean_only=True)
