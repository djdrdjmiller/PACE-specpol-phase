"""265 K OCI all-liquid reference. Source conventions and limits: Documentation/algorithm_details.md.

The ocean MS LUT is read one spectral band at a time (~88 MiB/band).
Single scattering is constructed at native COT/CER nodes before interpolation.
No Mie calculations or temperature interpolation are performed here.
"""

from dataclasses import dataclass, asdict
from pathlib import Path
import hashlib
import numpy as np
from scipy.interpolate import RegularGridInterpolator

VERSION = "0.1.0"
SOURCE_COMMIT = "5d057d85310e594a2efc4ccf5ff01dfd43954c14"


def sha256(path):
    h = hashlib.sha256()
    with Path(path).expanduser().open("rb") as f:
        for b in iter(lambda: f.read(2**20), b""):
            h.update(b)
    return h.hexdigest()


def hdf_read(path, names):
    from pyhdf.SD import SD, SDC

    f = SD(str(Path(path).expanduser()), SDC.READ)
    try:
        result = {}
        for name in names:
            s = f.select(name)
            try:
                result[name] = np.asarray(s.get(), dtype="float64")
            finally:
                s.endaccess()
        return result
    finally:
        f.end()


def brackets(axis, x):
    x = np.asarray(x, dtype=float)
    lo = np.clip(np.searchsorted(axis, x, side="right") - 1, 0, len(axis) - 2)
    hi = lo + 1
    w = (x - axis[lo]) / (axis[hi] - axis[lo])
    return lo, hi, w


def scattering_angle(sza, vza, lut_raa):
    """CHIMAERA ScatAngle; lut_raa = 180 - folded solar/sensor azimuth difference."""
    mu0, mu = np.cos(np.deg2rad(sza)), np.cos(np.deg2rad(vza))
    c = -mu0 * mu + np.sqrt(np.maximum(0, (1 - mu0**2) * (1 - mu**2))) * np.cos(
        np.deg2rad(lut_raa)
    )
    return np.rad2deg(np.arccos(np.clip(c, -1, 1)))


@dataclass(frozen=True)
class LiquidConfig:
    ms_path: str
    phase_path: str
    transmittance_path: str
    cer_source: str = "oci"
    oci_reference_band: int = 2260
    # Ice CER has a different physical interpretation; retain with explicit flags.
    oci_phase_policy: str = "all"  # all | liquid; liquid includes phase code 2 only
    ocean_only: bool = False
    batch_samples: int = 8192

    def validate(self):
        if self.cer_source not in ("oci", "harp2"):
            raise ValueError("cer_source: oci or harp2")
        if self.oci_reference_band not in (2130, 2260):
            raise ValueError("reference band: 2130 or 2260")
        if self.oci_phase_policy not in ("all", "liquid"):
            raise ValueError("phase policy: all or liquid")
        if self.batch_samples < 1:
            raise ValueError("batch_samples must be positive")
        for p in (self.ms_path, self.phase_path, self.transmittance_path):
            if not Path(p).expanduser().is_file():
                raise FileNotFoundError(Path(p).expanduser())

    def metadata(self):
        self.validate()
        return dict(
            version=VERSION,
            config=asdict(self),
            hashes={
                k: sha256(getattr(self, k))
                for k in ("ms_path", "phase_path", "transmittance_path")
            },
            source_commit=SOURCE_COMMIT,
            temperature_K=265,
            surface="single supplied ocean wind-speed LUT",
            interpolation="native mu/mu0/azimuth; SS at tau/re nodes, then linear tau/re",
            pressure_interpolation="nearest pressure; bilinear PW/mu_effective",
        )


class LiquidLUT:
    def __init__(self, config):
        config.validate()
        self.cfg = config
        keys = [
            "Wavelengths",
            "OpticalThickness",
            "ParticleRadius",
            "ReflectanceSolarZenith",
            "ReflectanceSensorZenith",
            "ReflectanceRelativeAzimuth",
            "ExtinctionCoefficient",
            "SingleScatterAlbedo",
            "TruncationFactor",
            "PhaseFuncNormConstant",
        ]
        self.d = hdf_read(config.ms_path, keys)
        self.p = hdf_read(
            config.phase_path,
            [
                "WaveLengths",
                "ParticleRadiusWater",
                "ScatAnglesWater",
                "WaterPhaseFuncVals",
            ],
        )
        for a, b in [
            ("Wavelengths", "WaveLengths"),
            ("ParticleRadius", "ParticleRadiusWater"),
        ]:
            if not np.allclose(self.d[a], self.p[b]):
                raise ValueError("MS and phase grids disagree")
        self.axes = tuple(
            self.d[k]
            for k in [
                "ReflectanceSensorZenith",
                "ReflectanceSolarZenith",
                "ReflectanceRelativeAzimuth",
                "OpticalThickness",
                "ParticleRadius",
            ]
        )
        for x in (*self.axes, self.p["ScatAnglesWater"]):
            if not np.all(np.diff(x) > 0):
                raise ValueError("LUT axes must increase strictly")
        if not np.allclose(
            self.d["Wavelengths"], [0.645, 0.865, 1.25, 1.616, 2.13, 2.26], atol=0.002
        ):
            raise ValueError(
                "This reader is validated for the supplied six-band OCI LUT"
            )
        if not Path(config.ms_path).name.startswith("ocean_msr_"):
            raise ValueError(
                "Use an ocean MS LUT; Lambertian albedo correction is not implemented"
            )
        self._band = None
        self._interp = None

    def _load_band(self, index):
        if self._band == index:
            return
        # Release previous band BEFORE allocating another one.
        self._interp = None
        from pyhdf.SD import SD, SDC

        f = SD(str(Path(self.cfg.ms_path).expanduser()), SDC.READ)
        try:
            s = f.select("MultiScatBDReflectance")
            try:
                dims = tuple(s.info()[2])
                expected = tuple(len(x) for x in self.axes[:4]) + (6, len(self.axes[4]))
                if dims != expected:
                    raise ValueError(f"Unexpected MS shape: {dims}")
                count = list(dims)
                count[4] = 1
                data = s.get(start=(0, 0, 0, 0, index, 0), count=tuple(count)).squeeze(
                    axis=4
                )
            finally:
                s.endaccess()
        finally:
            f.end()
        self._interp = RegularGridInterpolator(
            self.axes, data, bounds_error=False, fill_value=np.nan
        )
        self._band = index

    def reflectance(self, sza, vza, raa, cot, cer, band_nm):
        """raa is the folded solar/sensor azimuth difference in the existing cache.

        Ocean below-cloud Rayleigh/aerosol SS terms are omitted for OCI SWIR by
        NASA's Single_Scattering_Calcs_Ocean (wave_limit=number_wavelengths-3).
        """
        index = {1615: 3, 2130: 4, 2260: 5}[band_nm]
        self._load_band(index)
        sza, vza, raa, cot, cer = np.broadcast_arrays(
            *[np.asarray(x, dtype=float) for x in (sza, vza, raa, cot, cer)]
        )
        shape = sza.shape
        sza, vza, raa, cot, cer = [x.ravel() for x in (sza, vza, raa, cot, cer)]
        mu0, mu = np.cos(np.deg2rad(sza)), np.cos(np.deg2rad(vza))
        az = 180 - raa
        ang = scattering_angle(sza, vza, az)
        out = np.full(len(cot), np.nan)
        for start in range(0, len(cot), self.cfg.batch_samples):
            sl = slice(start, start + self.cfg.batch_samples)
            ms = self._interp(
                np.column_stack((mu[sl], mu0[sl], az[sl], cot[sl], cer[sl]))
            )
            tl, th, tw = brackets(self.axes[3], cot[sl])
            rl, rh, rw = brackets(self.axes[4], cer[sl])
            al, ah, aw = brackets(self.p["ScatAnglesWater"], ang[sl])
            ss = np.zeros(ms.shape)
            for ti, wt in [(tl, 1 - tw), (th, tw)]:
                for ri, wr in [(rl, 1 - rw), (rh, rw)]:
                    phase = (1 - aw) * self.p["WaterPhaseFuncVals"][
                        ri, index, al
                    ] + aw * self.p["WaterPhaseFuncVals"][ri, index, ah]
                    phase = phase / self.d["PhaseFuncNormConstant"][ri, index]
                    omega = self.d["SingleScatterAlbedo"][ri, index]
                    f = self.d["TruncationFactor"][ri, index]
                    scale = 1 - f * omega
                    tau = (
                        self.axes[3][ti]
                        * self.d["ExtinctionCoefficient"][ri, index]
                        / self.d["ExtinctionCoefficient"][ri, 0]
                    )
                    exponent = tau * scale * (1 / mu[sl] + 1 / mu0[sl])
                    # NASA shortcut exp(-x)=0 for x>10, preserved here.
                    expval = np.exp(-np.minimum(exponent, 700))
                    expval[exponent > 10] = 0
                    ss += (
                        wt
                        * wr
                        * omega
                        * (phase / scale)
                        * (1 - expval)
                        / (4 * (mu[sl] + mu0[sl]))
                    )
            out[sl] = ms + ss
        ok = (
            (sza >= 0)
            & (sza < 90)
            & (vza >= 0)
            & (vza < 90)
            & (raa >= 0)
            & (raa <= 180)
        )
        ok &= (
            (cot >= self.axes[3][0])
            & (cot <= self.axes[3][-1])
            & (cer >= self.axes[4][0])
            & (cer <= self.axes[4][-1])
        )
        out[~ok] = np.nan
        return out.reshape(shape)


class TransmissionLUT:
    def __init__(self, path):
        self.d = hdf_read(
            path,
            [
                "Cloud_Top_Pressure",
                "View_Angle_Cosine",
                "Precipitable_Water",
                "Wavelengths",
                "Transmittance",
                "counts",
            ],
        )
        if not np.allclose(
            self.d["Wavelengths"],
            [0.64, 0.86, 0.94, 1.2, 1.38, 1.6, 2.1, 2.23],
            atol=0.005,
        ):
            raise ValueError("Unexpected transmission band ordering")
        self.p = self.d["Cloud_Top_Pressure"]
        self.u = self.d["View_Angle_Cosine"]
        self.w = self.d["Precipitable_Water"]
        if not np.allclose(self.p, np.arange(100, 1001, 100)):
            raise ValueError("Unexpected pressure grid")

    def evaluate(self, sza, vza, pressure_hpa, water_g_cm2, band_nm):
        """GetTransmittanceData convention, with strict rejection of missing corners.

        Deliberate difference: do NOT extrapolate across unpopulated PW bins or
        silently clamp pressure. Return NaN instead; diagnostic counts retain loss.
        """
        b = {1615: 5, 2130: 6, 2260: 7}[band_nm]
        sza, vza, p, w = np.broadcast_arrays(
            *[np.asarray(x, float) for x in (sza, vza, pressure_hpa, water_g_cm2)]
        )
        mu0 = np.cos(np.deg2rad(sza))
        mu = np.cos(np.deg2rad(vza))
        u = mu0 * mu / (mu0 + mu)
        pi = np.clip(
            np.floor((np.nan_to_num(p, nan=100) - 100) / 100 + 0.5).astype(int), 0, 9
        )
        ul, uh, uw = brackets(self.u, u)
        wl, wh, ww = brackets(self.w, w)
        out = np.zeros(p.shape)
        valid = (
            np.isfinite(p + w + u)
            & (p >= 100)
            & (p <= 1000)
            & (w >= 0)
            & (w <= self.w[-1])
            & (u >= self.u[0])
            & (u <= self.u[-1])
        )
        for ui, wu in [(ul, 1 - uw), (uh, uw)]:
            for wi, wt in [(wl, 1 - ww), (wh, ww)]:
                v = self.d["Transmittance"][b, ui, wi, pi]
                n = self.d["counts"][b, ui, wi, pi]
                weight = wu * wt
                needed = weight > 1e-12
                valid &= (~needed) | ((n > 0) & np.isfinite(v) & (v > 0) & (v <= 1))
                out += np.where(needed, v * weight, 0)
        return np.where(valid, out, np.nan)


def correct_samples(ds, config, lut=None, transmission=None):
    """Input includes oci_cer/cot_21/22, oci_phase_21/22, ctp, above_cloud_water.
    Preserves every sample and records why each mode excludes samples.
    """
    lut = lut or LiquidLUT(config)
    transmission = transmission or TransmissionLUT(config.transmittance_path)
    suffix = {2130: "21", 2260: "22"}[config.oci_reference_band]
    required = [
        "sza",
        "vza",
        "raa",
        "raw_ratio",
        "rho1615",
        "rho2260",
        "ctp",
        "above_cloud_water",
        f"oci_cot_{suffix}",
        f"oci_cer_{suffix}",
        f"oci_phase_{suffix}",
        "oci_ocean",
    ]
    for name in required:
        if name not in ds:
            raise ValueError(f"Missing companion field: {name}")
    cot = ds[f"oci_cot_{suffix}"].values
    cer = (
        ds[f"oci_cer_{suffix}"].values
        if config.cer_source == "oci"
        else ds["re"].values
    )
    phase = ds[f"oci_phase_{suffix}"].values
    sza, vza, raa = [ds[k].values for k in ["sza", "vza", "raa"]]
    ctp, pw = [ds[k].values for k in ["ctp", "above_cloud_water"]]
    r16 = lut.reflectance(sza, vza, raa, cot, cer, 1615)
    r22 = lut.reflectance(sza, vza, raa, cot, cer, 2260)
    t16 = transmission.evaluate(sza, vza, ctp, pw, 1615)
    t22 = transmission.evaluate(sza, vza, ctp, pw, 2260)
    status = np.zeros(ds.sizes["sample"], dtype="uint16")
    status[
        ~(np.isfinite(cer) & (cer >= lut.axes[4][0]) & (cer <= lut.axes[4][-1]))
    ] |= 1
    status[
        ~(np.isfinite(cot) & (cot >= lut.axes[3][0]) & (cot <= lut.axes[3][-1]))
    ] |= 2
    status[~(np.isfinite(t16) & np.isfinite(t22))] |= 4
    status[~(np.isfinite(r16) & np.isfinite(r22) & (r16 > 0) & (r22 > 0))] |= 8
    if config.ocean_only:
        status[ds.oci_ocean.values != 1] |= 16
    if config.oci_phase_policy == "liquid":
        status[phase != 2] |= 32
    else:
        status[~np.isin(phase, [2, 3, 4])] |= 32
    if "oci_match_valid" in ds:
        status[ds.oci_match_valid.values != 1] |= 64
    out = ds.copy()

    def put(k, v):
        out[k] = ("sample", np.asarray(v, dtype="float32"))

    with np.errstate(divide="ignore", invalid="ignore"):
        liquid = r22 / r16
        corrected = ds.raw_ratio.values * t16 / t22
        norm = corrected / liquid
        baseline = liquid * t22 / t16  # Existing dashboard divides raw_ratio by this.
    status[~np.isfinite(norm)] |= 8
    for k, v in {
        "reference_cer": cer,
        "reference_cot": cot,
        "liquid_rho1615": r16,
        "liquid_rho2260": r22,
        "transmission1615": t16,
        "transmission2260": t22,
        "liquid_ratio": liquid,
        "corrected_ratio": corrected,
        "normalized_ratio": np.where(status == 0, norm, np.nan),
        "effective_reference": np.where(status == 0, baseline, np.nan),
    }.items():
        put(k, v)
    out["reference_status"] = ("sample", status)
    out.reference_status.attrs.update(
        flag_masks=np.array([1, 2, 4, 8, 16, 32, 64], dtype="uint16"),
        flag_meanings="invalid_cer invalid_cot invalid_transmission invalid_reflectance_or_ratio non_ocean phase_policy unmatched_oci",
    )
    out["ice_cer_as_liquid"] = (
        "sample",
        ((phase == 3) & (config.cer_source == "oci")).astype("int8"),
    )
    return out


@dataclass(frozen=True)
class CachedReference:
    """Adapter for the existing build_index/dashboard, without altering raw_ratio."""

    mode: str = "normalized"
    population: str = "own"  # own | common (common_valid supplied by comparison cache)

    def metadata(self):
        if self.mode not in ("normalized", "corrected", "raw"):
            raise ValueError("invalid mode")
        if self.population not in ("own", "common"):
            raise ValueError("invalid population")
        return {
            "mode": self.mode,
            "population": self.population,
            "version": VERSION,
            "description": "265 K liquid reference; provenance in source manifest",
        }

    def evaluator(self):
        self.metadata()

        def evaluate(ds):
            if self.mode == "normalized":
                v = ds.effective_reference.values.astype(float).copy()
            elif self.mode == "corrected":
                v = (ds.transmission2260 / ds.transmission1615).values
            else:
                v = np.ones(ds.sizes["sample"])
            if self.population == "common":
                v = np.where(ds.common_valid.values == 1, v, np.nan)
            return v

        return evaluate
