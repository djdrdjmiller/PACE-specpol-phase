# PACE spectropolarimetric cloud phase

Research and development algorithm combining PACE HARP2 polarimetry and OCI
spectral reflectance observations to explore liquid, ice, and liquid-topped
mixed-phase (LTMP) cloud candidates. The method combines polarimetric evidence
for liquid at the optical cloud top with spectral departures from a simulated
all-liquid cloud. The spectral component is intended to detect the influence of
ice: differences between ice and liquid water in the wavelength dependence of
their complex refractive indices affect absorption and scattering across the
1615 and 2260 nm bands and can shift the observed reflectance ratio away from
its all-liquid expectation. It is an exploratory algorithm;
the phase thresholds have not been independently calibrated for OCI.

## Algorithm

1. Pair refined, non-NRT **HARP2 GPC V4.0** with **OCI L1C V3** by granule
   timestamp, then verify geolocation and nadir timing. Read HARP2 liquid index
   (LI); do not recompute it. Select both OCI SWIR bands from the same view and
   convert radiances to TOA reflectances.
2. Form the observed ratio `R_TOA = rho2260_TOA / rho1615_TOA`.
3. Match time-compatible native **OCI CLD V3.1** microphysics to each L1C sample.
   Retrieve the original ancillary inputs named by that cloud product and
   estimate above-cloud water. Evaluate the supplied transmission table to obtain
   the effective two-way transmissions at each band.
4. Correct atmospheric attenuation: `R_corrected = R_TOA * T1615 / T2260`.
   The effective-path transmission is already downward × upward; do not square it.
5. Reconstruct both liquid-cloud reflectances from the supplied **265 K** optical
   properties and ocean LUT, using observation geometry and selected CER/COT.
   Form `R_liquid = rho2260_liquid / rho1615_liquid`, then
   `R_normalized = R_corrected / R_liquid`.
6. Classify each eligible observation using LI and the selected ratio. Retain
   continuous evidence, exclusion flags, and provenance alongside phase labels.

The default reference uses OCI 2260-derived effective radius (CER) and its
corresponding optical thickness (COT). Alternatives use OCI 2130 CER/COT, or
HARP2 CER with either OCI COT. Missing HARP2 CER is never filled with OCI CER.
Reference comparisons can use each variant's valid samples or their common
intersection. Raw, corrected, and normalized views of one selected variant use
the same reference-valid population.

The current provisional rule is:

| Per-observation evidence | Candidate class |
|---|---|
| Ratio < cutoff | Liquid, regardless of LI |
| Ratio ≥ cutoff and LI < LI cutoff | Ice |
| Ratio ≥ cutoff and LI ≥ LI cutoff | LTMP |

Defaults are LI = **0.3** and ratio = **1.27**. The ratio threshold originated in
an RSP demonstration and is not an OCI calibration. Its interpretation differs
between raw, corrected, and normalized spaces. A normalized ratio of one means
agreement with the modeled liquid reference; it is not proof of an all-liquid
column. Neither LI nor normalized ratio is an ice or liquid mass fraction.

The current processing unit is the **matched L1C sample**, approximately 5 km.
Native OCI microphysics are transferred from a nearest time-compatible pixel,
not averaged over that footprint. Future processing on OCI pixels with coarser
polarimeter information requires an explicit collocation design. `Operational/`
is reserved for that future development and currently contains no algorithm.

## Scientific limits

The selected 3 m/s ocean surface model is deliberately applied globally,
including land, snow, and ice. Near-surface profiles requiring operational
repair and unpopulated transmission corners are excluded. Cloud parallax,
advection, and 3D radiative effects are not corrected. Using an OCI ice-retrieval
radius in a liquid LUT is recorded as a numerical counterfactual, not interpreted
as a retrieved liquid radius. OCI2260 microphysics also depend partly on the
ratio's numerator band. These choices constrain interpretation of departures
from the liquid reference.

See [algorithm conventions and equations](Documentation/algorithm_details.md),
[input data and LUTs](Documentation/data_and_luts.md), and
[validation scope](Validation/README.md). The tests verify implementation
behavior; they do not establish retrieval skill.

## Install and run

Use Python 3.11 or newer. From the repository root, in the environment used by
your notebook kernel:

```bash
python -m pip install -e ".[validation,test]"
cp config/paths.example.toml config/paths.local.toml
```

Edit the local file to point to your LUTs and experiment-specific caches. For
the June 10 Arctic workflow use `config/arcsix.example.toml` as the template for
`config/arcsix.local.toml`. Select it explicitly before starting Jupyter:

```bash
export PACE_CONFIG="/absolute/path/to/pace-specpol-phase/config/arcsix.local.toml"
jupyter lab
```

Alternatively, set `os.environ["PACE_CONFIG"]` in the notebook before
`load_paths()`. Relative data paths resolve against the TOML file's directory.
The loader never moves data. Do not use a September cache for the Arctic run.
HDF4 LUTs require `pyhdf` and a working HDF4 installation; `h5py` alone cannot
read them. Earthdata authentication is needed for remote discovery/downloads.
Perform large OCI extraction near the data in AWS us-west-2 when possible.

Open an [analysis notebook](Validation/README.md). Importable processing code
lives under `src/pace_specpol/`; visualization lives in its `validation/vis_*`
modules. Existing notebooks that import the former flat modules should use the
migrated versions; see [migration notes](Documentation/migration.md).

Run offline tests with `python -m pytest -q`. Scientific data and LUT binaries
are excluded from Git. See the
[GitHub repository](https://github.com/djdrdjmiller/PACE-specpol-phase)
for its license.
