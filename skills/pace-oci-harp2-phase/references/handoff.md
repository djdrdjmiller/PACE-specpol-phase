# Research handoff — 2026-09-22

## Repository organization update

The source is now an installable `pace_specpol` package under `src/`. Analysis
notebooks live in `Validation/notebooks/`, with explicit `analysis_` prefixes.
Shared validation and visualization code lives in `src/pace_specpol/validation/`;
`Operational/` is reserved and contains only a Git placeholder. Current setup
is in the root README and `Documentation/migration.md`. Historical version and
cache conventions below remain intentional. Local data, LUTs, and caches are
external inputs configured through `PACE_CONFIG`; never hard-code a developer's
machine paths. Existing external notebooks need the documented import updates.

The user subsequently supplied Arctic corrected/normalized dashboard images
from their AWS run. Actual independently reviewed satellite-aircraft matchup
validation remains unestablished by those images. Historical statements below
about local Arctic testing describe the earlier local verification scope.

## Purpose and evolution

The original three-day global HARP2 liquid-index figure used calibration V3, which the user reports introduced a cross-swath LI decrease through georegistration/polarimetric-calibration problems. The aim is to reproduce September 3–5, 2024 with refined GPC V4.0, then combine LI with OCI spectral evidence for liquid-topped mixed-phase clouds. The date range is configurable. Preserve the user's distinction between detecting liquid at the optical cloud top and detecting ice influence within the spectrally sensed cloud.

Cloudbow LI comes from HARP2 GPC, not an LI recomputation. OCI supplies SWIR reflectances near 1615 and 2260 nm, rather than RSP's 1590/2260 pair. The initial RSP/ACTIVATE demonstration used an all-liquid-normalized ratio threshold 1.27. This number and LI=0.3 are starting choices, not validated OCI boundaries. The current map is a 3-category heuristic, not the 4-category classification shown in the collaborator's ARCSIX figure. No independent algorithm for its extra “mixed-phase” class has been implemented.

## Package map

All paths below are relative to the analysis package root.

| File | Responsibility |
|---|---|
| `src/pace_specpol/matching.py` | Refined HARP2/OCI L1C discovery and sample matching/cache |
| `src/pace_specpol/classification.py` | Per-observation provisional phase rule |
| `src/pace_specpol/reference.py` | Raw, scalar, and external ratio-reference adapters |
| `src/pace_specpol/validation/aggregation.py` | Sparse index, counts, fractions, and optional geographic aggregation |
| `src/pace_specpol/validation/vis_dashboard.py` | Histogram/maps/dashboard and reference selector |
| `src/pace_specpol/validation/vis_samples.py` and `vis_airborne.py` | Sample-center and airborne context plots |
| `src/pace_specpol/liquid_reference.py` | 265 K ocean MS/SS reconstruction, transmission LUT, correction status and normalized ratios |
| `src/pace_specpol/oci_companion.py` | Native OCI CLD discovery, time-compatible nearest pixel, original ancillary downloads, above-cloud water |
| `src/pace_specpol/workflow.py` | Resumable companion/reference stages; four reference variants and own/common population adapters |
| `src/pace_specpol/validation/arcsix.py` | Airborne discovery, actual HSRL reader, geometry/time matching, per-satellite summaries, LI tests, reviewed-label ratio fitting, ICARTT reader |
| `Validation/notebooks/analysis_paired_observations.ipynb` | Original paired-cache workflow |
| `Validation/notebooks/analysis_liquid_reference.ipynb` | Reference extension starting from a completed paired cache |
| `Validation/notebooks/analysis_arcsix_single_granule.ipynb` | Fresh June 10 pair extraction, reference analysis, appended airborne cells |
| `Validation/examples/airborne_validation_cells.py` | Copyable additions for an existing Arctic notebook; expects its stage/path variables |
| `Documentation/algorithm_details.md` | Detailed equations, conventions, dependencies, NASA source links, resource/implementation limits |
| `Validation/records/LIQUID_VALIDATION.json` | Actual September single-pair test record; not an Arctic result |
| `Documentation/arcsix_validation.md` / `Validation/records/AIRBORNE_VALIDATION.json` | Airborne code scope and verified file evidence |

Tests live in `tests/`, organized by matching, classification, reference reconstruction,
ancillary inputs, workflow resume, aggregation, airborne validation, and path handling.
Run `python -m pytest -q` from the repository root after installing the package and
test dependencies. Do not assume a developer's temporary virtualenv exists on AWS.

## Products, timestamps, locations

- CMR: `PACE_HARP2_L2_CLOUD_GPC`, version `4.0`; exact non-NRT names `PACE_HARP2.<stamp>.L2.CLOUD_GPC.V4_0.nc`.
- CMR: `PACE_OCI_L1C_SCI`, version `3`; names `PACE_OCI.<stamp>.L1C.V3.5km.nc`.
- Microphysics: refined `PACE_OCI_L2_CLOUD`, version `3.1`.
- Arctic figure target was verified in CMR: **20240610T154205**, 2024-06-10 15:42:05–15:47:04 UTC. The notebook selects exactly one paired L1C/GPC granule. Adjacent native OCI L2 files may still be needed for valid observation-time matching.
- User AWS LUT path example: `~/Data/LUT/OCI_Cloud/PACE_265_LUTS_2026`; expand `~` before HDF4 reads. Local development used an external directory named `PACE_265K_LUTS_2026`; the naming difference is real. Paths are configurable.
- LUTs are external inputs, not distributed in the package. Examples:
  - `LIQUID/ocean_msr_water_wspeed_3_v6.PACE.1.1.5.2026144071240.hdf`
  - `IceAndWaterPhaseFunctionData_v6.PACE.1.1.5.2026142144440.hdf`
  - `Transmittance_OCI.hdf`
- HDF4 requires `pyhdf`, not h5py. Airborne HSRL files are HDF5.

## LUT reconstruction details that must survive refactoring

Microphysical wavelengths are 0.645, 0.865, 1.25, 1.616, 2.13, 2.26 µm. CER spans 2–30 µm; COT spans 0–158.78 on nonuniform grids. The supplied MS array has shape (30,33,37,35,6,18): sensor cosine, solar cosine, relative azimuth, COT, wavelength, CER. Do not treat the cosine axes as angles in degrees.

MS is multilinearly interpolated on native axes. Cached folded solar/sensor azimuth separation becomes `lut_raa=180-raa`. With mu0 and mu, scattering angle follows `acos(-mu0*mu + sqrt(1-mu0²)*sqrt(1-mu²)*cos(lut_raa))`.

SS uses the phase function interpolated at actual scattering angle and divided by `PhaseFuncNormConstant`. Reconstruct at neighboring CER/COT nodes before interpolating those dimensions. Scale COT by extinction relative to the first wavelength; apply the supplied delta-truncation quantities. The implemented form is `omega*PF/(1-f*omega) * (1-exp(-tau_prime*(1/mu0+1/mu))) / (4*(mu0+mu))`, with `tau_prime=tau_lambda*(1-f*omega)`. Source exp(-x)=0 for x>10 is preserved. The inspected OCI ocean source omits additional below-cloud Rayleigh/aerosol SS for these SWIR bands. `stddev` files are not used as mean reflectance.

Reference variants:
- OCI2260: `cer_22`, `cot_22`; default.
- OCI2130: `cer_21`, `cot_21`.
- HARP2260/HARP2130: HARP2 `re`, with the corresponding OCI COT.
- An OCI ice-retrieval CER used numerically in a liquid LUT is explicitly a counterfactual. Keep `ice_cer_as_liquid`; do not call it a liquid CER retrieval.
- Do not clamp out-of-domain CER. HARP2 CER availability selects a different population; compare common-valid samples when attributing changes to reference mode.

## Transmittance and original ancillary handling

Transmittance shape is (8 wavelengths,20 cosines,53 water values,10 pressures). Labels are legacy 0.64,0.86,0.94,1.2,1.38,1.6,2.1,2.23 µm. NASA's explicit OCI remapping selects zero-based 5/6/7 for 1615/2130/2260; do not use nearest-label spectral interpolation or “correct” labels by assumption.

Evaluate with `mu_eff=mu0*mu/(mu0+mu)`: one LUT evaluation gives the combined downward/upward attenuation. Pressure is nearest 100 hPa, ties upward; interpolate water and effective cosine bilinearly. Above-cloud precipitable water is g/cm² (equivalent to cm liquid-water depth). Do not feed total-column water directly or confuse it with kg/m².

Read `anc_profile1/2/3`, `met1/2/3` from each OCI CLD file's processing input parameters. The inspected September files use GMAO MERRA2 PROFILE and MET inputs from OB.DAAC. Preserve original names/hashes; do not silently replace with forecast GEOS5 or a differently processed reanalysis. QV is kg/kg, surface pressure Pa. Convert q to g/kg mixing ratio, interpolate in log pressure to the source 101-level grid, and integrate above CTP with the documented 980.616 trapezoid factor.

Current limitations intentionally retained for the pilot:
- Near-surface profiles needing operational level repair/extrapolation are rejected.
- Missing profile values and unpopulated transmission corners are rejected instead of applying all NASA fallback logic.
- Native OCI microphysics/ancillary inputs are from one nearest time-compatible ~1 km pixel (3 km,10 s defaults), transferred to a ~5 km L1C sample. This is not footprint-averaged forward modeling.
- No new cloud-parallax correction.

NASA source commit inspected: `5d057d85310e594a2efc4ccf5ff01dfd43954c14`; inspected OCI L2 was produced with a different earlier l2gen release. Links and source filenames are in `Documentation/algorithm_details.md`. Do not claim bit-exact reproduction.

## Cache and timing pitfalls

NASA's OCI L1C offset is **nadir time minus observation time**. Extractor version 1.1 uses `observation=nadir-offset`; old 1.0 used the wrong sign. The companion reader repairs older completed caches using stored time minus twice the offset. Samples lost through the old date-boundary filter cannot be restored without extraction. Do not resume v1.1 raw extraction into a v1.0 directory.

Preserve exact timestamp pairing and grid/time checks. A `HARP2/OCI nadir-time mismatch` is not a reason to silently disable validation. Completed per-pair files can be reused after interruptions, but changed configurations require new directories.

AWS CryoCloud is used for data proximity. Distinguish persistent disk, managed download cache, and RAM. Only one ~88 MiB MS band should be resident; process individual granules. Default ancillary cache cap is 512 MiB, affecting only managed downloads—not derived caches, indexes, user files, or pip's cache. Do not imply the system automatically frees user data when disk fills. Kernel lifetime/idle policy is platform-dependent, not guaranteed by closing the laptop.

The L3 0.1-degree grid is a display/aggregation choice; it is not the input footprint. Cell width varies with latitude, and repeat sampling changes per-cell counts. Average sample ratios after correction/normalization as implemented; do not replace mean ratios with ratios of means without an explicit change. Phase is classified per observation before modal L3 display. Counts/fractions must remain available because modal phase can conceal rare LTMP/ice observations.

## Evidence from actual runs

September single-pair pilot: 11,498 samples, all assigned a native OCI CLD match after correcting the offset sign and allowing neighboring native granules. Above-cloud water finite for 10,609. Of 10,762 with finite CTP, 153 lacked computed water (1.42%); that combines profile-screening causes and is only an upper bound on surface-repair exclusions in that pilot.

With globally applied 3 m/s ocean LUT, valid counts were OCI2260 7,765; OCI2130 8,755; HARP2260 7,567; HARP2130 8,256; common 7,356. Indexing/dashboard callbacks and PNG/NetCDF export were exercised. This is plumbing/consistency evidence, not scientific calibration or a global validation.

The Arctic notebook was configured and regional display tested using existing data, but the actual Arctic satellite retrieval/matchup run has not been completed locally. Never imply otherwise from the notebook filename or a successful syntax check.

## Airborne evidence and validation labels

Actual June 10 file: `ARCSIX-HSRL-CloudAndSurface_G3_20240610_R1_L1.h5`, in `ARCSIX_AircraftRemoteSensing_LaRC-G3_HALO_Data`. Local inspection confirmed 31,683 profiles and an 801-altitude cloud mask, 15 m bins. Time is UTC decimal hours. `Nav/gps_lat/lon` describes the aircraft, not a separately ray-traced cloud intersection. Cloud-top and mask altitudes are metres. Readme: no uncertainty fields in this release.

Cloud phase code 1 is “Ice/Water Mix OR Ice/HOI Mix.” Codes 0/2 are water/ice dominant, 3 HOI,4 aerosol; NaN is unclassified. Lidar record ends 15:24:12.572 UTC, before the ~15:42 satellite granule. A same-day overlay is not proof of close temporal collocation. Empty matches are a legitimate result; do not silently enlarge tolerances.

The property-retrieval `Cloud Top Height` field may target lower water under ice; the validation reader instead uses the highest classified cloud-mask bin below the aircraft and above the surface, retaining the property height separately. Top-window defaults are 90 m, >=3 cloud bins, >=75% agreement, no mixed/HOI bins; these are analysis choices. Per-satellite summary requires >=3 profiles and >=80% consistent labels. Preserve ambiguous/unknown fractions and sampling/temporal offsets. Top water can support LI validation but does not establish all-liquid cloud interiors. Ice bins below the top may belong to a separate layer; they are evidence, not automatic LTMP truth. No advection correction is implemented yet.

P-3 FCDP R1 header verified: concentration #/L, LWC g/m³; 3 seconds already added to raw timestamps. Do not reapply. Navigation is separate. Optional ICARTT reader supports FFI1001 scaling/missing/detection-limit flags. No automatic imaging-probe phase retrieval, navigation-field guessing, or 3D in situ-to-cloud association has been implemented. Assess altitude and probe quality before assigning labels.

Threshold workflow: first inspect along-track/context plots, then use independently reviewed evidence labels. Current optional ratio fitting is liquid vs LTMP among water-topped/high-LI eligible samples. It does not optimize the full three-category phase map. Group train/test by whole cloud/flight segment, preferably separate flight days. The skill must not convert exploratory fitted cutoffs into campaign-wide defaults.

## Pending science and implementation priorities

1. Run the Arctic satellite/reference cache and quantify actual aircraft overlap at each time/distance window. If inadequate, choose a genuinely better collocated overpass/day rather than broadening until a result appears.
2. Complete P-3 navigation field mapping, image-based ice evidence, visibility/altitude review, and segment labels. Refine time/advection and cloud-height parallax/footprint handling as warranted.
3. Port and test operational surface-profile repair if population bias warrants it. User accepted the present exclusion for the exploratory test; this is not evidence that it is negligible globally.
4. Evaluate 3D radiative effects on normalized 2260/1615: unequal absorption/path-length distributions, illumination/shadowing, horizontal photon transport, unresolved heterogeneity, vertical weighting, ice/snow surface contribution, and spatial alignment. Their sign is geometry/state dependent. The OCI2260-derived CER/COT also responds to the same numerator band, so normalization can couple to retrieval biases. Compare OCI2130 and HARP2-CER on common samples; they are sensitivity tests, not fully independent measurements of all terms.
5. User requested discussion of these 3D effects after finishing validation cells. `Documentation/notes/3d_radiative_effects.md` now provides a source-supported starting point, fractional-ratio sensitivity equations, and actual SSA values from the supplied liquid LUT. Continue that discussion; no 3D correction has been implemented.
6. Assess surface-model and temperature assumptions only as requested; current accepted baseline is fixed 265 K and a globally applied 3 m/s ocean model.

User communication: provide concise scientific explanations, actionable notebook cells, honest scope of verification, and exact output links. Avoid repeating lengthy preambles or asking again about choices already settled. Keep this reference updated as actual results replace current hypotheses.

## Display update — 2026-09-23
`dashboard` now accepts `layout="slide"` (16:9), `map_framing="granule"` (regional occupied-cell framing), and `graticules=True`. The Arctic notebook opts into these; defaults remain classic/extent for existing callers. Revert with `layout="classic", map_framing="extent"`. This is display-only; scientific caches need no rebuilding. See `Documentation/dashboard.md`. No PowerPoint or interactive slide deck was created.

## Presentation workflow — 2026-09-24

`Validation/notebooks/presentation_arcsix_single_granule.ipynb` is the presentation
companion to the unchanged Arctic analysis notebook. It verifies the completed
single-granule file provenance chain and reads prepared references; no scientific
cache resume or LUT rebuilding occurs. New modules `validation/vis_truecolor.py`
and `validation/vis_presentation.py` extract/cache visible L1C TOA RGB and export
aligned threshold states. Default thresholds 1.10/1.27/1.45 are illustrative.
All full frames and map crops share the dashboard's axes and fixed framing.
RGB is a full-scene, gamma-enhanced 645/555/469 nm approximate true-color view at
5 km L1C resolution, not atmospheric correction or independent phase evidence.
Actual wavelengths, unknown QC, view selection, geometry, and display parameters
are recorded. First RGB read needs exact-granule Earthdata/S3 access (or a local
original L1C); subsequent reads use the extracted visible cache. See
`Documentation/dashboard.md` for settings, caveats, and slide assembly.
