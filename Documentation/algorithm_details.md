# 265 K liquid-reference extension

Start with `Validation/notebooks/analysis_liquid_reference.ipynb` after installing the package and configuring paths as described in the root README.
The original HARP2/OCI paired cache is the input. Its reflectances are never overwritten.
This is an exploratory reference calculation, not a reproduction of the complete OCI L2 retrieval.

## Choices

Default: OCI `cer_22` with its matching `cot_22`. The comparison cache prepares:

- `oci_2260`: OCI 2260 nm CER and COT.
- `oci_2130`: OCI 2130 nm CER and COT.
- `harp2_2260`: HARP2 cloudbow CER with OCI `cot_22`.
- `harp2_2130`: HARP2 cloudbow CER with OCI `cot_21`.

All modes use the supplied 265 K optical-property tables. The water_comb `.dat` file is the bulk complex refractive index used to generate those tables; it is not a radius grid and is not an additional runtime input. No temperature interpolation or Mie calculation is performed.

Each mode can show its own valid population or the intersection valid in ALL configured variants. Remove variants from the dictionary before creating a new cache to change that comparison population. The dashboard offers raw TOA, atmospherically corrected, and liquid-normalized ratios on the selected population. Default ratio threshold remains 1.27 and remains provisional. Continuous density uses Turbo; LI and ratio maps retain RdBu_r/TwoSlopeNorm.

OCI phase codes 2 (water), 3 (ice), and 4 (unknown, assumed water) are accepted by default. **Using ice-retrieval CER in a liquid LUT is a same-numerical-radius counterfactual, not a retrieved liquid radius.** `ice_cer_as_liquid` records this. CER outside the supplied liquid grid is rejected, never clamped. Set `oci_phase_policy='liquid'` to retain code 2 only. Missing HARP2 CER is not replaced by OCI CER.

## What is implemented

### Reflectance

`LiquidLUT` reads one of the supplied `ocean_msr_water_wspeed_*` files, plus the matching phase-function file. Wind speed is fixed by the chosen file. The notebook uses 3 m/s globally (`ocean_only=False`), including land and snow/ice, as an exploratory surface-model assumption requested for this comparison. It does not represent their actual surface reflectance. A Lambertian surface-albedo correction has NOT been implemented. Set `ocean_only=True` to restrict to OCI pixels using the ocean nonabsorbing band (code 2).

The multiple-scattering component is multilinearly interpolated in native sensor cosine, solar cosine, relative azimuth, COT and CER. Single scattering is constructed at bracketing COT/CER nodes and interpolated after reconstruction. The phase function is linearly interpolated at the actual scattering angle, divided by `PhaseFuncNormConstant`, and used with extinction ratios relative to band 1, single-scattering albedo, and truncation factor. The exp(-x)=0 shortcut for x>10 is preserved. For these three OCI SWIR bands, the NASA ocean SS routine omits its additional below-cloud Rayleigh/aerosol SS terms; this implementation follows that branch.

The cache's folded solar/sensor azimuth difference is converted to the LUT convention with `180 - raa`. Both wavelengths are simulated separately before their ratio is formed. No generic interpolation of a precomputed ratio is used.

### Atmospheric correction and provenance

Actual OCI CLD V3.1 metadata inspected for September 3, 2024 specifies:

- `GMAO_MERRA2.20240903T000000.PROFILE.nc`
- `GMAO_MERRA2.20240903T030000.PROFILE.nc`
- `GMAO_MERRA2.20240903T000000.MET.nc`
- `GMAO_MERRA2.20240903T010000.MET.nc`

The code reads `anc_profile1/2/3` and `met1/2/3` from EACH matched CLD granule. It retrieves those exact OB.DAAC filenames, not substitute GEOS forecasts or independently selected MERRA collections. PROFILE QV is kg/kg; MET PS is Pa. File hashes, native cloud product, software version, and relevant input parameters are saved with every companion.

QV and pressure are interpolated bilinearly in latitude/longitude (periodic longitude), then linearly in time. QV is converted to mixing ratio in g/kg; humidity is interpolated linearly in log pressure to the CHIMAERA 101-level grid. Above-cloud water is integrated with its source trapezoidal rule and 980.616 factor, producing g/cm². The upper-profile water floor is 0.003 g/kg.

**Profiles whose cloud-top integration needs surface-adjusted levels are rejected.** The operational source includes additional surface-level repair/extrapolation logic; we do not silently approximate that logic. Similarly, missing spatial/profile values are rejected rather than replaced with neighboring values. Therefore the exact input datasets are reproduced, but the full operational ancillary repair behavior and acceptance population are not. Rejected above-cloud-water inputs appear in the transmission exclusion count.

Transmittance uses `mu_eff = mu0*mu/(mu0+mu)`, nearest 100-hPa pressure level (ties upward), and bilinear interpolation in water amount and effective cosine. It is already the two-way transmission; no squaring or separate-path multiplication is applied. The table's fixed band indices are 5 (1615), 6 (2130), 7 (2260), zero-based; legacy wavelength labels are recorded by the reader, not treated as interpolatable spectral coordinates. Missing/unpopulated corners, pressure outside 100–1000 hPa and out-of-range water are rejected. This is stricter than NASA's fallback behavior, which can walk through water bins or extrapolate.

The final calculation is:

`normalized = raw_2260_over_1615 * T1615 / T2260 / (liquid_rho2260 / liquid_rho1615)`

The source observation reflectances remain TOA. Corrected and normalized values are separate fields.

### Spatial and time matching

OCI CLD is native approximately 1 km, whereas the original paired cache is L1C approximately 5 km. The new code transfers the nearest native cloud retrieval compatible with the L1C observation time, searching the eight closest native pixels per overlapping CLD granule. Default limits: 3 km and 10 seconds. It computes ancillary profiles at that native pixel and scan time. **This is a center-point approximation, not an average of the native microphysical retrievals or reconstructed reflectances over the full L1C footprint.** Geometry for the liquid simulation remains the cached L1C geometry. No new cloud-parallax correction is applied. Distance, time difference, native row/column and source-granule index are retained for auditing.

NASA's L1C writer sets `view_time_offsets = nadir_time - mean_native_scan_time`. The older paired script used the opposite sign when constructing observation time. The new companion reader repairs legacy cached timestamps (stored time minus twice the offset) and searches adjacent CLD granules. The paired extractor is updated to version 1.1 and the correct sign for future extractions. Existing completed caches remain usable here; do not rerun `cache_pairs` with 1.1 into a 1.0 directory. The prior date-boundary mask may have excluded some valid samples near the start/end of the requested dates; these cannot be recovered from an existing cache. A future boundary-only or fresh extraction would be needed to restore those samples.

## Running, restarting and resource use

1. Install requirements in the existing notebook environment.
2. Set the existing paired-cache path and your LUT directory.
3. Run the default three-pair pilot and inspect coverage/status summaries.
4. Set `MAX_PAIRS=None` and use the full-run directories to process the complete manifest.
5. Run the reference stage once; use the dashboard selector afterward.

A failed run can be restarted by rerunning the same cell with the same configuration. Each completed file is saved atomically; source hashes and configuration fingerprints are checked. Changed scientific settings require a new output directory. Companion downloads and corrected reference stages are separate, so changing CER mode does not repeat satellite/meteorological reads.

Only one MS band is resident (~88 MiB float32 for the supplied table). Native OCI arrays and interpolated ancillary fields are processed one granule at a time. No three-day collection is loaded into RAM. The managed ancillary download cache defaults to 512 MiB, with 512 MiB additional free-disk headroom; it evicts only files listed in its own managed-files.json. Preexisting unmanaged files are never deleted and are not counted against that managed-cache cap. Companion/reference outputs and sparse dashboard indexes are separate permanent disk consumers. Inspect disk usage after the pilot before a full run. Downloaded profiles in the live test were ~53 MiB each and MET files ~3 MiB each.

The dashboard creates one sparse index per variant/population/ratio-space selection. First use requires a local pass through correction files. Subsequent threshold changes need no remote data. Export filenames should use a separate output directory for each experiment to avoid replacing earlier figures.

## Source references and validation boundary

The implementation was checked against NASA OCSSW source commit `5d057d85310e594a2efc4ccf5ff01dfd43954c14`. The inspected OCI CLD product was generated by `l2gen 9.11.0-3f35fce06`; no bit-for-bit equivalence with that release is claimed.

- [MOD06 C6 LUT overview](https://atmosphere-imager.gsfc.nasa.gov/sites/default/files/ModAtmo/C6_LUT_document_final.pdf)
- [interpolate_libraries.f90](https://git.smce.nasa.gov/oel/ocssw/-/blob/5d057d85310e594a2efc4ccf5ff01dfd43954c14/src/libcloud/interpolate_libraries.f90): phase normalization, geometry, tau scaling, ocean SS.
- [atmospheric_correction.f90](https://git.smce.nasa.gov/oel/ocssw/-/blob/5d057d85310e594a2efc4ccf5ff01dfd43954c14/src/libcloud/atmospheric_correction.f90): effective two-way cosine, pressure selection, PW/mu interpolation.
- [specific_ancillary.f90](https://git.smce.nasa.gov/oel/ocssw/-/blob/5d057d85310e594a2efc4ccf5ff01dfd43954c14/src/libcloud/specific_ancillary.f90): explicit OCI spectral-band remapping.
- [ancillary_module.f90](https://git.smce.nasa.gov/oel/ocssw/-/blob/5d057d85310e594a2efc4ccf5ff01dfd43954c14/src/libcloud/ancillary_module.f90): above-cloud water integration.
- [profile_management.c](https://git.smce.nasa.gov/oel/ocssw/-/blob/5d057d85310e594a2efc4ccf5ff01dfd43954c14/src/libcloud/profile_management.c): 101-level grid and log-pressure interpolation.
- [anc_acq.c](https://git.smce.nasa.gov/oel/ocssw/-/blob/5d057d85310e594a2efc4ccf5ff01dfd43954c14/src/l2gen/anc_acq.c), [get_cmp.c](https://git.smce.nasa.gov/oel/ocssw/-/blob/5d057d85310e594a2efc4ccf5ff01dfd43954c14/src/l2gen/get_cmp.c): original profile handling and mixing-ratio conversion.
- [l1c_latlongrid.cpp](https://git.smce.nasa.gov/oel/ocssw/-/blob/5d057d85310e594a2efc4ccf5ff01dfd43954c14/src/l1cgen/l1c_latlongrid.cpp): view-offset sign.

See `Validation/records/LIQUID_VALIDATION.json` for the actual live-test results. No complete three-day scientific validation or threshold calibration has been performed.

Pilot profile screening: of 11,498 matched samples, 736 lacked finite cloud-top pressure and another 153 lacked computed above-cloud water despite finite cloud-top pressure (1.42% of those with finite pressure). The latter combines all profile-screening causes; it is an upper bound on near-surface exclusions, not an isolated count or global estimate. Such exclusions can selectively remove low clouds.
