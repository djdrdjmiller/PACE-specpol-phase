# Validation — September 21, 2026

## Catalog and real observations

NASA CMR metadata returned 437 HARP2 CLOUD_GPC V4.0 non-NRT granules and 437 OCI
L1C V3 granules for September 3–5, 2024. All 437 filename timestamps paired.
See VERIFIED_CATALOG.json. Runtime discovery checks again and fails on missing
partners unless the user explicitly enables partial coverage.

The tested real pair was:

- PACE_HARP2.20240903T001241.L2.CLOUD_GPC.V4_0.nc
- PACE_OCI.20240903T001241.L1C.V3.5km.nc

OCI was read by authenticated HTTPS range access; a two-band subset retaining
its geometry, timing and metadata was used for repeat tests. HARP2 was read
from a real downloaded source file. No synthetic observations were substituted
for these real-file checks.

The pair has exactly equal grid latitudes/longitudes. Selected band centers are
1618.034 and 2258.429 nm; solar irradiances are 235.0 and 73.97 W m^-2 um^-1.
The OCI file is about 951 MB and stores these spectral variables contiguously.
The initial subset read took about 41 seconds in the development environment;
this is not an AWS performance benchmark.

The real pair contained 62,132 finite LI samples, of which **11,498 matched
samples** passed the notebook's default screening. All matched samples had
unknown OCI spectral QC because both selected bands' QC fields were entirely
fill. Requiring zero QC rejected all of them, as expected. Selecting a later
UTC date rejected the whole granule.

The reflectance ratio agreed with the independent radiance ratio multiplied by
the ratio of band-specific solar irradiances. Tested conversion retains low
and negative LI values and never applies a liquid-only retrieval mask.

## Automated and rendering checks

Seven standalone offline regression tests pass (`test_phase_logic.py`):

- Exact threshold equality, low-LI/low-ratio region, NaN classification and ties.
- Sparse threshold counts vs direct observation-by-observation classification,
  including underflow/overflow, scalar normalization, off-grid threshold errors,
  and NetCDF/index round trips.
- LUT interpolation and exclusion outside its domain; invalid reference rejection.
- Synthetic multi-view fixture: bands existing only in different OCI views are
  rejected; valid duplicate views are counted once; date/quality/grid checks work.
- Dateline wrapping and geographic endpoint binning.

Additional local checks exercised widget slider and minimum-sample callbacks,
Cartopy rendering, the export-button callback, and an index built from the real
pair. The real-pair two-panel image was visually inspected. The preview is a
single-granule raw-ratio test, not a full-period phase retrieval. Notebook JSON,
cell syntax and schema were validated.

## Not tested or claimed

- Full-period (437-pair) extraction and performance.
- Authenticated direct S3 execution in the user's AWS environment; real-file tests
  used HTTPS range access instead.
- Browser-side widget-manager behavior in the user's Jupyter installation.
- A supplied radiative-transfer LUT or validated OCI-specific normalized threshold.
- Production scientific validation of LTMP identification, saturation QA,
  atmospheric correction, or cloud-top parallax correction.

The seven regression tests can be rerun with:

```bash
python -m pip install pytest
python -m pytest -q test_phase_logic.py
```

## Guarded OCI midnight-epoch repair

Verified that OCI granules 20240905T000311, 000811, 001311, and 001811
have global coverage on September 5 but CF time units referencing September 4
with seconds reset near zero. The preceding 20240904T235811 granule correctly
crosses midnight. The new reader repairs only a uniform +/-86400-second offset
when the corrected times also agree with OCI global coverage and HARP2 times
within the original one-second tolerance. It records the offset and source units
in new pair reports. Other mismatches still fail.

Real-file validation of 20240905T000311 passed: +86400-second correction, zero
residual timing difference, 205005 aligned grid cells. This was a metadata/grid
validation, not full SWIR extraction. Seven regression tests passed, including
rejection of non-day and coverage-inconsistent errors. Existing successful
pair caches remain compatible; their extraction fingerprint is unchanged.

## September 22 dashboard update
Eight offline tests passed, including exact continuous sample means (negative LI and out-of-histogram ratios), reference normalization, persistence and minimum-count export masking. Synthetic four-panel rendering, slider callbacks, sample-count masking, and unaggregated sample plotting passed with cached coastlines. Full three-day AWS data were not rerun.

Diverging color update: eight regression tests passed. Synthetic callback checks verified midpoint=0.5 at thresholds 0.5, 2.0, 3.5 and 1.27; LI midpoint=0.3; ratio arrays unchanged as cutoff moves; empty minimum-count mask handled; PNG rendered and inspected.
