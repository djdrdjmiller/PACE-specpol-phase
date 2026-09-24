# Presentation workflow verification — 2026-09-24

Scope: new Arctic presentation notebook, visible OCI L1C reader/cache, and aligned
figure exports. The analysis notebook is byte-for-byte unchanged from commit
`84128ce`; the scientific matching, correction, normalization, and classification
modules are unchanged. The dashboard's default interactive behavior is retained;
a new `show=False` option supports batch rendering using the same Figure.

- All 30 repository tests passed locally, including four new tests covering
  visible radiance-to-TOA conversion; same-view RGB selection; missing-band and
  strict-QC rejection; RGB cache reuse/settings changes; completed scientific
  cache file provenance; finite-distance RGB gap masking; threshold-dependent
  exports; constant LI map; and exact full-frame/map-crop agreement.
- The pre-existing input-preflight test now skips the CryoCloud `%pip` setup cell
  and accepts the notebook's explicit config argument. It still verifies missing
  LUTs fail before Earthdata login and no data/output directories are created.
- New notebook validated against nbformat; every code cell compiled through the
  IPython transformer. Stored outputs and execution counts are empty.
- Every new notebook code cell except package installation executed locally with
  synthetic single-granule caches and a synthetic visible-band L1C file using the
  production variable schema. The data-access call was redirected to that local
  fixture; no Earthdata authentication or S3 access was exercised.
- Eight PNGs were generated: three threshold dashboards, three corresponding
  phase-map crops, a true-color context dashboard, and its RGB crop. Full images
  measured 3200 × 1800 pixels; map crops measured 1354 × 531. Geographic framing,
  titles, legends, and colorbars were visually inspected using local coastlines.
- New helper/test import and undefined-name lint checks and Git whitespace checks
  passed. No LUTs, source datasets, synthetic fixtures, or figure binaries are
  included in the commit.

The actual Arctic data reside in the user's CryoCloud environment. Its first
visible-band read and resulting true-color scene remain to be run and inspected
there. Synthetic rendering checks establish software behavior, not Arctic cloud
phase performance, spectral calibration, or independent validation. Illustrative
thresholds 1.10/1.27/1.45 do not replace a scientific calibration.
