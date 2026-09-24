# Reorganization verification

This record covers the structural refactor into the `pace_specpol` package. It
does not establish scientific phase discrimination or airborne matchup skill.

## Preservation

The original flat package was retained untouched outside this repository and
archived before migration. All 36 regular source/record/preview/bytecode files
in that archive matched the retained originals byte-for-byte after verification
(macOS archive metadata entries were excluded from that comparison).

## Automated tests and package checks

- Original package: **22 tests passed** in the verification environment.
- Reorganized package: **24 tests passed**, retaining the original scientific
  assertions and adding configuration-location and missing-config checks.
- Editable package installation succeeded in an isolated verification environment.
- All three analysis notebooks passed schema validation and code-cell syntax
  checks. Their initial setup cells and imports were executed without remote
  discovery or extraction. Published copies have empty outputs/execution counts.
- Undefined-name checks passed for source modules and tests.
- Python source and tests were formatted; unused imports were removed.

## Real local September sample

Recomputed both old and new implementations on the same **11,498** companion
samples using the supplied real 265 K HDF4 LUTs. Full returned datasets, including
correction/reference fields and status flags, were identical for all variants:

| Reference | Valid samples | Old/new comparison |
|---|---:|---|
| OCI2260 | 7,765 | Identical |
| OCI2130 | 8,755 | Identical |
| HARP2260 | 7,567 | Identical |
| HARP2130 | 8,256 | Identical |

For OCI2260, rebuilt old/new indices in separate verification directories for
raw, corrected, and normalized ratios, each on own and common populations.
All metadata/fingerprints, histograms, and sums matched. Counts and exported
dataset contents matched at thresholds 0.5, 1.10, 1.27, 2.0, and 3.5, excluding
the intentionally variable export creation timestamp. Own population contained
7,765 samples and common population 7,356. An existing schema-2 index also loaded
and produced identical phase counts. Existing source caches were not rewritten.

## Visualization

Classic and slide dashboard canvases were pixel-identical between the preserved
and reorganized versions at thresholds 1.27 and 1.10 in the same environment.
Raw/corrected/normalized selector callbacks retained the expected 7,765 samples.
PNG and NetCDF export callbacks produced their expected files. The slide preview
was visually inspected. This used the actual September pilot's global map,
not a new Arctic retrieval. No separate regional framing or airborne field
campaign was executed as part of this migration.

## Repository review and limits

The local repository has no remote, commits, or pushes. Candidate publication
files were reviewed for scientific binaries, large files, absolute developer
paths, credential patterns, and signed URLs. None were found by those checks.
Ignore rules were exercised against representative LUT, data, results, local
configuration, bytecode, and credential filenames. This is a targeted check,
not a claim of an exhaustive security audit.

The verification environment used Python 3.11 from a task-local environment
with the existing PACE scientific packages and locally installed pytest,
pyhdf, and widgets. Its Cartopy/PROJ test setup needed an explicit PROJ data
directory because it reused a Conda environment through a virtual environment;
no machine-specific path was added to the package.

Full notebook execution with fresh satellite downloads and the user's AWS
Arctic caches was not repeated. No scientific defaults, LUT conventions,
cache versions, or threshold rules were changed. License selection and any
third-party redistribution review remain separate from this refactor.

## Follow-up: Arctic input check

Moved local LUT configuration and file validation before Earthdata login,
discovery, and extraction in the Arctic notebook. Missing files now produce an
actionable error identifying the path and `lut_root` configuration. Two regression
cases execute notebook cells with network calls intercepted: missing tables stop
before login or output creation; present files reach login without loading LUTs.
The suite now has **26 passing tests**. Notebook schema and code-cell syntax
also pass. This changes startup ordering and error guidance, not the science.
