# External data and lookup tables

No mission observations, ancillary downloads, or LUT binaries are distributed
with this repository. Configure their locations in an ignored local TOML file.
The loader expands `~`; relative values resolve against that configuration file.

The Arctic notebook checks that all three required LUT files exist before
Earthdata login, discovery, or extraction. Missing inputs stop the run with the
missing path and configuration guidance. This is a file-presence check, not a
full integrity check; the actual HDF contents are read during reference processing.

Required LUT files beneath `lut_root` are:

- `LIQUID/ocean_msr_water_wspeed_3_v6.PACE.1.1.5.2026144071240.hdf`
- `IceAndWaterPhaseFunctionData_v6.PACE.1.1.5.2026142144440.hdf`
- `Transmittance_OCI.hdf`

These are the supplied 265 K tables, not generic replacements available through
this package. Obtain the authorized input set from its provider and verify file
hashes against the applicable experiment's provenance. The water-combination
`.dat` file is refractive-index provenance, not an additional runtime LUT.

Satellite discovery uses HARP2 GPC V4.0, OCI L1C V3, and native OCI CLD V3.1.
Ancillary PROFILE/MET filenames come from each native cloud product. Preserve
their original names and hashes. See [algorithm details](algorithm_details.md)
for the matching, units, and screening conventions.

`.gitignore` excludes research data formats, generated cache/output directories,
and local settings. Small metadata records under `Validation/records/` describe
inputs and verification; they are not substitutes for the observations. Do not
commit Earthdata credentials or signed download links. Review the staged file
list before publishing; ignore rules do not untrack already committed files.

The eventual software license will not change the rights or attribution terms
of external observations, tables, or adapted third-party code.
