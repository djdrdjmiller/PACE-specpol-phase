# Interactive PACE cloud-top phase exploration

Open **pace_mixed_phase.ipynb** in Jupyter. Keep `pace_mixed_phase.py` and
`requirements.txt` beside it. The original LI-only notebook is not required.

1. Install the notebook's dependencies and restart its kernel if necessary.
2. Set dates and extraction settings. Default: September 3–5, 2024 inclusive UTC.
3. Run discovery/extraction in AWS us-west-2 using Earthdata Login.
4. Choose raw exploration, your scalar simulated reference, or your reference LUT.
5. Build the local threshold index, then use the slider to update the phase map.
6. Click **Save current PNG + NetCDF** to export the two-panel figure and native-grid data.

The notebook pairs **HARP2 GPC V4.0 non-NRT** with **OCI L1C V3.0**, checking actual
coordinates and nadir times. It converts OCI radiance to TOA reflectance and uses
nominal 2260/1615 nm bands from the same OCI view. These are not RSP's identical
bands. Matching is on the shared surface grid; cloud parallax is not corrected.

No simulated all-liquid baseline is invented. The default UI states **RAW ratio;
unnormalized**. Supply a reference as documented in REFERENCE_INPUT.md before
interpreting normalized-ratio phase candidates. The initial 1.27 ratio threshold
comes from the user's ACTIVATE example and is provisional for PACE.

Phase regions follow the supplied diagram: ratio below cutoff is liquid; above
or equal to cutoff is ice if LI <0.3 and LTMP if LI >=0.3. This includes the
low-LI/low-ratio quadrant in liquid by the spectral rule. Observations are
classified before aggregation. The three-color L3-like map shows the most
frequent class in each grid cell. Ties, insufficient counts and no data are white.
Phase fractions and counts are exported as well; they are sample-population
statistics, not cloud-volume or water fractions.

## Reuse and performance

The extraction cache contains compact matched samples. Reference changes and LI
cutoff changes rebuild the sparse index from those local files. Ratio-slider
changes use the index alone. Reopen an existing index with `mp.PhaseIndex.load`.
Slider updates occur on release to avoid queuing expensive renders.

Default grid spacing is 0.1°. The on-screen Robinson preview samples native cells
without blending class colors. NetCDF export preserves the full native grid.
Use 0.5° or 1° if index memory is constrained. At 0.1°, six global int64 arrays
alone occupy about 311 MB; export and sparse-index overhead require additional
memory. Monthly caches can be substantial; start with a shorter interval.

OCI L1C files can approach 1 GB. Their contiguous spectral storage means selected
band access can still read much of a file. Run extraction in-region; do not treat
it as a lightweight metadata operation. Two workers/four files per batch are the
initial defaults. Source files are not fully retained locally; matched samples
are. Earthdata sessions are refreshed between batches after 40 minutes.

If a cache's extraction settings change, choose a new cache directory. A failed
extraction leaves its manifest incomplete and resumes completed pair files on
rerun. Missing partners and grid mismatches fail explicitly. Do not use an
incomplete manifest as a complete-period result.

## Known limitations

The real tested OCI L1C `qc` field is entirely fill. Default ignore policy means
unknown spectral QC, not good QC. All screening choices and QC availability are
recorded. There is no implemented independent saturation test or atmospheric-gas
correction. No full-period validated LTMP retrieval is claimed.

See VALIDATION.md for tests and VERIFIED_CATALOG.json for the discovery audit.

## Midnight time-reader fix

The reader includes a guarded correction for observed OCI V3 CF time epochs
that differ by one day from the global coverage date. It requires independent
agreement with global coverage and HARP2 row times, retains the one-second
tolerance, and records the adjustment. Previously completed caches remain
compatible. Replace the original module, restart the kernel, then rerun setup
and extraction with the same configuration and cache directory.

## Dashboard update (September 22, 2026)
Three maps now show dominant phase, arithmetic mean liquid index, and arithmetic mean per-sample ratio. Continuous plots use Turbo. Export retains both mean fields. Same eligible matched population for all fields; no atmospheric correction is added.

Replace pace_mixed_phase.py and reload/restart, then rerun build_index on the existing completed pair cache. Schema 2 creates a new index automatically (omit any old explicit index_dir). No NASA extraction is repeated. Two float64 sum arrays add about 104 MB at 0.1 degree resolution, excluding temporary/display allocations.

Use index.sampling_summary() for exact accepted sample counts. A 0.1-degree cell covers approximately 123.6*cos(latitude) square km, or 4.57*cos(latitude) native 5.2 km cell areas per full coverage. This is an area estimate, not an actual observation count. Repeat passes add samples; screening removes them.

plot_cached_samples provides optional per-sample center plots without geographic averaging; overlaps overplot in manifest order and marker size is not footprint size. The fast dashboard retains nearest-cell display sampling; exported grids retain every cell.

### Diverging map color update
Dashboard mean LI now uses RdBu_r/TwoSlopeNorm centered at the LI threshold (0.3), with li_clim=(-0.5,4). Mean ratio uses RdBu_r centered at the moving ratio cutoff and ratio_clim=(0.5,2). Requested limits expand only as needed to bracket the cutoff. Under/over colors retain outliers; the histogram remains Turbo and its x bounds are unchanged. Native-cell ratio summaries quantify the upper tail. This display-only update needs no index rebuild; reload the module and recreate the dashboard using the existing index.

## 265 K liquid reference (new)

Use **pace_liquid_reference.ipynb** for the all-liquid normalization workflow.
Read **LIQUID_REFERENCE.md** for source conventions, matching, resource use and limitations.
It reuses completed paired caches and adds original OCI CLD/MERRA-2 ancillary inputs,
OCI 2260/2130 CER-COT modes, HARP2 CER alternatives, a common-population comparison,
and raw/corrected/normalized dashboard selections. Requirements now include pyhdf.

The new companion reader repairs legacy L1C view-offset signs without redownloading
reflectances. The base extractor is now version 1.1 for future extractions; do not
resume a 1.0 extraction directory with 1.1. Completed old caches are accepted as
inputs to the new notebook. Precise date-boundary completeness may require a fresh
or boundary-only extraction; see LIQUID_REFERENCE.md.

## ARCSIX Arctic pilot
Open `pace_arcsix_single_granule.ipynb` and run in order on AWS. It selects only `20240610T154205` (June 10, 2024, 15:42:05–15:47:04 UTC), prepares fresh paired data, and runs the same four liquid-reference variants. The regional polar map is a display crop; the histogram retains the full granule. No flight track is loaded. This Arctic retrieval has not been run locally.

## Airborne validation and agent handoff
The Arctic notebook now includes the airborne-validation cells. To preserve your existing notebook, copy from `arcsix_validation_cells.ipynb` or `arcsix_validation_cells.py` and add `pace_arcsix_validation.py` beside it. Read `ARCSIX_VALIDATION.md` for verified June 10 timing, lidar-code ambiguity, optional in situ field mapping, and remaining validation work.

The portable agent skill is `skills/pace-oci-harp2-phase/SKILL.md`, with detailed scientific/project context in its `references/handoff.md`. Tell a future agent to read that entrypoint and locate this package. It is bundled documentation, not automatically installed in an agent environment and not model-weight training. `3D_RADIATIVE_EFFECTS.md` starts the requested discussion of wavelength-dependent 3D effects.
