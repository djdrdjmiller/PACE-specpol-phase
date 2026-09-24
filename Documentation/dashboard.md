# Widescreen analysis dashboard

Install the package in your notebook environment and restart the kernel after code updates. Display changes do not require extraction or LUT calculations to be repeated. Import the viewer with:

```python
from pace_specpol.validation.vis_dashboard import reference_dashboard
```

Add these keyword arguments to your existing `reference_dashboard(...)` call:

```python
layout="slide",
map_framing="granule",
graticules=True,
```

Keep the existing Arctic `map_extent=(-85,-10,58,87)` and all scientific settings. The updated Arctic notebook already contains the options. Use `layout="classic", map_framing="extent"` to restore the old panel layout and framing. The classic version still supports `graticules=False` if you prefer no graticules.

## Display changes

- Slide preset: 13⅓ × 7.5 inches, exactly 16:9. PNG export is 3200 × 1800 pixels at 240 dpi without tight bounding-box cropping.
- Four equally sized panel slots, with shorter titles and reduced gaps. The geographic projection is never stretched.
- Colorbars have dedicated axes. The lower bars, including their extensions, fit the actual projected map width rather than the outer subplot slot.
- Lower colorbars have a vertical midpoint marker and a label above: LI threshold or ratio threshold. The ratio label, ticks, and normalization update with the slider. LI remains fixed by the index configuration.
- Sparse light graticules on all maps; small geographic labels sit inside the map edges in slide mode. Graticules have lower drawing order than the data. Coastlines remain above the data for geographic context.
- Granule framing uses occupied cells within the requested regional extent, adds vertical padding, and preserves or expands the horizontal span to match panel aspect. It does not throw away data, change counts, or restrict the histogram. It is intended for regional single-granule use and requires `map_extent`.
- Export names include layout and framing, so classic and slide images do not overwrite one another. Scientific NetCDF data are identical across display presets at the same thresholds.

The current phase slider changes phase counts and ratio color scaling; it does not change the mean LI or mean ratio values. Reframing is fixed for a given dashboard instance, so maps do not jump when adjusting thresholds.

## Verification

Rendered and inspected both initial and changed-threshold slide layouts using a clearly marked synthetic Arctic placement of the September pilot values. The actual Arctic dataset is still in the user's AWS environment. Checked slide dimensions, updating colorbar threshold, and existing scientific tests. No synthetic map is provided as a scientific result.

## PowerPoint options (not implemented here)

Microsoft defines widescreen as 13.333 × 7.5 inches. A PNG fits that slide ratio directly; leave room within the slide if adding a separate PowerPoint title. [Microsoft slide-size documentation](https://support.microsoft.com/en-us/powerpoint/change-the-size-of-your-powerpoint-slides).

For a few preset thresholds, export one figure per value and place each on an identically arranged slide. Small linked buttons can jump among those slides, creating a selectable comparison without running Python or requiring network access. [Microsoft action-button documentation](https://support.microsoft.com/en-us/powerpoint/add-commands-to-your-presentation-with-action-buttons).

A pasted figure cannot retain its ipywidgets connection. Hosting a separate web dashboard or using a compatible web add-in is a different deployment path with platform/network dependencies; it is not required for a preset slide comparison. No PowerPoint file or interactive presentation was created in this change.

## Preview/spacing refinement — 2026-09-23
The slide preview now renders the full figure canvas into a responsive PNG inside the notebook, bypassing inline tight-bounding-box cropping. It scales down to fit the available notebook width without trimming the right edge. The histogram colorbar has more separation from both its x-axis label and the lower-map title. Explicit map-title positions avoid automatic title shifts from geographic-label layout. The saved PNG remains 3200 × 1800; use the dashboard Save button for the presentation-resolution image.

## Arctic presentation notebook — 2026-09-24

Run `Validation/notebooks/presentation_arcsix_single_granule.ipynb` after completing
its analysis counterpart. It starts from the existing single-granule caches and
checks the raw → companion → reference file-hash chain before indexing. It never
rewrites scientific cache manifests or calls the LUT correction stages. The
analysis notebook is preserved, including its current CryoCloud setup.

Select the existing Arctic local TOML (`config/arcsix.local.toml`, or explicitly
`PACE_CONFIG`). Repository discovery works from inside the clone; a kernel started
elsewhere needs the notebook's `REPO_OVERRIDE`. No local path file is committed or
generated. After pulling new Python code, restart the kernel. After a full
CryoCloud image reset, run the notebook's `%pip` install cell again.

Default exports use the normalized OCI 2260 / own population at ratio thresholds
1.10, 1.27, and 1.45, with LI=0.3. These are illustrative, not calibrated values.
The mean fields and histogram population remain fixed, while phase counts and
ratio color centering change. Settings are collected in one notebook cell.

`vis_presentation.export_presentation` reuses the existing dashboard Figure at
240 dpi (3200 × 1800). All frames use one projection, extent, and pixel crop.
It produces three dashboard PNGs, three map-only phase PNGs, a true-color context
dashboard and map-only RGB PNG, and `presentation_manifest.json` containing
settings, provenance, dimensions, crop coordinates, and output hashes. An
existing nonempty export folder requires a new `EXPORT_LABEL` or an explicit
`OVERWRITE=True`. The normal interactive save button remains separate.

### True-color meaning and data access

`vis_truecolor` reads the exact paired OCI L1C V3 product, selecting the closest
channels to 645/555/469 nm and recording their actual wavelengths. Radiance is
converted using the same TOA convention as the SWIR extractor:
`rho = pi * I * distance² / (F0 * cos(sza))`. The source geometry and irradiance
fields follow the [NASA L1C format specification](https://oceancolor.gsfc.nasa.gov/files/NASA_TM2024219027v12_PACE_Level_1C_Format.pdf).
This is an approximate gamma-enhanced visible composite, without atmospheric
correction, per-channel auto balancing, or a colorimetric transform.

All three channels must be valid in one view; the smallest absolute nadir time
offset wins. The paired cache supplies geometry limits and OCI QC policy. If
that policy ignores QC, the notebook reports unknown-QC counts rather than
claiming those pixels passed a quality flag. HARP2 cloud/LI and liquid-reference
masks are not applied, so true color includes additional scene context. RGB view
selection can differ from SWIR where spectral availability differs.

The first cloud run authenticates with Earthdata, verifies the exact source name
and S3 URL against CMR, and range-reads the three bands. Contiguous spectral
storage can still incur substantial underlying I/O. It saves a local visible
reflectance cache; later gamma/threshold edits do not require more satellite
reads. Outside AWS, pass an already downloaded exact L1C file through
`LOCAL_OCI_SOURCE`. The cached RGB and exported figures stay under configured
work/export directories, outside tracked source files.

The default common-channel display stretch is `clip(rho, 0, 1) ** (1/2.2)`.
Nearest-center sampling onto the dashboard projection has a 4 km distance cutoff;
it neither fills distant gaps nor increases the native 5 km resolution. White
RGB pixels indicate no displayed RGB coverage; white phase pixels can instead
reflect scientific exclusions or tied phase counts. Both maps share geolocation
and framing, not a common validity mask or independent phase evidence.

### Slide assembly

Use identically sized and positioned images on consecutive 16:9 slides with a
short Fade transition. For a clickable demonstration, place three consistently
positioned threshold buttons on each slide and link to the corresponding slides.
The map-only RGB and phase crops have identical pixel sizes and can be swapped
without manual geographic alignment. Full context and phase dashboards can also
be swapped; the histogram/lower panels on the context frame use the explicitly
recorded `CONTEXT_THRESHOLD`. No PowerPoint file is generated by this notebook.
