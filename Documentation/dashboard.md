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
