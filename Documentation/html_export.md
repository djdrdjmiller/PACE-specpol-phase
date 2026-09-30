# HTML Dashboard Export

## Overview

The HTML export feature creates standalone, shareable dashboard files from ARCSIX phase analysis results. These HTML files embed all visualizations and include an interactive threshold slider, making them ideal for sharing with colleagues who don't have Python environments set up.

## Quick Start

```python
from pace_specpol.validation.vis_html import export_html_dashboard

# Export with default settings
export_html_dashboard(
    index,
    "dashboard.html",
    thresholds=(1.10, 1.27, 1.45),
    default_threshold=1.27,
)
```

See `Validation/notebooks/export_html_dashboard.ipynb` for a complete example.

## Features

### Interactive threshold slider
Users can adjust the ratio threshold using a slider control. The dashboard pre-renders all specified threshold values, allowing instant switching without recomputation.

### Embedded visualizations
All dashboard panels are embedded as base64-encoded PNG images:
- Joint histogram with phase regions
- Modal phase map
- Mean liquid index map
- Mean ratio map

### Metadata display
The HTML includes comprehensive metadata:
- Sample totals and statistics
- Cell statistics (mean, median, percentiles)
- Configuration parameters
- Reference mode and variant information

### Self-contained
The HTML file includes:
- All images (base64 encoded)
- All JavaScript (inline)
- All CSS (inline)
- No external dependencies

## Parameters

### Required
- `index`: Pre-built aggregation index from `build_index()`
- `output_path`: Path where HTML file will be saved

### Threshold configuration
- `thresholds`: Tuple of threshold values to pre-render (default: `(1.10, 1.27, 1.45)`)
- `default_threshold`: Initial threshold to display (default: `1.27`)
- `min_count`: Minimum sample count per cell (default: `1`)

### Visual settings
- `dpi`: Figure resolution (default: `150`)
  - 150 DPI: Good for web viewing, ~2-3 MB per threshold
  - 240 DPI: High resolution for presentations, ~5-7 MB per threshold
- `preview_width`: Map preview width in pixels (default: `900`)
- `ratio_clim`: Ratio colorbar limits (default: `(0.5, 2.0)`)
- `li_clim`: Liquid index colorbar limits (default: `(-0.5, 4.0)`)
- `map_extent`: Regional extent for Arctic maps (default: `None`)
  - Format: `(west, east, south, north)`
  - Example: `(-85, -10, 58, 87)` for Arctic region
- `graticules`: Show graticule lines (default: `True`)

### Customization
- `title`: Custom dashboard title (default: `"ARCSIX Phase Dashboard"`)
- `**dashboard_options`: Additional options passed to `dashboard()`

## File Size Considerations

Typical file sizes:
- 3 thresholds at 150 DPI: ~5-8 MB
- 5 thresholds at 150 DPI: ~8-12 MB
- 3 thresholds at 240 DPI: ~12-18 MB

File size scales approximately linearly with:
- Number of thresholds
- DPI setting
- Map preview width

## Sharing Options

### Email or messaging
Attach the HTML file directly. Recipients download and open in their browser.

### Cloud storage
Upload to:
- Google Drive (download to view; "Open with Google Docs" won't work)
- Dropbox (share a direct download link)
- OneDrive, Box, etc.

### Web server
Upload to an internal or public web server for permanent URL access. The file is completely self-contained.

### Version control
HTML files can be committed to Git repositories, but large file sizes may make this impractical for repositories with many dashboards. Consider using Git LFS or keeping HTML exports separate from source code.

## Browser Compatibility

The HTML dashboards work in all modern browsers:
- Chrome/Edge (Chromium) 90+
- Firefox 88+
- Safari 14+
- Any browser supporting HTML5, CSS3, and ES6 JavaScript

No plugins or extensions required.

## Workflow

Typical workflow for creating shareable dashboards:

1. **Run analysis notebook** (`analysis_arcsix_single_granule.ipynb`)
   - Extract satellite data
   - Generate reference comparisons
   - Build aggregation index

2. **Export HTML dashboard** (`export_html_dashboard.ipynb`)
   - Load cached index
   - Configure threshold values
   - Generate HTML file

3. **Share with colleagues**
   - Email, cloud storage, or web hosting
   - Recipients open in browser
   - No Python environment needed

## Advanced Usage

### Multiple configurations

Export dashboards for different reference modes:

```python
configs = [
    {"mode": "normalized", "label": "normalized"},
    {"mode": "corrected", "label": "corrected"},
    {"mode": "raw", "label": "raw"},
]

for config in configs:
    reference = VariantReference(
        variant="oci_2260",
        population="own",
        mode=config["mode"]
    )
    index = build_index(companion_dir, ref_dir, index_config, reference)
    
    export_html_dashboard(
        index,
        f"dashboard_{config['label']}.html",
        title=f"ARCSIX - {config['label'].title()}"
    )
```

### Custom threshold ranges

For detailed phase boundary analysis:

```python
# Fine-grained threshold exploration
thresholds = [1.0 + i * 0.05 for i in range(21)]  # 1.0 to 2.0 by 0.05

export_html_dashboard(
    index,
    "detailed_thresholds.html",
    thresholds=thresholds,
    default_threshold=1.27,
)
```

### Regional focus

For specific geographic regions:

```python
# Arctic region
export_html_dashboard(
    index,
    "arctic_dashboard.html",
    map_extent=(-85, -10, 58, 87),
    title="Arctic ARCSIX Analysis"
)

# Antarctic region
export_html_dashboard(
    index,
    "antarctic_dashboard.html",
    map_extent=(-180, 180, -90, -60),
    title="Antarctic Analysis"
)
```

## Comparison with PNG Export

| Feature | HTML Export | PNG Export (presentation) |
|---------|-------------|---------------------------|
| File format | Single HTML | Multiple PNG files |
| Interactivity | Threshold slider | Static images |
| File size | 5-15 MB | 2-5 MB per PNG |
| Sharing | Single file | Multiple files + manifest |
| Viewing | Any browser | Image viewer or slides |
| Threshold switching | Instant | Requires new file |
| Best for | Exploration & sharing | Presentations & papers |

## Technical Details

### Image encoding
Dashboard panels are rendered as PNG images and encoded using base64. This increases file size by ~33% compared to raw binary PNG, but ensures single-file portability.

### JavaScript interactions
The threshold slider uses vanilla JavaScript (no frameworks). Threshold changes update the displayed image instantly by switching the base64 source.

### Responsive design
The dashboard uses responsive CSS that adapts to different screen sizes. Maps maintain aspect ratios, and the layout adjusts for mobile devices.

### Rendering backend
HTML export uses matplotlib's `Agg` backend (non-interactive) to avoid GUI dependencies. This works in headless environments like Jupyter servers.

## Troubleshooting

### Large file size
- Reduce DPI (try 120 or 100)
- Reduce number of thresholds
- Reduce preview_width

### Slow generation
- Reduce number of thresholds
- Reduce DPI
- Check available memory

### Missing data
Ensure the aggregation index was built successfully:
```python
print(f"Indexed samples: {index.metadata['totals']['indexed_samples']:,}")
```

### Display issues in browser
- Try a different modern browser
- Check browser console for JavaScript errors
- Verify HTML file wasn't corrupted during transfer

## Related Documentation

- [Dashboard usage](dashboard.md) - Interactive dashboard for analysis
- [ARCSIX validation](arcsix_validation.md) - ARCSIX campaign context
- [Algorithm details](algorithm_details.md) - Phase retrieval methodology
