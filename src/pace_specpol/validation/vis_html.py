"""HTML dashboard export for sharing ARCSIX analysis results."""

from datetime import datetime, timezone
from pathlib import Path
import json
import base64
from io import BytesIO
import numpy as np

from pace_specpol.oci_companion import atomic_json
from .vis_dashboard import dashboard


def export_html_dashboard(
    index,
    output_path,
    *,
    thresholds=(1.10, 1.27, 1.45),
    default_threshold=1.27,
    min_count=1,
    dpi=150,
    preview_width=900,
    ratio_clim=(0.5, 2.0),
    li_clim=(-0.5, 4.0),
    map_extent=None,
    graticules=True,
    title=None,
    **dashboard_options,
):
    """Export an interactive HTML dashboard with threshold slider.
    
    Creates a standalone HTML file that can be shared and viewed in any browser.
    The dashboard includes:
    - Interactive threshold slider
    - Joint histogram with phase regions
    - Modal phase map
    - Mean liquid index map
    - Mean ratio map
    - Summary statistics
    
    Parameters
    ----------
    index : IndexConfig
        Pre-built aggregation index from build_index()
    output_path : str or Path
        Output HTML file path
    thresholds : tuple of float
        Threshold values to pre-render (default: 1.10, 1.27, 1.45)
    default_threshold : float
        Initial threshold to display (default: 1.27)
    min_count : int
        Minimum sample count per cell (default: 1)
    dpi : int
        Figure resolution (default: 150)
    preview_width : int
        Map preview width in pixels (default: 900)
    ratio_clim : tuple
        Ratio colorbar limits (default: 0.5, 2.0)
    li_clim : tuple
        Liquid index colorbar limits (default: -0.5, 4.0)
    map_extent : tuple or None
        Regional extent (west, east, south, north) for Arctic maps
    graticules : bool
        Show graticule lines (default: True)
    title : str or None
        Custom dashboard title (default: auto-generated)
    **dashboard_options
        Additional options passed to dashboard()
    """
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    
    thresholds = sorted(set(float(t) for t in thresholds))
    if not thresholds:
        raise ValueError("Provide at least one threshold")
    if default_threshold not in thresholds:
        thresholds.append(default_threshold)
        thresholds.sort()
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Generate dashboard figures for each threshold
    frames = {}
    print(f"Generating {len(thresholds)} dashboard frames...")
    
    for i, threshold in enumerate(thresholds):
        print(f"  Rendering threshold {threshold:.2f} ({i+1}/{len(thresholds)})...")
        
        # Create dashboard figure
        ui = dashboard(
            index,
            threshold=threshold,
            min_count=min_count,
            preview_width=preview_width,
            ratio_clim=ratio_clim,
            li_clim=li_clim,
            map_extent=map_extent,
            layout="slide",
            graticules=graticules,
            show=False,
            **dashboard_options
        )
        
        fig = ui["figure"]
        
        # Render to PNG in memory
        buffer = BytesIO()
        with plt.rc_context({"savefig.bbox": None}):
            fig.savefig(buffer, format="png", dpi=dpi, bbox_inches=None, facecolor="white")
        
        # Encode as base64
        buffer.seek(0)
        img_data = base64.b64encode(buffer.read()).decode('utf-8')
        frames[threshold] = img_data
        
        plt.close(fig)
    
    # Collect metadata
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "index_fingerprint": index.metadata["fingerprint"],
        "reference": index.metadata["reference"],
        "data_scope": str(index.metadata.get("data_scope", "")),
        "thresholds": thresholds,
        "default_threshold": default_threshold,
        "min_count": min_count,
        "liquid_index_threshold": index.config.li_threshold,
        "resolution_deg": index.config.resolution,
        "ratio_clim": ratio_clim,
        "li_clim": li_clim,
        "total_samples": index.metadata["totals"]["indexed_samples"],
        "histogram_outside": index.metadata["totals"]["histogram_outside_samples"],
        "reference_invalid": index.metadata["totals"]["reference_invalid_samples"],
    }
    
    # Get sampling statistics
    sampling = index.sampling_summary()
    stats = {
        "mean_per_cell": f"{sampling['mean_samples_per_occupied_cell']:.2f}",
        "median": f"{sampling['median_samples']:.0f}",
        "p90": f"{sampling['p90_samples']:.0f}",
        "maximum": sampling['maximum_samples'],
        "single_sample_fraction": f"{sampling['single_sample_cell_fraction']:.1%}",
    }
    
    # Get mode label
    raw = index.metadata["reference"]["mode"] == "raw"
    corrected = index.metadata["reference"]["mode"] == "corrected"
    mode_label = (
        "RAW ratio; unnormalized — phase candidates only"
        if raw
        else "Atmospherically corrected ratio — provisional phase candidates"
        if corrected
        else "Normalized ratio — provisional phase candidates"
    )
    if "variant" in index.metadata["reference"]:
        mode_label += (
            f" | {index.metadata['reference']['variant']}"
            f" | {index.metadata['reference'].get('population', 'own')} population"
        )
    
    if title is None:
        title = "ARCSIX Phase Dashboard"
    
    # Generate HTML
    html_content = _generate_html(
        frames=frames,
        default_threshold=default_threshold,
        metadata=metadata,
        stats=stats,
        mode_label=mode_label,
        title=title,
    )
    
    # Write HTML file
    output_path.write_text(html_content, encoding='utf-8')
    print(f"\nHTML dashboard saved to: {output_path}")
    print(f"File size: {output_path.stat().st_size / 1024 / 1024:.1f} MB")
    
    return output_path


def _generate_html(frames, default_threshold, metadata, stats, mode_label, title):
    """Generate the HTML content with embedded images and JavaScript."""
    
    thresholds_json = json.dumps(metadata["thresholds"])
    
    # Build frame data for JavaScript
    frame_data = {}
    for threshold, img_data in frames.items():
        frame_data[f"{threshold:.2f}"] = img_data
    frame_data_json = json.dumps(frame_data)
    
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
            color: #333;
        }}
        
        .container {{
            max-width: 1800px;
            margin: 0 auto;
            background: white;
            border-radius: 16px;
            box-shadow: 0 20px 60px rgba(0, 0, 0, 0.3);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, #1e3c72 0%, #2a5298 100%);
            color: white;
            padding: 30px 40px;
            border-bottom: 4px solid #667eea;
        }}
        
        .header h1 {{
            font-size: 2.2em;
            margin-bottom: 10px;
            font-weight: 600;
        }}
        
        .header .subtitle {{
            font-size: 1.1em;
            opacity: 0.95;
            line-height: 1.5;
        }}
        
        .controls {{
            background: #f8f9fa;
            padding: 25px 40px;
            border-bottom: 1px solid #dee2e6;
        }}
        
        .control-group {{
            display: flex;
            align-items: center;
            gap: 20px;
            margin-bottom: 15px;
            flex-wrap: wrap;
        }}
        
        .control-group label {{
            font-weight: 600;
            font-size: 1.05em;
            color: #495057;
        }}
        
        .slider-container {{
            flex: 1;
            min-width: 300px;
            max-width: 600px;
        }}
        
        input[type="range"] {{
            width: 100%;
            height: 8px;
            border-radius: 5px;
            background: #d3d3d3;
            outline: none;
            -webkit-appearance: none;
        }}
        
        input[type="range"]::-webkit-slider-thumb {{
            -webkit-appearance: none;
            appearance: none;
            width: 24px;
            height: 24px;
            border-radius: 50%;
            background: #667eea;
            cursor: pointer;
            box-shadow: 0 2px 6px rgba(0, 0, 0, 0.2);
            transition: all 0.2s;
        }}
        
        input[type="range"]::-webkit-slider-thumb:hover {{
            background: #764ba2;
            transform: scale(1.1);
        }}
        
        input[type="range"]::-moz-range-thumb {{
            width: 24px;
            height: 24px;
            border-radius: 50%;
            background: #667eea;
            cursor: pointer;
            border: none;
            box-shadow: 0 2px 6px rgba(0, 0, 0, 0.2);
            transition: all 0.2s;
        }}
        
        input[type="range"]::-moz-range-thumb:hover {{
            background: #764ba2;
            transform: scale(1.1);
        }}
        
        .threshold-value {{
            font-size: 1.3em;
            font-weight: 700;
            color: #667eea;
            min-width: 80px;
            text-align: center;
            background: white;
            padding: 8px 16px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }}
        
        .dashboard-image {{
            width: 100%;
            height: auto;
            display: block;
            background: white;
        }}
        
        .metadata {{
            padding: 30px 40px;
            background: #f8f9fa;
            border-top: 1px solid #dee2e6;
        }}
        
        .metadata h2 {{
            font-size: 1.5em;
            margin-bottom: 20px;
            color: #1e3c72;
        }}
        
        .info-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
        }}
        
        .info-card {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
        }}
        
        .info-card h3 {{
            font-size: 1.1em;
            color: #667eea;
            margin-bottom: 12px;
            font-weight: 600;
        }}
        
        .info-card p {{
            line-height: 1.6;
            color: #495057;
            margin-bottom: 8px;
        }}
        
        .info-card strong {{
            color: #1e3c72;
        }}
        
        .loading {{
            text-align: center;
            padding: 60px;
            font-size: 1.2em;
            color: #6c757d;
        }}
        
        .footer {{
            padding: 20px 40px;
            text-align: center;
            background: #e9ecef;
            color: #6c757d;
            font-size: 0.9em;
        }}
        
        @media (max-width: 768px) {{
            body {{
                padding: 10px;
            }}
            
            .header {{
                padding: 20px;
            }}
            
            .header h1 {{
                font-size: 1.6em;
            }}
            
            .controls {{
                padding: 20px;
            }}
            
            .control-group {{
                flex-direction: column;
                align-items: stretch;
            }}
            
            .metadata {{
                padding: 20px;
            }}
            
            .info-grid {{
                grid-template-columns: 1fr;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>{title}</h1>
            <div class="subtitle">{mode_label}</div>
            <div class="subtitle" style="margin-top: 8px; font-size: 0.95em;">{metadata['data_scope']}</div>
        </div>
        
        <div class="controls">
            <div class="control-group">
                <label for="threshold-slider">Ratio Threshold:</label>
                <div class="slider-container">
                    <input type="range" 
                           id="threshold-slider" 
                           min="0" 
                           max="{len(metadata['thresholds']) - 1}" 
                           value="{metadata['thresholds'].index(default_threshold)}" 
                           step="1">
                </div>
                <div class="threshold-value" id="threshold-value">{default_threshold:.2f}</div>
            </div>
            <div style="font-size: 0.95em; color: #6c757d; margin-top: 5px;">
                Adjust the ratio threshold to see how phase classification changes. 
                Mean fields remain constant; only phase counts and ratio color centering update.
            </div>
        </div>
        
        <div id="dashboard-container">
            <img id="dashboard-img" class="dashboard-image" alt="Phase Dashboard">
        </div>
        
        <div class="metadata">
            <h2>Metadata & Statistics</h2>
            <div class="info-grid">
                <div class="info-card">
                    <h3>Sample Totals</h3>
                    <p><strong>Indexed samples:</strong> {metadata['total_samples']:,}</p>
                    <p><strong>Outside histogram:</strong> {metadata['histogram_outside']:,}</p>
                    <p><strong>Reference invalid:</strong> {metadata['reference_invalid']:,}</p>
                    <p><strong>Min count per cell:</strong> {metadata['min_count']}</p>
                </div>
                
                <div class="info-card">
                    <h3>Cell Statistics</h3>
                    <p><strong>Mean per occupied cell:</strong> {stats['mean_per_cell']}</p>
                    <p><strong>Median:</strong> {stats['median']}</p>
                    <p><strong>90th percentile:</strong> {stats['p90']}</p>
                    <p><strong>Maximum:</strong> {stats['maximum']}</p>
                    <p><strong>Single-sample cells:</strong> {stats['single_sample_fraction']}</p>
                </div>
                
                <div class="info-card">
                    <h3>Configuration</h3>
                    <p><strong>Resolution:</strong> {metadata['resolution_deg']:.4f}°</p>
                    <p><strong>Liquid index threshold:</strong> {metadata['liquid_index_threshold']}</p>
                    <p><strong>Ratio limits:</strong> {metadata['ratio_clim'][0]:.2f} – {metadata['ratio_clim'][1]:.2f}</p>
                    <p><strong>LI limits:</strong> {metadata['li_clim'][0]:.2f} – {metadata['li_clim'][1]:.2f}</p>
                </div>
                
                <div class="info-card">
                    <h3>Reference Mode</h3>
                    <p><strong>Mode:</strong> {metadata['reference']['mode']}</p>
                    <p><strong>Variant:</strong> {metadata['reference'].get('variant', 'N/A')}</p>
                    <p><strong>Population:</strong> {metadata['reference'].get('population', 'N/A')}</p>
                    <p style="margin-top: 12px; font-size: 0.9em; color: #6c757d;">
                        <strong>Index fingerprint:</strong><br>{metadata['index_fingerprint'][:16]}...
                    </p>
                </div>
            </div>
        </div>
        
        <div class="footer">
            Generated on {metadata['created_utc'][:19].replace('T', ' ')} UTC | 
            PACE Spectropolarimetric Phase Retrieval
        </div>
    </div>
    
    <script>
        // Embedded frame data
        const thresholds = {thresholds_json};
        const frameData = {frame_data_json};
        
        // Get elements
        const slider = document.getElementById('threshold-slider');
        const valueDisplay = document.getElementById('threshold-value');
        const dashboardImg = document.getElementById('dashboard-img');
        
        // Update dashboard when slider changes
        function updateDashboard() {{
            const index = parseInt(slider.value);
            const threshold = thresholds[index];
            const thresholdKey = threshold.toFixed(2);
            
            valueDisplay.textContent = thresholdKey;
            
            if (frameData[thresholdKey]) {{
                dashboardImg.src = 'data:image/png;base64,' + frameData[thresholdKey];
            }}
        }}
        
        // Initialize with default threshold
        updateDashboard();
        
        // Listen for slider changes
        slider.addEventListener('input', updateDashboard);
    </script>
</body>
</html>
"""
    
    return html
