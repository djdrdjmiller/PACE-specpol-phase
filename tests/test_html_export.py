"""Test HTML dashboard export functionality."""

import pytest
from pathlib import Path
import tempfile


def test_html_export_basic(sample_index):
    """Test basic HTML export with default parameters."""
    from pace_specpol.validation.vis_html import export_html_dashboard
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "dashboard.html"
        
        result = export_html_dashboard(
            sample_index,
            output_path,
            thresholds=(1.1, 1.27, 1.45),
            default_threshold=1.27,
        )
        
        assert result.exists()
        assert result.stat().st_size > 0
        
        # Check HTML content
        content = result.read_text(encoding='utf-8')
        assert '<!DOCTYPE html>' in content
        assert 'ARCSIX Phase Dashboard' in content
        assert 'data:image/png;base64,' in content
        assert 'threshold-slider' in content


def test_html_export_custom_title(sample_index):
    """Test HTML export with custom title."""
    from pace_specpol.validation.vis_html import export_html_dashboard
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "custom_dashboard.html"
        custom_title = "My Custom ARCSIX Analysis"
        
        result = export_html_dashboard(
            sample_index,
            output_path,
            title=custom_title,
            thresholds=(1.27,),
        )
        
        content = result.read_text(encoding='utf-8')
        assert custom_title in content


def test_html_export_multiple_thresholds(sample_index):
    """Test HTML export with multiple threshold values."""
    from pace_specpol.validation.vis_html import export_html_dashboard
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "multi_threshold.html"
        
        thresholds = (1.0, 1.1, 1.2, 1.3, 1.4, 1.5)
        
        result = export_html_dashboard(
            sample_index,
            output_path,
            thresholds=thresholds,
            default_threshold=1.3,
        )
        
        content = result.read_text(encoding='utf-8')
        
        # Check that all thresholds are in the JavaScript data
        for threshold in thresholds:
            assert f'"{threshold:.2f}"' in content or f'{threshold:.2f}' in content


def test_html_export_invalid_threshold(sample_index):
    """Test that invalid threshold configuration raises error."""
    from pace_specpol.validation.vis_html import export_html_dashboard
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "invalid.html"
        
        # Empty thresholds should raise error
        with pytest.raises(ValueError, match="at least one threshold"):
            export_html_dashboard(
                sample_index,
                output_path,
                thresholds=(),
            )


def test_html_export_creates_directory(sample_index):
    """Test that HTML export creates parent directory if needed."""
    from pace_specpol.validation.vis_html import export_html_dashboard
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "subdir" / "nested" / "dashboard.html"
        
        assert not output_path.parent.exists()
        
        result = export_html_dashboard(
            sample_index,
            output_path,
            thresholds=(1.27,),
        )
        
        assert result.exists()
        assert output_path.parent.exists()


@pytest.fixture
def sample_index():
    """Create a sample index for testing."""
    from pace_specpol.validation.aggregation import IndexConfig, SparseIndex
    import numpy as np
    
    # Create a minimal test index
    config = IndexConfig(
        resolution=0.1,
        ratio_step=0.01,
        li_threshold=0.3,
        ratio_edges=np.linspace(0.5, 2.0, 31),
        li_edges=np.linspace(-0.5, 4.0, 31),
    )
    
    # Create sparse arrays with some sample data
    shape = (3600, 7200)  # 0.1 degree global grid
    histogram = np.zeros((30, 30), dtype=np.int32)
    histogram[10:20, 10:20] = np.random.poisson(100, (10, 10))
    
    # Create some sample cell data
    n_thresholds = len(config.ratio_edges) - 1
    n_cells = 100
    
    low = np.random.poisson(50, (n_thresholds, n_cells)).astype(np.int32)
    high = np.random.poisson(50, (n_thresholds, n_cells)).astype(np.int32)
    
    metadata = {
        "fingerprint": "test_fingerprint_12345678",
        "reference": {
            "mode": "normalized",
            "variant": "oci_2260",
            "population": "own"
        },
        "data_scope": "Test Arctic scene | 2024-06-10",
        "totals": {
            "indexed_samples": 10000,
            "histogram_outside_samples": 500,
            "reference_invalid_samples": 200,
        }
    }
    
    index = SparseIndex(
        config=config,
        histogram=histogram,
        low=low,
        high=high,
        metadata=metadata,
    )
    
    return index
