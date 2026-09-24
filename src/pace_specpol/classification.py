"""Per-observation provisional phase classification."""

import numpy as np


def sample_phase(ratio, li, ratio_threshold=1.27, li_threshold=0.3):
    """Reference-figure regions: 1 liquid, 2 ice, 3 LTMP, 0 invalid."""
    ratio, li = np.broadcast_arrays(ratio, li)
    result = np.zeros(ratio.shape, dtype="int8")
    ok = np.isfinite(ratio) & np.isfinite(li)
    result[ok & (ratio < ratio_threshold)] = 1
    result[ok & (ratio >= ratio_threshold) & (li < li_threshold)] = 2
    result[ok & (ratio >= ratio_threshold) & (li >= li_threshold)] = 3
    return result
