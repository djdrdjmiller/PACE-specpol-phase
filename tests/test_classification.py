import numpy as np
from pace_specpol.classification import sample_phase
from pace_specpol.validation.aggregation import PhaseIndex


def test_regions_and_ties():
    np.testing.assert_array_equal(
        sample_phase([1, 1.27, 1.27, 2, np.nan], [0, 0.29, 0.3, 3, 2]), [1, 2, 3, 3, 0]
    )
    assert PhaseIndex.dominant(
        np.array([[2, 1, 0], [2, 0, 0], [1, 0, 0]])
    ).tolist() == [0, 1, 0]
