import platform
import numpy as np
import pytest
import molecularnodes as mn
from .utils import NumpySnapshotExtension


@pytest.mark.skipif(
    platform.system() != "Darwin" or platform.machine() != "arm64",
    reason="Test only runs on macOS ARM",
)
def test_random_rgb(snapshot_custom):
    n = 100
    colors = np.array(list(map(mn.color.random_rgb, range(n))))
    assert snapshot_custom == colors


def test_colos_from_atomic_numbers(snapshot_custom: NumpySnapshotExtension):
    length = len(mn.color.iupac_colors_rgb)
    colors = mn.color.colors_from_elements(np.array(list(range(length))))
    assert snapshot_custom == colors
