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


@pytest.mark.parametrize(
    "color, expected",
    [
        ("#FFFFFF", (1.0, 1.0, 1.0)),
        ("000000", (0.0, 0.0, 0.0)),
        ("#0f8", (0.0, 1.0, 0.24620133)),
        ("#00FF8880", (0.0, 1.0, 0.24620133, 128 / 255)),
        ("#808080", (0.21586050, 0.21586050, 0.21586050)),
        # CSS (MolViewSpec, ChimeraX) wins for bare names: green is #008000
        ("green", (0.0, 0.21586050, 0.0)),
        ("Dark Slate Gray", mn.color.srgb_to_linear((47 / 255, 79 / 255, 79 / 255))),
        ("pymol:green", (0.0, 1.0, 0.0)),
        ("vmd:gray", mn.color.srgb_to_linear((0.35, 0.35, 0.35))),
        # names only PyMOL or VMD have need no prefix
        ("slate", mn.color.srgb_to_linear((0.5, 0.5, 1.0))),
        ("carbon", mn.color.srgb_to_linear((0.2, 1.0, 0.2))),
        ("iceblue", mn.color.srgb_to_linear((0.5, 0.5, 0.75))),
        ("tab:blue", mn.color.srgb_to_linear((31 / 255, 119 / 255, 180 / 255))),
    ],
)
def test_parse_color(color, expected):
    assert np.allclose(mn.color.parse_color(color), expected)


@pytest.mark.parametrize(
    "color", ["#12345", "#GGGGGG", "", "#-12345", "notacolor", "pymol:notacolor"]
)
def test_parse_color_invalid(color):
    with pytest.raises(ValueError):
        mn.color.parse_color(color)
