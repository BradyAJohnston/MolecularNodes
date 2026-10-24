from typing import cast
import bpy
import numpy as np
import pytest
from bpy.types import Camera, SpaceView3D
import molecularnodes as mn
from molecularnodes.utils import frame_mapper
from .constants import codes, data_dir


def _viewport_spaces():
    return [
        cast(SpaceView3D, space)
        for screen in bpy.data.screens
        for area in screen.areas
        if area.type == "VIEW_3D"
        for space in area.spaces
        if space.type == "VIEW_3D"
    ]


def test_view_distance_increases():
    context = bpy.context
    scene = context.scene
    assert scene
    camera = cast(Camera, scene.camera.data)
    assert camera.clip_end == pytest.approx(100.0)
    for space in _viewport_spaces():
        assert space.clip_end == pytest.approx(1000.0)
    bpy.ops.mn.import_molecule(code=codes[0], cache_dir=str(data_dir))
    assert camera.clip_end == pytest.approx(mn.utils._INCREASED_CLIP_END)
    for space in _viewport_spaces():
        assert space.clip_end == pytest.approx(mn.utils._INCREASED_CLIP_END)


def test_correct_1d():
    assert np.allclose(
        mn.utils.correct_periodic_1d(np.array((0.9, 0.1)), np.array((0.1, 0.9)), 1.0),
        np.array((1.1, -0.1)),
    )


def test_correct_1d_is_pure():
    # the correction must not touch its inputs — in-place mutation poisoned the
    # frame position cache, unwrapping trajectories further out every crossing
    value1 = np.array((0.9, 0.1))
    value2 = np.array((0.1, 0.9))
    mn.utils.correct_periodic_1d(value1, value2, 1.0)
    assert np.allclose(value2, (0.1, 0.9))


def test_correct_1d_multiple_wraps():
    # a jump of several box lengths is folded back in a single correction
    assert np.allclose(
        mn.utils.correct_periodic_1d(np.array((0.0,)), np.array((2.3,)), 1.0),
        np.array((0.3,)),
    )


def test_minimum_image_orthorhombic_matches_old():
    # orthorhombic boxes give the same positions as the orthorhombic-only
    # correction, and MDAnalysis's minimize_vectors agrees with both
    from MDAnalysis.lib.distances import minimize_vectors

    rng = np.random.default_rng(0)
    box = np.array([17.6, 8.2, 10.9, 90.0, 90.0, 90.0])
    pos_1 = (rng.random((1000, 3)) * box[:3]).astype(np.float32)
    pos_2 = ((pos_1 + rng.normal(0, 1.5, pos_1.shape)) % box[:3]).astype(np.float32)
    old = mn.utils.correct_periodic_positions(pos_1, pos_2, box[:3])
    new = mn.utils.minimum_image_positions(pos_1, pos_2, box)
    assert new.dtype == pos_2.dtype
    assert np.allclose(new, old, atol=1e-5)
    assert np.allclose(pos_1 + minimize_vectors(pos_2 - pos_1, box), old, atol=1e-5)


def test_minimum_image_triclinic():
    # a step out through the slanted face of a rhombic dodecahedron box is
    # brought back as the short step, where per-axis wrapping would not be
    from MDAnalysis.lib.mdamath import triclinic_vectors

    box = np.array([100.0, 100.0, 100.0, 60.0, 60.0, 90.0])
    a, b, c = triclinic_vectors(box)
    start = 0.5 * a + 0.5 * b + 0.98 * c
    step = np.array([0.0, 0.0, 2.0])
    pos_1 = np.array([start], dtype=np.float32)
    pos_2 = np.array([start + step - c], dtype=np.float32)
    corrected = mn.utils.minimum_image_positions(pos_1, pos_2, box)
    assert corrected.dtype == np.float32
    assert np.allclose(corrected, start + step, atol=1e-4)
    # inputs are not modified
    assert np.allclose(pos_2, start + step - c)


@pytest.mark.parametrize(
    "dimensions,expected",
    [
        (None, False),
        ([0.0, 0.0, 0.0, 90.0, 90.0, 90.0], False),
        ([1.0, 1.0, 1.0, 90.0, 90.0, 90.0], False),
        ([10.0, 10.0, 10.0, 90.0, 90.0, 200.0], False),
        ([10.0, 10.0, np.nan, 90.0, 90.0, 90.0], False),
        ([10.0, 12.0, 14.0, 90.0, 90.0, 90.0], True),
        ([98.0, 98.0, 130.0, 90.0, 90.0, 120.0], True),
        ([100.0, 100.0, 100.0, 60.0, 60.0, 90.0], True),
    ],
)
def test_is_periodic_box(dimensions, expected):
    assert mn.utils.is_periodic_box(dimensions) == expected


def test_frame_mapper_basic():
    assert frame_mapper(10) == 10
    assert frame_mapper(0) == 0
    assert frame_mapper(-2) == 0


def test_frame_mapper_with_offset():
    assert frame_mapper(10, offset=2) == 8
    assert frame_mapper(2, offset=2) == 0
    assert frame_mapper(0, offset=5) == 0


def test_frame_mapper_with_subframes():
    assert frame_mapper(10, subframes=1) == 5
    assert frame_mapper(9, subframes=2) == 3


def test_frame_mapper_with_mapping():
    mapping = np.array([0, 0, 0, 1, 2])
    assert frame_mapper(1, mapping=mapping) == 0
    assert frame_mapper(3, mapping=mapping) == 1
    assert frame_mapper(3, mapping=list(mapping)) == 1


def test_frame_mapper_with_mapping_and_subframes():
    mapping = np.array([0, 0, 0, 1, 2])
    assert frame_mapper(5, subframes=1, mapping=mapping) == 0


def test_resolve_file_path(tmp_path):
    from molecularnodes.blender.utils import resolve_file_path

    existing = tmp_path / "structure.pdb"
    existing.write_text("")
    assert resolve_file_path(str(existing)) == existing
    assert resolve_file_path(existing, "Topology") == existing

    with pytest.raises(ValueError, match="Topology path is empty"):
        resolve_file_path("   ", "Topology")
    with pytest.raises(ValueError, match="File not found"):
        resolve_file_path(str(tmp_path / "missing.pdb"))
    # a directory is not an acceptable file
    with pytest.raises(ValueError, match="Map file not found"):
        resolve_file_path(str(tmp_path), "Map file")
