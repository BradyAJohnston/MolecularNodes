from typing import cast
import bpy
import numpy as np
import pytest
from bpy.types import Camera, SpaceView3D
import molecularnodes as mn
from molecularnodes.utils import frame_mapper
from .constants import codes


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
    bpy.ops.mn.import_molecule(code=codes[0])
    assert camera.clip_end == pytest.approx(mn.utils._INCREASED_CLIP_END)
    for space in _viewport_spaces():
        assert space.clip_end == pytest.approx(mn.utils._INCREASED_CLIP_END)


def test_correct_1d():
    assert np.allclose(
        mn.utils.correct_periodic_1d(np.array((0.9, 0.1)), np.array((0.1, 0.9)), 1.0),
        np.array((1.1, -0.1)),
    )


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
