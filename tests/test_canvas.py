import bpy
import pytest
from molecularnodes.scene.base import Canvas


@pytest.fixture
def canvas():
    return Canvas()


def test_resolution(canvas):
    # Test getting and setting resolution
    canvas.resolution = (1920, 1080)
    assert canvas.resolution == (1920, 1080)


def test_animation_settings(canvas):
    # Test FPS
    canvas.fps = 30
    assert canvas.fps == 30

    # Test frame range
    canvas.frame_start = 1
    canvas.frame_end = 250
    assert canvas.frame_start == 1
    assert canvas.frame_end == 250


def test_frame_object(canvas):
    # Create test object
    bpy.ops.mesh.primitive_cube_add()
    test_obj = bpy.context.active_object

    # Test framing
    canvas.frame_object(test_obj)

    # Cleanup
    bpy.data.objects.remove(test_obj, do_unlink=True)
