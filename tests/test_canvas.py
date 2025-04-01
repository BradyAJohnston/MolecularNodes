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


def test_render_engine(canvas):
    # Test setting different render engines
    engines = ["CYCLES", "EEVEE", "WORKBENCH"]
    for engine in engines:
        canvas.render_engine = engine
        if engine == "EEVEE":
            # Account for EEVEE_NEXT in newer Blender versions
            assert "EEVEE" in canvas.render_engine
        else:
            assert engine in canvas.render_engine


def test_samples(canvas):
    # Test Cycles samples
    canvas.samples_cycles = 128
    assert canvas.samples_cycles == 128

    # Test EEVEE samples
    canvas.samples_eevee = 64
    assert canvas.samples_eevee == 64


def test_cycles_device(canvas):
    # Test setting render device
    canvas.cycles_device = "CPU"
    assert canvas.cycles_device == "CPU"

    with pytest.raises(ValueError):
        canvas.cycles_device = "INVALID"


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
