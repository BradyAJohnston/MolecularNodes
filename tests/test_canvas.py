import bpy
import MDAnalysis as mda
import pytest
import molecularnodes as mn
from .constants import data_dir


@pytest.fixture
def canvas():
    return mn.scene.base.Canvas()


@pytest.fixture()
def universe():
    topo = data_dir / "md_ppr/box.gro"
    traj = data_dir / "md_ppr/first_5_frames.xtc"
    return mda.Universe(topo, traj)


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


def test_frame_view(canvas, universe):
    t1 = mn.Trajectory(universe)

    # save views for later framing
    # view of resid 1
    v1 = t1.get_view(selection="resid 1")
    # view of resid 1 at trajectory frame 3
    v13 = t1.get_view(selection="resid 1", frame=3)
    # view of resid 2
    v2 = t1.get_view(selection="resid 2")
    # view of whole trajectory
    v0 = t1.get_view()

    camera = bpy.context.scene.camera
    initial_location = camera.location.copy()

    # frame v1 (resid 1 of trajectory)
    canvas.frame_view(v1)
    l1 = camera.location.copy()
    assert l1 != initial_location
    # frame v13 (resid 1 of trajectory at frame 3)
    canvas.frame_view(v13)
    l13 = camera.location.copy()
    assert l13 != l1
    # frame v2 (resid 2 of trajectory)
    canvas.frame_view(v2)
    l2 = camera.location.copy()
    assert l2 != l1 and l2 != l13
    # frame whole trajectory
    canvas.frame_view(v0)
    l0 = camera.location.copy()
    assert l0 != initial_location

    # test different viewpoints
    canvas.frame_view(v1, viewpoint="front")
    r1f = camera.rotation_euler.copy()
    canvas.frame_view(v1, viewpoint="back")
    r1b = camera.rotation_euler.copy()
    assert r1f != r1b

    # test framing multiple views
    # frame v1 and v2 (resid 1 and resid 2 of trajectory)
    canvas.frame_view(v1 + v2)
    l12 = camera.location.copy()
    assert l12 != l1 and l12 != l2

    # test camera lens (zoom in and out)
    # location should not change until re-framed
    canvas.camera.lens = 35
    l3 = camera.location.copy()
    assert l3 == l12
    canvas.camera.lens = 85
    l4 = camera.location.copy()
    assert l4 == l3 and l4 == l12
