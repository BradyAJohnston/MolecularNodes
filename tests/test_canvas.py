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


def test_flat_settings(canvas):
    canvas.samples = 12
    assert canvas.samples == 12
    canvas.frame = 7
    assert canvas.frame == 7
    canvas.frame_range = (3, 42)
    assert canvas.frame_range == (3, 42)
    canvas.render_scale = 50
    assert canvas.render_scale == 50
    canvas.exposure = 0.25
    assert canvas.exposure == pytest.approx(0.25)
    canvas.gamma = 1.2
    assert canvas.gamma == pytest.approx(1.2)


def test_passes(canvas):
    canvas.passes = ["combined", "z", "mist"]
    assert set(canvas.passes) == {"combined", "z", "mist"}
    with pytest.raises(ValueError):
        canvas.passes = ["not_a_pass"]


def test_world(canvas):
    canvas.world.background = (0.1, 0.2, 0.3, 1.0)
    assert tuple(canvas.background) == pytest.approx((0.1, 0.2, 0.3, 1.0))
    canvas.world.hdri_strength = 1.5
    assert canvas.world.hdri_strength == pytest.approx(1.5)


def test_load_blend():
    # Test loading a .blend file
    file = data_dir / "blendfiles/suzanne.blend"
    suzanne_name = "MOLECULAR NODES Suzanne"
    assert not bpy.data.objects.get(suzanne_name)
    canvas = mn.Canvas(template=file)
    assert bpy.data.objects.get(suzanne_name)
    canvas.scene_reset(None)
    assert not bpy.data.objects.get(suzanne_name)
    canvas.load(file)
    assert bpy.data.objects.get(suzanne_name)
    with pytest.raises(ValueError):
        mn.Canvas(template="x")


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
    t1 = mn.Molecule(universe)

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


def test_compositor_setup(canvas):
    # the default compositor composites the annotations image on top of the
    # render with an Alpha Over node feeding the group output
    node_tree = canvas.scene.compositing_node_group
    assert node_tree is not None
    bl_idnames = [n.bl_idname for n in node_tree.nodes]
    assert "CompositorNodeRLayers" in bl_idnames
    assert "CompositorNodeAlphaOver" in bl_idnames
    assert "CompositorNodeImage" in bl_idnames

    alpha_over = next(
        n for n in node_tree.nodes if n.bl_idname == "CompositorNodeAlphaOver"
    )
    output_node = node_tree.nodes["Group Output"]
    assert alpha_over.outputs["Image"].links[0].to_node == output_node

    # the annotations image is the foreground composited on top
    image_node = next(
        n for n in node_tree.nodes if n.bl_idname == "CompositorNodeImage"
    )
    assert image_node.image.name == mn.scene.compositor.annotations_image


def test_compositor_reset_and_annotations(canvas):
    from nodebpy import compositor as c

    # reset clears everything (including annotations) back to render -> output
    with canvas.compositor.reset() as (image, output):
        image >> c.Glare(type="Bloom") >> output
    bl_idnames = [n.bl_idname for n in canvas.compositor.tree.nodes]
    assert "CompositorNodeGlare" in bl_idnames
    assert "CompositorNodeAlphaOver" not in bl_idnames  # annotations were cleared

    # annotations are opt-in after a reset
    canvas.compositor.add_annotations()
    bl_idnames = [n.bl_idname for n in canvas.compositor.tree.nodes]
    assert "CompositorNodeAlphaOver" in bl_idnames
