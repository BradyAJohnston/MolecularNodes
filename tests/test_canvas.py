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
    assert not bpy.data.objects.get("Suzanne")
    canvas = mn.Canvas(template=file)
    assert bpy.data.objects.get("Suzanne")
    canvas.scene_reset(None)
    assert not bpy.data.objects.get("Suzanne")
    assert bpy.data.objects.get("Cube")
    canvas.load(file)
    assert bpy.data.objects.get("Suzanne")
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


def test_look_at_object(canvas):
    # Create test object
    bpy.ops.mesh.primitive_cube_add()
    test_obj = bpy.context.active_object

    # Test looking at a plain Blender object
    canvas.look_at(test_obj)

    # Cleanup
    bpy.data.objects.remove(test_obj, do_unlink=True)


def test_look_at_views(canvas, universe):
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
    canvas.look_at(v1)
    l1 = camera.location.copy()
    assert l1 != initial_location
    # frame v13 (resid 1 of trajectory at frame 3)
    canvas.look_at(v13)
    l13 = camera.location.copy()
    assert l13 != l1
    # frame v2 (resid 2 of trajectory)
    canvas.look_at(v2)
    l2 = camera.location.copy()
    assert l2 != l1 and l2 != l13
    # frame whole trajectory
    canvas.look_at(v0)
    l0 = camera.location.copy()
    assert l0 != initial_location

    # looking at the entity itself dispatches through its current view
    canvas.look_at(t1)
    assert (camera.location - l0).length < 1e-6

    # test different viewpoints
    canvas.look_at(v1, viewpoint="front")
    r1f = camera.rotation_euler.copy()
    canvas.look_at(v1, viewpoint="back")
    r1b = camera.rotation_euler.copy()
    assert r1f != r1b

    # a custom viewpoint is an XYZ Euler rotation in degrees
    canvas.look_at(v1, viewpoint=(90, 0, -45))
    assert canvas.camera.rotation == pytest.approx((90, 0, -45), abs=1e-3)

    # test framing multiple views
    # frame v1 and v2 (resid 1 and resid 2 of trajectory)
    canvas.look_at(v1 + v2)
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


def test_scene_reset_rebinds(canvas):
    # the world and compositor builders are bound to the scene's node trees,
    # so they must be rebound after the scene is replaced
    canvas.world.background = (1.0, 0.0, 0.0, 1.0)
    _ = canvas.compositor.tree
    canvas.scene_reset()
    # builders work again and the annotation compositor is re-created
    canvas.world.background = (0.0, 1.0, 0.0, 1.0)
    assert tuple(canvas.background) == pytest.approx((0.0, 1.0, 0.0, 1.0))
    assert canvas.scene.compositing_node_group is not None
    bl_idnames = [n.bl_idname for n in canvas.compositor.tree.nodes]
    assert "CompositorNodeAlphaOver" in bl_idnames


def test_engine_syncs_with_scene(canvas):
    # an engine switch made outside of the Canvas is reflected, without the
    # constructor defaults overwriting the scene's settings
    canvas.scene.render.engine = "CYCLES"
    canvas.scene.cycles.samples = 7
    assert isinstance(canvas.engine, mn.scene.Cycles)
    assert canvas.samples == 7
    canvas.scene.render.engine = "BLENDER_EEVEE"
    assert isinstance(canvas.engine, mn.scene.EEVEE)


def test_engine_settings_applied_on_activation(canvas):
    # constructing an engine stores its settings without touching the scene;
    # they are applied when the engine is activated on the Canvas
    engine = mn.scene.Cycles(samples=1234, device="CPU")
    assert canvas.scene.cycles.samples != 1234
    canvas.engine = engine
    assert canvas.scene.render.engine == "CYCLES"
    assert canvas.samples == 1234
    assert canvas.engine.device == "CPU"


@pytest.fixture
def render_canvas():
    # cheapest possible renders: tiny resolution, Cycles with a single sample
    # (EEVEE requires a GPU context that isn't available headless)
    canvas = mn.Canvas(resolution=(32, 32))
    canvas.engine = "CYCLES"
    canvas.samples = 1
    return canvas


def test_snapshot(render_canvas, tmp_path):
    out = tmp_path / "snapshot.png"
    render_canvas.frame = 3
    image = render_canvas.snapshot(out, frame=7, render_scale=50)
    # saved to the given path as a PNG
    assert out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    # the current frame is restored after rendering a different frame
    assert render_canvas.frame == 3
    # render_scale is restored
    assert render_canvas.render_scale == 100
    # returns a displayable image when IPython is available
    if mn.scene.base.Image is not None:
        assert isinstance(image, mn.scene.base.Image)


def test_snapshot_non_displayable_format(render_canvas, tmp_path):
    out = tmp_path / "snapshot.tif"
    # TIFF renders and saves fine but can't be embedded in a notebook
    image = render_canvas.snapshot(out, file_format="TIFF")
    assert out.stat().st_size > 0
    assert image is None


def test_animation_mp4(render_canvas, tmp_path):
    out = tmp_path / "animation.mp4"
    original_fps = render_canvas.fps
    video = render_canvas.animation(out, frame_start=1, frame_end=3, fps=12)
    # an ftyp box marks an MP4 container
    assert out.read_bytes()[4:8] == b"ftyp"
    # the fps override is restored
    assert render_canvas.fps == original_fps
    if mn.scene.base.Video is not None:
        assert isinstance(video, mn.scene.base.Video)


def test_animation_gif(render_canvas, tmp_path):
    PILImage = pytest.importorskip("PIL.Image")
    out = tmp_path / "animation.gif"
    # animate the exposure so the frames differ, otherwise pillow merges
    # identical consecutive frames into one
    scene = render_canvas.scene
    scene.view_settings.exposure = -10
    scene.keyframe_insert(data_path="view_settings.exposure", frame=1)
    scene.view_settings.exposure = 10
    scene.keyframe_insert(data_path="view_settings.exposure", frame=3)
    image = render_canvas.animation(out, frame_start=1, frame_end=3)
    assert out.read_bytes()[:6] in (b"GIF87a", b"GIF89a")
    with PILImage.open(out) as gif:
        assert gif.n_frames == 3
    if mn.scene.base.Image is not None:
        assert isinstance(image, mn.scene.base.Image)


def test_animation_validation(canvas):
    with pytest.raises(ValueError, match="cannot be less than"):
        canvas.animation(frame_start=10, frame_end=5)
    with pytest.raises(ValueError, match="must be either"):
        canvas.animation(frame_start=1, frame_end=2, format="AVI")


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
