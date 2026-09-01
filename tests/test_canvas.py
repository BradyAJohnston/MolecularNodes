import bpy
import MDAnalysis as mda
import numpy as np
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
    canvas.load_preset(None)
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


def test_load_preset_rebinds(canvas):
    # the world and compositor builders are bound to the scene's node trees,
    # so they must be rebound after the scene is replaced
    canvas.world.background = (1.0, 0.0, 0.0, 1.0)
    _ = canvas.compositor.tree
    canvas.load_preset()
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
    # (EEVEE requires a GPU that GHA runners don't have access to). CPU
    # explicitly - both for rendering and for the render-time compositor,
    # which defaults to GPU in Blender 5.x - as the CI runners have no GPU
    canvas = mn.Canvas(resolution=(32, 32))
    canvas.engine = mn.scene.Cycles(samples=1, device="CPU")
    canvas.compositor.device = "CPU"
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


def test_compositor_settings(canvas):
    comp = canvas.compositor
    # values are case-insensitive and write through to the scene
    comp.device = "cpu"
    assert comp.device == "CPU"
    assert canvas.scene.render.compositor_device == "CPU"
    comp.precision = "full"
    assert comp.precision == "FULL"
    comp.denoise_device = "cpu"
    assert comp.denoise_device == "CPU"
    comp.denoise_preview_quality = "fast"
    assert comp.denoise_preview_quality == "FAST"
    comp.denoise_final_quality = "balanced"
    assert comp.denoise_final_quality == "BALANCED"
    # invalid values are rejected by Blender's enum validation
    with pytest.raises(TypeError):
        comp.device = "TPU"


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


def test_clear_removes_content_but_keeps_lighting(canvas):
    "clear() removes content objects, not just the Molecular Nodes entities."
    cube = bpy.data.objects.new("UserCube", bpy.data.meshes.new("UserCubeMesh"))
    bpy.context.collection.objects.link(cube)
    mn.Molecule.fetch("4ozs").add_style("cartoon")

    canvas.clear()
    remaining = {obj.name: obj.type for obj in bpy.data.objects}
    # the molecule, the user's cube and the backdrop are content and go
    assert "4ozs" not in remaining
    assert "UserCube" not in remaining
    assert "Backdrop" not in remaining
    # the camera and lights are how the scene is rendered, and stay
    assert set(remaining.values()) == {"CAMERA", "LIGHT"}
    assert canvas.scene.camera is not None
    assert len(mn.session.get_session().entities) == 0


def test_clear_preserves_lighting_in_the_render(canvas):
    "Removing the lights would silently flatten every subsequent render."
    lights_before = [o.name for o in bpy.data.objects if o.type == "LIGHT"]
    mn.Molecule.fetch("4ozs").add_style("cartoon")
    canvas.clear()
    assert [o.name for o in bpy.data.objects if o.type == "LIGHT"] == lights_before


def test_clear_keeps_render_settings_and_world(canvas):
    canvas.engine = "CYCLES"
    canvas.resolution = (400, 300)
    canvas.samples = 8
    world = canvas.scene.world
    mn.Molecule.fetch("4ozs").add_style("cartoon")

    canvas.clear()
    assert canvas.resolution == (400, 300)
    assert canvas.samples == 8
    assert canvas.scene.render.engine == "CYCLES"
    assert canvas.scene.world == world
    assert canvas.scene.compositing_node_group is not None


def test_clear_does_not_leak_datablocks(canvas):
    """A style leaves >100 node groups behind, which used to accumulate.

    Only a recursive purge collects them, as they hang off the object's
    modifier tree rather than being directly orphaned.
    """
    mn.Molecule.fetch("4ozs").add_style("cartoon")
    canvas.clear()
    baseline = (len(bpy.data.node_groups), len(bpy.data.meshes))

    for _ in range(3):
        mn.Molecule.fetch("4ozs").add_style("cartoon")
        canvas.clear()
        assert (len(bpy.data.node_groups), len(bpy.data.meshes)) == baseline


def test_load_preset_restores_the_preset_scene(canvas):
    canvas.clear()
    assert "Backdrop" not in {obj.name for obj in bpy.data.objects}

    canvas.load_preset()
    names = {obj.name for obj in bpy.data.objects}
    assert {"Backdrop", "Camera", "Key Light", "Rim Light"} <= names


def test_load_preset_prunes_dangling_entities(canvas):
    "Objects do not survive the scene swap, so their entities must not either."
    mn.Molecule.fetch("4ozs").add_style("cartoon")
    assert len(mn.session.get_session().entities) == 1

    canvas.load_preset()
    assert len(mn.session.get_session().entities) == 0


def test_clear_survives_an_externally_deleted_object(canvas):
    "Deleting a molecule from the outliner used to make clear() raise."
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    bpy.data.objects.remove(mol.object, do_unlink=True)
    assert len(mn.session.get_session().entities) == 1

    canvas.clear()
    assert len(mn.session.get_session().entities) == 0


def test_clear_leaves_other_scenes_alone(canvas):
    "clear() empties this canvas's scene, not every scene in the file."
    other = bpy.data.scenes.new("OtherScene")
    cube = bpy.data.objects.new("OtherCube", bpy.data.meshes.new("OtherCubeMesh"))
    other.collection.objects.link(cube)
    mn.Molecule.fetch("4ozs").add_style("cartoon")

    canvas.clear()
    assert "4ozs" not in canvas.scene.objects
    assert "OtherCube" in other.objects
    bpy.data.scenes.remove(other)


def test_canvas_loads_the_preset_into_an_empty_scene():
    "The out-of-the-box lighting still arrives on a first Canvas()."
    canvas = mn.Canvas()
    names = {obj.name for obj in canvas.scene.objects}
    assert {"Backdrop", "Camera", "Key Light", "Rim Light"} <= names


def test_canvas_rerun_does_not_wipe_the_scene(canvas):
    """Re-running `mn.Canvas()` used to destroy every molecule silently.

    A notebook cell doing exactly this re-runs on its own in marimo.
    """
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    session = mn.session.get_session()
    assert session.n_items == 1

    again = mn.Canvas()
    assert session.n_items == 1
    assert mol.name == "4ozs"  # the handle is still live, not LinkedObjectError
    assert again.scene == canvas.scene


def test_canvas_rerun_keeps_a_custom_compositor(canvas):
    "Binding must not rebuild the compositor over the top of one already set up."
    from nodebpy import compositor as c

    with canvas.compositor.reset() as (image, output):
        image >> c.Glare() >> output
    mn.Molecule.fetch("4ozs").add_style("cartoon")

    again = mn.Canvas()
    assert "CompositorNodeGlare" in [n.bl_idname for n in again.compositor.tree.nodes]


def test_canvas_with_explicit_template_always_loads(canvas):
    "Asking for a template is asking for the scene it describes."
    mn.Molecule.fetch("4ozs").add_style("cartoon")
    assert mn.session.get_session().n_items == 1

    mn.Canvas(template="Molecular Nodes")
    assert mn.session.get_session().n_items == 0


def test_canvas_template_none_binds_to_the_current_scene(canvas):
    mn.Molecule.fetch("4ozs").add_style("cartoon")
    again = mn.Canvas(template=None)
    assert mn.session.get_session().n_items == 1
    assert again.scene == canvas.scene


def test_canvas_ignores_dangling_entities_when_deciding(canvas):
    "A molecule deleted from the outliner is not work worth preserving."
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    bpy.data.objects.remove(mol.object, do_unlink=True)

    again = mn.Canvas()
    assert "Backdrop" in again.scene.objects
    assert mn.session.get_session().n_items == 0


def _camera_distance(canvas, points):
    return float(
        np.linalg.norm(np.asarray(points).mean(axis=0) - canvas.camera.camera.location)
    )


def test_look_at_accepts_any_number_of_points(canvas):
    """A target with more than 8 points used to silently empty the frame.

    `look_at` was annotated `list[tuple]` with no arity check, and anything but
    a bounding box put the subject a few pixels wide in the middle of nothing.
    """
    mol = mn.Molecule.fetch("4ozs").add_style("spheres", sphere_geometry="Mesh")
    points = np.asarray(mol.get_view())
    assert len(points) > 8

    canvas.look_at(points)
    offsets = points - np.asarray(canvas.camera.camera.location)
    depth = offsets @ canvas.camera.basis[2]
    assert (depth > 0).all(), "the subject ended up behind the camera"


def test_look_at_frames_what_is_drawn_not_the_whole_entity(canvas):
    "Styling one chain of four should frame that chain, not all four."
    mol = mn.Molecule.fetch("8H1B").add_style("cartoon", selection="chainID A")
    canvas.look_at(mol, viewpoint="front")
    on_chain = _camera_distance(canvas, mol.get_view("chainID A"))

    mol.add_style("cartoon")  # now the whole thing is drawn
    canvas.look_at(mol, viewpoint="front")
    on_everything = _camera_distance(canvas, mol.get_view("chainID A"))

    assert on_chain < on_everything


def test_get_view_is_not_stale_after_adding_a_style(canvas):
    """`bound_box` is evaluated geometry, and used to be read before the
    depsgraph had caught up - so the view depended on whether anything happened
    to have triggered an update."""
    mol = mn.Molecule.fetch("8H1B")
    mol.add_style("cartoon", selection="chainID A")

    before = np.asarray(mol.get_view())
    bpy.context.view_layer.update()
    after = np.asarray(mol.get_view())
    assert np.allclose(before.min(axis=0), after.min(axis=0))
    assert np.allclose(before.max(axis=0), after.max(axis=0))


def test_look_at_margin_moves_the_camera_back(canvas):
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    points = mol.get_view()

    distances = []
    for margin in (-0.1, 0.0, 0.2, 0.5):
        canvas.look_at(mol, viewpoint="front", margin=margin)
        distances.append(_camera_distance(canvas, points))
    assert distances == sorted(distances)


def test_look_at_frames_consistently_across_viewpoints(canvas):
    """Framing an axis-aligned box projected differently from each direction.

    Solving on the points themselves keeps the subject filling the frame from
    any angle - it can't be tight from one direction and loose from another.
    """
    mol = mn.Molecule.fetch("8H1B").add_style("cartoon")
    points = np.asarray(mol.get_view())

    filled = []
    for viewpoint in ("default", "front", "top", "left"):
        canvas.look_at(mol, viewpoint=viewpoint, margin=0.0)
        offsets = points - np.asarray(canvas.camera.camera.location)
        basis = canvas.camera.basis
        depth = offsets @ basis[2]
        left, right, bottom, top = canvas.camera.frame_bounds(canvas.scene)
        # the fraction of the frame the subject reaches on its widest axis
        filled.append(
            max(
                np.max(np.abs(offsets @ basis[0]) / depth) / max(right, -left),
                np.max(np.abs(offsets @ basis[1]) / depth) / max(top, -bottom),
            )
        )
    assert all(fill == pytest.approx(1.0, abs=1e-6) for fill in filled)


def test_look_at_rejects_targets_it_cannot_frame(canvas):
    with pytest.raises(ValueError, match=r"shape \(N, 3\)"):
        canvas.look_at([(1, 2), (3, 4)])
    with pytest.raises(ValueError, match=r"shape \(N, 3\)"):
        canvas.look_at([])


def test_look_at_extends_the_far_clip_to_reach_the_subject(canvas):
    "A subject beyond the far clip renders as nothing at all."
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    canvas.camera.clip_end = 0.1

    canvas.look_at(mol)
    points = np.asarray(mol.get_view())
    depth = (points - np.asarray(canvas.camera.camera.location)) @ canvas.camera.basis[
        2
    ]
    assert canvas.camera.clip_end >= depth.max()


@pytest.mark.parametrize(
    "style,kwargs",
    [
        ("cartoon", {}),
        ("spheres", {}),
        ("spheres", {"sphere_geometry": "Mesh"}),
        ("spheres", {"sphere_geometry": "Instance"}),
        ("ball_and_stick", {}),
        ("surface", {}),
    ],
)
def test_framing_covers_every_geometry_component(canvas, style, kwargs):
    """A style can render components an object never exposes through `data`.

    Spheres evaluate to an empty mesh with the real point cloud alongside it,
    and instanced styles put the geometry in an instances component - reading
    `obj.data` alone framed an empty scene.
    """
    mol = mn.Molecule.fetch("4ozs").add_style(style, **kwargs)
    points = np.asarray(mol.get_view())
    assert len(points) > 8

    # bounds Blender reports for the object cover every component it draws
    bounds = np.array([corner[:] for corner in mol.object.bound_box])
    extent = points.max(axis=0) - points.min(axis=0)
    # framing is allowed to be a little generous, never short
    assert (extent >= (bounds.max(axis=0) - bounds.min(axis=0)) - 1e-4).all()


def test_spheres_are_framed_by_their_surface_not_their_centres(canvas):
    "Framing the centres cuts the outermost spheres in half at the frame edge."
    mol = mn.Molecule.fetch("4ozs").add_style("spheres")
    points = np.asarray(mol.get_view())
    centres = mol._world_positions(mol.atoms)

    extent = points.max(axis=0) - points.min(axis=0)
    centre_extent = centres.max(axis=0) - centres.min(axis=0)
    assert (extent > centre_extent).all()


def test_look_at_leaves_breathing_room_by_default(canvas):
    "The default margin keeps the subject off the edge of the frame."
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    points = np.asarray(mol.get_view())

    canvas.look_at(mol, viewpoint="front")
    offsets = points - np.asarray(canvas.camera.camera.location)
    basis = canvas.camera.basis
    depth = offsets @ basis[2]
    left, right, bottom, top = canvas.camera.frame_bounds(canvas.scene)
    filled = max(
        np.max(np.abs(offsets @ basis[0]) / depth) / max(right, -left),
        np.max(np.abs(offsets @ basis[1]) / depth) / max(top, -bottom),
    )
    assert filled == pytest.approx(0.95, abs=1e-6)
