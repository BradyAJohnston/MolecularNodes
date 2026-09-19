import math
import bpy
import numpy as np
import pytest
import molecularnodes as mn
from molecularnodes.nodes import geometry as mg
from molecularnodes.scene.timeline import (
    Channel,
    Easing,
    Orbit,
    Tween,
    resolve_channels,
)
from .constants import data_dir


@pytest.fixture
def canvas():
    return mn.Canvas()


@pytest.fixture
def cartoon_mol(canvas):
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    with mol.tree as tree:
        cartoon = mg.StyleCartoon()
        tree.atoms >> cartoon >> tree.join
    canvas.look_at(mol, viewpoint="front")
    return mol, cartoon


@pytest.fixture
def traj(canvas):
    mol = mn.Molecule.load(
        data_dir / "md_ppr/box.gro", data_dir / "md_ppr/first_5_frames.xtc"
    )
    mol.add_style("ribbon")
    return mol


def value_at(canvas, socket, frame):
    canvas.frame = frame
    return socket.default_value


# -- easing -----------------------------------------------------------------


@pytest.mark.parametrize("name", ["linear", "smooth", "sine", "ease_in", "ease_out"])
def test_easing_functions_span_zero_to_one(name):
    rate = Easing.function(name)
    assert rate(0.0) == pytest.approx(0.0)
    assert rate(1.0) == pytest.approx(1.0)
    samples = [rate(t) for t in np.linspace(0, 1, 11)]
    assert samples == sorted(samples)


def test_easing_rejects_unknown_names():
    with pytest.raises(ValueError, match="Unknown easing"):
        Easing.validate("bouncy")


# -- time -------------------------------------------------------------------


def test_cursor_and_frames(canvas):
    t = canvas.timeline(fps=24, start=1)
    assert t.now == 0 and t.frame == 1
    t.wait(1.5)
    assert t.frame == 37
    with pytest.raises(ValueError):
        t.wait(-1)
    with pytest.raises(ValueError):
        t.play()


def test_context_fits_the_frame_range(canvas, cartoon_mol):
    _, cartoon = cartoon_mol
    canvas.frame = 10
    with canvas.timeline(fps=10, start=1) as t:
        t.play(t.tween(cartoon.i.loop_radius, 1.0), run_time=2)
        t.wait(1)
    assert canvas.frame_range == (1, 31)
    # rewound to the start so the scene shows the opening state
    assert canvas.frame == 1


# -- tweens on node sockets -------------------------------------------------


def test_tween_keys_a_node_input(canvas, cartoon_mol):
    _, cartoon = cartoon_mol
    socket = cartoon.i.loop_radius.socket
    socket.default_value = 0.2
    with canvas.timeline(fps=10, start=1) as t:
        t.play(t.tween(cartoon.i.loop_radius, 1.2), run_time=1, easing="linear")
    assert value_at(canvas, socket, 1) == pytest.approx(0.2)
    assert value_at(canvas, socket, 6) == pytest.approx(0.7)
    assert value_at(canvas, socket, 11) == pytest.approx(1.2)
    # held either side of the clip
    assert value_at(canvas, socket, 30) == pytest.approx(1.2)


def test_tween_continues_from_the_previous_clip(canvas, cartoon_mol):
    _, cartoon = cartoon_mol
    socket = cartoon.i.loop_radius.socket
    socket.default_value = 0.0
    with canvas.timeline(fps=10, start=1) as t:
        t.play(t.tween(cartoon.i.loop_radius, 1.0), run_time=1, easing="linear")
        t.wait(1)
        t.play(t.tween(cartoon.i.loop_radius, 0.0), run_time=1, easing="linear")
    # the value holds through the wait ...
    assert value_at(canvas, socket, 16) == pytest.approx(1.0)
    # ... and the second tween starts from where the first left off
    assert value_at(canvas, socket, 26) == pytest.approx(0.5)
    assert value_at(canvas, socket, 31) == pytest.approx(0.0)


def test_easing_is_written_onto_the_keyframes(canvas, cartoon_mol):
    _, cartoon = cartoon_mol
    channel = Channel(cartoon.i.loop_radius.socket, "default_value")
    with canvas.timeline(fps=10) as t:
        t.play(t.tween(cartoon.i.loop_radius, 1.0), easing="sine")
    fcurve = channel.fcurve()
    assert fcurve is not None
    first, last = fcurve.keyframe_points[0], fcurve.keyframe_points[-1]
    assert first.interpolation == "SINE" and first.easing == "EASE_IN_OUT"
    assert (first.co.x, last.co.x) == (1, 11)


def test_smooth_easing_is_flat_at_both_ends(canvas, cartoon_mol):
    _, cartoon = cartoon_mol
    socket = cartoon.i.loop_radius.socket
    socket.default_value = 0.0
    with canvas.timeline(fps=20) as t:
        t.play(t.tween(cartoon.i.loop_radius, 1.0), run_time=1, easing="smooth")
    early = value_at(canvas, socket, 2) - value_at(canvas, socket, 1)
    middle = value_at(canvas, socket, 11) - value_at(canvas, socket, 10)
    assert early < middle / 4


def test_set_steps_without_advancing(canvas, cartoon_mol):
    _, cartoon = cartoon_mol
    socket = cartoon.i.loop_radius.socket
    socket.default_value = 0.2
    with canvas.timeline(fps=10) as t:
        t.wait(1)
        t.set(cartoon.i.loop_radius, 2.0)
        assert t.frame == 11
    assert value_at(canvas, socket, 10) == pytest.approx(0.2)
    assert value_at(canvas, socket, 11) == pytest.approx(2.0)


def test_vector_values_key_every_component(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    pairs = resolve_channels((mol.object, "location"), (1.0, 2.0, 3.0))
    assert [c.index for c, _ in pairs] == [0, 1, 2]
    with pytest.raises(ValueError, match="components"):
        resolve_channels((mol.object, "location"), (1.0, 2.0))
    with pytest.raises(ValueError, match="scalar"):
        resolve_channels((mol.object, "pass_index"), (1.0, 2.0))


def test_entity_properties_are_channels(canvas, traj):
    with canvas.timeline(fps=10) as t:
        t.play(t.tween((traj, "subframes"), 4), run_time=1, easing="linear")
    canvas.frame = 11
    assert traj.subframes == 4


def test_annotation_parameters_are_channels(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    label = mol.annotations.add_label_2d(text="hello", location=(0.5, 0.5))
    label.text_size = 10
    with canvas.timeline(fps=10) as t:
        t.play(t.tween((label, "text_size"), 30), run_time=1, easing="linear")
    canvas.frame = 6
    assert label.text_size == 20
    canvas.frame = 11
    assert label.text_size == 30
    with pytest.raises(ValueError, match="no input"):
        resolve_channels((label, "nonsense"), 1.0)


def test_unsupported_targets_are_rejected(canvas):
    with pytest.raises(TypeError):
        resolve_channels("lens", 1.0)
    with pytest.raises(TypeError):
        resolve_channels(("lens", "x"), 1.0)


# -- conflicts ----------------------------------------------------------------


def test_two_clips_on_one_channel_conflict(canvas, cartoon_mol):
    _, cartoon = cartoon_mol
    t = canvas.timeline(fps=10)
    with pytest.raises(ValueError, match="already animates"):
        t.play(
            t.tween(cartoon.i.loop_radius, 1.0),
            t.tween(cartoon.i.loop_radius, 2.0),
        )


def test_clips_back_to_back_do_not_conflict(canvas, cartoon_mol):
    _, cartoon = cartoon_mol
    with canvas.timeline(fps=10) as t:
        t.play(t.tween(cartoon.i.loop_radius, 1.0))
        t.play(t.tween(cartoon.i.loop_radius, 2.0))
        # layered at an explicit time that is free
        t.play(t.tween(cartoon.i.loop_radius, 3.0), at=5.0)
    assert len(t.entries) == 3


# -- camera -------------------------------------------------------------------


def camera_pose(canvas, frame):
    canvas.frame = frame
    bpy.context.view_layer.update()
    cam = canvas.camera.camera
    return np.array(cam.matrix_world.translation), np.array(cam.rotation_euler)


def test_look_at_clip_lands_on_the_look_at_pose(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    canvas.look_at(mol, viewpoint="front")
    start_location, _ = camera_pose(canvas, 1)
    with canvas.timeline(fps=10) as t:
        t.play(canvas.camera.look_at(mol, viewpoint="top"), run_time=1)
    end_location, end_rotation = camera_pose(canvas, 11)
    # the clip did not move the camera when it was created, only when played
    assert not np.allclose(start_location, end_location)
    # the destination is exactly what the immediate look_at gives
    canvas.look_at(mol, viewpoint="top")
    assert np.allclose(end_location, canvas.camera.camera.location, atol=1e-5)
    assert np.allclose(end_rotation, canvas.camera.camera.rotation_euler, atol=1e-6)


def test_orbit_keeps_its_distance_from_the_pivot(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    with canvas.timeline(fps=10) as t:
        t.play(canvas.camera.orbit(180, about=mol), run_time=2, easing="linear")
    pivot = mn.scene.timeline.target_centre(mol)
    start, _ = camera_pose(canvas, 1)
    radius = np.linalg.norm(start - pivot)
    for frame in range(1, 22, 5):
        location, _ = camera_pose(canvas, frame)
        assert np.linalg.norm(location - pivot) == pytest.approx(radius, rel=1e-5)
    end, _ = camera_pose(canvas, 21)
    # half a turn about z puts the camera on the far side of the pivot
    assert np.allclose(end[:2], 2 * np.array(pivot[:2]) - start[:2], atol=1e-5)
    assert end[2] == pytest.approx(start[2])
    # the camera still looks at the pivot
    forward = canvas.camera.basis[2]
    to_pivot = (np.array(pivot) - end) / np.linalg.norm(np.array(pivot) - end)
    assert np.dot(forward, to_pivot) > 0.99


def test_orbit_is_keyed_every_frame_and_continuous(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    with canvas.timeline(fps=10) as t:
        t.play(canvas.camera.orbit(360, about=mol), run_time=1, easing="linear")
    fcurve = Channel(canvas.camera.camera, "rotation_euler", 2).fcurve()
    assert len(fcurve.keyframe_points) == 11
    angles = [p.co.y for p in fcurve.keyframe_points]
    steps = np.diff(angles)
    # no wrap-around jump in the euler curve over a full turn
    assert np.allclose(steps, steps[0], atol=1e-6)
    assert abs(angles[-1] - angles[0]) == pytest.approx(2 * math.pi, abs=1e-5)


def test_clips_start_from_where_the_previous_clip_left_the_camera(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    with canvas.timeline(fps=10) as t:
        t.play(canvas.camera.look_at(mol, viewpoint="top"), run_time=1)
        t.wait(1)
        # the camera is on frame 1 of the scene while this is authored, but the
        # orbit must start from the top view the look_at ended on
        t.play(canvas.camera.orbit(90, about=mol), run_time=1)
    before, _ = camera_pose(canvas, 11)
    start, _ = camera_pose(canvas, 21)
    assert np.allclose(before, start, atol=1e-5)
    # and the focus pull reads the camera as it is at its own start frame
    with canvas.timeline(fps=10, start=31) as t:
        t.play(canvas.camera.dolly(1.0), run_time=1)
        t.play(canvas.camera.focus(mol), run_time=1)
    empty = canvas.camera.camera_data.dof.focus_object
    canvas.frame = 51
    assert np.allclose(empty.location, mn.scene.timeline.target_centre(mol), atol=1e-6)


def test_orbit_axis_options(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    for axis in ("x", "up", (0, 1, 0)):
        clip = canvas.camera.orbit(90, axis=axis, about=mol)
        clip.prepare()
    with pytest.raises(ValueError, match="orbit axis"):
        Orbit(canvas.camera, 90, axis="diagonal", about=mol).prepare()


def test_dolly_and_zoom(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    start, _ = camera_pose(canvas, 1)
    forward = canvas.camera.basis[2].copy()
    with canvas.timeline(fps=10) as t:
        t.play(canvas.camera.dolly(2.0), canvas.camera.zoom(85), run_time=1)
    end, _ = camera_pose(canvas, 11)
    assert np.allclose(end - start, forward * 2.0, atol=1e-5)
    assert canvas.camera.lens == pytest.approx(85)


def test_move_to(canvas):
    with canvas.timeline(fps=10) as t:
        t.play(canvas.camera.move_to(location=(1, 2, 3), rotation=(90, 0, 45)))
    location, rotation = camera_pose(canvas, 11)
    assert np.allclose(location, (1, 2, 3), atol=1e-6)
    assert np.allclose(np.degrees(rotation), (90, 0, 45), atol=1e-4)
    with pytest.raises(ValueError):
        canvas.camera.move_to()


def test_focus_uses_a_tracked_empty(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    data = canvas.camera.camera_data
    assert not data.dof.use_dof
    view = mol.get_view("resid 1-10")
    with canvas.timeline(fps=10) as t:
        t.play(canvas.camera.focus(view, fstop=1.4), run_time=1)
    canvas.frame = 11
    assert data.dof.use_dof
    empty = data.dof.focus_object
    assert empty is not None and empty.hide_render
    centre = mn.scene.timeline.target_centre(view)
    assert np.allclose(empty.location, centre, atol=1e-6)
    assert data.dof.aperture_fstop == pytest.approx(1.4)
    # pulling focus again re-uses the same empty
    with canvas.timeline(fps=10, start=11) as t:
        t.play(canvas.camera.focus(mol, fstop=4))
    assert data.dof.focus_object is empty


# -- trajectory playback --------------------------------------------------------


def test_frames_clip_plays_the_trajectory(canvas, traj):
    assert traj.update_with_scene
    with canvas.timeline(fps=10) as t:
        t.play(t.frames(traj, 0, 4), run_time=2)
    assert not traj.update_with_scene
    canvas.frame = 1
    assert traj.uframe == 0
    canvas.frame = 11
    assert traj.uframe == 2
    canvas.frame = 21
    assert traj.uframe == 4
    # held on the last frame afterwards
    canvas.frame = 40
    assert traj.uframe == 4


def test_frames_clip_honours_subframes(canvas, traj):
    traj.subframes = 1
    with canvas.timeline(fps=10) as t:
        t.play(t.frames(traj, 0, 4), run_time=1)
    canvas.frame = 11
    assert traj.frame == 8
    assert traj.uframe == 4


# -- storyboard ------------------------------------------------------------------


def test_repr_reads_as_a_storyboard(canvas, cartoon_mol):
    mol, cartoon = cartoon_mol
    with canvas.timeline(fps=10) as t:
        t.play(canvas.camera.orbit(90, about=mol), run_time=1)
        t.play(t.tween(cartoon.i.loop_radius, 1.0), run_time=0.5)
    text = repr(t)
    assert "camera.orbit 90deg" in text
    assert "tween" in text
    assert "frames 1-16" in text
    assert repr(Tween(cartoon.i.loop_radius, 1.0)).startswith("<tween")


def test_render_a_storyboard(tmp_path, cartoon_mol):
    canvas = mn.Canvas(resolution=(32, 32))
    canvas.engine = mn.scene.Cycles(samples=1, device="CPU")
    canvas.compositor.device = "CPU"
    mol, cartoon = cartoon_mol
    with canvas.timeline(fps=2) as t:
        t.play(canvas.camera.orbit(45, about=mol), run_time=1)
    out = tmp_path / "story.mp4"
    t.render(out)
    assert out.read_bytes()[4:8] == b"ftyp"
    assert canvas.frame_range == (1, 3)
