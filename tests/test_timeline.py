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


def keyframes(channel):
    """The keyframe points of a channel's F-curve, which must exist."""
    fcurve = channel.fcurve()
    assert fcurve is not None
    return fcurve.keyframe_points


def update_view_layer():
    view_layer = bpy.context.view_layer
    assert view_layer is not None
    view_layer.update()


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


def test_a_rejected_clip_writes_nothing(canvas, cartoon_mol):
    _, cartoon = cartoon_mol
    channel = Channel(cartoon.i.loop_radius.socket, "default_value")
    t = canvas.timeline(fps=10)
    with pytest.raises(ValueError, match="already animates"):
        t.play(
            t.tween(cartoon.i.loop_radius, 1.0),
            t.tween(cartoon.i.loop_radius, 2.0),
        )
    # only the first tween's two keys, at the value it asked for
    points = keyframes(channel)
    assert [p.co.x for p in points] == [1, 11]
    assert points[1].co.y == pytest.approx(1.0)


def test_clips_on_one_property_play_in_time_order(canvas, cartoon_mol):
    _, cartoon = cartoon_mol
    channel = Channel(cartoon.i.loop_radius.socket, "default_value")
    t = canvas.timeline(fps=10)
    t.play(t.tween(cartoon.i.loop_radius, 1.0), at=2.0)
    # slotting a clip in before existing keys would break the hold between them
    with pytest.raises(ValueError, match="time order"):
        t.play(t.tween(cartoon.i.loop_radius, 3.0), at=0.0)
    assert [p.co.x for p in keyframes(channel)] == [21, 31]


# -- camera -------------------------------------------------------------------


def camera_pose(canvas, frame):
    """World location and rotation of the camera at a frame."""
    canvas.frame = frame
    update_view_layer()
    matrix = canvas.camera.camera.matrix_world
    return np.array(matrix.translation), np.array(matrix.to_euler("XYZ"))


def evaluated_matrix(canvas, frame=None):
    """The camera's world matrix as Blender evaluates it, constraint included."""
    if frame is not None:
        canvas.frame = frame
    update_view_layer()
    return canvas.camera.camera.matrix_world.copy()


def test_only_timeline_camera_moves_rig_the_camera(canvas, cartoon_mol):
    mol, cartoon = cartoon_mol
    cam = canvas.camera
    # an immediate look_at and a timeline that leaves the camera alone do not rig
    canvas.look_at(mol, viewpoint="left")
    with canvas.timeline(fps=10) as t:
        t.play(t.tween(cartoon.i.loop_radius, 1.0))
    assert cam.camera.parent is None and not cam.camera.constraints
    assert cam.pivot is None and cam.target is None
    before = evaluated_matrix(canvas)
    with canvas.timeline(fps=10, start=11) as t:
        t.play(cam.orbit(30), run_time=1)
    # one empty parent to orbit about, one to look at and one to focus on
    pivot, target, focus = cam.pivot, cam.target, cam.focus_point
    assert cam.camera.parent is pivot and pivot.type == "EMPTY"
    assert pivot.parent is None
    (constraint,) = cam.camera.constraints
    assert constraint.type == "DAMPED_TRACK" and constraint.target is target
    assert target.type == "EMPTY" and target.hide_render
    assert focus.type == "EMPTY" and focus.hide_render and focus is not target
    assert cam.camera_data.dof.focus_object is focus
    # both sit on the view axis, so rigging did not move the view
    assert np.allclose(evaluated_matrix(canvas, 11), before, atol=1e-5)
    assert np.allclose(pivot.location, target.location)
    assert np.allclose(focus.location, target.location)
    forward = np.array(target.location) - np.array(cam.location)
    assert np.dot(forward / np.linalg.norm(forward), cam.basis[2]) > 1 - 1e-6


def test_a_parented_camera_is_not_rigged(canvas):
    parent = bpy.data.objects.new("rig_parent", None)
    canvas.scene.collection.objects.link(parent)
    canvas.camera.camera.parent = parent
    t = canvas.timeline(fps=10)
    with pytest.raises(ValueError, match="parented"):
        t.play(canvas.camera.dolly(1.0))
    canvas.camera.camera.parent = None


def test_matrix_world_follows_the_look_at_constraint(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    cam = canvas.camera
    with canvas.timeline(fps=10) as t:
        t.play(cam.dolly(0.5))
    # the target moves off the pivot, so the constraint turns the camera
    cam.target.location.x += 1.0
    assert not np.allclose(cam.pivot.location, cam.target.location)
    assert np.allclose(cam.matrix_world, evaluated_matrix(canvas), atol=1e-5)


def test_look_at_on_a_rigged_camera(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    cam = canvas.camera
    with canvas.timeline(fps=10) as t:
        t.play(cam.dolly(0.5))
    view = mol.get_view("resid 1-10")
    canvas.look_at(view)
    # pivot and target together at the subject's depth on the view axis
    centre = np.array(mn.scene.timeline.target_centre(view))
    offset = centre - np.array(cam.location)
    assert np.allclose(cam.pivot.location, cam.target.location)
    depth = np.linalg.norm(np.array(cam.target.location) - np.array(cam.location))
    assert depth == pytest.approx(np.dot(offset, cam.basis[2]))
    # the rotation lives on the pivot and the camera looks down its -Z
    canvas.look_at(mol, viewpoint=(90, 0, -45))
    assert np.allclose(np.degrees(cam.pivot.rotation_euler), (90, 0, -45), atol=1e-4)
    assert np.allclose(cam.camera.rotation_euler, 0.0)
    assert cam.rotation == pytest.approx((90, 0, -45), abs=1e-3)
    assert np.allclose(cam.matrix_world, evaluated_matrix(canvas), atol=1e-5)


def test_clear_removes_the_rig_and_keeps_the_view(canvas, cartoon_mol):
    cam = canvas.camera
    with canvas.timeline(fps=10) as t:
        t.play(cam.orbit(90), cam.zoom(80), run_time=1)
    pose = evaluated_matrix(canvas, 6)
    names = [cam.pivot.name, cam.target.name, cam.focus_point.name]
    canvas.clear()
    assert cam.camera.parent is None and not cam.camera.constraints
    assert not any(name in bpy.data.objects for name in names)
    # the keys were relative to the rig, so they go with it
    assert cam.camera.animation_data is None or cam.camera.animation_data.action is None
    assert np.allclose(evaluated_matrix(canvas, 20), pose, atol=1e-5)


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
    assert np.allclose(end_location, canvas.camera.location, atol=1e-5)
    assert np.allclose(np.degrees(end_rotation), canvas.camera.rotation, atol=1e-4)


def test_orbit_keeps_its_distance_from_the_pivot(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    with canvas.timeline(fps=10) as t:
        t.play(canvas.camera.orbit(180), run_time=2, easing="linear")
    pivot = canvas.camera.pivot.location.copy()
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


def test_turntable_orbit_is_two_keys_on_the_pivot(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    with canvas.timeline(fps=10) as t:
        t.play(canvas.camera.orbit(360, about=mol), run_time=1, easing="linear")
    pivot = canvas.camera.pivot
    points = keyframes(Channel(pivot, "rotation_euler", 2))
    assert [p.co.x for p in points] == [1, 11]
    angles = [p.co.y for p in points]
    assert angles[1] - angles[0] == pytest.approx(2 * math.pi)
    # the other pivot channels and the camera itself are untouched
    assert Channel(pivot, "rotation_euler", 0).fcurve() is None
    assert Channel(canvas.camera.camera, "location", 0).fcurve() is None
    # the camera's right axis is the pivot's X, so a tumble is two keys too
    with canvas.timeline(fps=10, start=11) as t:
        t.play(canvas.camera.orbit(30, axis="right"), run_time=1)
    assert len(keyframes(Channel(pivot, "rotation_euler", 0))) == 2


def test_orbit_about_another_axis_is_sampled_and_continuous(canvas, cartoon_mol):
    with canvas.timeline(fps=10) as t:
        t.play(
            canvas.camera.orbit(360, axis=(1, 1, 0)),
            run_time=1,
            easing="linear",
        )
    pivot = canvas.camera.pivot
    start, _ = camera_pose(canvas, 1)
    centre = np.array(pivot.location)
    radius = np.linalg.norm(start - centre)
    for index in range(3):
        points = keyframes(Channel(pivot, "rotation_euler", index))
        assert len(points) == 11
        # no wrap-around jump in the euler curve over a full turn
        steps = np.abs(np.diff([p.co.y for p in points]))
        assert steps.max() < math.pi
    for frame in (1, 4, 7, 11):
        location, _ = camera_pose(canvas, frame)
        assert np.linalg.norm(location - centre) == pytest.approx(radius, rel=1e-5)
    end, _ = camera_pose(canvas, 11)
    assert np.allclose(end, start, atol=1e-4)


def test_dolly_during_an_orbit_does_not_conflict(canvas, cartoon_mol):
    with canvas.timeline(fps=10) as t:
        t.play(canvas.camera.orbit(90), canvas.camera.dolly(1.0), run_time=1)
    centre = np.array(canvas.camera.pivot.location)
    start, _ = camera_pose(canvas, 1)
    radius = np.linalg.norm(start - centre)
    end, _ = camera_pose(canvas, 11)
    # the orbit turned the pivot and the dolly pulled the camera in on it
    assert np.linalg.norm(end - centre) == pytest.approx(radius - 1.0, abs=1e-4)
    assert not np.allclose(end[:2], start[:2], atol=0.1)


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
    empty = canvas.camera.dof.focus_object
    assert empty is not None
    canvas.frame = 51
    centre = mn.scene.timeline.target_centre(mol)
    assert np.allclose(empty.location.to_tuple(), centre.to_tuple(), atol=1e-6)


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
    assert np.allclose(empty.location.to_tuple(), centre.to_tuple(), atol=1e-6)
    assert data.dof.aperture_fstop == pytest.approx(1.4)
    # pulling focus again re-uses the same empty
    with canvas.timeline(fps=10, start=11) as t:
        t.play(canvas.camera.focus(mol, fstop=4))
    assert data.dof.focus_object is empty


def test_look_at_after_a_full_orbit_does_not_spin(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    with canvas.timeline(fps=24) as t:
        t.play(canvas.camera.orbit(360, easing="linear"), run_time=1)
        t.play(canvas.camera.look_at(mol, viewpoint="front"), run_time=1)
    # the orbit ends a whole turn round, which is the front view already
    z = Channel(canvas.camera.pivot, "rotation_euler", 2)
    assert z.value_at(49) == pytest.approx(z.value_at(25), abs=1e-5)


def test_orbit_about_moves_the_pivot_and_target_onto_the_subject(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    view = mol.get_view("resid 1-10")
    with canvas.timeline(fps=10) as t:
        t.play(canvas.camera.orbit(90, about=view), run_time=1)
    canvas.frame = 11
    centre = np.array(mn.scene.timeline.target_centre(view))
    cam = canvas.camera
    assert np.allclose(cam.pivot.location, centre, atol=1e-6)
    assert np.allclose(cam.target.location, centre, atol=1e-6)
    location, _ = camera_pose(canvas, 11)
    to_centre = (centre - location) / np.linalg.norm(centre - location)
    assert np.dot(cam.basis[2], to_centre) == pytest.approx(1.0)


def test_focus_leaves_the_frames_before_it_alone(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    dof = canvas.camera.camera_data.dof
    assert not dof.use_dof
    with canvas.timeline(fps=10) as t:
        t.wait(2)
        t.play(canvas.camera.focus(mol, fstop=1.4), run_time=1)
    canvas.frame = 1
    assert not dof.use_dof
    canvas.frame = 21
    # switched on as the pull starts, opening up from a sharp aperture
    assert dof.use_dof
    assert dof.aperture_fstop == pytest.approx(mn.scene.timeline._SHARP_FSTOP)
    canvas.frame = 31
    assert dof.aperture_fstop == pytest.approx(1.4)


def test_a_focus_pull_holds_the_framing(canvas, cartoon_mol):
    mol, _ = cartoon_mol
    cam = canvas.camera
    helix, far = mol.get_view("resid 5-19"), mol.get_view("resid 100-120")
    with canvas.timeline(fps=10) as t:
        # framing and focus share a play: they key different empties
        t.play(cam.look_at(helix, margin=0.3), cam.focus(helix), run_time=1)
        t.play(cam.focus(far, fstop=1.4), run_time=1)
    framed = evaluated_matrix(canvas, 11)
    assert np.allclose(evaluated_matrix(canvas, 21), framed, atol=1e-6)
    centre = np.array(mn.scene.timeline.target_centre(far))
    assert np.allclose(cam.focus_point.location, centre, atol=1e-6)


# -- trajectory playback --------------------------------------------------------


def test_styles_are_reachable_by_index_and_name(canvas, traj):
    traj.add_style("cartoon", name="thin")
    assert len(traj.styles) == 2
    assert [style.name for style in traj.styles] == ["Style Ribbon", "Style Cartoon"]
    assert traj.styles[1].name == traj.styles["Style Cartoon"].name == "Style Cartoon"
    # a label set through add_style(name=) works too
    assert traj.styles["thin"].name == "Style Cartoon"
    with pytest.raises(KeyError, match="no style 'Style Surface'"):
        traj.styles["Style Surface"]
    with pytest.raises(IndexError):
        traj.styles[2]
    # the wrapped node's inputs animate from the simple API
    cartoon = traj.styles["Style Cartoon"]
    with canvas.timeline(fps=10) as t:
        t.play(t.tween(cartoon.i.loop_radius, 1.5), run_time=1)
    socket = cartoon.i.loop_radius.socket
    assert value_at(canvas, socket, 11) == pytest.approx(1.5)


def test_frames_clip_plays_the_trajectory(canvas, traj):
    assert traj.update_with_scene
    with canvas.timeline(fps=10) as t:
        t.play(traj.play(0, 4), run_time=2)
    assert not traj.update_with_scene
    canvas.frame = 1
    assert traj.universe.trajectory.frame == 0
    canvas.frame = 11
    assert traj.universe.trajectory.frame == 2
    canvas.frame = 21
    assert traj.universe.trajectory.frame == 4
    # held on the last frame afterwards
    canvas.frame = 40
    assert traj.universe.trajectory.frame == 4


def test_frames_clip_ignores_subframes(canvas, traj):
    # the manual frame is a trajectory frame, so subframes do not scale it
    traj.subframes = 1
    with canvas.timeline(fps=10) as t:
        t.play(traj.play(0, 4), run_time=1)
    canvas.frame = 11
    assert traj.frame == 4
    assert traj.universe.trajectory.frame == 4


# -- style clips -----------------------------------------------------------------


@pytest.fixture
def two_styles(canvas):
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    mol.add_style("spheres", color="#ff0000").add_style("cartoon")
    return mol


def point_count(mol):
    from molecularnodes.blender import utils as blender_utils

    obj = mol.object
    points = blender_utils.evaluated_points(obj)
    # an object that renders nothing falls back to its 8 bounding box corners
    return 0 if len(points) == 8 else len(points)


def mean_color(mol):
    obj = mol.object.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh = obj.to_mesh()
    try:
        attribute = mesh.attributes["Color"]
        values = np.zeros(len(attribute.data) * 4)
        attribute.data.foreach_get("color", values)
        return values.reshape(-1, 4).mean(axis=0)
    finally:
        obj.to_mesh_clear()


def fade_nodes(mol):
    return [n for n in mol.modifier_node_tree.nodes if n.label == "Timeline Fade"]


def test_fade_adds_one_node_per_style_and_reuses_it(canvas, two_styles):
    mol = two_styles
    before = point_count(mol)
    with canvas.timeline(fps=10) as t:
        clip = t.fade(mol.styles[0], 0.0)
        # splicing the node in changes nothing until it is keyed
        assert point_count(mol) == before
        t.play(clip, run_time=1)
        t.play(t.fade(mol.styles[0], 1.0), run_time=1)
    assert len(fade_nodes(mol)) == 1
    assert "fade Style Spheres -> 0" in repr(t)
    fade = fade_nodes(mol)[0].inputs["Fade"]
    assert value_at(canvas, fade, 1) == pytest.approx(1.0)
    assert value_at(canvas, fade, 11) == pytest.approx(0.0)
    assert point_count(mol) < before  # the spheres are gone, the cartoon stays
    assert value_at(canvas, fade, 21) == pytest.approx(1.0)
    assert point_count(mol) == before


def test_fade_an_entity_fades_every_style(canvas, two_styles):
    mol = two_styles
    with canvas.timeline(fps=10) as t:
        t.play(t.fade(mol, 0.0), run_time=1)
    assert len(fade_nodes(mol)) == 2
    canvas.frame = 11
    assert point_count(mol) == 0


def test_set_visible_cuts_and_shares_the_fade_node(canvas, two_styles):
    mol = two_styles
    cartoon = mol.styles["Style Cartoon"]
    with canvas.timeline(fps=10) as t:
        t.wait(1)
        t.set_visible(cartoon, False)
        t.play(t.fade(cartoon, 1.0), run_time=1)
    assert len(fade_nodes(mol)) == 1
    fade = fade_nodes(mol)[0].inputs["Fade"]
    assert value_at(canvas, fade, 10) == pytest.approx(1.0)
    assert value_at(canvas, fade, 11) == pytest.approx(0.0)
    assert value_at(canvas, fade, 21) == pytest.approx(1.0)
    assert "hide Style Cartoon" in repr(t)


def test_fade_rejects_things_that_are_not_styles(canvas):
    t = canvas.timeline()
    with pytest.raises(TypeError, match="style"):
        t.fade(42, 0.0)


def test_recolor_blends_into_the_new_color(canvas, two_styles):
    mol = two_styles
    # the cartoon is the mesh; the spheres are a point cloud beside it
    cartoon = mol.styles["Style Cartoon"]
    blue = np.array([0.0, 0.0, 1.0, 1.0])
    with canvas.timeline(fps=10) as t:
        start = mean_color(mol)
        t.play(t.recolor(cartoon, "#0000ff"), run_time=1)
    distances = []
    for frame in (1, 6, 11):
        canvas.frame = frame
        distances.append(np.linalg.norm(mean_color(mol) - blue))
    canvas.frame = 1
    assert np.allclose(mean_color(mol), start, atol=1e-4)
    assert distances[0] > distances[1] > distances[2]
    assert distances[2] == pytest.approx(0.0, abs=1e-3)
    assert "recolor Style Cartoon -> #0000ff" in repr(t)


def test_recolor_stacks_and_takes_a_selection(canvas, two_styles):
    mol = two_styles
    spheres = mol.styles["Style Spheres"]
    with canvas.timeline(fps=10) as t:
        # "chain A" is rewritten to "chainID A" until MDAnalysis 2.12
        t.play(t.recolor(spheres, (0, 0, 1), selection="chain A"), run_time=1)
        t.play(t.recolor(spheres, lambda: mg.ColorResidueClass()), run_time=1)
    nodes = mol.modifier_node_tree.nodes
    recolors = [n for n in nodes if n.label == "Timeline Recolor"]
    assert len(recolors) == 2
    # the later recolour sits directly in front of the style
    style = nodes["Style Spheres"]
    assert style.inputs["Atoms"].links[0].from_node == recolors[1]
    assert recolors[1].inputs["Atoms"].links[0].from_node == recolors[0]
    assert recolors[0].inputs["Selection"].is_linked
    assert not recolors[1].inputs["Selection"].is_linked
    with pytest.raises(ValueError, match="selection"):
        t.recolor(spheres, "red", selection="not a selection (")


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
