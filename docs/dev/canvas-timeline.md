# A timeline layer on `Canvas`

Design notes for the prototype in `molecularnodes/scene/timeline.py`, following
item 5.1 of `protein-motion-comparison.md`. Written 2026-09-19 against Blender 5.2.

## 1. What problem this solves

A Molecular Nodes render script today has two ways to make a movie:

- `canvas.animation()` plays back whatever the scene already animates: trajectory
  frames through the entity handlers, node trees that read Scene Time, and any
  keyframes placed by hand.
- `canvas.record()` hands the loop to Python: mutate the scene, render a frame,
  repeat.

Neither lets a script *say* "over two seconds, swing the camera round and
thicken the ribbon, then pull focus onto the active site". The first has no
authoring API for camera and parameter motion; the second throws away
Blender's timeline, so the result cannot be scrubbed, edited or rendered with
motion blur, and easing is the author's problem.

ProteinMotion's appeal is exactly that a script reads as a storyboard. The
comparison doc's verdict was that the gap is an API layer, not capability, and
that Blender's animation system should remain the engine.

## 2. Design principles

1. **Compile to keyframes; never run at render time.** A clip only inserts
   keyframes on Blender properties. There is no Python interpolation during
   playback, no frame-change handler owned by the timeline, and no snapshot
   and replay machinery. Seeking is exact because a property has one F-curve,
   and everything the timeline produces is visible and editable in the Graph
   Editor. ProteinMotion needs a compiled timeline because it has no Blender;
   we do not.

2. **Geometry already animates itself.** Trajectory positions follow the scene
   frame, DSSP and selections update per frame, and node trees read Scene
   Time. The timeline therefore does not touch geometry. Its channels are the
   things that today need hand-placed keys: the camera, node inputs, entity
   playback properties, annotation parameters, material and world values.

3. **Easing lives on the keys.** A tween writes two keys and sets the
   interpolation on the first (Blender attaches the segment's interpolation to
   its left key). Only motion that has no two-key form is sampled per frame:
   an orbit, whose path is an arc. Even then the samples are ordinary linear
   keys, so a scrub between frames is still exact.

4. **One channel, one owner.** Two clips writing the same property over
   overlapping frames raise at authoring time. Blender would silently let the
   later `keyframe_insert` win; ProteinMotion raises for the same reason.

5. **A clip is a description.** `canvas.camera.orbit(90)` moves nothing. It is
   played by a timeline, which hands it a frame range and an easing, and the
   camera pose it starts from is whatever earlier clips left at that frame.
   That is what makes clips composable and re-orderable.

## 3. The model

```
Timeline            cursor in seconds, fps from the canvas, start frame
  .play(*clips, run_time=, easing=, at=)   clips run together; cursor advances past the longest
  .wait(seconds)                           hold
  .set(target, value)                      instant step at the cursor
  .tween(target, value) / .frames(entity, a, b)   clip factories
  .finish() / with-block exit              fits scene frame range, rewinds
  .render(path)                            Canvas.animation over the storyboard

Clip                run_time, easing, apply(timeline, frame_start, frame_end, easing) -> channels
  Tween             two keys per channel, Blender easing
  Step              hold-then-jump (CONSTANT)
  Sampled           one linear key per frame, Python easing on progress
  Frames            entity.frame keyed from universe frame a to b, entity detached from the scene frame
  CameraMove        Tween whose end pose is solved at play time: LookAt, MoveTo, Dolly
  Orbit             Sampled rotation of the camera matrix about a pivot
  Focus             DOF via a focus empty whose location is tweened, plus f-stop

Channel             (owner struct, property, index): the unit of keyframing
```

`Channel` is deliberately thin. It wraps `keyframe_insert`, finds its own
F-curve (slotted actions, Blender 4.4+), evaluates itself at a frame so that a
tween can start from where the previous clip left the property, and casts to
the property's type so integer and boolean properties can be tweened too.

### Targets

`resolve_channels(target, value)` is the single place that turns user-facing
objects into channels:

| target | resolves to |
| --- | --- |
| nodebpy socket `style.i.loop_radius` or a `bpy` node socket | the socket's `default_value` on the node tree |
| `(camera.data, "lens")`, `(obj, "location")`, `(scene.view_settings, "exposure")` | any RNA property on any struct |
| `(entity, "subframes")` | the entity's `mn` property group |
| `(annotation, "text_size")` | the annotation's entry in `object.mn_annotations`, or its typed inputs group |
| a sequence value on an array property | one channel per component |

This means material sockets, world sockets and compositor sockets are covered
without special cases; they are node sockets.

### Camera

The camera stays a free object with keys on `location` and `rotation_euler`.
The alternative, a pivot-empty rig with the camera parented to it, was
considered and rejected for the prototype: it changes what `Camera.frame_points`
writes (local versus world location), it makes `look_at` after an orbit depend
on rig state, and a script re-run is the primary edit path anyway. The trade is
that an orbit is a dense curve in the Graph Editor rather than one editable key
pair. Section 6 revisits this.

`LookAt` reuses `Camera.set_viewpoint` and `Camera.frame_points` in a dry run
(pose captured, camera restored), so the destination is exactly what an
immediate `canvas.look_at` gives. It keys `clip_end` as well when the framing
would otherwise clip the subject.

`Focus` uses Blender's own depth of field with a focus object rather than a
distance, so every later camera move keeps the subject in focus without extra
keys. The preset scene already ships a `focal_point` empty; the clip adopts it.

### Entity playback

`Frames` sets `update_with_scene = False` and keys `mn.frame` from universe
frame `a` to `b`, respecting the entity's subframes (the property is stepped
through them like the scene frame). Before the clip the trajectory holds `a`;
after it, `b`. Playing at any speed, pausing, or scrubbing backwards falls out
of ordinary F-curve behaviour.

## 4. What the prototype does not do

- **Track a moving selection.** `Focus` places the empty at the target's
  centre when played; if the trajectory then moves the residues, the focus
  stays where they were. Following a selection needs a per-frame update of
  the empty, which can be done by a constraint to a helper object with a
  Geometry Nodes centroid, or by the existing entity frame handler. Worth
  doing; a `follow=True` flag is the obvious shape.
- **Fade or opacity clips.** `Set Color` alpha reaches the `Default` material
  but not `Flat` or `AmbientOcclusion`; a `Set Opacity` node (comparison doc
  5.6) is the prerequisite. A tween on a material socket already works for
  materials that expose one.
- **Write-on labels.** Needs a `progress` property on text annotations.
- **Retiming or removing clips.** The timeline appends keys; it does not own
  them. Re-running the script after `canvas.clear()` or on a fresh preset is
  the way to change a storyboard. Deleting a clip's keys is a small addition
  if it turns out to matter.
- **NLA.** Each storyboard writes into the objects' active actions. Wrapping a
  storyboard as an NLA strip per object would let several be layered and
  offset in the UI, and is the natural next step if scripts start composing
  storyboards.

## 5. Node-tree interaction (what could be automated)

The comparison doc suggested the timeline "could result in automatic creation
of nodes, storing / updating of attributes". Having built the prototype, the
line to draw is:

- **Keys on existing node inputs, yes.** This is the whole point: the timeline
  is how a script animates a style parameter. It works today for every input
  of every style node, and for `Set Color`, `Animate Value` and the symmetry
  nodes' `Factor`.
- **Creating nodes, no.** A clip that inserted an `Animate Value` node into the
  entity tree would make the timeline responsible for tree layout and for
  the node's lifetime, and would leave two ways of animating the same thing.
  The stagger node from the comparison doc (4.4) belongs in the tree, driven
  by `Frame Start` / `Frame End` inputs that the timeline can key like any
  other socket. That is the interface between the two: nodes read the scene
  frame and expose inputs, the timeline keys inputs.
- **Attributes, via nodes.** Animating a per-atom attribute (a morph target,
  a fade mask) means a `Store Named Attribute` in the tree with a keyable
  factor, not the timeline writing attribute arrays per frame. Python-side
  per-frame writes are exactly the `record()` loop this layer exists to avoid.

## 6. Open questions

1. **Camera rig or free camera?** A rig (pivot empty, camera child, keys on the
   empty's rotation and the camera's local Z) gives one key pair per orbit and
   a familiar turntable setup for GUI users, at the cost of changing the
   `Camera` contract. If the GUI-editability of orbits matters more than the
   script being the source of truth, switch. The clip classes would not
   change; only `Orbit.apply` and `Camera.frame_points` would.
2. **Where do clip factories live?** Camera moves sit on `Camera` because they
   read as `canvas.camera.orbit(...)`. Entity clips are on the timeline
   (`t.frames(traj, 0, 100)`) to keep `Molecule` untouched; `traj.play(0, 100)`
   would read better and is a one-line addition once the shape settles.
   Annotation clips (`label.write()`) would follow the same pattern.
3. **`add_style` should return the style node.** Every parameter tween in the
   examples needs a handle on the node, which today means building the tree
   with `with mol.tree`. Returning the node (or a thin style object) from
   `add_style` would make `t.tween(style.i.quality, 4)` available from the
   simple API.
4. **Seconds or frames?** The cursor is in seconds so a storyboard survives an
   fps change. Frame-based `at=` and `run_time=` overrides may be wanted for
   syncing to trajectory frames; the conversion is trivial to expose.
5. **Conflict granularity.** Overlap detection is per channel and per frame
   range, so a tween and an orbit on the camera conflict only if both touch
   `rotation_euler`. A `dolly` during an `orbit` conflicts on `location`;
   composing them would need the orbit to sample the dolly's contribution,
   i.e. a clip stack per channel. Not needed until someone asks.

## 7. Verification

`tests/test_timeline.py` (30 tests) covers easing endpoints and monotonicity,
cursor arithmetic, node-input tweens evaluated through `frame_set`, holds
between clips, interpolation written onto keys, steps, vector and integer
channels, entity and annotation channels, conflict detection, `look_at`
matching the immediate call, orbit radius and pointing, continuity of the
euler curve over a full turn, dolly and zoom, focus via the empty, trajectory
playback with and without subframes, and a three-frame Cycles render of a
storyboard to MP4.

## 8. Picking this up

State as of 2026-09-19, PR "Canvas timeline layer (design discussion + prototype)".

- **Code**: `molecularnodes/scene/timeline.py` (clips, channels, timeline),
  clip factories at the bottom of `molecularnodes/scene/camera.py`,
  `Canvas.timeline()` in `molecularnodes/scene/base.py`, exports in
  `molecularnodes/scene/__init__.py`.
- **Tests**: `uv run pytest tests/test_timeline.py -q` (31 tests, ~10 s, one
  tiny Cycles render). `tests/test_canvas.py` and `tests/test_framing.py` were
  run alongside and pass.
- **Demo**: `uv run python docs/dev/canvas-timeline-demo.py <out_dir>` renders
  the storyboard in section 3 to `storyboard.mp4`, a contact sheet and a
  `.blend` holding the keyframes (EEVEE, 960x540, ~2 min). Open the blend to
  see what the timeline wrote in the Graph Editor.
- **Environment**: `uv sync --all-extras` in a fresh checkout or worktree,
  otherwise `bpy` is missing. Set `BLENDER_USER_EXTENSIONS` before importing
  `bpy` in ad-hoc scripts (the demo does).
- **Known quirk seen while testing, not from this work**: saving a `.blend`
  after loading a trajectory from a relative path raises in the session's
  `save_post` handler (`_remap_trajectory_paths`). Use absolute paths, or
  ignore; the file still saves.
- **Next steps**, in the order they unblock each other: settle section 6
  questions 1 to 3 on the PR; then `add_style` returning its node (small,
  makes `t.tween(style.i.x, v)` reachable from the simple API); then
  `follow=True` on `Focus`; then a `Set Opacity` node so fades are possible.
