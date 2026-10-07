"""
A timeline layer on :class:`~molecularnodes.Canvas`.

A render script written against this module reads as a storyboard: a sequence
of *clips* - camera moves, value changes, trajectory playback - each played for
a length of time with an easing. The timeline does not run those clips itself.
It compiles them into keyframes on the Blender properties they touch, so the
result is an ordinary animated scene: the timeline in the Blender UI seeks
exactly, ``Canvas.animation`` renders it, and every curve can be inspected or
edited in the Graph Editor afterwards.

Design notes
------------
* **Blender's animation system is the engine.** A clip only ever inserts
  keyframes. There is no Python-side interpolation at render time, no frame
  handler and no snapshot/replay machinery: seeking to a frame is exact because
  a property has one F-curve.
* **One channel, one owner.** Two clips writing the same property over
  overlapping frames is an error at authoring time, not a silent last-write-wins.
* **Geometry already animates itself.** Trajectory positions follow the scene
  frame through the entity handlers and node trees read the scene time, so this
  layer is mostly about the camera, node inputs, entity playback and annotation
  parameters - the things that today need hand-placed keyframes.
* **Easing lives on the keyframes.** Tweens set the interpolation of the
  Blender keys so the curve is editable. The camera is a pivot rig so that an
  orbit is two keys on the pivot's rotation; only an orbit about an axis that
  is not one of the pivot's Euler components is sampled per frame.
"""

from __future__ import annotations
import math
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, Sequence
import bpy
import numpy as np
from mathutils import Euler, Matrix, Vector
from .. import framing
from ..blender import utils as blender_utils
from ..entities.base import MolecularEntity
from ..session import get_session

if TYPE_CHECKING:
    from .base import Canvas
    from .camera import Camera


# ---------------------------------------------------------------------------
# Easing
# ---------------------------------------------------------------------------


class Easing:
    """
    Rate functions a clip can be played with.

    Each name maps both to a Blender keyframe interpolation (used by tweens, so
    the curve stays editable in the Graph Editor) and to a Python function of
    ``t`` in ``[0, 1]`` (used by clips that are sampled per frame).
    """

    LINEAR = "linear"
    SMOOTH = "smooth"
    SINE = "sine"
    EASE_IN = "ease_in"
    EASE_OUT = "ease_out"

    #: Blender ``(interpolation, easing)`` per rate function.
    _BLENDER = {
        LINEAR: ("LINEAR", "AUTO"),
        SMOOTH: ("BEZIER", "AUTO"),
        SINE: ("SINE", "EASE_IN_OUT"),
        EASE_IN: ("QUAD", "EASE_IN"),
        EASE_OUT: ("QUAD", "EASE_OUT"),
    }

    _FUNCTIONS: dict[str, Callable[[float], float]] = {
        LINEAR: lambda t: t,
        # the quintic smoothstep, flat at both ends like Blender's auto-clamped
        # Bezier handles between two keys
        SMOOTH: lambda t: t * t * t * (t * (t * 6.0 - 15.0) + 10.0),
        SINE: lambda t: (1.0 - math.cos(math.pi * t)) / 2.0,
        EASE_IN: lambda t: t * t,
        EASE_OUT: lambda t: 1.0 - (1.0 - t) * (1.0 - t),
    }

    @classmethod
    def validate(cls, easing: str) -> str:
        key = str(easing).strip().lower()
        if key not in cls._BLENDER:
            raise ValueError(
                f"Unknown easing '{easing}'. Choose from {sorted(cls._BLENDER)}."
            )
        return key

    @classmethod
    def blender(cls, easing: str) -> tuple[str, str]:
        return cls._BLENDER[cls.validate(easing)]

    @classmethod
    def function(cls, easing: str) -> Callable[[float], float]:
        return cls._FUNCTIONS[cls.validate(easing)]


# ---------------------------------------------------------------------------
# Channels: the Blender properties a clip writes to
# ---------------------------------------------------------------------------


def _find_fcurve(
    id_data: bpy.types.ID, data_path: str, index: int
) -> bpy.types.FCurve | None:
    """
    The F-curve animating ``data_path[index]`` on an ID, if there is one.

    Actions are slotted (Blender 4.4+): the curves for this ID live in the
    channel bag of the slot the ID's animation data is assigned to.
    """
    anim = id_data.animation_data
    if anim is None or anim.action is None:
        return None
    slot = anim.action_slot
    if slot is None:
        return None
    for layer in anim.action.layers:
        for strip in layer.strips:
            bag = strip.channelbag(slot)
            if bag is None:
                continue
            fcurve = bag.fcurves.find(data_path, index=max(index, 0))
            if fcurve is not None:
                return fcurve
    return None


@dataclass(frozen=True)
class Channel:
    """
    One animatable Blender property: a struct, the property's name on it and,
    for array properties, the index.

    Anything ``keyframe_insert`` accepts can be a channel - a node socket's
    ``default_value``, an object's ``location[2]``, a camera's
    ``dof.focus_distance``, an entity's ``mn.frame``.
    """

    owner: bpy.types.bpy_struct
    attr: str
    index: int = -1

    @property
    def id_data(self) -> bpy.types.ID:
        return self.owner.id_data

    @property
    def data_path(self) -> str:
        return self.owner.path_from_id(self.attr)

    @property
    def key(self) -> tuple[str, str, int]:
        """Identity of the channel for conflict detection."""
        return (repr(self.id_data), self.data_path, self.index)

    def __repr__(self) -> str:
        suffix = "" if self.index < 0 else f"[{self.index}]"
        return f"{self.id_data.name}:{self.data_path}{suffix}"

    def get(self) -> float:
        value = getattr(self.owner, self.attr)
        return float(value if self.index < 0 else value[self.index])

    def _cast(self, value: float):
        # F-curves are floats; integer and boolean properties want their own type
        kind = self.owner.bl_rna.properties[self.attr].type
        if kind == "INT":
            return int(round(value))
        if kind == "BOOLEAN":
            return bool(value)
        return float(value)

    def set(self, value: float) -> None:
        value = self._cast(value)
        if self.index < 0:
            setattr(self.owner, self.attr, value)
        else:
            current = getattr(self.owner, self.attr)
            current[self.index] = value
            setattr(self.owner, self.attr, current)

    def fcurve(self) -> bpy.types.FCurve | None:
        return _find_fcurve(self.id_data, self.data_path, self.index)

    def value_at(self, frame: int) -> float:
        """The value the property has at ``frame`` - from its F-curve when it
        has one, otherwise the value it holds now."""
        fcurve = self.fcurve()
        if fcurve is None:
            return self.get()
        return float(fcurve.evaluate(frame))

    def insert(
        self, frame: int, value: float, interpolation: tuple[str, str] | None = None
    ) -> None:
        """
        Key the property to ``value`` at ``frame``.

        ``interpolation`` is the Blender ``(interpolation, easing)`` pair applied
        to *this* key, which in Blender governs the segment leading from it to
        the next key.
        """
        self.set(value)
        self.owner.keyframe_insert(self.attr, index=self.index, frame=frame)
        if interpolation is None:
            return
        fcurve = self.fcurve()
        if fcurve is None:  # pragma: no cover - keyframe_insert succeeded
            return
        for point in fcurve.keyframe_points:
            if int(round(point.co.x)) == frame:
                point.interpolation, point.easing = interpolation
                if point.interpolation == "BEZIER":
                    point.handle_left_type = "AUTO_CLAMPED"
                    point.handle_right_type = "AUTO_CLAMPED"
                break
        fcurve.update()


def _array_length(owner: bpy.types.bpy_struct, attr: str) -> int:
    prop = owner.bl_rna.properties.get(attr)
    if prop is None:
        raise ValueError(f"{owner!r} has no property '{attr}'.")
    return getattr(prop, "array_length", 0)


def _annotation_struct(interface, attr: str) -> bpy.types.bpy_struct:
    """
    The property group holding ``attr`` for an annotation interface.

    Common parameters (``visible``, ``text_size``, ...) live on the annotation's
    entry in ``object.mn_annotations``; the annotation type's own inputs live in
    a nested group named after the entity and annotation type.
    """
    instance = interface._instance
    # molecule annotations keep their entity as `trajectory`
    entity = getattr(instance, "entity", None) or getattr(instance, "trajectory", None)
    if entity is None:
        raise TypeError("Annotation is not attached to an entity.")
    prop = entity.object.mn_annotations[interface._uuid]
    if attr in prop.bl_rna.properties:
        return prop
    entity_type = entity._get_annotation_entity_type()
    inputs = getattr(prop, f"{entity_type}_{prop.type}", None)
    if inputs is not None and attr in inputs.bl_rna.properties:
        return inputs
    raise ValueError(f"Annotation '{interface.name}' has no input '{attr}'.")


def resolve_channels(target, value) -> list[tuple[Channel, float]]:
    """
    Turn a user-facing target and value into ``(channel, value)`` pairs.

    Parameters
    ----------
    target
        One of:

        - a nodebpy socket (``style.i.loop_radius``) or a Blender node socket,
          which animates its ``default_value``
        - a ``(owner, "attr")`` pair, where ``owner`` is any Blender struct
          (``(camera.data, "lens")``, ``(obj, "location")``), a Molecular Nodes
          entity (its ``mn`` properties: ``(traj, "subframes")``), an
          annotation interface (``(label, "text_size")``) or a :class:`Camera`
        - a :class:`Channel`
    value
        A number, or a sequence of numbers for an array property, in which case
        one channel per component is returned.
    """
    if isinstance(target, Channel):
        return [(target, float(value))]

    if hasattr(target, "socket") and isinstance(target.socket, bpy.types.NodeSocket):
        target = target.socket
    if isinstance(target, bpy.types.NodeSocket):
        owner, attr = target, "default_value"
    elif isinstance(target, tuple) and len(target) == 2:
        owner, attr = target
        from .camera import Camera

        if isinstance(owner, MolecularEntity):
            owner = owner.object.mn
        elif isinstance(owner, Camera):
            owner = (
                owner.camera_data
                if attr in owner.camera_data.bl_rna.properties
                else owner.camera
            )
        elif hasattr(owner, "_instance") and hasattr(owner, "_uuid"):
            owner = _annotation_struct(owner, attr)
        if not isinstance(owner, bpy.types.bpy_struct):
            raise TypeError(f"Cannot animate properties of {type(owner).__name__}.")
    else:
        raise TypeError(
            "A target must be a node socket, a (struct, 'property') pair or a "
            f"Channel, not {type(target).__name__}."
        )

    length = _array_length(owner, attr)
    if length == 0:
        if isinstance(value, (Sequence, np.ndarray)) and not isinstance(value, str):
            raise ValueError(f"'{attr}' is a scalar property; got {value!r}.")
        return [(Channel(owner, attr), float(value))]
    values = np.asarray(value, dtype=np.float64).ravel()
    if len(values) != length:
        raise ValueError(f"'{attr}' has {length} components; got {len(values)} values.")
    return [(Channel(owner, attr, i), float(v)) for i, v in enumerate(values)]


def target_points(target) -> np.ndarray:
    """
    World-space points for a framing or pivot target.

    An entity or object is read as the geometry it renders; ``None`` means every
    entity in the session; anything else is taken as ``(N, 3)`` positions.
    """
    if target is None:
        points = [
            blender_utils.evaluated_points(entity.object)
            for entity in get_session().entities.values()
            if entity.object.name in bpy.context.scene.objects
        ]
        if not points:
            raise ValueError("The scene holds no entities to use as a target.")
        return np.concatenate(points)
    if isinstance(target, MolecularEntity):
        target = target.object
    if isinstance(target, bpy.types.Object):
        return blender_utils.evaluated_points(target)
    return framing.as_points(target)


def target_centre(target) -> Vector:
    centre, _ = framing.enclosing_sphere(target_points(target))
    return Vector(centre)


# ---------------------------------------------------------------------------
# Clips
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Key:
    """One keyframe a clip asks for: a channel, a frame, a value, and the
    ``(interpolation, easing)`` of the segment leading on from it."""

    channel: Channel
    frame: int
    value: float
    interpolation: tuple[str, str] | None = None


def _tween_keys(
    pairs: Sequence[tuple[Channel, float]], frame_start: int, frame_end: int, easing
) -> list[Key]:
    """Keys easing each channel from its value at ``frame_start`` to a value at
    ``frame_end``, with the easing on the first key."""
    interpolation = Easing.blender(easing)
    keys = []
    for channel, value in pairs:
        if frame_end > frame_start:
            keys.append(
                Key(channel, frame_start, channel.value_at(frame_start), interpolation)
            )
        keys.append(Key(channel, frame_end, value, interpolation))
    return keys


def _step_keys(pairs: Sequence[tuple[Channel, float]], frame: int) -> list[Key]:
    """Keys holding each channel's value up to the frame before ``frame`` and
    jumping to a new value on it."""
    hold = ("CONSTANT", "AUTO")
    keys = []
    for channel, value in pairs:
        keys.append(Key(channel, frame - 1, channel.value_at(frame - 1), hold))
        keys.append(Key(channel, frame, value, hold))
    return keys


class Clip(ABC):
    """
    Something that can be played on a timeline.

    A clip is a description - nothing happens when it is created. When it is
    played, :meth:`keys` is given the frame range and easing and returns the
    keyframes to write. The timeline checks them against what earlier clips
    wrote before any of them are inserted, so a clip that would fight another
    over a property leaves the scene untouched.
    """

    run_time: float = 1.0
    easing: str = Easing.SMOOTH
    label: str = "clip"

    def __init__(self, run_time: float | None = None, easing: str | None = None):
        if run_time is not None:
            if run_time < 0:
                raise ValueError(f"run_time must not be negative, got {run_time}.")
            self.run_time = run_time
        if easing is not None:
            self.easing = Easing.validate(easing)

    @abstractmethod
    def keys(
        self, timeline: "Timeline", frame_start: int, frame_end: int, easing: str
    ) -> list[Key]: ...

    def commit(self, timeline: "Timeline") -> None:
        """Changes other than keyframes, made once the keys are accepted."""

    def __repr__(self) -> str:
        return f"<{self.label} {self.run_time:g}s {self.easing}>"


class Tween(Clip):
    """
    Ease one or more properties from what they are to new values.

    Two keys per channel: the value the property has at the start frame (read
    off its F-curve, so a tween continues from wherever an earlier clip left
    the property) and the target value at the end frame, with the easing set on
    the first key so the curve between them is Blender's own.
    """

    label = "tween"

    def __init__(
        self,
        target=None,
        value=None,
        run_time: float | None = None,
        easing: str | None = None,
        *,
        pairs: Sequence[tuple[Channel, float]] | None = None,
    ):
        super().__init__(run_time, easing)
        if pairs is None:
            if target is None:
                raise ValueError("A tween needs a target and a value.")
            pairs = resolve_channels(target, value)
        self._pairs = list(pairs)
        self._name = self._pairs[0][0].data_path if self._pairs else "tween"

    def channels(self) -> list[Channel]:
        return [channel for channel, _ in self._pairs]

    def keys(self, timeline, frame_start, frame_end, easing) -> list[Key]:
        return _tween_keys(self._pairs, frame_start, frame_end, easing)

    def __repr__(self) -> str:
        return f"<tween {self._name} {self.run_time:g}s {self.easing}>"


class Step(Tween):
    """
    Jump properties to new values at a frame, holding the previous value up to
    the frame before. The instant form of a tween; ``run_time`` is always 0.
    """

    label = "set"
    run_time = 0.0

    def keys(self, timeline, frame_start, frame_end, easing) -> list[Key]:
        return _step_keys(self._pairs, frame_start)


class Sampled(Clip):
    """
    A clip keyed on every frame it spans.

    For motion that is not linear in the animated properties - a camera arc,
    say - two keys with an easing would cut the corner. Instead the clip is
    asked for its state at each frame, with the easing already applied to the
    ``0..1`` progress, and keyed linearly so that seeking between frames stays
    exact.
    """

    label = "sampled"

    def channels(self) -> list[Channel]:
        raise NotImplementedError

    def prepare(self) -> None:
        """Capture the starting state, right before the frames are written."""

    def sample(self, u: float) -> Sequence[float]:
        """Values for :meth:`channels` at eased progress ``u`` in ``[0, 1]``."""
        raise NotImplementedError

    def keys(self, timeline, frame_start, frame_end, easing) -> list[Key]:
        rate = Easing.function(easing)
        self.prepare()
        channels = self.channels()
        linear = Easing.blender(Easing.LINEAR)
        span = max(frame_end - frame_start, 1)
        keys = []
        for frame in range(frame_start, frame_end + 1):
            u = rate(min((frame - frame_start) / span, 1.0))
            for channel, value in zip(channels, self.sample(u)):
                keys.append(Key(channel, frame, value, linear))
        return keys


class Frames(Clip):
    """
    Play a stretch of an entity's trajectory over the clip's run time.

    Detaches the entity from the scene frame and keys its own ``frame``
    property instead, so a trajectory can be played at any speed, paused, or
    scrubbed backwards from a script. Universe frames are given; subframes on
    the entity are honoured (the property steps through them like the scene
    frame does).
    """

    label = "frames"
    easing = Easing.LINEAR

    def __init__(
        self,
        entity: MolecularEntity,
        start: int,
        end: int,
        run_time: float | None = None,
        easing: str | None = None,
    ):
        super().__init__(run_time, easing)
        self._entity = entity
        self._start, self._end = int(start), int(end)

    def keys(self, timeline, frame_start, frame_end, easing) -> list[Key]:
        step = int(getattr(self._entity, "subframes", 0)) + 1
        channel = Channel(self._entity.object.mn, "frame")
        interpolation = Easing.blender(easing)
        return [
            Key(channel, frame_start, self._start * step, interpolation),
            Key(channel, frame_end, self._end * step, interpolation),
        ]

    def commit(self, timeline) -> None:
        self._entity.update_with_scene = False

    def __repr__(self) -> str:
        return (
            f"<frames {self._entity.name} {self._start}->{self._end} "
            f"{self.run_time:g}s {self.easing}>"
        )


# ---------------------------------------------------------------------------
# Camera clips
# ---------------------------------------------------------------------------
#
# Camera clips rig the camera first (see `Camera.rig`): a pivot empty it is
# parented to and orbits about, a target empty it tracks, both sitting on the
# view axis at the subject, and a focus empty for depth of field. A framing
# move keys the pivot, target and camera, an orbit keys the pivot's rotation,
# a dolly keys the camera's own location and a focus pull keys the focus
# empty, so they compose.


def _rig(camera: "Camera") -> None:
    """Rig the camera if it is not already, with the pivot at the depth of the
    scene's entities, or at the focus distance if there are none."""
    if camera.is_rigged:
        return
    distance = camera.camera_data.dof.focus_distance
    try:
        depth = (target_centre(None) - camera.location).dot(Vector(camera.basis[2]))
    except ValueError:
        depth = 0.0
    if depth > camera.clip_start:
        distance = depth
    camera.rig(distance)


def _transform_pairs(
    obj: bpy.types.Object,
    location: Vector | None = None,
    rotation: Euler | None = None,
    compatible: Euler | None = None,
) -> list[tuple[Channel, float]]:
    """
    ``(channel, value)`` pairs putting an object at a local location and
    rotation. The rotation is made compatible with ``compatible`` (the
    object's current rotation by default), so the curve does not take the long
    way round.
    """
    pairs: list[tuple[Channel, float]] = []
    if location is not None:
        pairs += [(Channel(obj, "location", i), float(location[i])) for i in range(3)]
    if rotation is not None:
        rotation = rotation.copy()
        rotation.make_compatible(
            obj.rotation_euler if compatible is None else compatible
        )
        pairs += [
            (Channel(obj, "rotation_euler", i), float(rotation[i])) for i in range(3)
        ]
    return pairs


class CameraMove(Tween):
    """
    Ease the camera rig from where it is to a state solved when the clip is
    played.

    Subclasses provide :meth:`end_state`, evaluated with the rig as earlier
    clips leave it at the start of this one. Only the channels the move changes
    are keyed, so a move that leaves the pivot alone leaves it free for an
    orbit played at the same time.
    """

    label = "camera"

    def __init__(self, camera: "Camera", run_time=None, easing=None):
        Clip.__init__(self, run_time, easing)
        self._camera = camera
        self._pairs = []

    def end_state(self) -> list[tuple[Channel, float]]:
        raise NotImplementedError

    def keys(self, timeline, frame_start, frame_end, easing) -> list[Key]:
        _rig(self._camera)
        self._pairs = [
            (channel, value)
            for channel, value in self.end_state()
            if abs(channel.value_at(frame_start) - value) > 1e-9
        ]
        return super().keys(timeline, frame_start, frame_end, easing)

    def __repr__(self) -> str:
        return f"<camera.{self.label} {self.run_time:g}s {self.easing}>"


class _PoseMove(CameraMove):
    """
    A camera move whose destination is what an immediate change to the camera
    gives: :meth:`pose` is run on the rig, its transforms (and the far clip
    and orthographic scale, which framing can change) are read off, and the rig
    is put back as it was.
    """

    def pose(self) -> None:
        raise NotImplementedError

    def end_state(self) -> list[tuple[Channel, float]]:
        cam = self._camera
        obj, pivot, target = cam.camera, cam.pivot, cam.target
        data = cam.camera_data
        saved = (
            pivot.location.copy(),
            pivot.rotation_euler.copy(),
            obj.location.copy(),
            obj.rotation_euler.copy(),
            target.location.copy(),
            data.clip_end,
            data.ortho_scale,
        )
        try:
            self.pose()
            # compatible with the rotations the move starts from, not the
            # ones the dry run just wrote
            pairs = _transform_pairs(
                pivot, pivot.location, pivot.rotation_euler, compatible=saved[1]
            )
            pairs += _transform_pairs(
                obj, obj.location, obj.rotation_euler, compatible=saved[3]
            )
            pairs += _transform_pairs(target, target.location)
            if data.clip_end > saved[5]:
                pairs.append((Channel(data, "clip_end"), data.clip_end))
            if data.type == "ORTHO":
                pairs.append((Channel(data, "ortho_scale"), data.ortho_scale))
            return pairs
        finally:
            pivot.location, pivot.rotation_euler = saved[0], saved[1]
            obj.location, obj.rotation_euler = saved[2], saved[3]
            target.location = saved[4]
            data.clip_end, data.ortho_scale = saved[5], saved[6]


class LookAt(_PoseMove):
    """Ease the rig to the framing ``Canvas.look_at`` would jump to."""

    label = "look_at"

    def __init__(
        self,
        camera,
        target,
        viewpoint=None,
        margin: float = 0.05,
        run_time=None,
        easing=None,
    ):
        super().__init__(camera, run_time, easing)
        self._target, self._viewpoint, self._margin = target, viewpoint, margin

    def pose(self) -> None:
        cam = self._camera
        if self._viewpoint is not None:
            cam.set_viewpoint(self._viewpoint)
        cam.frame_points(target_points(self._target), margin=self._margin)


class MoveTo(_PoseMove):
    """Ease the camera to an explicit world location and/or rotation (degrees)."""

    label = "move_to"

    def __init__(
        self, camera, location=None, rotation=None, run_time=None, easing=None
    ):
        super().__init__(camera, run_time, easing)
        if location is None and rotation is None:
            raise ValueError("move_to needs a location, a rotation, or both.")
        self._location, self._rotation = location, rotation

    def pose(self) -> None:
        cam = self._camera
        if self._rotation is not None:
            cam.rotation = self._rotation
        if self._location is not None:
            cam.location = self._location


class Dolly(CameraMove):
    """Move the camera along its view axis; positive is towards the target."""

    label = "dolly"

    def __init__(self, camera, distance: float, run_time=None, easing=None):
        super().__init__(camera, run_time, easing)
        self._distance = float(distance)

    def end_state(self) -> list[tuple[Channel, float]]:
        cam = self._camera
        world = cam.location + Vector(cam.basis[2]) * self._distance
        return _transform_pairs(
            cam.camera, location=cam.parent_matrix.inverted() @ world
        )


class Orbit(Sampled):
    """
    Turn the camera's pivot about an axis, keeping the camera pointed the same
    way relative to the subject.

    Turning about the world ``z`` axis or the camera's ``right`` axis changes
    one component of the pivot's XYZ Euler rotation, so it is written as a
    pair of keys with the easing on them - one editable curve. Any other axis
    is sampled per frame. With ``about`` given, the pivot and the target ease
    to that centre over the same frames.
    """

    label = "orbit"
    run_time = 2.0

    def __init__(
        self,
        camera: "Camera",
        angle: float,
        axis="z",
        about=None,
        run_time=None,
        easing=None,
    ):
        super().__init__(run_time, easing)
        self._camera = camera
        self._angle = math.radians(angle)
        self._axis = axis
        self._about = about
        self._start: Euler | None = None
        self._axis_vector: Vector | None = None
        self._previous: Euler | None = None

    def channels(self) -> list[Channel]:
        return [Channel(self._camera.pivot, "rotation_euler", i) for i in range(3)]

    def prepare(self) -> None:
        _rig(self._camera)
        pivot = self._camera.pivot
        self._start = pivot.rotation_euler.copy()
        self._previous = pivot.rotation_euler.copy()
        axis = self._axis
        if isinstance(axis, str):
            named = {
                "x": Vector((1, 0, 0)),
                "y": Vector((0, 1, 0)),
                "z": Vector((0, 0, 1)),
                # the camera's own axes, for tumbling relative to the view
                "up": Vector(self._camera.basis[1]),
                "right": Vector(self._camera.basis[0]),
            }
            if axis.lower() not in named:
                raise ValueError(
                    f"Unknown orbit axis '{axis}'; use {sorted(named)} or a vector."
                )
            axis = named[axis.lower()]
        self._axis_vector = Vector(axis).normalized()

    def euler_component(self) -> tuple[int, float] | None:
        """
        The index of the pivot's Euler component that turns about the orbit
        axis, and the sign to turn it with, if there is one.

        For an XYZ Euler ``(x, y, z)`` the z component turns about world Z, the
        x component about the pivot's own X axis and the y component about
        world Z's rotation of Y; a matching axis can be keyed with two keys.
        """
        assert self._start is not None and self._axis_vector is not None
        x, y, z = self._start
        candidates = {
            2: Vector((0.0, 0.0, 1.0)),
            0: self._start.to_matrix() @ Vector((1.0, 0.0, 0.0)),
            1: Matrix.Rotation(z, 3, "Z") @ Vector((0.0, 1.0, 0.0)),
        }
        for index, axis in candidates.items():
            dot = axis.dot(self._axis_vector)
            if abs(dot) > 1.0 - 1e-6:
                return index, math.copysign(1.0, dot)
        return None

    def sample(self, u: float) -> Sequence[float]:
        assert self._start is not None and self._previous is not None
        matrix = Matrix.Rotation(self._angle * u, 3, self._axis_vector)
        euler = (matrix @ self._start.to_matrix()).to_euler("XYZ", self._previous)
        # keep successive eulers continuous so the curve never wraps
        self._previous = euler
        return list(euler)

    def keys(self, timeline, frame_start, frame_end, easing) -> list[Key]:
        self.prepare()
        cam = self._camera
        keys: list[Key] = []
        if self._about is not None:
            centre = target_centre(self._about)
            pairs = [
                (Channel(obj, "location", i), float(centre[i]))
                for obj in (cam.pivot, cam.target)
                for i in range(3)
                if abs(obj.location[i] - centre[i]) > 1e-9
            ]
            keys += _tween_keys(pairs, frame_start, frame_end, easing)
        component = self.euler_component()
        if component is None:
            return keys + super().keys(timeline, frame_start, frame_end, easing)
        index, sign = component
        channel = Channel(cam.pivot, "rotation_euler", index)
        end = channel.value_at(frame_start) + sign * self._angle
        return keys + _tween_keys([(channel, end)], frame_start, frame_end, easing)

    def __repr__(self) -> str:
        return f"<camera.orbit {math.degrees(self._angle):g}deg {self.run_time:g}s {self.easing}>"


#: The f-stop a focus pull opens up from when depth of field was off: stopped
#: down far enough that nothing visibly blurs.
_SHARP_FSTOP = 128.0


class Focus(Clip):
    """
    Pull focus onto a target with depth of field.

    Moves the rig's focus empty to the target's centre without moving the
    view, so the focus distance follows every later camera move for free. The f-stop is
    eased too. If depth of field is off where the clip starts, it is keyed on
    at the first frame and the aperture opens up from :data:`_SHARP_FSTOP`, so
    the blur eases in instead of popping and nothing before the clip changes.
    """

    label = "focus"

    def __init__(
        self, camera: "Camera", target, fstop: float = 2.8, run_time=None, easing=None
    ):
        super().__init__(run_time, easing)
        self._camera = camera
        self._target = target
        self._fstop = float(fstop)

    def keys(self, timeline, frame_start, frame_end, easing) -> list[Key]:
        cam = self._camera
        _rig(cam)
        dof = cam.camera_data.dof
        centre = target_centre(self._target)
        pairs = [(Channel(cam.focus_point, "location", i), centre[i]) for i in range(3)]
        keys = _tween_keys(pairs, frame_start, frame_end, easing)
        fstop = Channel(dof, "aperture_fstop")
        use_dof = Channel(dof, "use_dof")
        if use_dof.value_at(frame_start):
            return keys + _tween_keys(
                [(fstop, self._fstop)], frame_start, frame_end, easing
            )
        keys += _step_keys([(use_dof, 1.0)], frame_start)
        interpolation = Easing.blender(easing)
        keys.append(Key(fstop, frame_start, _SHARP_FSTOP, interpolation))
        keys.append(Key(fstop, frame_end, self._fstop, interpolation))
        return keys

    def __repr__(self) -> str:
        return f"<camera.focus f/{self._fstop:g} {self.run_time:g}s {self.easing}>"


# ---------------------------------------------------------------------------
# The timeline
# ---------------------------------------------------------------------------


@dataclass
class Entry:
    clip: Clip
    frame_start: int
    frame_end: int

    @property
    def frames(self) -> int:
        return self.frame_end - self.frame_start


class Timeline:
    """
    Compile a storyboard of clips into keyframes on the scene.

    Created with [](`~mn.Canvas.timeline`). Time is kept in seconds and
    converted to frames at the canvas's frame rate; a cursor advances as clips
    are played, so a script reads top to bottom as the animation plays.

    Examples
    --------
    ::

        cartoon = mol.styles["Style Cartoon"]
        with canvas.timeline(fps=30) as t:
            t.play(canvas.camera.look_at(mol, viewpoint="front"), run_time=1)
            t.play(canvas.camera.orbit(90), t.tween(cartoon.i.loop_radius, 1.2), run_time=2)
            t.wait(0.5)
            t.play(canvas.camera.focus(mol.get_view("resid 40-60"), fstop=2), run_time=1)
        canvas.animation("story.mp4")

    A trajectory plays through its own clip, ``traj.play(0, 100)``, which
    detaches it from the scene frame and keys its frame property instead.

    Leaving the ``with`` block sets the scene's frame range to what was played,
    so [](`~mn.Canvas.animation`) renders exactly the storyboard.
    """

    def __init__(
        self, canvas: "Canvas", fps: float | None = None, start: int | None = None
    ):
        self.canvas = canvas
        if fps is not None:
            canvas.fps = fps
        self.start = int(canvas.frame_start if start is None else start)
        self._cursor = 0.0
        self._entries: list[Entry] = []
        self._written: dict[tuple, list[tuple[int, int, Clip]]] = {}

    # -- time -------------------------------------------------------------

    @property
    def fps(self) -> float:
        return self.canvas.fps

    @property
    def now(self) -> float:
        """The cursor, in seconds from the start of the timeline."""
        return self._cursor

    @property
    def frame(self) -> int:
        """The scene frame the cursor sits on."""
        return self.to_frame(self._cursor)

    @property
    def frame_end(self) -> int:
        """The last frame any clip wrote, or the start frame if none did."""
        if not self._entries:
            return self.frame
        return max(self.frame, max(entry.frame_end for entry in self._entries))

    def to_frame(self, seconds: float) -> int:
        return self.start + int(round(seconds * self.fps))

    @property
    def entries(self) -> list[Entry]:
        return list(self._entries)

    # -- playing ----------------------------------------------------------

    def play(
        self,
        *clips: Clip,
        run_time: float | None = None,
        easing: str | None = None,
        at: float | None = None,
    ) -> "Timeline":
        """
        Play clips together, then advance the cursor past the longest of them.

        Parameters
        ----------
        *clips : Clip
            Clips to play at the same time.
        run_time : float, optional
            Length in seconds for every clip given, overriding each clip's own.
        easing : str, optional
            Easing for every clip given, overriding each clip's own.
        at : float, optional
            Start the clips at this time instead of the cursor, leaving the
            cursor where it is - for layering a clip over an earlier passage
            on properties that are not animated later on.

        Raises
        ------
        ValueError
            If a clip writes a property that another clip already animates
            over overlapping frames, or that is already keyed later than the
            clip ends: clips on one property are played in time order. Nothing
            the clip would have written is kept.
        """
        if not clips:
            raise ValueError("play() needs at least one clip.")
        if easing is not None:
            easing = Easing.validate(easing)
        t0 = self._cursor if at is None else float(at)
        frame_start = self.to_frame(t0)
        longest = 0.0
        for clip in clips:
            if not isinstance(clip, Clip):
                raise TypeError(f"{clip!r} is not a clip.")
            length = clip.run_time if run_time is None else run_time
            frame_end = self.to_frame(t0 + length)
            # put the scene at the clip's first frame, so the clip sees the
            # camera, geometry and values as earlier clips leave them there
            self.canvas.scene.frame_set(frame_start)
            keys = clip.keys(self, frame_start, frame_end, easing or clip.easing)
            self._claim(clip, keys, frame_start, frame_end)
            clip.commit(self)
            for key in keys:
                key.channel.insert(key.frame, key.value, key.interpolation)
            self._entries.append(Entry(clip, frame_start, frame_end))
            longest = max(longest, length)
        if at is None:
            self._cursor = t0 + longest
        return self

    def wait(self, seconds: float) -> "Timeline":
        """Hold everything as it is for a while."""
        if seconds < 0:
            raise ValueError("wait() takes a non-negative number of seconds.")
        self._cursor += seconds
        return self

    def set(self, target, value) -> "Timeline":
        """Jump a property to a value at the cursor, without advancing it."""
        return self.play(Step(target, value))

    # -- clip factories ---------------------------------------------------

    def tween(
        self, target, value, run_time: float | None = None, easing: str | None = None
    ) -> Tween:
        """A clip easing ``target`` to ``value``; see :func:`resolve_channels`."""
        return Tween(target, value, run_time, easing)

    # -- bookkeeping ------------------------------------------------------

    def _claim(
        self, clip: Clip, keys: Sequence[Key], frame_start: int, frame_end: int
    ) -> None:
        """
        Check the keys a clip wants to write against what is already there,
        then record its channels. Raises before anything is recorded.
        """
        last: dict[tuple, tuple[Channel, int]] = {}
        for key in keys:
            _, frame = last.get(key.channel.key, (key.channel, key.frame))
            last[key.channel.key] = (key.channel, max(frame, key.frame))
        for channel_key, (channel, frame) in last.items():
            for a, b, other in self._written.get(channel_key, []):
                if a < frame_end and frame_start < b:
                    raise ValueError(
                        f"{clip!r} animates {channel!r} over frames {frame_start}-{frame_end}, "
                        f"which {other!r} already animates over {a}-{b}."
                    )
            fcurve = channel.fcurve()
            later = fcurve and [
                p.co.x for p in fcurve.keyframe_points if p.co.x > frame + 0.5
            ]
            if later:
                # the keys after this clip were written for the value the
                # property had then, so slotting a clip in before them would
                # leave the property drifting back instead of holding
                raise ValueError(
                    f"{clip!r} animates {channel!r} up to frame {frame}, but it is "
                    f"already keyed at frame {int(later[0])}. Play the clips on one "
                    "property in time order."
                )
        for channel_key in last:
            self._written.setdefault(channel_key, []).append(
                (frame_start, frame_end, clip)
            )

    def finish(self) -> None:
        """Fit the scene's frame range to the storyboard and rewind to its start."""
        self.canvas.frame_range = (self.start, self.frame_end)
        self.canvas.frame = self.start

    def __enter__(self) -> "Timeline":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        if exc_type is None:
            self.finish()

    def render(self, path=None, **kwargs):
        """Render the storyboard with [](`~mn.Canvas.animation`)."""
        self.finish()
        return self.canvas.animation(
            path, frame_start=self.start, frame_end=self.frame_end, **kwargs
        )

    def __repr__(self) -> str:
        lines = [
            f"<Timeline {self.fps:g} fps, frames {self.start}-{self.frame_end}, cursor {self.now:g}s>"
        ]
        for entry in self._entries:
            t0 = (entry.frame_start - self.start) / self.fps
            t1 = (entry.frame_end - self.start) / self.fps
            lines.append(
                f"  {t0:6.2f}s - {t1:6.2f}s  [{entry.frame_start:4d}-{entry.frame_end:4d}]  {entry.clip!r}"
            )
        return "\n".join(lines)


__all__ = [
    "Channel",
    "Clip",
    "Dolly",
    "Easing",
    "Entry",
    "Focus",
    "Frames",
    "Key",
    "LookAt",
    "MoveTo",
    "Orbit",
    "Sampled",
    "Step",
    "Timeline",
    "Tween",
    "resolve_channels",
    "target_centre",
    "target_points",
]
