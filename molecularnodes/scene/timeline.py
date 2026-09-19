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
  Blender keys so the curve is editable; only clips that cannot be expressed as
  two keys (an orbit, whose path is an arc) are sampled per frame.
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


class Clip(ABC):
    """
    Something that can be played on a timeline.

    A clip is a description - nothing happens when it is created. When it is
    played, :meth:`apply` is given the frame range and easing to compile into
    keyframes, and returns the channels it wrote so the timeline can detect two
    clips fighting over one property.
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
    def apply(
        self, timeline: "Timeline", frame_start: int, frame_end: int, easing: str
    ) -> list[Channel]: ...

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

    def apply(self, timeline, frame_start, frame_end, easing) -> list[Channel]:
        interpolation = Easing.blender(easing)
        for channel, value in self._pairs:
            if frame_end > frame_start:
                channel.insert(
                    frame_start, channel.value_at(frame_start), interpolation
                )
            channel.insert(frame_end, value, interpolation)
        return self.channels()

    def __repr__(self) -> str:
        return f"<tween {self._name} {self.run_time:g}s {self.easing}>"


class Step(Tween):
    """
    Jump properties to new values at a frame, holding the previous value up to
    the frame before. The instant form of a tween; ``run_time`` is always 0.
    """

    label = "set"
    run_time = 0.0

    def apply(self, timeline, frame_start, frame_end, easing) -> list[Channel]:
        hold = ("CONSTANT", "AUTO")
        for channel, value in self._pairs:
            channel.insert(frame_start - 1, channel.value_at(frame_start - 1), hold)
            channel.insert(frame_start, value, hold)
        return self.channels()


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

    def apply(self, timeline, frame_start, frame_end, easing) -> list[Channel]:
        rate = Easing.function(easing)
        self.prepare()
        channels = self.channels()
        linear = Easing.blender(Easing.LINEAR)
        span = max(frame_end - frame_start, 1)
        for frame in range(frame_start, frame_end + 1):
            u = rate(min((frame - frame_start) / span, 1.0))
            for channel, value in zip(channels, self.sample(u)):
                channel.insert(frame, value, linear)
        return channels


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

    def apply(self, timeline, frame_start, frame_end, easing) -> list[Channel]:
        entity = self._entity
        entity.update_with_scene = False
        step = int(getattr(entity, "subframes", 0)) + 1
        channel = Channel(entity.object.mn, "frame")
        interpolation = Easing.blender(easing)
        channel.insert(frame_start, self._start * step, interpolation)
        channel.insert(frame_end, self._end * step, interpolation)
        return [channel]

    def __repr__(self) -> str:
        return (
            f"<frames {self._entity.name} {self._start}->{self._end} "
            f"{self.run_time:g}s {self.easing}>"
        )


# ---------------------------------------------------------------------------
# Camera clips
# ---------------------------------------------------------------------------


def _camera_channels(camera: bpy.types.Object) -> list[Channel]:
    return [Channel(camera, "location", i) for i in range(3)] + [
        Channel(camera, "rotation_euler", i) for i in range(3)
    ]


def _pose_values(matrix: Matrix, previous: Euler) -> list[float]:
    """Location and a rotation compatible with ``previous`` from a matrix."""
    euler = matrix.to_euler("XYZ", previous)
    return [*matrix.translation, *euler]


class CameraMove(Tween):
    """
    Ease the camera from where it is to a pose solved when the clip is played.

    Subclasses provide :meth:`end_pose`, which is evaluated with the camera in
    its state at the start of the clip - what an earlier clip left it as - and
    must leave the camera as it found it.
    """

    label = "camera"

    def __init__(self, camera: "Camera", run_time=None, easing=None):
        Clip.__init__(self, run_time, easing)
        self._camera = camera
        self._pairs = []

    def end_pose(self) -> Matrix:
        raise NotImplementedError

    def extra_pairs(self) -> list[tuple[Channel, float]]:
        """Channels beyond location and rotation (lens, clip range...)."""
        return []

    def apply(self, timeline, frame_start, frame_end, easing) -> list[Channel]:
        obj = self._camera.camera
        previous = obj.rotation_euler.copy()
        matrix = self.end_pose()
        values = _pose_values(matrix, previous)
        self._pairs = list(zip(_camera_channels(obj), values)) + self.extra_pairs()
        return super().apply(timeline, frame_start, frame_end, easing)

    def __repr__(self) -> str:
        return f"<camera.{self.label} {self.run_time:g}s {self.easing}>"


class LookAt(CameraMove):
    """Ease the camera to the framing ``Canvas.look_at`` would jump to."""

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
        self._clip_end: float | None = None

    def end_pose(self) -> Matrix:
        cam = self._camera
        obj = cam.camera
        saved = (obj.location.copy(), obj.rotation_euler.copy(), cam.clip_end)
        try:
            if self._viewpoint is not None:
                cam.set_viewpoint(self._viewpoint)
                # the basis is read from matrix_world, which lags the rotation
                bpy.context.view_layer.update()
            cam.frame_points(target_points(self._target), margin=self._margin)
            bpy.context.view_layer.update()
            self._clip_end = cam.clip_end
            return obj.matrix_world.copy()
        finally:
            obj.location, obj.rotation_euler, cam.clip_end = saved
            bpy.context.view_layer.update()

    def extra_pairs(self):
        if self._clip_end is None or self._clip_end <= self._camera.clip_end:
            return []
        return [(Channel(self._camera.camera_data, "clip_end"), self._clip_end)]


class MoveTo(CameraMove):
    """Ease the camera to an explicit location and/or rotation (degrees)."""

    label = "move_to"

    def __init__(
        self, camera, location=None, rotation=None, run_time=None, easing=None
    ):
        super().__init__(camera, run_time, easing)
        if location is None and rotation is None:
            raise ValueError("move_to needs a location, a rotation, or both.")
        self._location, self._rotation = location, rotation

    def end_pose(self) -> Matrix:
        obj = self._camera.camera
        location = obj.location if self._location is None else Vector(self._location)
        if self._rotation is None:
            rotation = obj.rotation_euler
        else:
            rotation = Euler([math.radians(a) for a in self._rotation], "XYZ")
        return Matrix.LocRotScale(location, rotation, None)


class Dolly(CameraMove):
    """Move the camera along its own view axis; positive is towards the subject."""

    label = "dolly"

    def __init__(self, camera, distance: float, run_time=None, easing=None):
        super().__init__(camera, run_time, easing)
        self._distance = float(distance)

    def end_pose(self) -> Matrix:
        obj = self._camera.camera
        forward = Vector(self._camera.basis[2])
        matrix = obj.matrix_world.copy()
        matrix.translation = obj.matrix_world.translation + forward * self._distance
        return matrix


class Orbit(Sampled):
    """
    Swing the camera around a pivot, keeping it pointed the same way relative
    to the subject.

    The camera's world matrix is rotated about the pivot by an angle that grows
    with the eased progress, so the distance to the pivot is constant and the
    framing holds. Sampled per frame as an arc has no two-key form.
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
        self._start: Matrix | None = None
        self._pivot: Vector | None = None
        self._axis_vector: Vector | None = None
        self._previous: Euler | None = None

    def channels(self) -> list[Channel]:
        return _camera_channels(self._camera.camera)

    def prepare(self) -> None:
        obj = self._camera.camera
        bpy.context.view_layer.update()
        self._start = obj.matrix_world.copy()
        self._previous = obj.rotation_euler.copy()
        self._pivot = target_centre(self._about)
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

    def sample(self, u: float) -> Sequence[float]:
        assert self._start is not None and self._pivot is not None
        rotation = Matrix.Rotation(self._angle * u, 4, self._axis_vector)
        matrix = (
            Matrix.Translation(self._pivot)
            @ rotation
            @ Matrix.Translation(-self._pivot)
            @ self._start
        )
        values = _pose_values(matrix, self._previous)
        # keep successive eulers continuous so the curve never wraps
        self._previous = Euler(values[3:], "XYZ")
        return values

    def __repr__(self) -> str:
        return f"<camera.orbit {math.degrees(self._angle):g}deg {self.run_time:g}s {self.easing}>"


_FOCUS_EMPTY = "MN Focus"


class Focus(Clip):
    """
    Pull focus onto a target with depth of field.

    Uses Blender's own depth of field: an empty is placed at the target's centre
    and set as the camera's focus object, so the focus distance follows every
    later camera move for free. The empty's location and the f-stop are eased,
    so the focus can be pulled from one subject to another.
    """

    label = "focus"

    def __init__(
        self, camera: "Camera", target, fstop: float = 2.8, run_time=None, easing=None
    ):
        super().__init__(run_time, easing)
        self._camera = camera
        self._target = target
        self._fstop = float(fstop)

    @staticmethod
    def focus_empty(scene: bpy.types.Scene, data: bpy.types.Camera) -> bpy.types.Object:
        """
        The object the camera focuses on: whatever it already has, else the
        preset's ``focal_point`` empty, else a new one.
        """
        empty = data.dof.focus_object
        if empty is None:
            empty = scene.objects.get("focal_point")
        if empty is None:
            empty = bpy.data.objects.new(_FOCUS_EMPTY, None)
            empty.empty_display_type = "SPHERE"
            empty.empty_display_size = 0.2
            scene.collection.objects.link(empty)
        empty.hide_render = True
        return empty

    def apply(self, timeline, frame_start, frame_end, easing) -> list[Channel]:
        data = self._camera.camera_data
        empty = self.focus_empty(timeline.canvas.scene, data)
        if data.dof.focus_object is not empty:
            # first focus pull: start from the camera's current focus distance
            # along its view axis, so the pull comes from where focus was
            origin = self._camera.camera.matrix_world.translation
            empty.location = (
                origin + Vector(self._camera.basis[2]) * data.dof.focus_distance
            )
            data.dof.focus_object = empty
        data.dof.use_dof = True
        centre = target_centre(self._target)
        pairs = [(Channel(empty, "location", i), centre[i]) for i in range(3)]
        pairs.append((Channel(data.dof, "aperture_fstop"), self._fstop))
        return Tween(pairs=pairs).apply(timeline, frame_start, frame_end, easing)

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

        with canvas.timeline(fps=30) as t:
            t.play(canvas.camera.look_at(mol, viewpoint="front"), run_time=1)
            t.play(canvas.camera.orbit(90), t.tween(cartoon.i.loop_radius, 1.2), run_time=2)
            t.wait(0.5)
            t.play(canvas.camera.focus(mol.get_view("resid 40-60"), fstop=2), run_time=1)
        canvas.animation("story.mp4")

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
            cursor where it is - for layering a clip over an earlier passage.

        Raises
        ------
        ValueError
            If a clip writes a property that another clip already animates
            over overlapping frames.
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
            written = clip.apply(self, frame_start, frame_end, easing or clip.easing)
            self._claim(clip, written, frame_start, frame_end)
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

    def frames(
        self,
        entity: MolecularEntity,
        start: int,
        end: int,
        run_time: float | None = None,
        easing: str | None = None,
    ) -> Frames:
        """A clip playing universe frames ``start`` to ``end`` of an entity."""
        return Frames(entity, start, end, run_time, easing)

    # -- bookkeeping ------------------------------------------------------

    def _claim(
        self, clip: Clip, channels: Sequence[Channel], frame_start: int, frame_end: int
    ) -> None:
        for channel in channels:
            claims = self._written.setdefault(channel.key, [])
            for a, b, other in claims:
                if a < frame_end and frame_start < b:
                    raise ValueError(
                        f"{clip!r} animates {channel!r} over frames {frame_start}-{frame_end}, "
                        f"which {other!r} already animates over {a}-{b}."
                    )
            claims.append((frame_start, frame_end, clip))

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
