from enum import StrEnum
from math import degrees, radians
from typing import Sequence
import bpy
import numpy as np
import numpy.typing as npt
from mathutils import Euler, Matrix, Vector
from .. import framing


class Viewpoint(StrEnum):
    DEFAULT = "default"
    FRONT = "front"
    BACK = "back"
    TOP = "top"
    BOTTOM = "bottom"
    LEFT = "left"
    RIGHT = "right"

    @classmethod
    def _missing_(cls, value: object) -> "Viewpoint | None":
        # called only when the exact-value lookup fails; allow the name or value
        # case-insensitively (e.g. "Front")
        if isinstance(value, str):
            key = value.strip().lower()
            for member in cls:
                if key in (member.value.lower(), member.name.lower()):
                    return member
        return None


_viewpoint_rotation_eulers = {
    # "default" is the camera rotation as per the template
    Viewpoint.DEFAULT: (radians(70.402), radians(0), radians(0)),
    Viewpoint.FRONT: (radians(90), radians(0), radians(0)),
    Viewpoint.BACK: (radians(90), radians(0), radians(-180)),
    Viewpoint.TOP: (radians(0), radians(0), radians(0)),
    Viewpoint.BOTTOM: (radians(-180), radians(0), radians(0)),
    Viewpoint.LEFT: (radians(-270), radians(0), radians(-90)),
    Viewpoint.RIGHT: (radians(-270), radians(0), radians(-270)),
}


#: Names of the rig the timeline builds for camera moves.
PIVOT_NAME = "camera_pivot"
TARGET_NAME = "camera_target"
FOCUS_NAME = "camera_focus"
LOOK_AT_CONSTRAINT = "MN Look At"
#: Custom property marking the rig's empties, so a user's own parent is never
#: mistaken for the pivot.
RIG_PROPERTY = "mn_camera_rig"


def world_matrix(obj: bpy.types.Object) -> Matrix:
    """
    The world matrix of an object composed from its transform channels.

    ``Object.matrix_world`` is only refreshed by a depsgraph update, so it lags
    a change to ``location`` or ``rotation_euler`` on the object or on any of
    its parents. Composing the parent chain's ``matrix_basis`` gives the
    matrix those channels currently describe. Constraints are not included.
    """
    if obj.parent is None:
        return obj.matrix_basis.copy()
    return world_matrix(obj.parent) @ obj.matrix_parent_inverse @ obj.matrix_basis


class Camera:
    """
    A class to handle camera settings in Blender.

    Camera moves played on a [](`~mn.Canvas.timeline`) rig the camera first
    (:meth:`rig`): it is parented to a pivot empty that it orbits about, a
    Damped Track constraint points it at a target empty, and a third empty is
    its depth of field focus object. Without a timeline the camera is left as
    a plain object. The ``location``, ``rotation`` and ``basis`` properties are
    world-space either way.
    """

    def __init__(self):
        # set defaults that match viewport virtual camera
        self.lens = 50
        self.clip_start = 0.01
        self.clip_end = 1000

    @property
    def camera(self) -> bpy.types.Object:
        """Get Camera object"""
        return bpy.context.scene.camera

    @property
    def camera_data(self) -> bpy.types.Camera:
        """Get Camera data"""
        return self.camera.data

    # -- rig ------------------------------------------------------------------

    @property
    def pivot(self) -> bpy.types.Object | None:
        """The empty the camera orbits about, if the camera is rigged."""
        parent = self.camera.parent
        if parent is not None and parent.get(RIG_PROPERTY) == "pivot":
            return parent
        return None

    @property
    def target(self) -> bpy.types.Object | None:
        """The empty the camera looks at and focuses on, if the camera is rigged."""
        constraint = self.camera.constraints.get(LOOK_AT_CONSTRAINT)
        return None if constraint is None else constraint.target

    @property
    def focus_point(self) -> bpy.types.Object | None:
        """The empty the camera focuses on, if the camera is rigged."""
        focus = self.camera_data.dof.focus_object
        if focus is not None and focus.get(RIG_PROPERTY) == "focus":
            return focus
        return None

    @property
    def is_rigged(self) -> bool:
        return self.pivot is not None and self.target is not None

    def rig(self, distance: float) -> None:
        """
        Rig the camera for timeline moves, without moving the view.

        The camera is parented to a pivot empty placed ``distance`` in front
        of it on its view axis, holding the viewing direction, and tracks a
        target empty at the same point. A third empty, also starting there,
        becomes the depth of field focus object. Turning the pivot orbits the
        camera, moving the target turns it, and moving the focus empty pulls
        focus without moving the view. Does nothing if the camera is already
        rigged.

        Raises
        ------
        ValueError
            If the camera already has a parent, or keyframes, which the rig
            would reinterpret.
        """
        if self.is_rigged:
            return
        camera = self.camera
        if camera.parent is not None:
            raise ValueError(
                f"The camera is parented to '{camera.parent.name}'; timeline camera "
                "moves need an unparented camera to rig."
            )
        if camera.animation_data is not None and camera.animation_data.action:
            raise ValueError(
                "The camera already has keyframes; timeline camera moves would "
                "reinterpret them. Clear them first."
            )
        location = camera.location.copy()
        rotation = camera.matrix_basis.to_euler("XYZ", camera.rotation_euler)

        scene = bpy.context.scene
        pivot = bpy.data.objects.new(PIVOT_NAME, None)
        pivot.empty_display_type = "PLAIN_AXES"
        pivot.empty_display_size = 0.5
        pivot[RIG_PROPERTY] = "pivot"
        target = bpy.data.objects.new(TARGET_NAME, None)
        target.empty_display_type = "SPHERE"
        target.empty_display_size = 0.2
        target.hide_render = True
        target[RIG_PROPERTY] = "target"
        focus = bpy.data.objects.new(FOCUS_NAME, None)
        focus.empty_display_type = "SPHERE"
        focus.empty_display_size = 0.1
        focus.hide_render = True
        focus[RIG_PROPERTY] = "focus"
        for obj in (pivot, target, focus):
            scene.collection.objects.link(obj)

        camera.parent = pivot
        camera.matrix_parent_inverse.identity()
        constraint = camera.constraints.new("DAMPED_TRACK")
        constraint.name = LOOK_AT_CONSTRAINT
        constraint.target = target
        constraint.track_axis = "TRACK_NEGATIVE_Z"
        self.camera_data.dof.focus_object = focus
        self._place(location, rotation, distance)
        focus.location = target.location

    def unrig(self) -> None:
        """
        Remove the rig, leaving a plain camera where the rigged one was.

        The animation the timeline wrote on the camera and its data goes with
        it, since the camera's keys are relative to the pivot.
        """
        if not self.is_rigged:
            return
        camera = self.camera
        empties = [self.pivot, self.target, self.focus_point]
        matrix = self.matrix_world
        camera.constraints.remove(camera.constraints[LOOK_AT_CONSTRAINT])
        camera.parent = None
        camera.animation_data_clear()
        self.camera_data.animation_data_clear()
        camera.matrix_basis = matrix
        for obj in empties:
            if obj is not None:
                bpy.data.objects.remove(obj, do_unlink=True)

    def _place(
        self, location: Sequence[float], rotation: Euler, distance: float
    ) -> None:
        """
        Pose the rig: the camera at a world ``location`` with a world
        ``rotation``, and the pivot and target together ``distance`` in front
        of it, so the camera already looks at the target.
        """
        pivot, target = self.pivot, self.target
        assert pivot is not None and target is not None
        rotation = Euler(rotation, "XYZ")
        rotation.make_compatible(pivot.rotation_euler)
        forward = rotation.to_matrix() @ Vector((0.0, 0.0, -1.0))
        centre = Vector(location) + forward * distance
        pivot.location = centre
        pivot.rotation_euler = rotation
        target.location = centre
        self.camera.location = (0.0, 0.0, distance)
        self.camera.rotation_euler = (0.0, 0.0, 0.0)

    @property
    def _distance(self) -> float:
        """How far the rigged camera is from its target."""
        target = self.target
        assert target is not None
        return (world_matrix(target).translation - self.location).length

    @property
    def parent_matrix(self) -> Matrix:
        """The world matrix the camera's own transform channels are relative to."""
        camera = self.camera
        if camera.parent is None:
            return Matrix.Identity(4)
        return world_matrix(camera.parent) @ camera.matrix_parent_inverse

    @property
    def matrix_world(self) -> Matrix:
        """
        The camera's world matrix, current with its transform channels and the
        rig's look-at constraint.
        """
        matrix = world_matrix(self.camera)
        target = self.target
        if target is None:
            return matrix
        location = matrix.translation
        direction = world_matrix(target).translation - location
        if direction.length < 1e-9:
            return matrix
        # the Damped Track constraint: the shortest turn that points -Z at the
        # target
        rotation = matrix.to_3x3().normalized()
        forward = rotation @ Vector((0.0, 0.0, -1.0))
        turn = forward.rotation_difference(direction.normalized()).to_matrix()
        return Matrix.LocRotScale(location, turn @ rotation, matrix.to_scale())

    @property
    def location(self) -> Vector:
        """Get the camera's world-space location"""
        return self.matrix_world.translation

    @location.setter
    def location(self, value: Sequence[float]) -> None:
        """
        Set the camera's world-space location, keeping the way it points. A
        rigged camera takes its pivot and target along with it.
        """
        if not self.is_rigged:
            self.camera.location = value
            return
        rotation = self.matrix_world.to_euler("XYZ", self.pivot.rotation_euler)
        self._place(value, rotation, self._distance)

    @property
    def lens(self) -> float:
        """Get Camera focal length"""
        return self.camera_data.lens

    @lens.setter
    def lens(self, value) -> None:
        """Set Camera focal length"""
        self.camera_data.lens = value

    @property
    def clip_start(self) -> float:
        """Get Camera near clipping distance"""
        return self.camera_data.clip_start

    @clip_start.setter
    def clip_start(self, value) -> None:
        """Set Camera near clipping distance"""
        self.camera_data.clip_start = value

    @property
    def clip_end(self) -> float:
        """Get Camera far clipping distance"""
        return self.camera_data.clip_end

    @clip_end.setter
    def clip_end(self, value) -> None:
        """Set Camera far clipping distance"""
        self.camera_data.clip_end = value

    @property
    def rotation(self) -> tuple[float, float, float]:
        """Get the camera's world rotation in degrees (XYZ)"""
        if not self.is_rigged:
            return tuple(degrees(angle) for angle in self.camera.rotation_euler)
        euler = self.matrix_world.to_euler("XYZ", self.pivot.rotation_euler)
        return tuple(degrees(angle) for angle in euler)

    @rotation.setter
    def rotation(self, angles: tuple[float, float, float]) -> None:
        """
        Set the camera's world rotation in degrees (XYZ).

        A rigged camera keeps its location; the rotation goes on the pivot, and
        the pivot and target move onto the new view axis at the same distance.
        """
        rotation = Euler(tuple(radians(angle) for angle in angles), "XYZ")
        if not self.is_rigged:
            self.camera.rotation_euler = rotation
            return
        self._place(self.location, rotation, self._distance)

    @property
    def basis(self) -> np.ndarray:
        """
        The camera's orthonormal basis as rows ``(right, up, forward)``.

        A Blender camera looks down its own ``-Z``, so ``forward`` is the
        negated third axis rather than the third axis itself.
        """
        matrix = self.matrix_world.to_3x3().normalized()
        return np.array(
            [
                matrix @ Vector((1.0, 0.0, 0.0)),
                matrix @ Vector((0.0, 1.0, 0.0)),
                matrix @ Vector((0.0, 0.0, -1.0)),
            ],
            dtype=np.float64,
        )

    def frame_bounds(
        self, scene: bpy.types.Scene | None = None
    ) -> tuple[float, float, float, float]:
        """
        The edges of what the camera sees, as ``(left, right, bottom, top)``.

        For a perspective camera these are ratios of offset to depth, so a point
        at depth ``d`` is in frame when ``left * d <= x <= right * d``. For an
        orthographic camera they are world-space offsets at its current scale.

        Taken from Blender's own view frame, so the sensor fit, render aspect
        ratio, pixel aspect and any lens shift are all accounted for - a shifted
        camera gives asymmetric bounds, which the framing solve handles.
        """
        if scene is None:
            scene = bpy.context.scene
        corners = self.camera_data.view_frame(scene=scene)
        xs = [corner.x for corner in corners]
        ys = [corner.y for corner in corners]
        if self.camera_data.type == "ORTHO":
            return (min(xs), max(xs), min(ys), max(ys))
        # every corner sits at the same depth, so one of them sets the scale
        depth = -corners[0].z
        return (min(xs) / depth, max(xs) / depth, min(ys) / depth, max(ys) / depth)

    def frame_points(
        self,
        points: npt.ArrayLike,
        margin: float = 0.05,
        scene: bpy.types.Scene | None = None,
    ) -> None:
        """
        Move the camera so that every one of these points is in frame.

        Solves for the closest position that still contains the points, without
        changing where the camera is pointing. See
        [](`molecularnodes.framing.fit_camera_to_points`) for the solve.

        A rigged camera's pivot and target move to the depth of the points'
        centre on the view axis, so an orbit that follows turns about what was
        framed.

        Parameters
        ----------
        points : array_like
            ``(N, 3)`` world-space positions to fit into the frame. Any number
            of points is fine.
        margin : float, default 0.05
            Fraction of the frame to leave empty around the subject. The small
            default keeps the subject off the edge of the frame; ``0`` fits it
            exactly, and a negative value crops in tighter.
        scene : bpy.types.Scene, optional
            Scene to read the render aspect ratio from. Defaults to the active
            scene.
        """
        if scene is None:
            scene = bpy.context.scene
        points = framing.as_points(points)
        bounds = self.frame_bounds(scene)
        basis = self.basis

        if self.camera_data.type == "ORTHO":
            location, scale_factor = framing.fit_orthographic_to_points(
                points, basis, bounds, margin
            )
            self.camera_data.ortho_scale *= scale_factor
        else:
            location = framing.fit_camera_to_points(points, basis, bounds, margin)

        if self.is_rigged:
            centre, _ = framing.enclosing_sphere(points)
            depth = float(np.dot(centre - location, basis[2]))
            rotation = self.matrix_world.to_euler("XYZ", self.pivot.rotation_euler)
            self._place(location, rotation, max(depth, self.clip_start))
        else:
            self.camera.location = location
        # a subject sitting beyond the far clip renders as nothing at all, so
        # make room for it rather than silently dropping it
        furthest = float(np.max((points - location) @ basis[2]))
        if furthest > self.clip_end:
            self.clip_end = furthest * 1.05

    # -- timeline clips -----------------------------------------------------
    # each returns a clip describing a camera move for `Canvas.timeline`; the
    # camera is not touched until the clip is played

    def look_at(
        self,
        target,
        viewpoint: Viewpoint | str | Sequence[float] | None = None,
        margin: float = 0.05,
        run_time: float | None = None,
        easing: str | None = None,
    ):
        """
        A clip easing the camera to the framing [](`~mn.Canvas.look_at`) would
        jump to. The pose is solved when the clip is played, against the
        geometry as it is at that frame.

        Parameters
        ----------
        target : MolecularEntity | bpy.types.Object | array_like
            What to frame, as for [](`~mn.Canvas.look_at`).
        viewpoint : Viewpoint | str | Sequence[float], optional
            Viewing direction to move to; the current one when left out.
        margin : float, default 0.05
            Fraction of the frame to leave empty around the target.
        run_time : float, optional
            Length of the move in seconds (default 1).
        easing : str, optional
            Rate function, one of ``"smooth"`` (default), ``"linear"``,
            ``"sine"``, ``"ease_in"``, ``"ease_out"``.

        Returns
        -------
        molecularnodes.scene.timeline.LookAt
        """
        from .timeline import LookAt

        return LookAt(self, target, viewpoint, margin, run_time, easing)

    def orbit(
        self,
        angle: float,
        axis: str | Sequence[float] = "z",
        about=None,
        run_time: float | None = None,
        easing: str | None = None,
    ):
        """
        A clip turning the camera's pivot by ``angle`` degrees.

        The camera keeps its distance from the pivot and its orientation
        relative to it, so a framed subject stays framed. A full turntable is
        ``orbit(360, easing="linear")``. Turning about the world ``z`` axis or
        the camera's own ``right`` axis is one pair of keys on the pivot's
        rotation, editable in the Graph Editor; any other axis is sampled per
        frame.

        Parameters
        ----------
        angle : float
            Degrees to rotate through; negative reverses the direction.
        axis : str | Sequence[float], default "z"
            World axis ``"x"``, ``"y"`` or ``"z"``, the camera's own ``"up"`` or
            ``"right"``, or any vector.
        about : MolecularEntity | bpy.types.Object | array_like, optional
            What to turn about, as for [](`~mn.Canvas.look_at`). The pivot
            and the target ease to its centre over the clip while turning, so
            the camera also turns to face it. Left out, the pivot stays where
            the last framing put it.
        run_time : float, optional
            Length in seconds (default 2).
        easing : str, optional
            Rate function (default ``"smooth"``).

        Returns
        -------
        molecularnodes.scene.timeline.Orbit
        """
        from .timeline import Orbit

        return Orbit(self, angle, axis, about, run_time, easing)

    def dolly(
        self, distance: float, run_time: float | None = None, easing: str | None = None
    ):
        """
        A clip moving the camera along its view axis by ``distance`` world
        units; positive moves towards the target. Keys the camera's own
        location, so it composes with an orbit of the pivot.

        Returns
        -------
        molecularnodes.scene.timeline.Dolly
        """
        from .timeline import Dolly

        return Dolly(self, distance, run_time, easing)

    def zoom(
        self, lens: float, run_time: float | None = None, easing: str | None = None
    ):
        """
        A clip easing the focal length to ``lens`` millimetres.

        Returns
        -------
        molecularnodes.scene.timeline.Tween
        """
        from .timeline import Tween

        return Tween((self.camera_data, "lens"), lens, run_time, easing)

    def move_to(
        self,
        location: Sequence[float] | None = None,
        rotation: Sequence[float] | None = None,
        run_time: float | None = None,
        easing: str | None = None,
    ):
        """
        A clip easing the camera to an explicit world ``location`` and/or XYZ
        Euler ``rotation`` in degrees, as setting :attr:`location` and
        :attr:`rotation` would: given only a rotation, the camera turns in
        place; given only a location, it moves without turning. The pivot and
        target travel with it.

        Returns
        -------
        molecularnodes.scene.timeline.MoveTo
        """
        from .timeline import MoveTo

        return MoveTo(self, location, rotation, run_time, easing)

    def focus(
        self,
        target,
        fstop: float = 2.8,
        run_time: float | None = None,
        easing: str | None = None,
    ):
        """
        A clip pulling focus onto ``target`` with depth of field.

        Moves the rig's focus empty to the centre of ``target`` without
        moving the view, so later camera moves keep it in focus. If depth of
        field was off it
        switches on as the clip starts, with the aperture opening up from fully
        stopped down. Playing it again on another target pulls focus across.

        Parameters
        ----------
        target : MolecularEntity | bpy.types.Object | array_like
            What to focus on, e.g. ``mol.get_view("resid 40-60")``.
        fstop : float, default 2.8
            Aperture; smaller is shallower.
        run_time : float, optional
            Length in seconds (default 1).
        easing : str, optional
            Rate function (default ``"smooth"``).

        Returns
        -------
        molecularnodes.scene.timeline.Focus
        """
        from .timeline import Focus

        return Focus(self, target, fstop, run_time, easing)

    def set_viewpoint(self, viewpoint: Viewpoint | str | Sequence[float]) -> None:
        """
        Set viewpoint to a preset or a custom Euler rotation.

        Parameters
        ----------
        viewpoint : Viewpoint | str | Sequence[float]
            Either a named viewpoint (a ``Viewpoint`` or its name, e.g. "front",
            "top") or a tuple/list of three Euler angles in degrees (XYZ),
            matching :attr:`rotation`.
        """
        # Viewpoint is a StrEnum, so named viewpoints (including bare strings) are
        # caught here; a Sequence[float] of Euler angles falls through
        if isinstance(viewpoint, str):
            euler = _viewpoint_rotation_eulers[Viewpoint(viewpoint)]
            self.rotation = tuple(degrees(angle) for angle in euler)
        else:
            self.rotation = viewpoint
