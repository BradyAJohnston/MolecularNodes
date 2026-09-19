from enum import StrEnum
from math import degrees, radians
from typing import Sequence
import bpy
import numpy as np
import numpy.typing as npt
from mathutils import Vector
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


class Camera:
    """
    A class to handle camera settings in Blender.

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
        """Get Camera rotation in degrees (XYZ)"""
        return tuple(degrees(angle) for angle in self.camera.rotation_euler)

    @rotation.setter
    def rotation(self, angles: tuple[float, float, float]) -> None:
        """Set Camera rotation in degrees (XYZ)"""
        self.camera.rotation_euler = tuple(radians(angle) for angle in angles)

    @property
    def basis(self) -> np.ndarray:
        """
        The camera's orthonormal basis as rows ``(right, up, forward)``.

        A Blender camera looks down its own ``-Z``, so ``forward`` is the
        negated third axis rather than the third axis itself.
        """
        matrix = self.camera.matrix_world.to_3x3().normalized()
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
        A clip swinging the camera around a pivot by ``angle`` degrees.

        The camera keeps its distance from the pivot and its orientation
        relative to it, so a framed subject stays framed. A full turntable is
        ``orbit(360, easing="linear")``.

        Parameters
        ----------
        angle : float
            Degrees to rotate through; negative reverses the direction.
        axis : str | Sequence[float], default "z"
            World axis ``"x"``, ``"y"`` or ``"z"``, the camera's own ``"up"`` or
            ``"right"``, or any vector.
        about : MolecularEntity | bpy.types.Object | array_like, optional
            The pivot, framed as for [](`~mn.Canvas.look_at`); its enclosing
            sphere's centre is used. Every entity in the scene when left out.
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
        units; positive moves towards the subject.

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
        A clip easing the camera to an explicit ``location`` and/or XYZ Euler
        ``rotation`` in degrees.

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

        Enables depth of field on the camera with an empty at the target's
        centre as the focus object, so later camera moves keep the subject in
        focus. Playing it again on another target pulls focus across.

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
            self.camera.rotation_euler = _viewpoint_rotation_eulers[
                Viewpoint(viewpoint)
            ]
        else:
            self.rotation = viewpoint
