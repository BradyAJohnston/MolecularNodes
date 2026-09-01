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
