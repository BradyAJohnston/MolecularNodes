from enum import StrEnum
from math import degrees, radians
from typing import Sequence
import bpy


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
    def camera(self) -> bpy.types.Camera:
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

    def set_viewpoint(self, viewpoint: Viewpoint | str | Sequence[float]) -> None:
        """
        Set viewpoint to a preset or a custom Euler rotation.

        Parameters
        ----------
        viewpoint : Viewpoint | str | Sequence[float]
            Either a named viewpoint (a ``Viewpoint`` or its name, e.g. "front",
            "top") or a tuple/list of three Euler angles in radians.
        """
        # Viewpoint is a StrEnum, so named viewpoints (including bare strings) are
        # caught here; a Sequence[float] of Euler angles falls through
        if isinstance(viewpoint, str):
            self.camera.rotation_euler = _viewpoint_rotation_eulers[Viewpoint(viewpoint)]
        else:
            self.camera.rotation_euler = viewpoint
