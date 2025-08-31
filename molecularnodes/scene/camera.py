from math import degrees, radians
from typing import Literal, Sequence, get_args
import bpy

Viewpoints = Literal["default", "front", "back", "top", "bottom", "left", "right"]

_viewpoint_rotation_eulers = {
    # "default" is the camera rotation as per the template
    "default": (radians(70.402), radians(0), radians(0)),
    "front": (radians(90), radians(0), radians(0)),
    "back": (radians(90), radians(0), radians(-180)),
    "top": (radians(0), radians(0), radians(0)),
    "bottom": (radians(-180), radians(0), radians(0)),
    "left": (radians(-270), radians(0), radians(-90)),
    "right": (radians(-270), radians(0), radians(-270)),
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

    def set_viewpoint(self, viewpoint: Viewpoints | Sequence[float]) -> None:
        """
        Set viewpoint to a preset or a custom Euler rotation.

        Parameters
        ----------
        viewpoint : Viewpoints | Sequence[float]
            Either a named viewpoint string (e.g. "front", "top") or a tuple/list of three Euler angles in radians.
        """
        if isinstance(viewpoint, str):
            if viewpoint not in get_args(Viewpoints):
                raise ValueError(f"{viewpoint} is not one of {get_args(Viewpoints)}")
            self.camera.rotation_euler = _viewpoint_rotation_eulers[viewpoint]
        else:
            self.camera.rotation_euler = viewpoint
