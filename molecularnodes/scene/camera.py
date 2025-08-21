import math
from typing import Literal, get_args
import bpy

viewpoints = Literal["default", "front", "back", "top", "bottom", "left", "right"]
viewpoint_rotation_eulers = {
    # "default" is the camera rotation as per the template
    "default": (math.radians(70.402), math.radians(0), math.radians(0)),
    "front": (math.radians(90), math.radians(0), math.radians(0)),
    "back": (math.radians(90), math.radians(0), math.radians(-180)),
    "top": (math.radians(0), math.radians(0), math.radians(0)),
    "bottom": (math.radians(-180), math.radians(0), math.radians(0)),
    "left": (math.radians(-270), math.radians(0), math.radians(-90)),
    "right": (math.radians(-270), math.radians(0), math.radians(-270)),
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
        """Camera object"""
        return bpy.context.scene.camera

    @property
    def camera_data(self) -> bpy.types.Camera:
        "Camera data"
        return self.camera.data

    @property
    def lens(self) -> float:
        """Camera focal length"""
        return self.camera_data.lens

    @lens.setter
    def lens(self, value) -> None:
        """Camera focal length"""
        self.camera_data.lens = value

    @property
    def clip_start(self) -> float:
        """Camera near clipping distance"""
        return self.camera_data.clip_start

    @clip_start.setter
    def clip_start(self, value) -> None:
        """Camera near clipping distance"""
        self.camera_data.clip_start = value

    @property
    def clip_end(self) -> float:
        """Camera far clipping distance"""
        return self.camera_data.clip_end

    @clip_end.setter
    def clip_end(self, value) -> None:
        """Camera far clipping distance"""
        self.camera_data.clip_end = value

    def set_viewpoint(self, viewpoint: viewpoints) -> None:
        """Set viewpoint to a preset"""
        if viewpoint not in get_args(viewpoints):
            raise ValueError(f"{viewpoint} is not one of {get_args(viewpoints)}")
        self.camera.rotation_euler = viewpoint_rotation_eulers[viewpoint]
