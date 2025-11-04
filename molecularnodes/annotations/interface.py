import bpy
import numpy as np


class AnnotationInterface:
    """
    Base class for Annotation Interface

    AnnotationManager creates a dynamic AnnotationInterface for every instance
    of an annotation. All the annotation specific inputs and the common
    annotation params are exposed though this interface for the API, GUI and the
    annotation drawing code.

    """

    # py annotations of common params for code complete
    # these are bound dynamically along with any annotation inputs
    visible: bool
    text_font: str
    text_color: np.ndarray | list | tuple
    text_size: int
    text_align: str
    text_rotation: float
    text_vspacing: float
    text_depth: bool
    text_falloff: float
    text_offset_x: float
    text_offset_y: float
    line_mode: str
    line_color: np.ndarray | list | tuple
    line_width: float
    line_arrow_size: float
    line_pointer_length: float
    mesh_wireframe: bool
    mesh_thickness: float
    mesh_color: np.ndarray | list | tuple
    mesh_material: str | bpy.types.Material
    mesh_shade_smooth: bool

    def __init__(self, instance):
        self._instance = instance
