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
    text_color: np.ndarray | list | tuple
    text_font: str
    text_size: int
    text_align: str
    text_rotation: float
    text_vspacing: float
    text_depth: bool
    text_falloff: float
    offset_x: float
    offset_y: float
    line_color: np.ndarray | list | tuple
    line_width: float
    arrow_size: int
    pointer_length: int

    def __init__(self, instance):
        self._instance = instance
