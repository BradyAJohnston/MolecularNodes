import numpy as np


class AnnotationInterface:
    """
    Base class for Annotation Interface

    AnnotationManager creates a dynamic AnnotationInterface for every instance
    of an annotation. All the annotation specific inputs and the common
    annotation params are exposed though this interface for the API, GUI and the
    annotation drawing code.

    """

    def __init__(self, instance):
        self._instance = instance
        # Common annotation params set here to aid code auto-complete
        # AnnotationManager resets these as dynamic properties that reference
        # the corresponding Annotation node (Geometry Node) inputs
        # Tests will ensure that these attributes match the Annotation node input names
        self.visible: bool = False
        self.text_color: np.ndarray | list | tuple = (1, 1, 1, 1)
        self.text_size: int = 16
        self.text_rotation: float = 0.0
        self.offset_x: float = 0.0
        self.offset_y: float = 0.0
        self.line_color: np.ndarray | list | tuple = (1, 1, 1, 1)
        self.line_width: float = 2.0
        self.arrow_size: int = 16
        self.pointer_length: int = 0
