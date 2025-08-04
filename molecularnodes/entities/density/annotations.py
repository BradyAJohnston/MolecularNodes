from ...annotations.base import BaseAnnotation
from ...annotations.manager import BaseAnnotationManager
from ..annotations import Label2D, Label3D


class DensityAnnotation(BaseAnnotation):
    """
    Base class for a Density Annotation

    All density annotations should derive from this base class and implement
    the 'draw' method. All derived classes will have access to the density
    instance (self.density) and all the annotation inputs and common params
    via self.interface.<property>

    An optional 'defaults' method can be provided to set default values
    to the annotation.

    """

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        # Auto-register any sub classes with the annotation manager
        DensityAnnotationManager.register(cls)

    def __init__(self, density):
        # Allow access to the density entity within the annotations
        self.density = density
        super().__init__()


class DensityAnnotationManager(BaseAnnotationManager):
    """
    Annotation Manager for Density Entity

    """

    _classes = {}  # Entity class specific annotation classes

    def __init__(self, entity):
        super().__init__(entity)
        self._interfaces = {}  # Entity instance specific annotation interfaces


class Label2D(DensityAnnotation, Label2D):
    """Common Label2D Annotation for all entities"""


class Label3D(DensityAnnotation, Label3D):
    """Common Label3D Annotation for all entities"""
