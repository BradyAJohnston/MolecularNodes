import os
import numpy as np
from ...annotations.base import BaseAnnotation
from ...annotations.manager import BaseAnnotationManager
from ..annotations import Label2D, Label3D
from ..base import EntityType


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

    _entity_type = EntityType.DENSITY
    _classes = {}  # Entity class specific annotation classes

    def __init__(self, entity):
        super().__init__(entity)
        self._interfaces = {}  # Entity instance specific annotation interfaces


class DensityInfo(DensityAnnotation):
    """
    Density Info Annotation

    Attributes
    ----------
    location: tuple[float, float]
        Normalized coordinates (0.0 - 1.0) of the postion in viewport / render

    show_filename: bool
        Whether or not to show the grid filename

    show_threshold: bool
        Whether or not to show the current threshold value

    show_origin: bool
        Whether or not to show the grid origin

    show_delta: bool
        Whether or not to show the grid delta

    show_shape: bool
        Whether or not to show the grid shape

    custom_text: str
        Any custom text to add at the end of the annotation

    """

    annotation_type = "density_info"

    location: tuple[float, float] = (0.025, 0.05)
    show_filename: bool = True
    show_threshold: bool = True
    show_origin: bool = True
    show_delta: bool = True
    show_shape: bool = True
    custom_text: str = ""

    def defaults(self) -> None:
        params = self.interface
        params.text_align = "left"

    def validate(self) -> bool:
        params = self.interface
        x, y = params.location
        if (not 0 <= x <= 1) or (not 0 <= y <= 1):
            raise ValueError("Normalized coordinates should lie between 0 and 1")
        return True

    def draw(self) -> None:
        params = self.interface
        grid = self.density.grid
        text = ""
        if params.show_filename:
            filename = os.path.basename(grid.metadata["filepath"])
            text = f"Filename: {filename}"
        if params.show_threshold:
            if self.density.node_group:
                iso_value = getattr(self.density.styles[0], "iso_value", None)
                if iso_value is None:
                    threshold = self.density.styles[0].threshold
                    text = text + f"|Threshold: {threshold:.2f}"
                else:
                    text = text + f"|ISO Value: {iso_value:.2f}"
        if params.show_origin:
            text = text + f"|Origin: {np.round(grid.origin, 2)}"
        if params.show_delta:
            text = text + f"|Delta: {np.round(grid.delta, 2)}"
        if params.show_shape:
            text = text + f"|Shape: {np.round(grid.grid.shape, 2)}"
        if params.custom_text != "":
            text = text + "|" + params.custom_text
        # Draw text at normalized coordinates wrt viewport / render
        self.draw_text_2d_norm(params.location, text)


class DensityGridAxes(DensityAnnotation):
    """
    Density Grid Axes Annotation

    Attributes
    ----------
    show_length: bool
        Whether or not to show the length of the grid axes

    units: str
        Units to use for length. Default: Å

    """

    annotation_type = "grid_axes"

    show_length: bool = True
    units: str = "Å"

    def defaults(self) -> None:
        params = self.interface
        params.text_depth = False
        params.arrow_size = 10
        params.text_color = (1.0, 1.0, 1.0, 0.5)
        params.line_color = (1.0, 1.0, 1.0, 0.5)

    def draw(self) -> None:
        params = self.interface
        grid = self.density.grid
        if grid.origin.size != 3:
            return
        origin = grid.origin.copy()
        if grid.metadata["center"]:
            origin = -np.array(grid.grid.shape) * 0.5 * grid.delta
        axes = ["X", "Y", "Z"]
        for i in range(3):
            length = grid.grid.shape[i] * grid.delta[i]
            mid_text = None
            if params.show_length:
                mid_text = f"{length:.2f} {params.units}"
            end = origin.copy()
            end[i] += length
            self.draw_line_3d(
                v1=origin, v2=end, mid_text=mid_text, v2_text=axes[i], v2_arrow=True
            )


class Label2D(DensityAnnotation, Label2D):
    """Common Label2D Annotation for all entities"""


class Label3D(DensityAnnotation, Label3D):
    """Common Label3D Annotation for all entities"""
