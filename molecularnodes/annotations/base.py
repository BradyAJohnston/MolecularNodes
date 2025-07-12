from abc import ABCMeta, abstractmethod
from math import radians
import blf
import bpy
from bpy_extras import view3d_utils
from mathutils import Vector
from .interface import AnnotationInterface


class BaseAnnotation(metaclass=ABCMeta):
    """
    Base class for an Annotation

    Any Entity that needs annotations support can derive from this base class
    for Entity specific annotations. The derived entity annotation  will have
    to implement '__init_subclass__' to register with the Entity's annotation
    manager and '__init__' to pass the entity to annotations.

    Entity annotations will have to implement the 'draw' method that specifies
    how to display the annotations An optional 'validate' method can be provided
    to validate annotation inputs An optional 'defaults' method can be provided
    to set default values to the annotation.

    Attributes
    ----------
    name: str
        Name (label) of the annotation
    interface: AnnotationInterface
        Dynamic interface of the annotation instance

    """

    name: str = None
    interface: AnnotationInterface

    def validate(self) -> bool:
        """
        Optional method to validate annotation inputs
        This is called during annotation creation and any time the inputs change
        either through the API or GUI. Can return False or raise an exception
        when validation fails. Returns True when all validations succeed

        """
        return True

    def defaults(self) -> None:
        """
        Optional method to set default annotation params
        This is called only once when the annotation instance is created

        """

    @abstractmethod
    def draw(self) -> None:
        """
        The main draw method for an annotation
        This is called multiple times in the 3D viewport draw handler

        """
        raise NotImplementedError("Subclasses must implement this method")

    def draw_text_3d(self, pos_3d: Vector, text: str) -> None:
        """
        Draw text at a given 3D position

        """
        pos_2d = self._get_2d_point(pos_3d)
        if pos_2d is None:
            return
        self.draw_text_2d(pos_2d, text)

    def draw_text_2d(self, pos_2d: Vector, text: str) -> None:
        """
        Draw text at a given 2D position of Viewport. (0, 0) is bottom left

        """
        if pos_2d is None:
            return
        font_id = 0
        scale = 1
        right_alignment_gap = 12
        params = self.interface
        pos_x, pos_y = pos_2d
        rgba = params.text_color
        # set the text size
        blf.size(font_id, params.text_size)
        # height of one line - for use in multiline text
        _, max_th = blf.dimensions(font_id, "Tp")  # uses high/low letters
        # split lines
        lines = text.split("|")
        n_lines = len(lines) - 1
        # draw all text lines
        for line in lines:
            line_width, _ = blf.dimensions(font_id, line)
            # x position based on alignment
            if params.text_align == "center":
                new_x = (pos_x * scale) - line_width / 2
            elif params.text_align == "right":
                new_x = (pos_x * scale) - line_width - right_alignment_gap
            else:
                new_x = pos_x * scale
                # rotation support only for left alignment
                blf.enable(font_id, blf.ROTATION)
                blf.rotation(font_id, radians(params.text_rotation))
            # calculate new y position
            new_y = (pos_y * scale) + (max_th * n_lines)
            # draw the text line
            blf.position(font_id, new_x, new_y, 0)
            blf.color(font_id, rgba[0], rgba[1], rgba[2], rgba[3])
            blf.draw(font_id, " " + line)
            n_lines -= 1  # for next iteration
            if params.text_align == "left":
                blf.disable(font_id, blf.ROTATION)

    def _set_viewport_region(
        self,
        region: bpy.types.Region | None,
        rv3d: bpy.types.RegionView3D | None,
    ) -> None:
        """Set the 3D viewport region and region data"""
        self._region = region
        self._rv3d = rv3d

    def _get_2d_point(self, pos_3d: Vector) -> Vector | None:
        """Get the 2D point in region corresponding to 3D position"""
        if self._rv3d is not None and self._region is not None:
            pos_3d_scaled = pos_3d * 0.01  # world scale for MN
            return view3d_utils.location_3d_to_region_2d(
                self._region, self._rv3d, pos_3d_scaled
            )
        else:
            return None
