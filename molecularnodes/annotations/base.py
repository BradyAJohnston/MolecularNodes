from abc import ABCMeta, abstractmethod
from math import cos, radians, sin
import blf
import bpy
import gpu
from bpy_extras import view3d_utils
from gpu_extras.batch import batch_for_shader
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
    viewport_width: int
        Width of the viewport region in pixels
    viewport_height: int
        Height of the viewport region in pixels

    """

    name: str = None
    interface: AnnotationInterface

    def __init__(self):
        self._world_scale = 0.01
        self._rad45 = radians(45)
        self._rad315 = radians(315)
        self._shader_line = (
            gpu.shader.from_builtin("POLYLINE_UNIFORM_COLOR")
            if not bpy.app.background
            else None
        )

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

    @property
    def viewport_width(self) -> int:
        """Get the viewport region width in pixels"""
        return self._region.width

    @property
    def viewport_height(self) -> int:
        """Get the viewport region height in pixels"""
        return self._region.height

    def draw_text_3d(self, pos_3d: Vector, text: str) -> None:
        """
        Draw text at a given 3D position

        Parameters
        ----------
        pos_3d: Vector
            Co-ordinates in 3D world space (x, y, z)

        text: str
            Text to display. '|' as multi-line separator

        """
        params = self.interface
        text_pos = pos_3d
        if params.pointer_length > 0:
            # draw a pointer to the 3D position
            nv = -Vector(pos_3d)
            nv.normalize()
            pointer_begin = pos_3d + (nv * params.pointer_length)
            text_pos = pointer_begin + Vector((0, 0, 0.5))  # offset text a bit
            self.draw_line_3d(pointer_begin, pos_3d, v2_arrow=True)

        pos_2d = self._get_2d_point(text_pos)
        if pos_2d is None:
            return
        self.draw_text_2d(pos_2d, text)

    def draw_text_2d_norm(self, pos_2d: Vector, text: str) -> None:
        """
        Draw text at a given 2D position (normalized co-ordinates) of Viewport.

        Parameters
        ----------
        pos_2d: Vector
            Normalized co-ordinates (0 - 1). (0, 0) is at bottom left

        text: str
            Text to display. '|' as multi-line separator

        """
        if pos_2d is None:
            return
        for comp in pos_2d:
            if not (0.0 <= comp <= 1.0):
                return
        pos_x, pos_y = pos_2d
        new_x = pos_x * self.viewport_width
        new_y = pos_y * self.viewport_height
        self._draw_text_2d((new_x, new_y), text)

    def draw_text_2d(self, pos_2d: Vector, text: str) -> None:
        """
        Draw text at a given 2D position (in pixels) of Viewport.

        Parameters
        ----------
        pos_2d: Vector
            Co-ordinates in pixels. (0, 0) is at bottom left

        text: str
            Text to display. '|' as multi-line separator

        """
        if pos_2d is None:
            return
        font_id = 0
        scale = 1
        right_alignment_gap = 12
        params = self.interface
        pos_x, pos_y = pos_2d
        # offset_x and offset_y
        pos_x += params.offset_x
        pos_y += params.offset_y
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

    def draw_line_3d(
        self,
        v1: Vector,
        v2: Vector,
        v1_text: str = None,
        v2_text: str = None,
        mid_text: str = None,
        v1_arrow: bool = False,
        v2_arrow: bool = False,
    ) -> None:
        """
        Draw a line between two points in 3D space

        Parameters
        ----------
        v1: Vector
            3D co-ordinates of the first point

        v2: Vector
            3D co-ordinates of the second point

        v1_text: str
            Optional text to display at v1

        v2_text: str
            Optional text to display at v2

        mid_text: str
            Optional text to display at the middle of the line

        v1_arrow: bool
            Whether to display an arrow at v1

        v2_arrow: bool
            Whether to display an arrow at v2

        """
        if v1 is None or v2 is None:
            return
        v1_2d = self._get_2d_point(v1)
        v2_2d = self._get_2d_point(v2)
        self.draw_line_2d(v1_2d, v2_2d, v1_text, v2_text, mid_text, v1_arrow, v2_arrow)

    def draw_line_2d(
        self,
        v1: Vector,
        v2: Vector,
        v1_text: str = None,
        v2_text: str = None,
        mid_text: str = None,
        v1_arrow: bool = False,
        v2_arrow: bool = False,
    ) -> None:
        """
        Draw a line between two points in 2D viewport space

        Parameters
        ----------
        v1: Vector
            2D co-ordinates of the first point

        v2: Vector
            2D co-ordinates of the second point

        v1_text: str
            Optional text to display at v1

        v2_text: str
            Optional text to display at v2

        mid_text: str
            Optional text to display at the middle of the line

        v1_arrow: bool
            Whether to display an arrow at v1

        v2_arrow: bool
            Whether to display an arrow at v2

        """
        if v1 is None or v2 is None:
            return
        self._draw_arrow_line_2d(v1, v2, v1_arrow, v2_arrow)
        if v1_text is not None:
            self.draw_text_2d(v1, v1_text)
        if v2_text is not None:
            self.draw_text_2d(v2, v2_text)
        if mid_text is not None:
            mid = (v1 + v2) / 2
            self.draw_text_2d(mid, mid_text)

    def distance(self, v1: Vector, v2: Vector) -> float:
        """
        Distance between two vectors

        Paramaters
        ----------
        v1: Vector
            A 3D or 2D vector or tuple

        v2: Vector
            A 3D or 2D vector or tuple

        Returns
        -------
        Distance between the two vectors

        """
        return (Vector(v2) - Vector(v1)).length

    def _draw_arrow_line_2d(
        self,
        v1: Vector,
        v2: Vector,
        v1_arrow: bool = False,
        v2_arrow: bool = False,
    ) -> None:
        """Internal: Draw a line between two 2D points with arrows"""
        if v1 is None or v2 is None:
            return
        # actual line
        self._draw_line_2d(v1, v2)
        if v1_arrow:
            # v1 arrow lines
            va, vb = self._get_arrow_end_points(v1, v2)
            self._draw_line_2d(v1, va)
            self._draw_line_2d(v1, vb)
        if v2_arrow:
            # v2 arrow lines
            va, vb = self._get_arrow_end_points(v2, v1)
            self._draw_line_2d(v2, va)
            self._draw_line_2d(v2, vb)

    def _get_arrow_end_points(self, v1, v2: Vector) -> tuple:
        """Internal: Get arrow end point positions"""
        params = self.interface
        v = self._interpolate_3d(
            (v1[0], v1[1], 0.0), (v2[0], v2[1], 0.0), params.arrow_size
        )
        vi = (v[0] - v1[0], v[1] - v1[1])
        va = (
            int(vi[0] * cos(self._rad45) - vi[1] * sin(self._rad45) + v1[0]),
            int(vi[1] * cos(self._rad45) + vi[0] * sin(self._rad45)) + v1[1],
        )
        vb = (
            int(vi[0] * cos(self._rad315) - vi[1] * sin(self._rad315) + v1[0]),
            int(vi[1] * cos(self._rad315) + vi[0] * sin(self._rad315) + v1[1]),
        )
        return (va, vb)

    def _draw_line_2d(self, v1: Vector, v2: Vector) -> None:
        """Internal: Draw a line between two 2D points"""
        if v1 is None or v2 is None:
            return
        params = self.interface
        rgba = params.line_color
        line_width = params.line_width
        viewport_size = (self.viewport_width, self.viewport_height)
        coords = [(v1[0], v1[1], 0), (v2[0], v2[1], 0)]
        gpu.state.blend_set("ALPHA")
        batch = batch_for_shader(self._shader_line, "LINES", {"pos": coords})
        try:
            self._shader_line.bind()
            self._shader_line.uniform_float("color", rgba)
            self._shader_line.uniform_float("lineWidth", line_width)
            self._shader_line.uniform_float("viewportSize", viewport_size)
            batch.draw(self._shader_line)
        except Exception as e:
            print(e)

    def _interpolate_3d(self, v1: Vector, v2: Vector, d1: float) -> Vector:
        """Internal: Interpolated 3D point between two 3D points at distance d1"""
        # calculate displacement vector
        v = Vector(v2) - Vector(v1)
        # calculate distance between points
        d0 = self.distance(v1, v2)
        # calculate interpolate factor (distance from origin / distance total)
        # if d1 > d0, the point is projected in 3D space
        if d0 > 0:
            x = d1 / d0
        else:
            x = d1
        return (v1[0] + (v[0] * x), v1[1] + (v[1] * x), v1[2] + (v[2] * x))

    def _set_viewport_region(
        self,
        region: bpy.types.Region | None,
        rv3d: bpy.types.RegionView3D | None,
    ) -> None:
        """Internal: Set the 3D viewport region and region data"""
        self._region = region
        self._rv3d = rv3d

    def _get_2d_point(self, pos_3d: Vector) -> Vector | None:
        """Internal: Get the 2D point in region corresponding to 3D position"""
        if self._rv3d is not None and self._region is not None:
            pos_3d_scaled = pos_3d * self._world_scale
            return view3d_utils.location_3d_to_region_2d(
                self._region, self._rv3d, pos_3d_scaled
            )
        else:
            return None
