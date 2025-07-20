from abc import ABCMeta, abstractmethod
from math import cos, radians, sin
import blf
import bpy
import gpu
from bpy_extras import view3d_utils
from gpu_extras.batch import batch_for_shader
from mathutils import Matrix, Vector
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
        when validation fails. Returns True when all validations succeed.

        Note: This method gets called when any inputs change, so updating values
        in here will lead to a recursive loop and should not be done.

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
        self._draw_text(pos_3d, text, is3d=True)

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
        self._draw_text((new_x, new_y), text, is3d=False)

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
        self._draw_text(pos_2d, text, is3d=False)

    def _draw_text(self, pos: Vector, text: str, is3d: bool = False) -> None:
        """Internal: Draw text 3D or 2D"""
        if pos is None:
            return
        pos_2d = pos
        if is3d:
            text_pos = pos
            params = self.interface
            # check if pointer is required
            if params.pointer_length > 0:
                # draw a pointer to the 3D position
                nv = -Vector(pos)
                nv.normalize()
                pointer_begin = pos + (nv * params.pointer_length)
                text_pos = pointer_begin + (nv * 0.5)  # offset text a bit
                self._draw_line(pointer_begin, pos, v2_arrow=True, is3d=True)
            pos_2d = self._get_2d_point(text_pos)
            if pos_2d is None:
                return
        # draw text at 2D position
        font_id = 0
        scale = 1
        right_alignment_gap = 12
        params = self.interface
        pos_x, pos_y = pos_2d
        # offset_x and offset_y
        pos_x += params.offset_x
        pos_y += params.offset_y
        rgba = params.text_color
        text_size = params.text_size
        # adjust the text size if depth enabled
        if is3d and params.text_depth:
            view_matrix = self._rv3d.view_matrix
            if self._rv3d.is_perspective:  # perspective
                dist = (view_matrix @ (Vector(pos) * self._world_scale)).length
            else:  # orthographic
                dist = -(view_matrix @ (Vector(pos) * self._world_scale)).z
            # adjust distance range based on falloff factor
            dist_range = self._dist_range * params.text_falloff
            r_factor = 0  # reduction factor
            if dist_range > 0:  # to avoid div by 0
                offset = dist - self._min_dist
                # clamp to within range
                if offset < 0:
                    offset = 0
                elif offset > dist_range:
                    offset = dist_range
                r_factor = 1.0 - (offset / dist_range)
            text_size *= r_factor
        # set the text size
        blf.size(font_id, text_size)
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
        self._draw_line(
            v1, v2, v1_text, v2_text, mid_text, v1_arrow, v2_arrow, is3d=True
        )

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
        self._draw_line(
            v1, v2, v1_text, v2_text, mid_text, v1_arrow, v2_arrow, is3d=False
        )

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

    def draw_circle_3d(
        self,
        center: Vector,
        radius: float,
        normal: Vector,
        angle: float = 360.0,
        start_dv: Vector = None,
        c_arrow: bool = False,
        cc_arrow: bool = False,
    ):
        """
        Draw a circle around a 3D point in the plane perpendicular to the
        given normal

        Parameters:
        -----------
        center: Vector
            A 3D position vector of the center

        radius: float
            The radius of the circle

        normal: Vector
            The normal vector of the plane on which the cirle is to be drawn

        angle: float, optional
            An angle less than 360 for partial circle (arc) - in degrees
            Default is 360 degrees

        start_dv: Vector, optional
            The direction vector along which to start the circle (arc)
            If not provided, a random point in the plane perpendicular to the
            normal is chosen

        c_arrow: bool, optional
            Whether to display clockwise arrow. Default is False

        cc_arrow: bool, optional
            Whether to display counter clockwise arrow. Default is False

        """
        # convert to vectors
        if not isinstance(center, Vector):
            center = Vector(center)
        if not isinstance(normal, Vector):
            normal = Vector(normal)
        # get a point in the circle plane to start the circle
        if start_dv is None:
            start_dv = self._get_a_normal_plane_point(normal)
        start_dv.normalize()
        start = center + (start_dv * radius)
        n_steps = 36  # number of individual line segments of the circle
        step = radians(angle) / n_steps
        # matrices to translate to center, rotate and translate back
        mat_trans1 = Matrix.Translation(-center)
        mat_rot = Matrix.Rotation(step, 4, normal)
        mat_trans2 = Matrix.Translation(center)
        p1 = start
        # draw individual line segments
        for i in range(n_steps):
            p2 = mat_trans2 @ mat_rot @ mat_trans1 @ p1
            if i == 0 and cc_arrow:
                self._draw_line(p1, p2, v1_arrow=True, is3d=True)
            elif i == n_steps - 1 and c_arrow:
                self._draw_line(p1, p2, v2_arrow=True, is3d=True)
            else:
                self._draw_line(p1, p2, is3d=True)
            p1 = p2.copy()

    def _get_a_normal_plane_point(self, normal: Vector):
        """Internal: Get a point in the plane perpendicular to the given normal"""
        # given there are infinite points, pick any standard non zero
        # ones whose dot product to given normal is 0
        v = Vector((0, normal[2], -normal[1]))
        if v.length != 0:
            return v
        v = Vector((normal[2], 0, -normal[0]))
        if v.length != 0:
            return v
        v = Vector((normal[1], -normal[0], 0))
        if v.length != 0:
            return v
        raise ValueError("No non-zero vector in normal plane")

    def _draw_line(
        self,
        v1: Vector,
        v2: Vector,
        v1_text: str = None,
        v2_text: str = None,
        mid_text: str = None,
        v1_arrow: bool = False,
        v2_arrow: bool = False,
        is3d: bool = False,
    ) -> None:
        """Internal: Draw line 3D or 2D"""
        if v1 is None or v2 is None:
            return
        self._draw_arrow_line(v1, v2, v1_arrow, v2_arrow, is3d=is3d)
        if v1_text is not None:
            self._draw_text(v1, v1_text, is3d=is3d)
        if v2_text is not None:
            self._draw_text(v2, v2_text, is3d=is3d)
        if mid_text is not None:
            mid = (v1 + v2) / 2
            self._draw_text(mid, mid_text, is3d=is3d)

    def _draw_arrow_line(
        self,
        v1: Vector,
        v2: Vector,
        v1_arrow: bool = False,
        v2_arrow: bool = False,
        is3d: bool = False,
    ) -> None:
        """Internal: Draw a line between two 2D points with arrows"""
        if v1 is None or v2 is None:
            return
        v1_2d = v1
        v2_2d = v2
        if is3d:
            v1_2d = self._get_2d_point(v1)
            v2_2d = self._get_2d_point(v2)
            if v1_2d is None or v2_2d is None:
                return
        # actual line
        self._draw_line_2d(v1_2d, v2_2d)
        if v1_arrow:
            # v1 arrow lines
            va, vb = self._get_arrow_end_points(v1_2d, v2_2d)
            self._draw_line_2d(v1_2d, va)
            self._draw_line_2d(v1_2d, vb)
        if v2_arrow:
            # v2 arrow lines
            va, vb = self._get_arrow_end_points(v2_2d, v1_2d)
            self._draw_line_2d(v2_2d, va)
            self._draw_line_2d(v2_2d, vb)

    def _get_arrow_end_points(self, v1: Vector, v2: Vector) -> tuple:
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
        min_dist: float,
        max_dist: float,
    ) -> None:
        """Internal: Set the 3D viewport region and region data"""
        self._region = region
        self._rv3d = rv3d
        self._min_dist = min_dist
        self._max_dist = max_dist
        self._dist_range = max_dist - min_dist

    def _get_2d_point(self, pos_3d: Vector) -> Vector | None:
        """Internal: Get the 2D point in region corresponding to 3D position"""
        if self._rv3d is not None and self._region is not None:
            pos_3d_scaled = pos_3d * self._world_scale
            return view3d_utils.location_3d_to_region_2d(
                self._region, self._rv3d, pos_3d_scaled
            )
        else:
            return None
