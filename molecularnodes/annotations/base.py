from abc import ABCMeta, abstractmethod
from math import cos, radians, sin, sqrt
from pathlib import Path
import blf
import bpy
import gpu
from bpy_extras import object_utils, view3d_utils
from gpu_extras.batch import batch_for_shader
from mathutils import Matrix, Vector
from PIL import Image, ImageDraw, ImageFont
from .interface import AnnotationInterface
from .utils import get_view_matrix, is_perspective_projection

FONT_INTER = Path(__file__).parent / "fonts/Inter.woff2"


class _get_params:
    def __init__(self, interface: AnnotationInterface, overrides: dict = None):
        self.interface = interface
        self.overrides = overrides

    def __getattr__(self, name):
        if self.overrides and name in self.overrides:
            return self.overrides[name]
        return getattr(self.interface, name)


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
        self._scene = None
        self._scale = 1.0
        self._render_mode = False
        self._dpi_scale = 2 * (72 / bpy.context.preferences.system.dpi)
        # distance params
        self._min_dist = 0
        self._max_dist = 0
        self._dist_range = 0
        # viewport params
        self._region = None
        self._rv3d = None
        # render params
        self._render_scale = 1.0
        self._image = None
        self._image_scale = 1

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
        if self._render_mode:
            return self._image.width / self._image_scale
        else:
            return self._region.width

    @property
    def viewport_height(self) -> int:
        """Get the viewport region height in pixels"""
        if self._render_mode:
            return self._image.height / self._image_scale
        else:
            return self._region.height

    def draw_text_3d(self, pos_3d: Vector, text: str, overrides: dict = None) -> None:
        """
        Draw text at a given 3D position

        Parameters
        ----------
        pos_3d: Vector
            Co-ordinates in 3D world space (x, y, z)

        text: str
            Text to display. '|' as multi-line separator

        overrides: dict, optional
            Optional dictionary to override common annotation params

        """
        self._draw_text(pos_3d, text, is3d=True, overrides=overrides)

    def draw_text_2d_norm(
        self, pos_2d: Vector, text: str, overrides: dict = None
    ) -> None:
        """
        Draw text at a given 2D position (normalized co-ordinates) of Viewport.

        Parameters
        ----------
        pos_2d: Vector
            Normalized co-ordinates (0 - 1). (0, 0) is at bottom left

        text: str
            Text to display. '|' as multi-line separator

        overrides: dict, optional
            Optional dictionary to override common annotation params

        """
        if pos_2d is None:
            return
        for comp in pos_2d:
            if not (0.0 <= comp <= 1.0):
                return
        pos_x, pos_y = pos_2d
        new_x = pos_x * self.viewport_width
        new_y = pos_y * self.viewport_height
        if not self._render_mode and self._rv3d.view_perspective == "CAMERA":
            # camera view mode in 3D viewport
            zoom_factor, camera_view_width, camera_view_height = (
                self._get_camera_view_info()
            )
            # offsets are based off the center of the viewport
            camera_offset_x = self._rv3d.view_camera_offset[0]
            camera_offset_y = self._rv3d.view_camera_offset[1]
            # calculate the origin (bottom left) of the camera view
            camera_view_x0 = (
                (self.viewport_width / 2)
                - (camera_view_width / 2)
                - (camera_offset_x * self.viewport_width * 2 * zoom_factor)
            )
            camera_view_y0 = (
                (self.viewport_height / 2)
                - (camera_view_height / 2)
                - (camera_offset_y * self.viewport_height * 2 * zoom_factor)
            )
            # calculate the actual position with respect to the camera view origin
            new_x = camera_view_x0 + (pos_x * camera_view_width)
            new_y = camera_view_y0 + (pos_y * camera_view_height)
        self._draw_text((new_x, new_y), text, is3d=False, overrides=overrides)

    def draw_text_2d(self, pos_2d: Vector, text: str, overrides: dict = None) -> None:
        """
        Draw text at a given 2D position (in pixels) of Viewport.

        Parameters
        ----------
        pos_2d: Vector
            Co-ordinates in pixels. (0, 0) is at bottom left

        text: str
            Text to display. '|' as multi-line separator

        overrides: dict, optional
            Optional dictionary to override common annotation params

        """
        self._draw_text(pos_2d, text, is3d=False, overrides=overrides)

    def _draw_text(
        self, pos: Vector, text: str, is3d: bool = False, overrides: dict = None
    ) -> None:
        """Internal: Draw text 3D or 2D"""
        if pos is None:
            return
        if not isinstance(pos, Vector):
            pos = Vector(pos)
        pos_2d = pos
        params = _get_params(self.interface, overrides)
        if is3d:
            text_pos = pos
            # check if pointer is required
            if params.pointer_length > 0:
                # draw a pointer to the 3D position
                nv = -pos
                nv.normalize()
                pointer_begin = pos + (nv * params.pointer_length)
                text_pos = pointer_begin + (nv * 0.5)  # offset text a bit
                self._draw_line(
                    pointer_begin, pos, v2_arrow=True, is3d=True, overrides=overrides
                )
            pos_2d = self._get_2d_point(text_pos)
            if pos_2d is None:
                return
        # draw text at 2D position
        right_alignment_gap = 12
        pos_x, pos_y = pos_2d
        # offset_x and offset_y
        pos_x += params.offset_x
        pos_y += params.offset_y
        rgba = params.text_color
        text_size = params.text_size
        # adjust the text size if depth enabled
        if is3d and params.text_depth:
            view_matrix = get_view_matrix(self)
            if is_perspective_projection(self):
                dist = (view_matrix @ (pos * self._world_scale)).length
            else:  # orthographic
                dist = -(view_matrix @ (pos * self._world_scale)).z
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
        # scale the text - applies for viewport, viewport camera view and renders
        text_size *= self._scale
        if round(text_size) == 0:
            return
        # set the text size
        if self._render_mode:
            load_default_font = True
            if params.text_font != "":
                try:
                    # try to load user defined font if specified
                    font = ImageFont.truetype(
                        params.text_font, text_size * self._image_scale
                    )
                    load_default_font = False
                except Exception:
                    pass  # loads default font
            if load_default_font:
                font = ImageFont.truetype(FONT_INTER, text_size * self._image_scale)
            image_draw = ImageDraw.Draw(self._image)
            bbox = image_draw.textbbox((0, 0), "Tp", font=font)
            max_th = bbox[3] - bbox[1]
            pos_scale = self._image_scale
        else:
            font_id = 0  # default font
            if params.text_font != "":
                # load user defined font if specified
                font_id = blf.load(params.text_font)
                if font_id == -1:  # use default if load fails
                    font_id = 0
            blf.size(font_id, text_size)
            # height of one line - for use in multiline text
            _, max_th = blf.dimensions(font_id, "Tp")  # uses high/low letters
            pos_scale = 1

        # split lines
        lines = text.split("|")
        n_lines = len(lines) - 1
        # draw all text lines
        for line in lines:
            if self._render_mode:
                bbox = image_draw.textbbox((0, 0), " " + line, font=font)
                line_width = bbox[2] - bbox[0]
            else:
                line_width, _ = blf.dimensions(font_id, line)
            # x position based on alignment
            if params.text_align == "center":
                new_x = (pos_x * pos_scale) - line_width / 2
            elif params.text_align == "right":
                new_x = (pos_x * pos_scale) - line_width - right_alignment_gap
            else:
                new_x = pos_x * pos_scale
                if not self._render_mode:
                    # rotation support only for left alignment and viewport
                    blf.enable(font_id, blf.ROTATION)
                    blf.rotation(font_id, radians(params.text_rotation))
            # calculate new y position
            new_y = (pos_y * pos_scale) + ((max_th * params.text_vspacing) * n_lines)
            # draw the text line
            if self._render_mode:
                image_draw.text(
                    (new_x, self._image.height - new_y - max_th),
                    text=" " + line,
                    font=font,
                    fill=(
                        round(rgba[0] * 255),
                        round(rgba[1] * 255),
                        round(rgba[2] * 255),
                        round(rgba[3] * 255),
                    ),
                )
            else:
                blf.position(font_id, new_x, new_y, 0)
                blf.color(font_id, rgba[0], rgba[1], rgba[2], rgba[3])
                blf.draw(font_id, " " + line)
            n_lines -= 1  # for next iteration
            if params.text_align == "left" and not self._render_mode:
                blf.disable(font_id, blf.ROTATION)
        if not self._render_mode:
            # unload any user defined font
            if font_id != 0:
                blf.unload(params.text_font)

    def draw_line_3d(
        self,
        v1: Vector,
        v2: Vector,
        v1_text: str = None,
        v2_text: str = None,
        mid_text: str = None,
        v1_arrow: bool = False,
        v2_arrow: bool = False,
        overrides: dict = None,
    ) -> None:
        """
        Draw a line between two points in 3D space

        Parameters
        ----------
        v1: Vector
            3D co-ordinates of the first point

        v2: Vector
            3D co-ordinates of the second point

        v1_text: str, optional
            Optional text to display at v1

        v2_text: str, optional
            Optional text to display at v2

        mid_text: str, optional
            Optional text to display at the middle of the line

        v1_arrow: bool, optional
            Whether to display an arrow at v1

        v2_arrow: bool, optional
            Whether to display an arrow at v2

        overrides: dict, optional
            Optional dictionary to override common annotation params

        """
        self._draw_line(
            v1,
            v2,
            v1_text,
            v2_text,
            mid_text,
            v1_arrow,
            v2_arrow,
            is3d=True,
            overrides=overrides,
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
        overrides: dict = None,
    ) -> None:
        """
        Draw a line between two points in 2D viewport space

        Parameters
        ----------
        v1: Vector
            2D co-ordinates of the first point

        v2: Vector
            2D co-ordinates of the second point

        v1_text: str, optional
            Optional text to display at v1

        v2_text: str, optional
            Optional text to display at v2

        mid_text: str, optional
            Optional text to display at the middle of the line

        v1_arrow: bool, optional
            Whether to display an arrow at v1

        v2_arrow: bool, optional
            Whether to display an arrow at v2

        overrides: dict, optional
            Optional dictionary to override common annotation params

        """
        self._draw_line(
            v1,
            v2,
            v1_text,
            v2_text,
            mid_text,
            v1_arrow,
            v2_arrow,
            is3d=False,
            overrides=overrides,
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
        overrides: dict = None,
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

        overrides: dict, optional
            Optional dictionary to override common annotation params

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
                self._draw_line(p1, p2, v1_arrow=True, is3d=True, overrides=overrides)
            elif i == n_steps - 1 and c_arrow:
                self._draw_line(p1, p2, v2_arrow=True, is3d=True, overrides=overrides)
            else:
                self._draw_line(p1, p2, is3d=True, overrides=overrides)
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
        overrides: dict = None,
    ) -> None:
        """Internal: Draw line 3D or 2D"""
        if v1 is None or v2 is None:
            return
        self._draw_arrow_line(
            v1, v2, v1_arrow, v2_arrow, is3d=is3d, overrides=overrides
        )
        if v1_text is not None:
            self._draw_text(v1, v1_text, is3d=is3d, overrides=overrides)
        if v2_text is not None:
            self._draw_text(v2, v2_text, is3d=is3d, overrides=overrides)
        if mid_text is not None:
            mid = (v1 + v2) / 2
            self._draw_text(mid, mid_text, is3d=is3d, overrides=overrides)

    def _draw_arrow_line(
        self,
        v1: Vector,
        v2: Vector,
        v1_arrow: bool = False,
        v2_arrow: bool = False,
        is3d: bool = False,
        overrides: dict = None,
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
        self._draw_line_2d(v1_2d, v2_2d, overrides=overrides)
        if v1_arrow:
            # v1 arrow lines
            va, vb = self._get_arrow_end_points(v1_2d, v2_2d, overrides=overrides)
            self._draw_line_2d(v1_2d, va, overrides=overrides)
            self._draw_line_2d(v1_2d, vb, overrides=overrides)
        if v2_arrow:
            # v2 arrow lines
            va, vb = self._get_arrow_end_points(v2_2d, v1_2d, overrides=overrides)
            self._draw_line_2d(v2_2d, va, overrides=overrides)
            self._draw_line_2d(v2_2d, vb, overrides=overrides)

    def _get_arrow_end_points(
        self, v1: Vector, v2: Vector, overrides: dict = None
    ) -> tuple:
        """Internal: Get arrow end point positions"""
        params = _get_params(self.interface, overrides)
        arrow_size = params.arrow_size * self._scale
        v = self._interpolate_3d((v1[0], v1[1], 0.0), (v2[0], v2[1], 0.0), arrow_size)
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

    def _draw_line_2d(self, v1: Vector, v2: Vector, overrides: dict = None) -> None:
        """Internal: Draw a line between two 2D points"""
        if v1 is None or v2 is None:
            return
        params = _get_params(self.interface, overrides)
        rgba = params.line_color
        line_width = params.line_width * self._scale
        if self._render_mode:
            v1 = tuple([v * self._image_scale for v in v1])
            v2 = tuple([v * self._image_scale for v in v2])
            image_height = self._image.height
            image_draw = ImageDraw.Draw(self._image)
            image_draw.line(
                [(v1[0], image_height - v1[1]), (v2[0], image_height - v2[1])],
                fill=(
                    round(rgba[0] * 255),
                    round(rgba[1] * 255),
                    round(rgba[2] * 255),
                    round(rgba[3] * 255),
                ),
                width=round(line_width),
            )
        else:
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

    def _set_distance_params(self, min_dist: float, max_dist: float) -> None:
        "Internal: Set the min, max distances and range"
        self._min_dist = min_dist
        self._max_dist = max_dist
        self._dist_range = max_dist - min_dist

    def _get_camera_view_info(self):
        # camera view mode in 3D viewport
        zoom = self._rv3d.view_camera_zoom
        # From: BKE_screen_view3d_zoom_to_fac in blender/blenkernel/intern/screen.cc
        zoom_factor = (((zoom / 50.0) + sqrt(2)) / 2) ** 2
        aspect_ratio = self._scene.render.resolution_x / self._scene.render.resolution_y
        # caculate actual camera width and height in view
        camera_view_width = self.viewport_width * zoom_factor
        if self.viewport_height > self.viewport_width:
            camera_view_width = self.viewport_height * zoom_factor
        camera_view_height = camera_view_width / aspect_ratio
        return zoom_factor, camera_view_width, camera_view_height

    def _set_viewport_params(
        self,
        scene: bpy.types.Scene,
        region: bpy.types.Region | None,
        rv3d: bpy.types.RegionView3D | None,
    ) -> None:
        """Internal: Set the 3D viewport region and region data"""
        self._scene = scene
        self._region = region
        self._rv3d = rv3d
        self._scale = bpy.context.preferences.system.ui_scale
        self._render_mode = False
        if self._rv3d.view_perspective == "CAMERA":
            # camera view mode in 3D viewport
            _, camera_view_width, _ = self._get_camera_view_info()
            self._scale *= (
                camera_view_width / self._scene.render.resolution_x
            ) * self._dpi_scale

    def _set_render_params(
        self,
        scene: bpy.types.Scene,
        render_scale: float,
        image: Image.Image,
        image_scale: float,
    ) -> None:
        """Internal: Set the render scene, image and scale"""
        self._scene = scene
        self._render_scale = render_scale
        self._image = image
        self._image_scale = image_scale
        self._scale = self._render_scale * self._dpi_scale
        self._render_mode = True

    def _get_2d_point(self, pos_3d: Vector) -> Vector | None:
        """Internal: Get the 2D point in region corresponding to 3D position"""
        if not isinstance(pos_3d, Vector):
            pos_3d = Vector(pos_3d)
        if self._render_mode:
            scene = self._scene
            pos_2d = object_utils.world_to_camera_view(
                scene, scene.camera, pos_3d * self._world_scale
            )
            x = round(pos_2d.x * (self._image.width / self._image_scale))
            y = round(pos_2d.y * (self._image.height / self._image_scale))
            return Vector((x, y))
        elif self._rv3d is not None and self._region is not None:
            return view3d_utils.location_3d_to_region_2d(
                self._region, self._rv3d, pos_3d * self._world_scale
            )
        else:
            return None
