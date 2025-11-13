import itertools
from abc import ABCMeta, abstractmethod
from math import cos, radians, sin, sqrt
from pathlib import Path
import blf
import bmesh
import bpy
import gpu
import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
from bpy_extras import object_utils, view3d_utils
from gpu_extras.batch import batch_for_shader
from mathutils import Matrix, Vector
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial import Voronoi
from ..blender.utils import new_bmesh
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
        if self.geometry:
            return
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
        if self.geometry:
            return
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
            if params.line_pointer_length > 0:
                # draw a pointer to the 3D position
                nv = -pos
                nv.normalize()
                pointer_begin = pos + (nv * params.line_pointer_length)
                text_pos = pointer_begin + (nv * 0.5)  # offset text a bit
                self._draw_line(
                    pointer_begin, pos, v2_arrow=True, is3d=True, overrides=overrides
                )
            pos_2d = self._get_2d_point(text_pos)
            if pos_2d is None:
                return
        if self.geometry:
            return
        # draw text at 2D position
        right_alignment_gap = 12
        pos_x, pos_y = pos_2d
        # text_offset_x and text_offset_y
        pos_x += params.text_offset_x
        pos_y += params.text_offset_y
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

        Parameters
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

        Parameters
        ----------
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
                self._draw_line(
                    p1,
                    p2,
                    v1_arrow=True,
                    arrow_plane_pt=start,
                    is3d=True,
                    overrides=overrides,
                )
            elif i == n_steps - 1 and c_arrow:
                self._draw_line(
                    p1,
                    p2,
                    v2_arrow=True,
                    arrow_plane_pt=start,
                    is3d=True,
                    overrides=overrides,
                )
            else:
                self._draw_line(p1, p2, is3d=True, overrides=overrides)
            p1 = p2.copy()

    def draw_sphere(
        self,
        location: Vector = (0, 0, 0),
        radius: float = 1.0,
        overrides: dict = None,
    ):
        """
        Draw a sphere

        Parameters
        ----------
        location: Vector
            A 3D position vector of the center

        radius: float
            Radius of the sphere

        overrides: dict, optional
            Optional dictionary to override common annotation params

        """
        if self.geometry is None:
            return
        with new_bmesh() as bm:
            bmesh.ops.create_icosphere(
                bm, radius=radius * self._world_scale, subdivisions=4
            )
            loc = Matrix.Translation(Vector(location) * self._world_scale)
            bm.transform(loc)
            self.draw_bmesh(bm, overrides=overrides)

    def draw_cone(
        self,
        location: Vector = (0, 0, 0),
        radius: float = 1.0,
        height: float = 1.0,
        axis: Vector = (0, 0, 1),
        cap_ends: bool = True,
        overrides: dict = None,
    ):
        """
        Draw a cone

        Parameters
        ----------
        location: Vector
            A 3D position vector of the base center

        radius: float
            Radius of the cone

        height: float
            Height of the cone

        axis: Vector
            Axis of the cone

        cap_ends: bool
            Whether to cap the base

        overrides: dict, optional
            Optional dictionary to override common annotation params

        """
        if self.geometry is None:
            return
        self._draw_cone(
            location,
            radius1=radius,
            radius2=0,
            height=height,
            axis=axis,
            cap_ends=cap_ends,
            overrides=overrides,
        )

    def draw_cylinder(
        self,
        location: Vector = (0, 0, 0),
        radius: float = 1.0,
        height: float = 1.0,
        axis: Vector = (0, 0, 1),
        cap_ends: bool = True,
        overrides: dict = None,
    ):
        """
        Draw a cylinder

        Parameters
        ----------
        location: Vector
            A 3D position vector of the base center

        radius: float
            Radius of the cylinder

        height: float
            Height of the cylinder

        axis: Vector
            Axis of the cylinder

        cap_ends: bool
            Whether to cap the ends of cylinder

        overrides: dict, optional
            Optional dictionary to override common annotation params

        """
        if self.geometry is None:
            return
        self._draw_cone(
            location,
            radius1=radius,
            radius2=radius,
            height=height,
            axis=axis,
            cap_ends=cap_ends,
            overrides=overrides,
        )

    def draw_triclinic_cell(
        self,
        a: float = 10.0,
        b: float = 10.0,
        c: float = 10.0,
        alpha: float = 90.0,
        beta: float = 90.0,
        gamma: float = 90.0,
        origin: Vector = (0, 0, 0),
        show_lattice: bool = False,
        overrides: dict = None,
    ):
        """
        Draw a triclinic box based on box vector lengths and angles

        Parameters
        ----------

        a: float
            Box vector a length

        b: float
            Box vector b length

        c: float
            Box vector c length

        alpha: float
            Angle between box vectors bc

        beta: float
            Angle between box vectors ac

        gamma: float
            Angle between box vectors ab

        origin: Vector
            Origin of the box

        show_lattice: bool
            Whether to show a 3x3x3 lattice

        overrides: dict, optional
            Optional dictionary to override common annotation params

        """
        if self.geometry is None:
            return
        dimensions = np.array([a, b, c, alpha, beta, gamma], dtype=np.float32)
        triclinic_vectors = mda.lib.mdamath.triclinic_vectors(dimensions)
        # convert to blender world scale
        box_vectors = triclinic_vectors * self._world_scale
        vo = Vector((0, 0, 0))
        vx = Vector(box_vectors[0])
        vxy = Vector(box_vectors[1])
        vz = Vector(box_vectors[2])
        vor = vx + vxy
        with new_bmesh() as bm:
            # create the four vertices in the xy plane
            v1 = bm.verts.new(vo)
            v2 = bm.verts.new(vx)
            v3 = bm.verts.new(vor)
            v4 = bm.verts.new(vxy)
            # create face in the xy plane
            bm.verts.ensure_lookup_table()
            face = bm.faces.new((v1, v2, v3, v4))
            # extrude face region for new verts
            ext_geom = bmesh.ops.extrude_face_region(bm, geom=[face])
            ext_verts = [
                v for v in ext_geom["geom"] if isinstance(v, bmesh.types.BMVert)
            ]
            # extrude verts along box vector c
            bmesh.ops.translate(bm, vec=vz, verts=ext_verts)
            # translate origin
            mat = Matrix.Translation(Vector(origin) * self._world_scale)
            bm.transform(mat)
            # update face normals
            bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
            # show a 3x3x3 lattice if enabled
            if show_lattice:
                box = bm.verts[:] + bm.edges[:] + bm.faces[:]
                pbc_img_coeffs = np.array(list(itertools.product([0, -1, 1], repeat=3)))
                lattice_points = np.matmul(pbc_img_coeffs, box_vectors)
                for lattice_point in lattice_points:
                    geom = bmesh.ops.duplicate(bm, geom=box)
                    dup_verts = [
                        v for v in geom["geom"] if isinstance(v, bmesh.types.BMVert)
                    ]
                    bmesh.ops.translate(bm, verts=dup_verts, vec=lattice_point)
            self.draw_bmesh(bm, overrides=overrides)

    def draw_wigner_seitz_cell(
        self,
        triclinic_vectors: npt.ArrayLike,
        center_to_origin: bool = False,
        show_lattice: bool = False,
        overrides: dict = None,
    ):
        """
        Draw a Wigner-Seitz cell from triclinic vectors

        Parameters
        ----------

        triclinic_vectors: npt.ArrayLike
            Vectors that represent the base triclinic cell

        center_to_origin: bool
            Move the center of the cell to origin (0, 0, 0)

        show_lattice: bool
            Whether to show a 3x3x3 lattice

        overrides: dict, optional
            Optional dictionary to override common annotation params

        """
        if self.geometry is None:
            return
        # convert to blender world scale
        box_vectors = triclinic_vectors * self._world_scale
        # From apply_compact_PBC of MDAnalysis
        pbc_img_coeffs = np.array(list(itertools.product([0, -1, 1], repeat=3)))
        pbc_img_vecs = np.matmul(pbc_img_coeffs, box_vectors)
        # add the center point
        box_center = np.sum(box_vectors, axis=0) / 2
        lattice_points = [box_center] + pbc_img_vecs
        # generate Voronoi
        voronoi = Voronoi(np.array(lattice_points))
        # The Wigner-Seitz cell for the center point (index 0)
        center_region = voronoi.point_region[0]
        center_region_vertices = voronoi.regions[center_region]
        # filter out invalid vertices
        valid_vertices = [
            voronoi.vertices[i] for i in center_region_vertices if i != -1
        ]
        # create bmesh
        with new_bmesh() as bm:
            # add vertices
            for vert in valid_vertices:
                bm.verts.new(vert)
            # use convex hull to create faces from vertices
            bm.verts.ensure_lookup_table()
            ch = bmesh.ops.convex_hull(bm, input=bm.verts)
            # delete interior verts
            bmesh.ops.delete(bm, geom=ch["geom_interior"], context="VERTS")
            # delete any duplicate vertices
            bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=1e-3)
            # convert triangles to n-gons for clean faces
            bmesh.ops.dissolve_limit(
                bm,
                angle_limit=1e-3,
                verts=bm.verts,
                edges=bm.edges,
            )
            # update face normals
            bm.normal_update()
            # transform center if needed
            if center_to_origin:
                mat = Matrix.Translation(-1 * box_center)
                bm.transform(mat)
            # show a 3x3x3 lattice if enabled
            if show_lattice:
                box = bm.verts[:] + bm.edges[:] + bm.faces[:]
                for lattice_point in lattice_points:
                    if np.array_equal(lattice_point, box_center):
                        continue
                    geom = bmesh.ops.duplicate(bm, geom=box)
                    dup_verts = [
                        v for v in geom["geom"] if isinstance(v, bmesh.types.BMVert)
                    ]
                    translation_vec = lattice_point - box_center
                    bmesh.ops.translate(bm, verts=dup_verts, vec=translation_vec)
            # draw bmesh
            self.draw_bmesh(bm, overrides=overrides)

    def draw_n_sided_pyramid(
        self,
        n: int = 6,
        radius: float = 10,
        height: float = 10,
        origin: Vector = (0, 0, 0),
        axis: Vector = (0, 0, 1),
        cap_ends: bool = True,
        overrides: dict = None,
    ):
        """
        Draw an n sided pyramid
        Eg: triangle, prism, square pyramid, pentagonal pyramid, etc

        Parameters
        ----------

        n: int
            Number of sides

        radius: float
            Radius of the pyramid

        height: float
            Height of the pyramid

        origin: Vector
            Center of the base of the pyramid

        axis: Vector
            Axis of the pyramid

        cap_ends: bool
            Whether to cap the ends of pyramid

        overrides: dict, optional
            Optional dictionary to override common annotation params

        """
        if self.geometry is None:
            return
        self._draw_cone(
            origin,
            radius1=radius,
            radius2=0,
            height=height,
            axis=axis,
            segments=n,
            cap_ends=cap_ends,
            overrides=overrides,
        )

    def draw_n_sided_cylinder(
        self,
        n: int = 6,
        radius: float = 10,
        height: float = 10,
        origin: Vector = (0, 0, 0),
        axis: Vector = (0, 0, 1),
        cap_ends: bool = True,
        overrides: dict = None,
    ):
        """
        Draw an n sided cylinder
        Eg: square, rectangle, triangular prism, cube, cuboid, hexagonal cell etc

        Parameters
        ----------

        n: int
            Number of sides

        radius: float
            Radius of the cylinder

        height: float
            Height of the cylinder

        origin: Vector
            Center of the base of the cylinder

        axis: Vector
            Axis of the cylinder

        cap_ends: bool
            Whether to cap the ends of cylinder

        overrides: dict, optional
            Optional dictionary to override common annotation params

        """
        if self.geometry is None:
            return
        self._draw_cone(
            origin,
            radius1=radius,
            radius2=radius,
            height=height,
            axis=axis,
            segments=n,
            cap_ends=cap_ends,
            overrides=overrides,
        )

    def draw_bmesh(
        self,
        bm: bmesh.types.BMesh,
        overrides: dict = None,
    ):
        """
        Draw a Blender bmesh

        Parameters
        ----------

        bm: bmesh.types.BMesh
            A bmesh object. A copy is made for internal use.
            Users will have to free the passed in object

        overrides: dict, optional
            Optional dictionary to override common annotation params

        """
        if self.geometry is None:
            return
        if not isinstance(bm, bmesh.types.BMesh):
            raise ValueError("Need a bmesh.types.BMesh object")
        geometry = self.geometry
        objects = geometry["objects"]
        # make a copy of the bmesh which will be freed after
        # updating the annotation object
        objects["meshes"].append(bm.copy())
        # add resolved params to be added as attributes
        params = _get_params(self.interface, overrides)
        objects["wireframe"].append(params.mesh_wireframe)
        objects["thickness"].append(params.mesh_thickness)
        objects["color"].append(params.mesh_color)
        objects["shade_smooth"].append(params.mesh_shade_smooth)
        self._add_material_to_geometry(objects, params.mesh_material)

    def _draw_cone(
        self,
        location: Vector,
        radius1: float,
        radius2: float,
        height: float,
        axis: Vector = (0, 0, 1),
        segments: int = 24,
        cap_ends: bool = True,
        overrides: dict = None,
    ):
        """Internal: Common method for cone, cylinder, pyramid bmesh"""
        if self.geometry is None:
            return
        radius1 = radius1 * self._world_scale
        radius2 = radius2 * self._world_scale
        height = height * self._world_scale
        with new_bmesh() as bm:
            bmesh.ops.create_cone(
                bm,
                cap_ends=cap_ends,
                segments=segments,
                radius1=radius1,
                radius2=radius2,
                depth=height,
            )
            base_offset = Matrix.Translation(Vector((0, 0, height / 2)))
            up = Vector((0, 0, 1))
            rot = up.rotation_difference(axis).to_matrix().to_4x4()
            loc = Matrix.Translation(Vector(location) * self._world_scale)
            mat = loc @ rot @ base_offset
            bm.transform(mat)
            self.draw_bmesh(bm, overrides=overrides)

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
        arrow_plane_pt: Vector = None,
        is3d: bool = False,
        overrides: dict = None,
    ) -> None:
        """Internal: Draw line 3D or 2D"""
        if v1 is None or v2 is None:
            return
        # convert to vectors from here on
        if not isinstance(v1, Vector):
            v1 = Vector(v1)
        if not isinstance(v2, Vector):
            v2 = Vector(v2)
        params = _get_params(self.interface, overrides)
        if is3d and self.geometry and params.line_mode != "overlay":
            # add line to geometry
            self._add_line_to_geometry(v1, v2, overrides=overrides)
            # add arrow ends to geometry
            if arrow_plane_pt is None:
                arrow_plane_pt = self._get_a_normal_plane_point(v2 - v1)
            if v1_arrow:
                self._add_arrow_end_to_geometry(
                    v1, v2, arrow_plane_pt, overrides=overrides
                )
            if v2_arrow:
                self._add_arrow_end_to_geometry(
                    v2, v1, arrow_plane_pt, overrides=overrides
                )

        self._draw_arrow_line(
            v1,
            v2,
            v1_arrow,
            v2_arrow,
            arrow_plane_pt=arrow_plane_pt,
            is3d=is3d,
            overrides=overrides,
        )
        if v1_text is not None:
            self._draw_text(v1, v1_text, is3d=is3d, overrides=overrides)
        if v2_text is not None:
            self._draw_text(v2, v2_text, is3d=is3d, overrides=overrides)
        if mid_text is not None:
            mid = (v1 + v2) / 2
            self._draw_text(mid, mid_text, is3d=is3d, overrides=overrides)

    def _add_line_to_geometry(
        self,
        v1: Vector,
        v2: Vector,
        overrides: dict = None,
    ):
        if self.geometry is None:
            return
        geometry = self.geometry
        lines = geometry["lines"]
        i = len(lines["vertices"])
        # add the ends of line as vertices
        lines["vertices"].append(v1 * self._world_scale)
        lines["vertices"].append(v2 * self._world_scale)
        # add an edge
        lines["edges"].append((i, i + 1))
        params = _get_params(self.interface, overrides)
        # add resolved params to be added as attributes
        lines["color"].append(params.mesh_color)
        lines["thickness"].append(params.mesh_thickness)
        self._add_material_to_geometry(lines, params.mesh_material)

    def _add_material_to_geometry(self, mesh_type, mesh_material):
        geometry = self.geometry
        if mesh_material:
            if isinstance(mesh_material, str):
                material = mesh_material
            else:
                material = mesh_material.name
            if material in geometry["materials"]:
                material_index = geometry["materials"][material]
            else:
                material_index = len(geometry["materials"])
                geometry["materials"][material] = material_index
            mesh_type["material_index"].append(material_index)
        else:
            mesh_type["material_index"].append(0)

    def _add_arrow_end_to_geometry(
        self,
        v1: Vector,
        v2: Vector,
        arrow_plane_pt: Vector,
        overrides: dict = None,
    ):
        if self.geometry is None:
            return
        # get arrow end points in 3d
        va, vb = self._get_arrow_end_points_3d(
            v1, v2, arrow_plane_pt, overrides=overrides
        )
        # add arrow lines to geometry
        self._add_line_to_geometry(v1, va, overrides=overrides)
        self._add_line_to_geometry(v1, vb, overrides=overrides)

    def _draw_arrow_line(
        self,
        v1: Vector,
        v2: Vector,
        v1_arrow: bool = False,
        v2_arrow: bool = False,
        arrow_plane_pt: Vector = None,
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
            if arrow_plane_pt is None:
                arrow_plane_pt = self._get_a_normal_plane_point(v2 - v1)
        # actual line
        self._draw_line_2d(v1_2d, v2_2d, overrides=overrides)
        # draw arrows
        if v1_arrow:
            self._draw_arrow(v1, v2, v1_2d, v2_2d, arrow_plane_pt, is3d, overrides)
        if v2_arrow:
            self._draw_arrow(v2, v1, v2_2d, v1_2d, arrow_plane_pt, is3d, overrides)

    def _draw_arrow(
        self,
        v1: Vector,
        v2: Vector,
        v1_2d: Vector,
        v2_2d: Vector,
        arrow_plane_pt: Vector = None,
        is3d: bool = False,
        overrides: dict = None,
    ) -> None:
        params = _get_params(self.interface, overrides)
        draw_3d_arrow_overlay = False
        if params.line_mode == "mesh_and_overlay":
            draw_3d_arrow_overlay = True
        # 3d or 2d arrow ends based on mode
        if is3d and draw_3d_arrow_overlay:
            va, vb = self._get_arrow_end_points_3d(
                v1, v2, arrow_plane_pt, overrides=overrides
            )
            # get the 2d points of the arrow ends for overlay
            va = self._get_2d_point(va)
            vb = self._get_2d_point(vb)
            if va is None or vb is None:
                return
        else:
            va, vb = self._get_arrow_end_points_2d(v1_2d, v2_2d, overrides=overrides)
        # draw the arrow ends
        self._draw_line_2d(v1_2d, va, overrides=overrides)
        self._draw_line_2d(v1_2d, vb, overrides=overrides)

    def _get_arrow_end_points_3d(
        self,
        v1: Vector,
        v2: Vector,
        arrow_plane_pt: Vector = None,
        overrides: dict = None,
    ) -> tuple:
        """Internal: Get arrow end point positions 3D"""
        if arrow_plane_pt is None:
            arrow_plane_pt = self._get_a_normal_plane_point(v2 - v1)
        # position vectors in arrow plane
        pv1 = v1 - arrow_plane_pt
        pv2 = v2 - arrow_plane_pt
        # calculate rotation axis
        ra = pv1.cross(pv2)
        ra.normalize()
        # use Rodrigues' Rotation Formula
        # https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        v = v2 - v1
        # +45 degrees direction vector
        dva = (
            v * cos(self._rad45)
            + (ra.cross(v) * sin(self._rad45))
            + (ra * ra.dot(v) * (1 - cos(self._rad45)))
        )
        dva.normalize()
        # -45 degrees direction vector
        dvb = (
            v * cos(-1 * self._rad45)
            + (ra.cross(v) * sin(-1 * self._rad45))
            + (ra * ra.dot(v) * (1 - cos(-1 * self._rad45)))
        )
        dvb.normalize()
        params = _get_params(self.interface, overrides)
        d = self.distance(v1, v2) * params.line_arrow_size
        return (v1 + (dva * d), v1 + (dvb * d))

    def _get_arrow_end_points_2d(
        self, v1: Vector, v2: Vector, overrides: dict = None
    ) -> tuple:
        """Internal: Get arrow end point positions 2D"""
        params = _get_params(self.interface, overrides)
        d = self.distance(v1, v2) * params.line_arrow_size
        v = self._interpolate_3d((v1[0], v1[1], 0.0), (v2[0], v2[1], 0.0), d)
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
        if self.geometry:
            return
        if v1 is None or v2 is None:
            return
        params = _get_params(self.interface, overrides)
        if params.line_mode == "mesh":
            return
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
