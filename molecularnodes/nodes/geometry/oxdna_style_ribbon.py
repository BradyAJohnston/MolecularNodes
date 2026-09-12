# Node-group asset "oxDNA Style Ribbon" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    ColorSocket,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    MaterialSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import (
    InputBoolean,
    InputColor,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputMaterial,
    InputMenu,
    InputVector,
)
from ._shared.set_instancer import SetInstancer
from ._shared.smooth_by_angle import SmoothByAngle
from .angstrom_to_world import AngstromToWorld
from .chain_id import ChainID
from .color import Color
from .color_res_name import ColorResName
from .fallback_geometry import FallbackGeometry
from .integer_distance import IntegerDistance
from .oxdna_vectors import OxDNAVectors
from .set_color import SetColor


class OxDNAStyleRibbon(AssetGeometryGroup):
    """
    oxDNA Style Ribbon

    Parameters
    ----------
    atoms : InputGeometry
        Atoms
    selection : InputBoolean
        Selection of atoms to apply this node to
    quality : InputInteger
        Quality
    backbone_shape : InputMenu | Literal["Arrows", "Curve"]
        Backbone Shape
    backbone_radius : InputFloat
        Backbone Radius
    ball_radius : InputFloat
        Ball Radius
    arrow_taper : InputFloat
        Arrow Taper
    base_shape : InputMenu | Literal["Sphere", "Cylinder", "None"]
        Base Shape
    base_geometry : InputGeometry
        Base Geometry
    base_scale : InputVector
        Base Scale
    stem_scale : InputVector
        Stem Scale
    base_colors : InputMenu | Literal["Uniform", "Base", "Color"]
        Base Colors
    bases : InputColor
        Bases
    a : InputColor
        A
    c : InputColor
        C
    g : InputColor
        G
    t_u : InputColor
        T / U
    strand_color : InputMenu | Literal["Uniform", "Strand", "Color"]
        Strand Color
    strands : InputColor
        Becomes the output value if it is chosen by the menu input
    strand_1 : InputColor
        Strand 1
    strand_2 : InputColor
        Strand 2
    strand_3 : InputColor
        Strand 3
    strand_4 : InputColor
        Strand 4
    shade_smooth : InputBoolean
        Shade Smooth
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.atoms : GeometrySocket
        Atoms
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.quality : IntegerSocket
        Quality
    i.backbone_shape : MenuSocket
        Backbone Shape
    i.backbone_radius : FloatSocket
        Backbone Radius
    i.ball_radius : FloatSocket
        Ball Radius
    i.arrow_taper : FloatSocket
        Arrow Taper
    i.base_shape : MenuSocket
        Base Shape
    i.base_geometry : GeometrySocket
        Base Geometry
    i.base_scale : VectorSocket
        Base Scale
    i.stem_scale : VectorSocket
        Stem Scale
    i.base_colors : MenuSocket
        Base Colors
    i.bases : ColorSocket
        Bases
    i.a : ColorSocket
        A
    i.c : ColorSocket
        C
    i.g : ColorSocket
        G
    i.t_u : ColorSocket
        T / U
    i.strand_color : MenuSocket
        Strand Color
    i.strands : ColorSocket
        Becomes the output value if it is chosen by the menu input
    i.strand_1 : ColorSocket
        Strand 1
    i.strand_2 : ColorSocket
        Strand 2
    i.strand_3 : ColorSocket
        Strand 3
    i.strand_4 : ColorSocket
        Strand 4
    i.shade_smooth : BooleanSocket
        Shade Smooth
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "oxDNA Style Ribbon"
    _asset_name = "oxDNA Style Ribbon"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "default_group_node_width": 160,
        "node_tool_idname": "geometry.mn_oxdna_style_ribbon",
    }

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atoms"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        quality: IntegerSocket
        """Quality"""
        backbone_shape: MenuSocket
        """Backbone Shape"""
        backbone_radius: FloatSocket
        """Backbone Radius"""
        ball_radius: FloatSocket
        """Ball Radius"""
        arrow_taper: FloatSocket
        """Arrow Taper"""
        base_shape: MenuSocket
        """Base Shape"""
        base_geometry: GeometrySocket
        """Base Geometry"""
        base_scale: VectorSocket
        """Base Scale"""
        stem_scale: VectorSocket
        """Stem Scale"""
        base_colors: MenuSocket
        """Base Colors"""
        bases: ColorSocket
        """Bases"""
        a: ColorSocket
        """A"""
        c: ColorSocket
        """C"""
        g: ColorSocket
        """G"""
        t_u: ColorSocket
        """T / U"""
        strand_color: MenuSocket
        """Strand Color"""
        strands: ColorSocket
        """Becomes the output value if it is chosen by the menu input"""
        strand_1: ColorSocket
        """Strand 1"""
        strand_2: ColorSocket
        """Strand 2"""
        strand_3: ColorSocket
        """Strand 3"""
        strand_4: ColorSocket
        """Strand 4"""
        shade_smooth: BooleanSocket
        """Shade Smooth"""
        material: MaterialSocket
        """Material to apply to the resulting geometry"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        quality: InputInteger = 2,
        backbone_shape: InputMenu | Literal["Arrows", "Curve"] = "Arrows",
        backbone_radius: InputFloat = 0.8,
        ball_radius: InputFloat = 2.0,
        arrow_taper: InputFloat = 0.3,
        base_shape: InputMenu | Literal["Sphere", "Cylinder", "None"] = "Sphere",
        base_geometry: InputGeometry = None,
        base_scale: InputVector = None,
        stem_scale: InputVector = None,
        base_colors: InputMenu | Literal["Uniform", "Base", "Color"] = "Uniform",
        bases: InputColor = None,
        a: InputColor = None,
        c: InputColor = None,
        g: InputColor = None,
        t_u: InputColor = None,
        strand_color: InputMenu | Literal["Uniform", "Strand", "Color"] = "Color",
        strands: InputColor = None,
        strand_1: InputColor = None,
        strand_2: InputColor = None,
        strand_3: InputColor = None,
        strand_4: InputColor = None,
        shade_smooth: InputBoolean = True,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Quality": quality,
                "Backbone Shape": backbone_shape,
                "Backbone Radius": backbone_radius,
                "Ball Radius": ball_radius,
                "Arrow Taper": arrow_taper,
                "Base Shape": base_shape,
                "Base Geometry": base_geometry,
                "Base Scale": base_scale,
                "Stem Scale": stem_scale,
                "Base Colors": base_colors,
                "Bases": bases,
                "A": a,
                "C": c,
                "G": g,
                "T / U": t_u,
                "Strand Color": strand_color,
                "Strands": strands,
                "Strand 1": strand_1,
                "Strand 2": strand_2,
                "Strand 3": strand_3,
                "Strand 4": strand_4,
                "Shade Smooth": shade_smooth,
                "Material": material,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry("Atoms")
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        quality = tree.inputs.integer("Quality", 2, min_value=1)
        with tree.inputs.panel("Backbone"):
            backbone_shape = tree.inputs.menu(
                "Backbone Shape", expanded=True, optional_label=True
            )
            backbone_radius = tree.inputs.float(
                "Backbone Radius",
                0.8,
                min_value=0.0,
                max_value=340_282_000_000_000_000_000_000_000_000_000_000_000.0,
            )
            ball_radius = tree.inputs.float("Ball Radius", 2.0, min_value=0.0)
            arrow_taper = tree.inputs.float(
                "Arrow Taper", 0.3, min_value=0.0, max_value=1.0, subtype="FACTOR"
            )
        with tree.inputs.panel("Bases"):
            base_shape = tree.inputs.menu(
                "Base Shape", expanded=True, optional_label=True
            )
            base_geometry = tree.inputs.geometry("Base Geometry")
            base_scale = tree.inputs.vector(
                "Base Scale", (0.1, 0.3, 0.2), min_value=0.0, subtype="XYZ"
            )
            stem_scale = tree.inputs.vector("Stem Scale", (0.7, 0.7, 0.7))
        with tree.inputs.panel("Colors"):
            base_colors = tree.inputs.menu(
                "Base Colors", expanded=True, optional_label=True
            )
            bases = tree.inputs.color("Bases", (0.0, 1.0, 1.0, 1.0))
            a = tree.inputs.color("A", (0.033104, 0.03310406, 1.0, 1.0))
            c_ = tree.inputs.color("C", (0.033104, 1.0, 0.033104, 1.0))
            g_ = tree.inputs.color("G", (1.0, 1.0, 0.033104, 1.0))
            t_u = tree.inputs.color("T / U", (1.0, 0.033104, 0.033104, 1.0))
            strand_color = tree.inputs.menu(
                "Strand Color", expanded=True, optional_label=True
            )
            strands = tree.inputs.color(
                "Strands",
                (0.8, 0.8, 0.8, 1.0),
                description="Becomes the output value if it is chosen by the menu input",
            )
            strand_1 = tree.inputs.color("Strand 1", (1.0, 0.0, 0.0, 1.0))
            strand_2 = tree.inputs.color("Strand 2", (0.0, 0.0, 1.0, 1.0))
            strand_3 = tree.inputs.color("Strand 3", (0.0, 1.0, 0.0, 1.0))
            strand_4 = tree.inputs.color("Strand 4", (1.0, 1.0, 0.0, 1.0))
        with tree.inputs.panel("Material"):
            shade_smooth = tree.inputs.boolean("Shade Smooth", True)
            material = tree.inputs.material(
                "Material", description="Material to apply to the resulting geometry"
            )
        geometry = tree.outputs.geometry("Geometry")

        separate_geometry = g.SeparateGeometry.point(atoms, selection)
        with g.Frame("Color strands if Auto-color is False"):
            integer_math = ChainID().o.chain_id.modulo(4)
            menu_switch = g.MenuSwitch.color(
                strand_color,
                {
                    "Uniform": (strands, "Single uniform color"),
                    "Strand": (
                        g.IndexSwitch.color(
                            integer_math, (strand_1, strand_2, strand_3, strand_4)
                        ),
                        "Set custom colors for the bases",
                    ),
                    "Color": (Color(), "Use the existing `Color` attribute"),
                },
            )
            group = SetColor(
                atoms=separate_geometry.o.selection,
                selection=integer_math,
                color=menu_switch.o.output,
            )
        group_1 = OxDNAVectors()
        capture = g.CaptureAttribute.point(geometry=group)
        rotation = capture.items.rotation("Rotation", group_1.o.rotation)
        with g.Frame("Colored bases"):
            menu_switch_1 = g.MenuSwitch.color(
                base_colors,
                {
                    "Uniform": bases,
                    "Base": ColorResName(
                        a=a, c=c_, g=g_, t=t_u, ra=a, rc=c_, rg=g_, ru=t_u
                    ),
                    "Color": Color(),
                },
            )
            instance_on_points = (
                SetColor(
                    atoms=SetInstancer(geometry=capture.o.geometry),
                    color=menu_switch_1.o.output,
                )
                >> g.SetPosition(offset=group_1.o.base_offset)
                >> g.InstanceOnPoints(
                    instance=FallbackGeometry(
                        geometry=base_geometry,
                        fallback=g.IcoSphere(subdivisions=quality),
                    ),
                    rotation=rotation.output,
                    scale=base_scale,
                )
            )
        group_2 = OxDNAVectors()
        with g.Frame():
            set_position = g.SetPosition(
                geometry=capture.o.geometry, offset=group_2.o.backbone_offset
            )
            with g.Frame("Backbone Stick"):
                with g.Frame(
                    "Each segment is it's own mesh, flipping the circular endpoints"
                ):
                    edge_vertices = g.EdgeVertices()
                    capture_1 = g.CaptureAttribute.edge(
                        geometry=set_position,
                        selection=IntegerDistance(
                            a=edge_vertices.o.vertex_index_1,
                            b=edge_vertices.o.vertex_index_2,
                        ).o.cutoff,
                    )
                    reverse_curve = (
                        capture_1.o.geometry
                        >> g.SplitEdges()
                        >> g.MeshToCurve()
                        >> g.ReverseCurve(selection=capture_1.o.selection)
                    )
                group_3 = AngstromToWorld(angstrom=backbone_radius)
                curve_circle = g.CurveCircle(resolution=quality * 4, radius=ball_radius)
                switch = g.EndpointSelection(start_size=0).o.selection.switch.float(
                    group_3, group_3.o.world * arrow_taper
                )
                set_spline_resolution = (
                    g.MeshToCurve(mesh=set_position)
                    >> g.SetCurveNormal(normal=group_2.o.base_normal, mode="Free")
                    >> g.SetSplineType.bezier()
                    >> g.SetSplineResolution(resolution=quality * 2)
                )
                curve_to_mesh = g.SetHandleType(
                    curve=set_spline_resolution
                ) >> g.CurveToMesh(profile_curve=curve_circle, scale=group_3)
                curve_to_mesh_1 = reverse_curve >> g.CurveToMesh(
                    profile_curve=curve_circle, scale=switch
                )
            with g.Frame("Backbone ball"):
                instance_on_points_1 = SetInstancer(
                    geometry=set_position
                ) >> g.InstanceOnPoints(
                    instance=g.IcoSphere(
                        radius=AngstromToWorld(angstrom=ball_radius),
                        subdivisions=quality,
                    )
                )
        with g.Frame("Base stem"):
            cylinder = g.Cylinder(
                vertices=quality * 5,
                side_segments=quality,
                radius=AngstromToWorld(angstrom=1.0),
                depth=1.0,
            )
            instance_on_points_2 = SetInstancer(
                geometry=set_position
            ) >> g.InstanceOnPoints(
                instance=g.TransformGeometry(
                    geometry=cylinder, translation=(0.0, 0.0, 0.5)
                ),
                rotation=rotation.output.rotate(
                    (math.pi / 9, 0.0, 0.0), rotation_space="LOCAL"
                ),
                scale=stem_scale,
            )
        menu_switch_2 = g.MenuSwitch.geometry(
            base_shape,
            {
                "Sphere": g.JoinGeometry(
                    geometry=(instance_on_points, instance_on_points_2)
                ),
                "Cylinder": instance_on_points_2,
                "None": None,
            },
        )
        menu_switch_3 = g.MenuSwitch.geometry(
            backbone_shape,
            {
                "Arrows": g.JoinGeometry(
                    geometry=(instance_on_points_1, curve_to_mesh_1)
                ),
                "Curve": curve_to_mesh,
            },
        )
        set_shade_smooth = g.SetShadeSmooth.face(
            g.JoinGeometry(geometry=(menu_switch_3, menu_switch_2)),
            shade_smooth=shade_smooth,
        )
        set_shade_smooth.node.warning_propagation = "ERRORS"
        group_4 = SmoothByAngle(mesh=set_shade_smooth, angle=math.pi / 3)
        group_4.node.warning_propagation = "ERRORS"
        group_4 >> g.SetMaterial(material=material) >> geometry

        backbone_shape.default_value = "Arrows"
        base_shape.default_value = "Sphere"
        base_colors.default_value = "Uniform"
        strand_color.default_value = "Color"


ASSET = OxDNAStyleRibbon

ASSET_METADATA = {
    "catalog_id": "0094c3e0-7885-427b-81b4-187a84dcff18",
}
