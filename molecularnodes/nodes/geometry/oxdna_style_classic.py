# Node-group asset "oxDNA Style Classic" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING, Literal
import bpy
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    ColorSocket,
    CustomGeometryGroup,
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
from ._shared.vector_in_angstroms import VectorInAngstroms
from .angstrom_to_world import AngstromToWorld
from .chain_id import ChainID
from .color import Color
from .color_res_name import ColorResName
from .edge_info import EdgeInfo
from .fallback_geometry import FallbackGeometry
from .integer_distance import IntegerDistance
from .oxdna_vectors import OxDNAVectors
from .set_color import SetColor
from .world_to_angstrom import WorldToAngstrom


class OxDNAAreIDs53(CustomGeometryGroup):
    _name = "oxDNA Are IDs 5'→3'"
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Extrapolate whether indices are assigned in 5'→3' order. Only works for relaxed DNA."
    }

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry("Atoms")
        result = tree.outputs.boolean("Result")

        chain_id = ChainID()
        index = g.Index()
        compare = g.Compare.integer.equal(
            chain_id.o.chain_id.point.at(index),
            chain_id.o.chain_id.point.at(index.o.index + 1),
        )
        vector_math = (
            g.Position().o.position.point.at(g.Index().o.index + 1) - g.Position()
        ).dot(OxDNAVectors().o.base_normal)
        (
            (g.AttributeStatistic.point.float(atoms, compare, vector_math).o.mean < 0.0)
            >> result
        )


class OxDNAStyleClassic(AssetGeometryGroup):
    """
    oxDNA Style Classic

    Parameters
    ----------
    atoms : InputGeometry
        Atoms
    selection : InputBoolean
        Selection of atoms to apply this node to
    quality : InputInteger
        Quality
    backbone_shape : InputMenu | Literal["Sticks", "Ribbon"]
        Backbone Shape
    backbone_radius : InputFloat
        Backbone Radius
    ball_radius : InputFloat
        Ball Radius
    backbone_taper : InputFloat
        Backbone Taper
    end_overhang : InputFloat
        End Overhang
    base_shape : InputMenu | Literal["Sphere", "Cylinder", "None"]
        Base Shape
    base_geometry : InputGeometry
        Base Geometry
    base_scale : InputVector
        Base Scale
    stem_geometry : InputGeometry
        Stem Geometry
    stem_scale : InputVector
        Stem Scale
    base_colors : InputMenu | Literal["Uniform", "Specific", "Strand"]
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
    strand_colors : InputMenu | Literal["Uniform", "Specific", "Auto"]
        Strand Colors
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
    i.backbone_taper : FloatSocket
        Backbone Taper
    i.end_overhang : FloatSocket
        End Overhang
    i.base_shape : MenuSocket
        Base Shape
    i.base_geometry : GeometrySocket
        Base Geometry
    i.base_scale : VectorSocket
        Base Scale
    i.stem_geometry : GeometrySocket
        Stem Geometry
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
    i.strand_colors : MenuSocket
        Strand Colors
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

    _name = "oxDNA Style Classic"
    _asset_name = "oxDNA Style Classic"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "default_group_node_width": 160,
        "node_tool_idname": "geometry.mn_oxdna_style_classic",
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
        backbone_taper: FloatSocket
        """Backbone Taper"""
        end_overhang: FloatSocket
        """End Overhang"""
        base_shape: MenuSocket
        """Base Shape"""
        base_geometry: GeometrySocket
        """Base Geometry"""
        base_scale: VectorSocket
        """Base Scale"""
        stem_geometry: GeometrySocket
        """Stem Geometry"""
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
        strand_colors: MenuSocket
        """Strand Colors"""
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
        backbone_shape: InputMenu | Literal["Sticks", "Ribbon"] = "Sticks",
        backbone_radius: InputFloat = 1.1,
        ball_radius: InputFloat = 2.3,
        backbone_taper: InputFloat = 0.0,
        end_overhang: InputFloat = 1.8,
        base_shape: InputMenu | Literal["Sphere", "Cylinder", "None"] = "Sphere",
        base_geometry: InputGeometry = None,
        base_scale: InputVector = None,
        stem_geometry: InputGeometry = None,
        stem_scale: InputVector = None,
        base_colors: InputMenu | Literal["Uniform", "Specific", "Strand"] = "Uniform",
        bases: InputColor = None,
        a: InputColor = None,
        c: InputColor = None,
        g: InputColor = None,
        t_u: InputColor = None,
        strand_colors: InputMenu | Literal["Uniform", "Specific", "Auto"] = "Auto",
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
                "Backbone Taper": backbone_taper,
                "End Overhang": end_overhang,
                "Base Shape": base_shape,
                "Base Geometry": base_geometry,
                "Base Scale": base_scale,
                "Stem Geometry": stem_geometry,
                "Stem Scale": stem_scale,
                "Base Colors": base_colors,
                "Bases": bases,
                "A": a,
                "C": c,
                "G": g,
                "T / U": t_u,
                "Strand Colors": strand_colors,
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
                1.1,
                min_value=0.0,
                max_value=340_282_000_000_000_000_000_000_000_000_000_000_000.0,
            )
            ball_radius = tree.inputs.float("Ball Radius", 2.3, min_value=0.0)
            backbone_taper = tree.inputs.float(
                "Backbone Taper", 0.0, min_value=0.0, max_value=1.0, subtype="FACTOR"
            )
            end_overhang = tree.inputs.float("End Overhang", 1.8)
        with tree.inputs.panel("Bases"):
            base_shape = tree.inputs.menu(
                "Base Shape", expanded=True, optional_label=True
            )
            base_geometry = tree.inputs.geometry("Base Geometry")
            base_scale = tree.inputs.vector(
                "Base Scale", (1.0, 2.3, 2.3), min_value=0.0, subtype="XYZ"
            )
            stem_geometry = tree.inputs.geometry("Stem Geometry")
            stem_scale = tree.inputs.vector(
                "Stem Scale", (1.1, 1.1, 1.0), min_value=0.0, subtype="XYZ"
            )
        with tree.inputs.panel("Base colors"):
            base_colors = tree.inputs.menu(
                "Base Colors", expanded=True, optional_label=True
            )
            bases = tree.inputs.color("Bases", (0.0, 1.0, 1.0, 1.0))
            a = tree.inputs.color("A", (0.033104, 0.03310406, 1.0, 1.0))
            c_ = tree.inputs.color("C", (0.033104, 1.0, 0.033104, 1.0))
            g_ = tree.inputs.color("G", (1.0, 1.0, 0.033104, 1.0))
            t_u = tree.inputs.color("T / U", (1.0, 0.033104, 0.033104, 1.0))
        with tree.inputs.panel("Strand colors"):
            strand_colors = tree.inputs.menu(
                "Strand Colors", expanded=True, optional_label=True
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
                "Material",
                bpy.data.materials.get("Default"),
                description="Material to apply to the resulting geometry",
            )
        geometry = tree.outputs.geometry("Geometry")

        separate_geometry = g.SeparateGeometry.point(atoms, selection)
        with g.Frame("Color strands if Auto-color is False"):
            index_switch = g.IndexSwitch.color(
                ChainID().o.chain_id.modulo(4), (strand_1, strand_2, strand_3, strand_4)
            )
            menu_switch = g.MenuSwitch.color(
                strand_colors,
                {
                    "Uniform": (strands, "Single uniform color"),
                    "Specific": (index_switch, "Set custom colors for the bases"),
                    "Auto": (Color(), "Use the existing `Color` attribute"),
                },
            )
            set_color = SetColor(
                atoms=separate_geometry.o.selection, color=menu_switch.o.output
            )
        oxdna_vectors = OxDNAVectors()
        vector_math = oxdna_vectors.o.stacking_offset - oxdna_vectors.o.backbone_offset
        axes_to_rotation = g.AxesToRotation(
            primary_axis=vector_math, secondary_axis=oxdna_vectors.o.base_normal
        )
        vector_math_1 = g.SampleIndex(
            geometry=atoms, value=vector_math, data_type="FLOAT_VECTOR"
        ).o.value.length()
        with g.Frame("Colored bases"):
            menu_switch_1 = g.MenuSwitch.color(
                base_colors,
                {
                    "Uniform": bases,
                    "Specific": ColorResName(
                        a=a, c=c_, g=g_, t=t_u, ra=a, rc=c_, rg=g_, ru=t_u
                    ),
                    "Strand": Color(),
                },
            )
            vector_math_2 = (
                oxdna_vectors.o.stacking_offset
                + vector_math.normalize() * ((stem_scale.z - 1.0) * vector_math_1)
            )
            instance_on_points = SetColor(
                atoms=SetInstancer(
                    geometry=g.SetPosition(geometry=set_color, offset=vector_math_2)
                ),
                color=menu_switch_1.o.output,
            ) >> g.InstanceOnPoints(
                instance=FallbackGeometry(
                    geometry=base_geometry, fallback=g.IcoSphere(subdivisions=quality)
                ),
                rotation=axes_to_rotation,
                scale=VectorInAngstroms(
                    vector=base_scale, normalize=False, angstrom=1.0
                ),
            )
        oxdna_vectors_1 = OxDNAVectors()
        set_position = g.SetPosition(
            geometry=set_color, offset=oxdna_vectors_1.o.backbone_offset
        )
        with g.Frame():
            with g.Frame("Backbone ball"):
                instance_on_points_1 = SetInstancer(
                    geometry=set_position
                ) >> g.InstanceOnPoints(
                    instance=g.IcoSphere(subdivisions=quality),
                    scale=AngstromToWorld(angstrom=ball_radius),
                )
            with g.Frame("Backbone Stick"):
                angstrom_to_world = AngstromToWorld(angstrom=backbone_radius)
                with g.Frame("Add overhang to strand ends"):
                    vector_math_3 = (
                        EdgeInfo(vertex_index=g.Index()).o.edge_vector.normalize()
                        * -1.0
                        * AngstromToWorld(angstrom=end_overhang)
                    )
                    set_position_1 = g.SetPosition(
                        geometry=set_position,
                        selection=g.Compare.integer.equal(
                            g.EdgesOfVertex(vertex_index=g.Index()).o.total, 1
                        ),
                        offset=vector_math_3,
                    )
                with g.Frame(
                    "Each edge to a curve pointing 5'->3'. Flip circular endpoints."
                ):
                    edge_vertices = g.EdgeVertices()
                    capture = g.CaptureAttribute.edge(
                        geometry=set_position,
                        selection=IntegerDistance(
                            a=edge_vertices.o.vertex_index_1,
                            b=edge_vertices.o.vertex_index_2,
                        ).o.cutoff,
                    )
                    reverse_curve = (
                        capture.o.geometry
                        >> g.SplitEdges()
                        >> g.MeshToCurve()
                        >> g.ReverseCurve(selection=capture.o.selection)
                    )
                    reverse_curve_1 = g.ReverseCurve(
                        curve=reverse_curve,
                        selection=~OxDNAAreIDs53(Atoms=reverse_curve).o.result,
                    )
                curve_circle = g.CurveCircle(resolution=quality * 4)
                switch = g.EndpointSelection(start_size=0).o.selection.switch.float(
                    angstrom_to_world,
                    angstrom_to_world.o.world - backbone_taper * angstrom_to_world,
                )
                set_spline_resolution = (
                    set_position_1
                    >> g.MeshToCurve()
                    >> g.SetCurveNormal(
                        normal=oxdna_vectors_1.o.base_normal, mode="Free"
                    )
                    >> g.SetSplineType.bezier()
                    >> g.SetSplineResolution(resolution=quality * 2)
                )
                curve_to_mesh = (
                    g.SetHandleType(curve=set_spline_resolution)
                    >> g.SetCurveNormal()
                    >> g.CurveToMesh(
                        profile_curve=curve_circle,
                        scale=angstrom_to_world,
                        fill_caps=True,
                    )
                )
                curve_to_mesh_1 = g.CurveToMesh(
                    curve=reverse_curve_1,
                    profile_curve=curve_circle,
                    scale=switch,
                    fill_caps=True,
                )
                menu_switch_2 = g.MenuSwitch.geometry(
                    backbone_shape,
                    {
                        "Sticks": g.JoinGeometry(
                            geometry=(curve_to_mesh_1, instance_on_points_1)
                        ),
                        "Ribbon": curve_to_mesh,
                    },
                )
        with g.Frame("Base stem"):
            world_to_angstrom = WorldToAngstrom(world=vector_math_1)
            transform_geometry = FallbackGeometry(
                geometry=stem_geometry,
                fallback=g.Cylinder(
                    vertices=quality * 5, side_segments=quality, depth=1.0
                ),
            ) >> g.TransformGeometry(
                translation=g.CombineXYZ(z=world_to_angstrom.o.angstrom * 0.5),
                scale=g.CombineXYZ(z=world_to_angstrom, x=1.0, y=1.0),
            )
            instance_on_points_2 = SetInstancer(
                geometry=set_position
            ) >> g.InstanceOnPoints(
                instance=transform_geometry,
                rotation=axes_to_rotation,
                scale=VectorInAngstroms(
                    vector=stem_scale, normalize=False, angstrom=1.0
                ),
            )
        menu_switch_3 = g.MenuSwitch.geometry(
            base_shape,
            {
                "Sphere": g.JoinGeometry(
                    geometry=(instance_on_points, instance_on_points_2)
                ),
                "Cylinder": instance_on_points_2,
                "None": None,
            },
        )
        set_shade_smooth = g.SetShadeSmooth.face(
            g.JoinGeometry(geometry=(menu_switch_3, menu_switch_2)),
            shade_smooth=shade_smooth,
        )
        set_shade_smooth.node.warning_propagation = "ERRORS"
        smooth_by_angle = SmoothByAngle(mesh=set_shade_smooth, angle=math.pi / 3)
        smooth_by_angle.node.warning_propagation = "ERRORS"
        smooth_by_angle >> g.SetMaterial(material=material) >> geometry

        backbone_shape.default_value = "Sticks"
        base_shape.default_value = "Sphere"
        base_colors.default_value = "Uniform"
        strand_colors.default_value = "Auto"


ASSET = OxDNAStyleClassic

ASSET_METADATA = {
    "catalog_id": "0094c3e0-7885-427b-81b4-187a84dcff18",
}

DATABLOCK_DEPENDENCIES = {
    "materials": ("Default",),
}
