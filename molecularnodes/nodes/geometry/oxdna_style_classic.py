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


class OxDNADirection(CustomGeometryGroup):
    _name = ".oxDNA Direction"
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Infer whether indices are assigned in 3'→5' order (True) or 5'→3' order (False) using base normals. Only accurate for relaxed systems."
    }

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms",
            description="Vertices and edges representing nucleotides and phosphodiester bonds, respectively",
        )
        result = tree.outputs.boolean(
            "Result",
            description="True if nucleotides are indexed in 3'→5' order, False if 5'→3' order",
        )

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
            (g.AttributeStatistic.point.float(atoms, compare, vector_math).o.mean > 0.0)
            >> result
        )


class OxDNAStyleClassic(AssetGeometryGroup):
    """
    oxDNA Style Classic

    Parameters
    ----------
    atoms : InputGeometry
        Vertices and edges representing nucleotides and phosphodiester bonds, respectively
    selection : InputBoolean
        Selection of atoms to apply this node to
    quality : InputInteger
        A lower value results in less geometry, while a higher value means better-looking but more dense geometry
    backbone_shape : InputMenu | Literal["Sticks", "Ribbon"]
        The visual style of the backbone
    backbone_radius : InputFloat
        Radius of the backbone sticks or ribbons
    ball_radius : InputFloat
        Radius of the backbone spheres
    backbone_taper : InputFloat
        Taper the backbone sticks so they point in the 5'→3' direction
    end_overhang : InputFloat
        Extend the backbone ribbon past the final nucleotides on each strand
    base_shape : InputMenu | Literal["Sphere", "Cylinder", "None"]
        Visual style of the bases
    base_geometry : InputGeometry
        Render the bases using custom geometry
    base_scale : InputVector
        Scale the bases along each axis
    stem_geometry : InputGeometry
        Render the base stems using custom geometry
    stem_scale : InputVector
        Scale the base stems along each axis
    base_colors : InputMenu | Literal["Uniform", "Specific", "Strand"]
        Method used to determine base colors
    bases : InputColor
        Color all bases uniformly
    a : InputColor
        Color adenines
    c : InputColor
        Color cytosines
    g : InputColor
        Color guanines
    t_u : InputColor
        Color thymines/uracils
    strand_colors : InputMenu | Literal["Uniform", "Specific", "Auto"]
        Method used to determine strand colors. `Auto` will color strands using a finite palette of distinct pastel hues.
    strands : InputColor
        Color all strands uniformly
    strand_1 : InputColor
        Color the first strand, and every 4th strand after that
    strand_2 : InputColor
        Color the second strand, and every 4th strand after that
    strand_3 : InputColor
        Color the third strand, and every 4th strand after that
    strand_4 : InputColor
        Color the fourth strand, and every 4th strand after that
    shade_smooth : InputBoolean
        Apply smooth shading to the created geometry
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.atoms : GeometrySocket
        Vertices and edges representing nucleotides and phosphodiester bonds, respectively
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.quality : IntegerSocket
        A lower value results in less geometry, while a higher value means better-looking but more dense geometry
    i.backbone_shape : MenuSocket
        The visual style of the backbone
    i.backbone_radius : FloatSocket
        Radius of the backbone sticks or ribbons
    i.ball_radius : FloatSocket
        Radius of the backbone spheres
    i.backbone_taper : FloatSocket
        Taper the backbone sticks so they point in the 5'→3' direction
    i.end_overhang : FloatSocket
        Extend the backbone ribbon past the final nucleotides on each strand
    i.base_shape : MenuSocket
        Visual style of the bases
    i.base_geometry : GeometrySocket
        Render the bases using custom geometry
    i.base_scale : VectorSocket
        Scale the bases along each axis
    i.stem_geometry : GeometrySocket
        Render the base stems using custom geometry
    i.stem_scale : VectorSocket
        Scale the base stems along each axis
    i.base_colors : MenuSocket
        Method used to determine base colors
    i.bases : ColorSocket
        Color all bases uniformly
    i.a : ColorSocket
        Color adenines
    i.c : ColorSocket
        Color cytosines
    i.g : ColorSocket
        Color guanines
    i.t_u : ColorSocket
        Color thymines/uracils
    i.strand_colors : MenuSocket
        Method used to determine strand colors. `Auto` will color strands using a finite palette of distinct pastel hues.
    i.strands : ColorSocket
        Color all strands uniformly
    i.strand_1 : ColorSocket
        Color the first strand, and every 4th strand after that
    i.strand_2 : ColorSocket
        Color the second strand, and every 4th strand after that
    i.strand_3 : ColorSocket
        Color the third strand, and every 4th strand after that
    i.strand_4 : ColorSocket
        Color the fourth strand, and every 4th strand after that
    i.shade_smooth : BooleanSocket
        Apply smooth shading to the created geometry
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        The generated geometry for the style node group
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
        """Vertices and edges representing nucleotides and phosphodiester bonds, respectively"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        quality: IntegerSocket
        """A lower value results in less geometry, while a higher value means better-looking but more dense geometry"""
        backbone_shape: MenuSocket
        """The visual style of the backbone"""
        backbone_radius: FloatSocket
        """Radius of the backbone sticks or ribbons"""
        ball_radius: FloatSocket
        """Radius of the backbone spheres"""
        backbone_taper: FloatSocket
        """Taper the backbone sticks so they point in the 5'→3' direction"""
        end_overhang: FloatSocket
        """Extend the backbone ribbon past the final nucleotides on each strand"""
        base_shape: MenuSocket
        """Visual style of the bases"""
        base_geometry: GeometrySocket
        """Render the bases using custom geometry"""
        base_scale: VectorSocket
        """Scale the bases along each axis"""
        stem_geometry: GeometrySocket
        """Render the base stems using custom geometry"""
        stem_scale: VectorSocket
        """Scale the base stems along each axis"""
        base_colors: MenuSocket
        """Method used to determine base colors"""
        bases: ColorSocket
        """Color all bases uniformly"""
        a: ColorSocket
        """Color adenines"""
        c: ColorSocket
        """Color cytosines"""
        g: ColorSocket
        """Color guanines"""
        t_u: ColorSocket
        """Color thymines/uracils"""
        strand_colors: MenuSocket
        """Method used to determine strand colors. `Auto` will color strands using a finite palette of distinct pastel hues."""
        strands: ColorSocket
        """Color all strands uniformly"""
        strand_1: ColorSocket
        """Color the first strand, and every 4th strand after that"""
        strand_2: ColorSocket
        """Color the second strand, and every 4th strand after that"""
        strand_3: ColorSocket
        """Color the third strand, and every 4th strand after that"""
        strand_4: ColorSocket
        """Color the fourth strand, and every 4th strand after that"""
        shade_smooth: BooleanSocket
        """Apply smooth shading to the created geometry"""
        material: MaterialSocket
        """Material to apply to the resulting geometry"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """The generated geometry for the style node group"""

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
        atoms = tree.inputs.geometry(
            "Atoms",
            description="Vertices and edges representing nucleotides and phosphodiester bonds, respectively",
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        quality = tree.inputs.integer(
            "Quality",
            2,
            description="A lower value results in less geometry, while a higher value means better-looking but more dense geometry",
            min_value=1,
        )
        with tree.inputs.panel("Backbone"):
            backbone_shape = tree.inputs.menu(
                "Backbone Shape",
                description="The visual style of the backbone",
                expanded=True,
                optional_label=True,
            )
            backbone_radius = tree.inputs.float(
                "Backbone Radius",
                1.1,
                description="Radius of the backbone sticks or ribbons",
                min_value=0.0,
                max_value=340_282_000_000_000_000_000_000_000_000_000_000_000.0,
            )
            ball_radius = tree.inputs.float(
                "Ball Radius",
                2.3,
                description="Radius of the backbone spheres",
                min_value=0.0,
            )
            backbone_taper = tree.inputs.float(
                "Backbone Taper",
                0.0,
                description="Taper the backbone sticks so they point in the 5'→3' direction",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
            end_overhang = tree.inputs.float(
                "End Overhang",
                1.8,
                description="Extend the backbone ribbon past the final nucleotides on each strand",
            )
        with tree.inputs.panel("Bases"):
            base_shape = tree.inputs.menu(
                "Base Shape",
                description="Visual style of the bases",
                expanded=True,
                optional_label=True,
            )
            base_geometry = tree.inputs.geometry(
                "Base Geometry", description="Render the bases using custom geometry"
            )
            base_scale = tree.inputs.vector(
                "Base Scale",
                (1.0, 2.3, 2.3),
                description="Scale the bases along each axis",
                min_value=0.0,
                subtype="XYZ",
            )
            stem_geometry = tree.inputs.geometry(
                "Stem Geometry",
                description="Render the base stems using custom geometry",
            )
            stem_scale = tree.inputs.vector(
                "Stem Scale",
                (1.1, 1.1, 1.0),
                description="Scale the base stems along each axis",
                min_value=0.0,
                subtype="XYZ",
            )
        with tree.inputs.panel("Base colors"):
            base_colors = tree.inputs.menu(
                "Base Colors",
                description="Method used to determine base colors",
                expanded=True,
                optional_label=True,
            )
            bases = tree.inputs.color(
                "Bases", (0.0, 1.0, 1.0, 1.0), description="Color all bases uniformly"
            )
            a = tree.inputs.color(
                "A", (0.033104, 0.03310406, 1.0, 1.0), description="Color adenines"
            )
            c_ = tree.inputs.color(
                "C", (0.033104, 1.0, 0.033104, 1.0), description="Color cytosines"
            )
            g_ = tree.inputs.color(
                "G", (1.0, 1.0, 0.033104, 1.0), description="Color guanines"
            )
            t_u = tree.inputs.color(
                "T / U",
                (1.0, 0.033104, 0.033104, 1.0),
                description="Color thymines/uracils",
            )
        with tree.inputs.panel("Strand colors"):
            strand_colors = tree.inputs.menu(
                "Strand Colors",
                description="Method used to determine strand colors. `Auto` will color strands using a finite palette of distinct pastel hues.",
                expanded=True,
                optional_label=True,
            )
            strands = tree.inputs.color(
                "Strands",
                (0.8, 0.8, 0.8, 1.0),
                description="Color all strands uniformly",
            )
            strand_1 = tree.inputs.color(
                "Strand 1",
                (1.0, 0.0, 0.0, 1.0),
                description="Color the first strand, and every 4th strand after that",
            )
            strand_2 = tree.inputs.color(
                "Strand 2",
                (0.0, 0.0, 1.0, 1.0),
                description="Color the second strand, and every 4th strand after that",
            )
            strand_3 = tree.inputs.color(
                "Strand 3",
                (0.0, 1.0, 0.0, 1.0),
                description="Color the third strand, and every 4th strand after that",
            )
            strand_4 = tree.inputs.color(
                "Strand 4",
                (1.0, 1.0, 0.0, 1.0),
                description="Color the fourth strand, and every 4th strand after that",
            )
        with tree.inputs.panel("Material"):
            shade_smooth = tree.inputs.boolean(
                "Shade Smooth",
                True,
                description="Apply smooth shading to the created geometry",
            )
            material = tree.inputs.material(
                "Material",
                bpy.data.materials.get("Default"),
                description="Material to apply to the resulting geometry",
            )
        geometry = tree.outputs.geometry(
            "Geometry", description="The generated geometry for the style node group"
        )

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
                    set_position_1 = g.SetPosition(
                        geometry=set_position,
                        selection=g.Compare.integer.equal(g.EdgesOfVertex().o.total, 1),
                        offset=EdgeInfo().o.edge_vector.normalize()
                        * -1.0
                        * AngstromToWorld(angstrom=end_overhang),
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
                        selection=OxDNADirection(Atoms=reverse_curve),
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
    "description": "Render oxDNA nucleotides with a ball-and-stick or ribbon style backbone",
    "catalog_id": "0094c3e0-7885-427b-81b4-187a84dcff18",
}

DATABLOCK_DEPENDENCIES = {
    "materials": ("Default",),
}
