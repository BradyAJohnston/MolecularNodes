# Node-group asset 'oxDNA Style Ball and Stick' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
import bpy
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
from ._shared.mn_units import MNUnits
from ._shared.set_instancer import SetInstancer
from .angstrom_to_world import AngstromToWorld
from .chain_id import ChainID
from .color import Color
from .color_res_name import ColorResName
from .integer_distance import IntegerDistance
from .oxdna_normal import OxDNANormal
from .oxdna_offset import OxDNAOffset
from .oxdna_vector import OxDNAVector
from .set_color import SetColor


class OxDNAStyleBallAndStick(AssetGeometryGroup):
    """
    oxDNA Style Ball and Stick

    Parameters
    ----------
    atoms : InputGeometry
        Atoms
    selection : InputBoolean
        Selection of atoms to apply this node to
    resolution : InputInteger
        Resolution
    ball_radius : InputFloat
        Ball Radius
    stick_radius : InputFloat
        Stick Radius
    stick_taper : InputFloat
        Stick Taper
    sphere : InputMenu | Literal["Points", "Instances"]
        Sphere
    stem_radius : InputFloat
        Stem Radius
    base_scale : InputVector
        Base Scale
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
    i.resolution : IntegerSocket
        Resolution
    i.ball_radius : FloatSocket
        Ball Radius
    i.stick_radius : FloatSocket
        Stick Radius
    i.stick_taper : FloatSocket
        Stick Taper
    i.sphere : MenuSocket
        Sphere
    i.stem_radius : FloatSocket
        Stem Radius
    i.base_scale : VectorSocket
        Base Scale
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

    _name = "oxDNA Style Ball and Stick"
    _asset_name = "oxDNA Style Ball and Stick"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
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
        resolution: IntegerSocket
        """Resolution"""
        ball_radius: FloatSocket
        """Ball Radius"""
        stick_radius: FloatSocket
        """Stick Radius"""
        stick_taper: FloatSocket
        """Stick Taper"""
        sphere: MenuSocket
        """Sphere"""
        stem_radius: FloatSocket
        """Stem Radius"""
        base_scale: VectorSocket
        """Base Scale"""
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
        resolution: InputInteger = 2,
        ball_radius: InputFloat = 1.0,
        stick_radius: InputFloat = 1.0,
        stick_taper: InputFloat = 0.3,
        sphere: InputMenu | Literal["Points", "Instances"] = "Instances",
        stem_radius: InputFloat = 1.0,
        base_scale: InputVector = None,
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
                "Resolution": resolution,
                "Ball Radius": ball_radius,
                "Stick Radius": stick_radius,
                "Stick Taper": stick_taper,
                "Sphere": sphere,
                "Stem Radius": stem_radius,
                "Base Scale": base_scale,
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

    def _build_group(self, tree):
        atoms = tree.inputs.geometry("Atoms")
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        resolution = tree.inputs.integer("Resolution", 2, min_value=1)
        with tree.inputs.panel("Backbone"):
            ball_radius = tree.inputs.float("Ball Radius", 1.0, min_value=0.0)
            stick_radius = tree.inputs.float(
                "Stick Radius",
                1.0,
                min_value=0.0,
                max_value=340_282_000_000_000_000_000_000_000_000_000_000_000.0,
            )
            stick_taper = tree.inputs.float(
                "Stick Taper", 0.3, min_value=0.0, max_value=1.0, subtype="FACTOR"
            )
        with tree.inputs.panel("Base and stem"):
            sphere = tree.inputs.menu("Sphere", expanded=True, optional_label=True)
            stem_radius = tree.inputs.float("Stem Radius", 1.0, min_value=0.0)
            base_scale = tree.inputs.vector(
                "Base Scale", (1.0, 1.0, 1.0), min_value=0.0, subtype="XYZ"
            )
        with tree.inputs.panel("Colors"):
            base_colors = tree.inputs.menu(
                "Base Colors", expanded=True, optional_label=True
            )
            bases = tree.inputs.color("Bases", (0.0, 524_941.2, 524_940.5, 1.0))
            a = tree.inputs.color("A", (0.0, 0.0, 1.0, 1.0))
            c_ = tree.inputs.color("C", (0.0, 1.0, 0.0, 1.0))
            g_ = tree.inputs.color("G", (1.0, 1.0, 0.0, 1.0))
            t_u = tree.inputs.color("T / U", (1.0, 0.0, 0.0, 1.0))
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
                "Material",
                bpy.data.materials.get("MN oxDNA Default"),
                description="Material to apply to the resulting geometry",
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
        with g.Frame("Recover oxDNA backbone, model proportions"):
            group_1 = OxDNAOffset()
            capture = g.CaptureAttribute.point(geometry=group)
            value = capture.items.vector("Value", group_1)
            set_position = g.SetPosition(
                geometry=capture.o.geometry, offset=value.output
            )
            vector_math = OxDNAVector().o.base_vector * 0.34
            vector_math_1 = vector_math - group_1
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
            transform_geometry = g.TransformGeometry(
                geometry=g.IcoSphere(subdivisions=resolution, radius=0.0414),
                scale=base_scale * (2.1, 5.1, 5.1),
            )
            instance_on_points = (
                SetColor(
                    atoms=SetInstancer(geometry=capture.o.geometry),
                    color=menu_switch_1.o.output,
                )
                >> g.SetPosition(offset=vector_math)
                >> g.InstanceOnPoints(
                    instance=transform_geometry,
                    rotation=g.AlignRotationToVector(vector=OxDNANormal(), axis="X"),
                )
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
            group_2 = AngstromToWorld(angstrom=stick_radius)
            switch = g.EndpointSelection(start_size=0).o.selection.switch.float(
                group_2, group_2.o.world * stick_taper
            )
            curve_to_mesh = reverse_curve >> g.CurveToMesh(
                profile_curve=g.CurveCircle(resolution=6), scale=switch
            )
        with g.Frame("Backbone ball"):
            math_1 = AngstromToWorld(angstrom=ball_radius).o.world * 2.0
            instance_on_points_1 = SetInstancer(
                geometry=set_position
            ) >> g.InstanceOnPoints(
                instance=g.IcoSphere(radius=math_1, subdivisions=resolution)
            )
            menu_switch_2 = g.MenuSwitch.geometry(
                sphere,
                {
                    "Points": g.MeshToPoints(mesh=set_position, radius=math_1),
                    "Instances": instance_on_points_1,
                },
            )
        with g.Frame("Base stem"):
            cylinder = g.Cylinder(
                vertices=resolution * 7,
                side_segments=resolution,
                radius=MNUnits(value=stem_radius).o.angstrom,
                depth=1.0,
            )
            instance_on_points_2 = SetInstancer(
                geometry=set_position
            ) >> g.InstanceOnPoints(
                instance=g.TransformGeometry(
                    geometry=cylinder, translation=(0.0, 0.0, 0.5)
                ),
                rotation=g.AlignRotationToVector(vector=vector_math_1),
                scale=(0.63, 0.63, 0.75),
            )
        set_shade_smooth = g.JoinGeometry(
            geometry=(
                curve_to_mesh,
                menu_switch_2,
                instance_on_points_2,
                instance_on_points,
            )
        ) >> g.SetShadeSmooth.face(shade_smooth=shade_smooth)
        set_shade_smooth.node.warning_propagation = "ERRORS"
        set_shade_smooth >> g.SetMaterial(material=material) >> geometry

        sphere.default_value = "Instances"
        base_colors.default_value = "Uniform"
        strand_color.default_value = "Color"


ASSET = OxDNAStyleBallAndStick

ASSET_METADATA = {
    "catalog_id": "0094c3e0-7885-427b-81b4-187a84dcff18",
}

DATABLOCK_DEPENDENCIES = {
    "materials": ("MN oxDNA Default",),
}
