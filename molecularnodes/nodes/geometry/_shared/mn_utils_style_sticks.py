# Node group '.MN_utils_style_sticks' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING, Literal
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    MaterialSocket,
    MenuSocket,
    SocketAccessor,
)
from nodebpy.types import (
    InputBoolean,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputMaterial,
    InputMenu,
)
from .is_extra_bonds import IsExtraBonds
from .map_bond_type import Map_bond_type
from .mn_units import MNUnits
from .sample_atomic_attributes_to_face_corner import SampleAtomicAttributesToFaceCorner


class MN_utils_style_sticks(CustomGeometryGroup):
    """
    .MN_utils_style_sticks

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this node to
    radius : InputFloat
        Radius of the bond mesh.
    resolution : InputInteger
        Resolution of the created bond cylinders.
    scale_extra_bond_radius : InputFloat
        Scale Extra Bond Radius
    extra_bond_offset : InputFloat
        Extra Bond Offset
    extra_bond_rotate : InputFloat
        Extra Bond Rotate
    menu : InputMenu | Literal["Single", "Double"]
        Menu
    shade_smooth : InputBoolean
        Apply smooth shading to the created geometry
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.radius : FloatSocket
        Radius of the bond mesh.
    i.resolution : IntegerSocket
        Resolution of the created bond cylinders.
    i.scale_extra_bond_radius : FloatSocket
        Scale Extra Bond Radius
    i.extra_bond_offset : FloatSocket
        Extra Bond Offset
    i.extra_bond_rotate : FloatSocket
        Extra Bond Rotate
    i.menu : MenuSocket
        Menu
    i.shade_smooth : BooleanSocket
        Apply smooth shading to the created geometry
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = ".MN_utils_style_sticks"
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry._mn_utils_style_sticks"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        radius: FloatSocket
        """Radius of the bond mesh."""
        resolution: IntegerSocket
        """Resolution of the created bond cylinders."""
        scale_extra_bond_radius: FloatSocket
        """Scale Extra Bond Radius"""
        extra_bond_offset: FloatSocket
        """Extra Bond Offset"""
        extra_bond_rotate: FloatSocket
        """Extra Bond Rotate"""
        menu: MenuSocket
        """Menu"""
        shade_smooth: BooleanSocket
        """Apply smooth shading to the created geometry"""
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
        radius: InputFloat = 0.3,
        resolution: InputInteger = 6,
        scale_extra_bond_radius: InputFloat = 0.65,
        extra_bond_offset: InputFloat = 0.4,
        extra_bond_rotate: InputFloat = 0.0,
        menu: InputMenu | Literal["Single", "Double"] = "Single",
        shade_smooth: InputBoolean = True,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Radius": radius,
                "Resolution": resolution,
                "Scale Extra Bond Radius": scale_extra_bond_radius,
                "Extra Bond Offset": extra_bond_offset,
                "Extra Bond Rotate": extra_bond_rotate,
                "Menu": menu,
                "Shade Smooth": shade_smooth,
                "Material": material,
            }
        )

    def _build_group(self, tree):
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        radius = tree.inputs.float(
            "Radius",
            0.3,
            description="Radius of the bond mesh.",
            min_value=0.0,
            max_value=1.0,
        )
        resolution = tree.inputs.integer(
            "Resolution",
            6,
            description="Resolution of the created bond cylinders.",
            min_value=3,
            max_value=512,
        )
        scale_extra_bond_radius = tree.inputs.float(
            "Scale Extra Bond Radius", 0.65, min_value=-10_000.0, max_value=10_000.0
        )
        extra_bond_offset = tree.inputs.float(
            "Extra Bond Offset", 0.4, min_value=-10_000.0, max_value=10_000.0
        )
        extra_bond_rotate = tree.inputs.float(
            "Extra Bond Rotate", 0.0, min_value=-10_000.0, max_value=10_000.0
        )
        menu = tree.inputs.menu("Menu", optional_label=True)
        with tree.inputs.panel("Material"):
            shade_smooth = tree.inputs.boolean(
                "Shade Smooth",
                True,
                description="Apply smooth shading to the created geometry",
            )
            material = tree.inputs.material(
                "Material", description="Material to apply to the resulting geometry"
            )
        geometry = tree.outputs.geometry("Geometry")

        separate_geometry = g.SeparateGeometry.point(atoms, selection)
        with g.Frame("Offset Vector"):
            with g.Frame("Vector for Extra Bond Offsetting"):
                edge_vertices = g.EdgeVertices()
                edge_vertices_1 = g.EdgeVertices()
                vector_math = edge_vertices.o.position_1 - edge_vertices.o.position_2
                evaluate_at_index = vector_math.edge.at(
                    g.EdgesOfVertex(
                        vertex_index=edge_vertices_1.o.vertex_index_1
                    ).o.edge_index
                )
                evaluate_at_index_1 = vector_math.edge.at(
                    g.EdgesOfVertex(
                        vertex_index=edge_vertices_1.o.vertex_index_2, sort_index=2
                    ).o.edge_index
                )
                vector_math_1 = (
                    evaluate_at_index.normalize()
                    .cross(evaluate_at_index_1.normalize())
                    .normalize()
                    * MNUnits(value=extra_bond_offset).o.angstrom
                )
            capture = g.CaptureAttribute.edge(geometry=separate_geometry.o.selection)
            vector = capture.items.vector("Vector", vector_math_1)
        split_edges = capture.o.geometry >> g.SplitEdges()
        with g.Frame("Rotate and offset Extra Bonds"):
            edge_vertices_2 = g.EdgeVertices()
            group = Map_bond_type()
            capture_1 = g.CaptureAttribute.point(
                geometry=split_edges, selection=group > 1
            )
            duplicate_elements = capture_1.o.geometry >> g.DuplicateElements.edge(
                amount=group
            )
            math_1 = g.Math.multiply_add(
                duplicate_elements.o.duplicate_index,
                math.tau / Map_bond_type(),
                g.Math.divide(math.pi, 2.0).o.value + extra_bond_rotate,
            )
            vector_rotate = g.VectorRotate(
                vector=vector.output,
                axis=edge_vertices_2.o.position_1 - edge_vertices_2.o.position_2,
                angle=math_1,
            )
            set_position = duplicate_elements >> g.SetPosition(
                selection=capture_1.o.selection, offset=vector_rotate.o.vector * 1.0
            )
        menu_switch = g.MenuSwitch.geometry(
            menu, {"Single": split_edges, "Double": set_position}
        )
        with g.Frame("Compute Radius"):
            switch = (
                menu_switch.o.double & IsExtraBonds().o.is_extra_bonds
            ).switch.float(1.0, scale_extra_bond_radius)
            capture_2 = g.CaptureAttribute.edge(geometry=menu_switch)
            output = capture_2.items.float("Output", switch)
        set_curve_radius = (
            capture_2.o.geometry
            >> g.MeshToCurve()
            >> g.SetCurveRadius(radius=output.output)
        )
        with g.Frame("Get correct index to sample from"):
            capture_3 = g.CaptureAttribute.curve(geometry=set_curve_radius)
            first_point = capture_3.items.integer(
                "First Point", g.PointsOfCurve().o.point_index
            )
            last_point = capture_3.items.integer(
                "Last Point", g.PointsOfCurve(sort_index=1).o.point_index
            )
            curve_index = capture_3.items.integer("Curve Index", g.Index())
            switch_1 = (
                resolution
                < g.AccumulateField.face.integer(
                    group_index=curve_index.output
                ).o.leading
            ).switch.integer(first_point.output, last_point.output)
        with g.Frame("curve to mesh"):
            curve_to_mesh = g.CurveToMesh(
                curve=g.SubdivideCurve(curve=capture_3.o.geometry),
                profile_curve=g.CurveCircle(
                    resolution=resolution, radius=MNUnits(value=1.0).o.angstrom
                ),
                scale=g.Radius().o.radius * radius,
            )
            group_1 = SampleAtomicAttributesToFaceCorner(
                geometry=curve_to_mesh,
                sample_atoms=capture_3.o.geometry,
                index=switch_1,
            )
        with g.Frame("Set up materials"):
            (
                group_1
                >> g.SetMaterial(material=material)
                >> g.SetShadeSmooth.face(shade_smooth=shade_smooth)
                >> geometry
            )

        menu.default_value = "Single"
