# Node-group asset "Style HBonds" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
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
    ObjectSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import (
    InputBoolean,
    InputColor,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputMaterial,
    InputMenu,
    InputObject,
)
from ._shared.mn_units import MNUnits
from .find_hbonds import FindHBonds


class StyleHBonds(AssetGeometryGroup):
    """
    Draw hydrogen bonds stored as edges on a data object, with positions sampled from the atoms

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that the hydrogen bond positions are sampled from
    selection : InputBoolean
        Selection of atoms, a hydrogen bond is shown only when both of its atoms are selected
    source : InputMenu | Literal["Data Object", "Find HBonds"]
        Where the hydrogen bonds come from: a data object calculated with MDAnalysis, or found by the Find HBonds node
    data_object : InputObject
        Object holding the hydrogen bonds as edges, with one vertex per atom
    quality : InputInteger
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    radius : InputFloat
        Radius of the hydrogen bonds in Angstroms
    dashes : InputInteger
        Number of dashes along each hydrogen bond, 1 draws a solid cylinder
    color : InputColor
        Color of the hydrogen bonds
    distance : InputFloat
        Maximum distance between the donor and acceptor in Angstroms, when found by nodes
    angle : InputFloat
        Minimum donor-hydrogen-acceptor angle, when found by nodes
    shade_smooth : InputBoolean
        Apply smooth shading to the created geometry
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that the hydrogen bond positions are sampled from
    i.selection : BooleanSocket
        Selection of atoms, a hydrogen bond is shown only when both of its atoms are selected
    i.source : MenuSocket
        Where the hydrogen bonds come from: a data object calculated with MDAnalysis, or found by the Find HBonds node
    i.data_object : ObjectSocket
        Object holding the hydrogen bonds as edges, with one vertex per atom
    i.quality : IntegerSocket
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    i.radius : FloatSocket
        Radius of the hydrogen bonds in Angstroms
    i.dashes : IntegerSocket
        Number of dashes along each hydrogen bond, 1 draws a solid cylinder
    i.color : ColorSocket
        Color of the hydrogen bonds
    i.distance : FloatSocket
        Maximum distance between the donor and acceptor in Angstroms, when found by nodes
    i.angle : FloatSocket
        Minimum donor-hydrogen-acceptor angle, when found by nodes
    i.shade_smooth : BooleanSocket
        Apply smooth shading to the created geometry
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        The generated geometry for the style node group
    """

    _name = "Style HBonds"
    _asset_name = "Style HBonds"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Draw hydrogen bonds stored as edges on a data object, with positions sampled from the atoms",
        "node_tool_idname": "geometry.style_hbonds",
        "show_modifier_manage_panel": False,
        "is_modifier": True,
    }

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that the hydrogen bond positions are sampled from"""
        selection: BooleanSocket
        """Selection of atoms, a hydrogen bond is shown only when both of its atoms are selected"""
        source: MenuSocket
        """Where the hydrogen bonds come from: a data object calculated with MDAnalysis, or found by the Find HBonds node"""
        data_object: ObjectSocket
        """Object holding the hydrogen bonds as edges, with one vertex per atom"""
        quality: IntegerSocket
        """A lower value results in less geometry, with a higher value meaning better looking but more dense geometry"""
        radius: FloatSocket
        """Radius of the hydrogen bonds in Angstroms"""
        dashes: IntegerSocket
        """Number of dashes along each hydrogen bond, 1 draws a solid cylinder"""
        color: ColorSocket
        """Color of the hydrogen bonds"""
        distance: FloatSocket
        """Maximum distance between the donor and acceptor in Angstroms, when found by nodes"""
        angle: FloatSocket
        """Minimum donor-hydrogen-acceptor angle, when found by nodes"""
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
        source: InputMenu | Literal["Data Object", "Find HBonds"] = "Data Object",
        data_object: InputObject = None,
        quality: InputInteger = 2,
        radius: InputFloat = 0.08,
        dashes: InputInteger = 4,
        color: InputColor = None,
        distance: InputFloat = 3.0,
        angle: InputFloat = 5 * math.pi / 6,
        shade_smooth: InputBoolean = True,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Source": source,
                "Data Object": data_object,
                "Quality": quality,
                "Radius": radius,
                "Dashes": dashes,
                "Color": color,
                "Distance": distance,
                "Angle": angle,
                "Shade Smooth": shade_smooth,
                "Material": material,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms",
            description="Atomic geometry that the hydrogen bond positions are sampled from",
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms, a hydrogen bond is shown only when both of its atoms are selected",
            hide_value=True,
        )
        source = tree.inputs.menu(
            "Source",
            "Data Object",
            description="Where the hydrogen bonds come from: a data object calculated with MDAnalysis, or found by the Find HBonds node",
            expanded=True,
        )
        data_object = tree.inputs.object(
            "Data Object",
            description="Object holding the hydrogen bonds as edges, with one vertex per atom",
            optional_label=True,
        )
        quality = tree.inputs.integer(
            "Quality",
            2,
            description="A lower value results in less geometry, with a higher value meaning better looking but more dense geometry",
            min_value=0,
            max_value=5,
        )
        radius = tree.inputs.float(
            "Radius",
            0.08,
            description="Radius of the hydrogen bonds in Angstroms",
            min_value=0.0,
            max_value=1.0,
        )
        dashes = tree.inputs.integer(
            "Dashes",
            4,
            description="Number of dashes along each hydrogen bond, 1 draws a solid cylinder",
            min_value=1,
            max_value=50,
        )
        color = tree.inputs.color(
            "Color", (0.8, 0.8, 0.8, 1.0), description="Color of the hydrogen bonds"
        )
        with tree.inputs.panel("Find HBonds", default_closed=True):
            distance = tree.inputs.float(
                "Distance",
                3.0,
                description="Maximum distance between the donor and acceptor in Angstroms, when found by nodes",
                min_value=0.0,
                max_value=10.0,
            )
            angle = tree.inputs.float(
                "Angle",
                5 * math.pi / 6,
                description="Minimum donor-hydrogen-acceptor angle, when found by nodes",
                min_value=0.0,
                max_value=math.pi,
                subtype="ANGLE",
            )
        with tree.inputs.panel("Material", default_closed=True):
            shade_smooth = tree.inputs.boolean(
                "Shade Smooth",
                True,
                description="Apply smooth shading to the created geometry",
            )
            material = tree.inputs.material(
                "Material",
                description="Material to apply to the resulting geometry",
                optional_label=True,
            )
        geometry = tree.outputs.geometry(
            "Geometry", description="The generated geometry for the style node group"
        )

        with g.Frame("Unit length dashed cylinder along Z"):
            math_1 = 1.0 / dashes
            cylinder = g.Cylinder(
                vertices=(quality * 4.0).max(3.0),
                depth=g.Compare.integer.equal(dashes, 1).o.result.switch.float(
                    math_1 * 0.5, math_1
                ),
            )
            realize_instances = (
                g.MeshLine(
                    count=dashes,
                    start_location=g.CombineXYZ(z=math_1 * 0.5 - 0.5),
                    offset=g.CombineXYZ(z=math_1),
                )
                >> g.InstanceOnPoints(instance=cylinder)
                >> g.RealizeInstances()
            )
        edge_vertices = g.EdgeVertices()
        with g.Frame("Find with nodes"):
            find_hbonds = FindHBonds(
                atoms=atoms, selection=selection, distance=distance, angle=angle
            )
            _string = g.String(
                "Points at the middle of each bond from a hydrogen to its acceptor, kept when both are selected, matching the edges from a data object."
            )
            capture = g.CaptureAttribute.point(geometry=atoms)
            vector = capture.items.vector(
                g.Position().o.position.point.at(find_hbonds.o.acceptor_index)
                - g.Position(),
                "Vector",
            )
            set_position = (
                capture.o.geometry
                >> g.SeparateGeometry.point(
                    selection=find_hbonds.o.is_hbonded
                    & selection.point.at(find_hbonds.o.acceptor_index)
                )
                >> g.MeshToPoints()
                >> g.SetPosition(offset=vector.output * 0.5)
            )
        with g.Frame("Sample from the atoms"):
            set_position_1 = g.SetPosition(
                geometry=g.ObjectInfo(data_object).o.geometry,
                position=g.SampleIndex.point.vector(atoms, g.Position(), g.Index()),
            )
            capture_1 = g.CaptureAttribute.point(geometry=set_position_1)
            selection_1 = capture_1.items.boolean(
                g.SampleIndex.point.boolean(atoms, selection, g.Index()), "Selection"
            )
            separate_geometry = capture_1.o.geometry >> g.SeparateGeometry.edge(
                selection=selection_1.output
            )
            _string_1 = g.String(
                "The data object has one vertex per atom in the same order, so positions and the selection are sampled from the atoms by index. The selection is captured on the points so the index is the atom index, then an edge is kept only when both of its atoms are selected."
            )
        mn_units = MNUnits(value=radius)
        menu_switch = g.MenuSwitch.boolean(
            source, {"Data Object": False, "Find HBonds": True}
        )
        capture_2 = g.CaptureAttribute.edge(geometry=separate_geometry.o.selection)
        vector_1 = capture_2.items.vector(
            edge_vertices.o.position_2 - edge_vertices.o.position_1, "Vector"
        )
        switch = menu_switch.o.output.switch.vector(vector_1.output, vector.output)
        (
            menu_switch.o.output.switch.geometry(
                capture_2.o.geometry >> g.MeshToPoints.edges(), set_position
            )
            >> g.InstanceOnPoints(
                instance=realize_instances,
                rotation=g.AlignRotationToVector(vector=switch),
                scale=g.CombineXYZ(
                    x=mn_units.o.angstrom, y=mn_units.o.angstrom, z=switch.length()
                ),
            )
            >> g.RealizeInstances()
            >> g.StoreNamedAttribute.point.color(name="Color", value=color)
            >> g.SetShadeSmooth.face(shade_smooth=shade_smooth)
            >> g.SetMaterial(material=material)
            >> geometry
        )


ASSET = StyleHBonds

ASSET_METADATA = {
    "description": "Draw hydrogen bonds stored as edges on a data object, with positions sampled from the atoms",
    "catalog_id": "541e6649-2ea6-4225-b1ee-5c0da6f5f1f6",
}
