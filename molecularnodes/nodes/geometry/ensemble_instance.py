# Node-group asset "Ensemble Instance" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    CollectionSocket,
    FloatSocket,
    GeometrySocket,
    MaterialSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import (
    InputBoolean,
    InputCollection,
    InputFloat,
    InputGeometry,
    InputMaterial,
    InputMenu,
)
from .chain_id import ChainID
from .slice_edge_instances import SliceEdgeInstances


class EnsembleInstance(AssetGeometryGroup):
    """
    Ensemble Instance

    Parameters
    ----------
    points : InputGeometry
        Points
    selection : InputBoolean
        Selection of atoms to apply this node to
    selection_type : InputMenu | Literal["Simple", "Precise"]
        Selection Type
    instances : InputCollection
        Instances
    fraction : InputFloat
        Fraction
    as_points : InputBoolean
        As Points
    point_radius : InputFloat
        Point Radius
    point_material : InputMaterial
        Point Material

    Inputs
    ------
    i.points : GeometrySocket
        Points
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.selection_type : MenuSocket
        Selection Type
    i.instances : CollectionSocket
        Instances
    i.fraction : FloatSocket
        Fraction
    i.as_points : BooleanSocket
        As Points
    i.point_radius : FloatSocket
        Point Radius
    i.point_material : MaterialSocket
        Point Material

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    """

    _name = "Ensemble Instance"
    _asset_name = "Ensemble Instance"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.ensemble_instance"}

    class _Inputs(SocketAccessor):
        points: GeometrySocket
        """Points"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        selection_type: MenuSocket
        """Selection Type"""
        instances: CollectionSocket
        """Instances"""
        fraction: FloatSocket
        """Fraction"""
        as_points: BooleanSocket
        """As Points"""
        point_radius: FloatSocket
        """Point Radius"""
        point_material: MaterialSocket
        """Point Material"""

    class _Outputs(SocketAccessor):
        instances: GeometrySocket
        """Instances"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        points: InputGeometry = None,
        selection: InputBoolean = True,
        selection_type: InputMenu | Literal["Simple", "Precise"] = "Simple",
        instances: InputCollection = None,
        fraction: InputFloat = 1.0,
        as_points: InputBoolean = True,
        point_radius: InputFloat = 0.1,
        point_material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Points": points,
                "Selection": selection,
                "Selection Type": selection_type,
                "Instances": instances,
                "Fraction": fraction,
                "As Points": as_points,
                "Point Radius": point_radius,
                "Point Material": point_material,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        points = tree.inputs.geometry("Points")
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        selection_type = tree.inputs.menu(
            "Selection Type", expanded=True, optional_label=True
        )
        instances = tree.inputs.collection("Instances", optional_label=True)
        fraction = tree.inputs.float(
            "Fraction", 1.0, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        as_points = tree.inputs.boolean("As Points", True)
        with tree.inputs.panel("Point", description="Points"):
            point_radius = tree.inputs.float(
                "Point Radius", 0.1, min_value=0.0, subtype="DISTANCE"
            )
            point_material = tree.inputs.material("Point Material", optional_label=True)
        instances_1 = tree.outputs.geometry("Instances")

        menu_switch = g.MenuSwitch.integer(selection_type, {"Simple": 0, "Precise": 1})
        boolean_math = selection & g.RandomValue.boolean(fraction)
        with g.Frame("Simple Selection"):
            separate_geometry = g.IndexSwitch.geometry(
                menu_switch.o.output,
                (g.SeparateGeometry.point(points, boolean_math).o.selection, points),
            ) >> g.SeparateGeometry.point(selection=as_points)
            instance_on_points = g.InstanceOnPoints(
                points=separate_geometry.o.inverted,
                instance=g.CollectionInfo(collection=instances, separate_children=True),
                instance_index=ChainID(),
                rotation=g.NamedAttribute.input_4x4_matrix("transform").o.attribute,
                pick_instance=True,
            )
            set_material = (
                separate_geometry
                >> g.MeshToPoints(radius=point_radius)
                >> g.SetMaterial(material=point_material)
            )
            join_geometry = g.JoinGeometry(geometry=(set_material, instance_on_points))
        with g.Frame("Precise Selection"):
            group = SliceEdgeInstances(
                instances=instance_on_points, selection=selection
            )
            set_material_1 = g.SetMaterial(
                geometry=g.JoinGeometry(
                    geometry=(group.o.realized_points, set_material)
                ),
                material=point_material,
            )
            join_geometry_1 = g.JoinGeometry(
                geometry=(group.o.instances, set_material_1)
            )
        (
            g.IndexSwitch.geometry(
                menu_switch.o.output, (join_geometry, join_geometry_1)
            )
            >> instances_1
        )

        selection_type.default_value = "Simple"


ASSET = EnsembleInstance

ASSET_METADATA = {
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
