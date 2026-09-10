# Node-group asset 'Evaluate Per Group' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ClosureSocket,
    GeometrySocket,
    IntegerSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputClosure, InputGeometry, InputInteger, InputMenu
from .chain_id import ChainID


class EvaluatePerGroup(AssetGeometryGroup):
    """
    Evaluate Per Group

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to split into two parts
    closure : InputClosure
        Closure
    group : InputMenu | Literal["chain_id", "Group ID"]
        Group
    group_id : InputInteger
        Group ID

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to split into two parts
    i.closure : ClosureSocket
        Closure
    i.group : MenuSocket
        Group
    i.group_id : IntegerSocket
        Group ID

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    o.instances : GeometrySocket
        Instances
    """

    _name = "Evaluate Per Group"
    _asset_name = "Evaluate Per Group"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to split into two parts"""
        closure: ClosureSocket
        """Closure"""
        group: MenuSocket
        """Group"""
        group_id: IntegerSocket
        """Group ID"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        instances: GeometrySocket
        """Instances"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        closure: InputClosure = None,
        group: InputMenu | Literal["chain_id", "Group ID"] = "chain_id",
        group_id: InputInteger = 0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Closure": closure,
                "Group": group,
                "Group ID": group_id,
            }
        )

    def _build_group(self, tree):
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry to split into two parts"
        )
        closure = tree.inputs.closure("Closure")
        group = tree.inputs.menu("Group", expanded=True, optional_label=True)
        group_id = tree.inputs.integer("Group ID", 0, hide_value=True)
        geometry_1 = tree.outputs.geometry("Geometry")
        instances = tree.outputs.geometry("Instances")

        with g.Frame("Count number of chains"):
            menu_switch = g.MenuSwitch.integer(
                group, {"chain_id": ChainID(), "Group ID": group_id}
            )
            separate_geometry = g.SeparateGeometry.point(
                geometry,
                ~g.AccumulateField.point.integer(
                    group_index=menu_switch.o.output
                ).o.trailing,
            )
            domain_size = g.DomainSize(geometry=separate_geometry.o.selection)
        repeat_zone = g.RepeatZone(domain_size.o.point_count)
        geometry_2 = repeat_zone.items.geometry("Geometry")
        instances_1 = repeat_zone.items.geometry("Instances")
        sample_index = g.SampleIndex(
            geometry=separate_geometry.o.selection,
            value=menu_switch.o.output,
            index=repeat_zone.iteration,
            data_type="INT",
        )
        separate_geometry_1 = g.SeparateGeometry.point(
            geometry, g.Compare.integer.equal(sample_index, menu_switch.o.output)
        )
        evaluate_closure = g.EvaluateClosure(closure)
        evaluate_closure.inputs.geometry("Geometry", separate_geometry_1.o.selection)
        evaluate_closure.inputs.integer("group_id", sample_index)
        geometry_3 = evaluate_closure.outputs.geometry("Geometry")
        store_named_attribute = geometry_3 >> g.StoreNamedAttribute.point.integer(
            name="group_id", value=repeat_zone.iteration
        )
        join_geometry = g.JoinGeometry(
            geometry=(instances_1.current, g.GeometryToInstance(store_named_attribute))
        )
        (
            g.JoinGeometry(geometry=(geometry_2.current, store_named_attribute))
            >> geometry_2.next
        )
        join_geometry >> instances_1.next

        geometry_2.result >> geometry_1
        instances_1.result >> instances

        group.default_value = "chain_id"


ASSET = EvaluatePerGroup

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
