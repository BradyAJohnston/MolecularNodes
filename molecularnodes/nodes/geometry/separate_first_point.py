# Node-group asset "Separate First Point" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputGeometry, InputInteger
from .group_parameter import GroupParameter


class SeparateFirstPoint(AssetGeometryGroup):
    """
    Separate the first point for each `Group ID` and return only those points. Optionally sort by the `Group ID` as well

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    sort : InputBoolean
        Sort by `Group ID` after separating
    group_id : InputInteger
        Group ID

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.sort : BooleanSocket
        Sort by `Group ID` after separating
    i.group_id : IntegerSocket
        Group ID

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Separate First Point"
    _asset_name = "Separate First Point"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Separate the first point for each `Group ID` and return only those points. Optionally sort by the `Group ID` as well"
    }

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        sort: BooleanSocket
        """Sort by `Group ID` after separating"""
        group_id: IntegerSocket
        """Group ID"""

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
        geometry: InputGeometry = None,
        sort: InputBoolean = True,
        group_id: InputInteger = 0,
    ):
        super().__init__(**{"Geometry": geometry, "Sort": sort, "Group ID": group_id})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry", hide_value=True)
        sort = tree.inputs.boolean(
            "Sort",
            True,
            description="Sort by `Group ID` after separating",
            structure_type="SINGLE",
            force_non_field=True,
        )
        group_id = tree.inputs.integer("Group ID", 0, hide_value=True)
        geometry_1 = tree.outputs.geometry("Geometry")

        separate_geometry = g.SeparateGeometry.point(
            geometry, GroupParameter(group_id=group_id).o.is_first
        )
        (
            sort.switch.geometry(
                separate_geometry.o.selection,
                g.SortElements.point(
                    separate_geometry.o.selection, sort_weight=group_id
                ),
            )
            >> geometry_1
        )


ASSET = SeparateFirstPoint

ASSET_METADATA = {
    "description": "Separate the first point for each `Group ID` and return only those points. Optionally sort by the `Group ID` as well",
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
