# Node-group asset "Relative Index" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from .group_parameter import GroupParameter
from .group_pick import GroupPick
from .group_pick_first import GroupPickFirst


class RelativeIndex(AssetGeometryGroup):
    """
    Get information about the points in a `Group ID` such as size and the start and end Indices

    Parameters
    ----------
    group_id : InputInteger
        Define the groups to get information about

    Inputs
    ------
    i.group_id : IntegerSocket
        Define the groups to get information about

    Outputs
    -------
    o.group_size : IntegerSocket
        The number of points in each `Group ID`
    o.relative_index : IntegerSocket
        The relative index of each point within each `Group ID`, starting from 0
    o.first_index : IntegerSocket
        The `Index` of the first point in each `Group ID`
    o.last_index : IntegerSocket
        The `Index` of the last point in each `Group ID`
    """

    _name = "Relative Index"
    _asset_name = "Relative Index"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Get information about the points in a `Group ID` such as size and the start and end Indices",
        "node_tool_idname": "geometry.relative_index",
    }

    class _Inputs(SocketAccessor):
        group_id: IntegerSocket
        """Define the groups to get information about"""

    class _Outputs(SocketAccessor):
        group_size: IntegerSocket
        """The number of points in each `Group ID`"""
        relative_index: IntegerSocket
        """The relative index of each point within each `Group ID`, starting from 0"""
        first_index: IntegerSocket
        """The `Index` of the first point in each `Group ID`"""
        last_index: IntegerSocket
        """The `Index` of the last point in each `Group ID`"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        group_id: InputInteger = 0,
    ):
        super().__init__(**{"Group ID": group_id})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="Define the groups to get information about",
            hide_value=True,
        )
        group_size = tree.outputs.integer(
            "Group Size", description="The number of points in each `Group ID`"
        )
        relative_index = tree.outputs.integer(
            "Relative Index",
            description="The relative index of each point within each `Group ID`, starting from 0",
        )
        first_index = tree.outputs.integer(
            "First Index",
            description="The `Index` of the first point in each `Group ID`",
        )
        last_index = tree.outputs.integer(
            "Last Index", description="The `Index` of the last point in each `Group ID`"
        )

        group = GroupParameter(group_id=group_id)
        _group_1 = GroupPick(pick=group.o.is_last, group_id=group_id)
        GroupPickFirst(pick=group.o.is_first, group_id=group_id).o.index >> first_index
        GroupPickFirst(pick=group.o.is_last, group_id=group_id).o.index >> last_index
        _group_2 = GroupPick(pick=group.o.is_first, group_id=group_id)

        group.o.group_size >> group_size
        group.o.relative_index >> relative_index


ASSET = RelativeIndex

ASSET_METADATA = {
    "description": "Get information about the points in a `Group ID` such as size and the start and end Indices",
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
