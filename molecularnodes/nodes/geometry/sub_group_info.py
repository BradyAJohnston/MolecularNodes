# Node-group asset "Sub Group Info" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from .group_info import GroupInfo
from .integer_run import IntegerRun
from .offset_integer import OffsetInteger


class SubGroupInfo(AssetGeometryGroup):
    """
    Sub Group Info

    Parameters
    ----------
    sub_group_id : InputInteger
        Sub Group ID
    group_id : InputInteger
        Group ID

    Inputs
    ------
    i.sub_group_id : IntegerSocket
        Sub Group ID
    i.group_id : IntegerSocket
        Group ID

    Outputs
    -------
    o.size : IntegerSocket
        Size of each comptued `Group ID`
    o.group_id : IntegerSocket
        A unique `Group ID` that increases whenever `Sub Group ID` or `Group ID` change
    o.index_of_first : IntegerSocket
        The absolute `Index` of the _first_ point in each computed `Group ID`
    o.index_of_last : IntegerSocket
        The absolute `Index` of the _last_ point in each computed `Group ID`
    o.index_in_group_id : IntegerSocket
        Index within each computed `Group ID`, starting at 0
    o.sub_group_id : IntegerSocket
        Sub Group ID
    o.sub_group_total : IntegerSocket
        Sub Group Total
    """

    _name = "Sub Group Info"
    _asset_name = "Sub Group Info"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        sub_group_id: IntegerSocket
        """Sub Group ID"""
        group_id: IntegerSocket
        """Group ID"""

    class _Outputs(SocketAccessor):
        size: IntegerSocket
        """Size of each comptued `Group ID`"""
        group_id: IntegerSocket
        """A unique `Group ID` that increases whenever `Sub Group ID` or `Group ID` change"""
        index_of_first: IntegerSocket
        """The absolute `Index` of the _first_ point in each computed `Group ID`"""
        index_of_last: IntegerSocket
        """The absolute `Index` of the _last_ point in each computed `Group ID`"""
        index_in_group_id: IntegerSocket
        """Index within each computed `Group ID`, starting at 0"""
        sub_group_id: IntegerSocket
        """Sub Group ID"""
        sub_group_total: IntegerSocket
        """Sub Group Total"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        sub_group_id: InputInteger = 0,
        group_id: InputInteger = 0,
    ):
        super().__init__(**{"Sub Group ID": sub_group_id, "Group ID": group_id})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        sub_group_id = tree.inputs.integer("Sub Group ID", 0, hide_value=True)
        group_id = tree.inputs.integer("Group ID", 0, hide_value=True)
        size = tree.outputs.integer(
            "Size", description="Size of each comptued `Group ID`"
        )
        group_id_1 = tree.outputs.integer(
            "Group ID",
            description="A unique `Group ID` that increases whenever `Sub Group ID` or `Group ID` change",
        )
        index_of_first = tree.outputs.integer(
            "Index of First",
            description="The absolute `Index` of the _first_ point in each computed `Group ID`",
            min_value=0,
        )
        index_of_last = tree.outputs.integer(
            "Index of Last",
            description="The absolute `Index` of the _last_ point in each computed `Group ID`",
            min_value=0,
        )
        index_in_group_id = tree.outputs.integer(
            "Index in Group ID",
            description="Index within each computed `Group ID`, starting at 0",
            min_value=0,
        )
        sub_group_id_1 = tree.outputs.integer("Sub Group ID")
        sub_group_total = tree.outputs.integer("Sub Group Total")

        group = IntegerRun(value=sub_group_id, group_id=group_id)
        group_1 = GroupInfo(group_id=group.o.group_id)
        compare = g.Compare.integer.not_equal(
            group.o.group_id, OffsetInteger(integer=group.o.group_id, offset=1)
        )
        accumulate_field = g.AccumulateField.point.integer(compare, group_id)

        group_1 >> size
        group.o.group_id >> group_id_1
        group_1.o.index_of_first >> index_of_first
        group_1.o.index_of_last >> index_of_last
        group_1.o.index_in_group >> index_in_group_id
        accumulate_field.o.trailing >> sub_group_id_1
        accumulate_field.o.total >> sub_group_total


ASSET = SubGroupInfo

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
