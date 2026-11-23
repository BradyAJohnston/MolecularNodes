# Node-group asset "Group Parameter" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from ._shared.index_to_factor import IndexToFactor
from .group_info import GroupInfo


class GroupParameter(AssetGeometryGroup):
    """
    Group Parameter

    Parameters
    ----------
    group_id : InputInteger
        The identifier specifying groupings of the points

    Inputs
    ------
    i.group_id : IntegerSocket
        The identifier specifying groupings of the points

    Outputs
    -------
    o.factor : FloatSocket
        A point's relative position in its group, 0 for the first point and 1 for the last
    o.index : IntegerSocket
        Index of the point within its group, starting at 0 for the first point
    o.count : IntegerSocket
        Number of points in the group
    o.is_first : BooleanSocket
        If the point is the first point in its group
    o.is_last : BooleanSocket
        If the point is the last point in its group
    o.index_of_first : IntegerSocket
        `Index` (in the whole geometry) of the first point in the group
    o.index_of_last : IntegerSocket
        `Index` (in the whole geometry) of the last point in the group
    """

    _name = "Group Parameter"
    _asset_name = "Group Parameter"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        group_id: IntegerSocket
        """The identifier specifying groupings of the points"""

    class _Outputs(SocketAccessor):
        factor: FloatSocket
        """A point's relative position in its group, 0 for the first point and 1 for the last"""
        index: IntegerSocket
        """Index of the point within its group, starting at 0 for the first point"""
        count: IntegerSocket
        """Number of points in the group"""
        is_first: BooleanSocket
        """If the point is the first point in its group"""
        is_last: BooleanSocket
        """If the point is the last point in its group"""
        index_of_first: IntegerSocket
        """`Index` (in the whole geometry) of the first point in the group"""
        index_of_last: IntegerSocket
        """`Index` (in the whole geometry) of the last point in the group"""

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
            description="The identifier specifying groupings of the points",
            hide_value=True,
        )
        factor = tree.outputs.float(
            "Factor",
            description="A point's relative position in its group, 0 for the first point and 1 for the last",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        index = tree.outputs.integer(
            "Index",
            description="Index of the point within its group, starting at 0 for the first point",
        )
        count = tree.outputs.integer(
            "Count", description="Number of points in the group"
        )
        is_first = tree.outputs.boolean(
            "Is First", description="If the point is the first point in its group"
        )
        is_last = tree.outputs.boolean(
            "Is Last", description="If the point is the last point in its group"
        )
        index_of_first = tree.outputs.integer(
            "Index of First",
            description="`Index` (in the whole geometry) of the first point in the group",
        )
        index_of_last = tree.outputs.integer(
            "Index of Last",
            description="`Index` (in the whole geometry) of the last point in the group",
        )

        group_info = GroupInfo(group_id=group_id)
        g.Compare.integer.equal(group_info.o.index_in_group, 0) >> is_first
        (
            IndexToFactor(index=group_info.o.index_in_group, size=group_info.o.size)
            >> factor
        )
        (
            g.Compare.integer.equal(group_info.o.index_in_group, group_info.o.size - 1)
            >> is_last
        )

        group_info.o.index_in_group >> index
        group_info >> count
        group_info.o.index_of_first >> index_of_first
        group_info.o.index_of_last >> index_of_last


ASSET = GroupParameter

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
