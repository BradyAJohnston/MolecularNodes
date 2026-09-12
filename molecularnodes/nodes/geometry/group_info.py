# Node-group asset "Group Info" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class GroupInfo(AssetGeometryGroup):
    """
    Group Info

    Parameters
    ----------
    group_id : InputInteger
        Group ID

    Inputs
    ------
    i.group_id : IntegerSocket
        Group ID

    Outputs
    -------
    o.size : IntegerSocket
        Size
    o.index_in_group : IntegerSocket
        Index in Group
    o.index_of_first : IntegerSocket
        Index of First
    o.index_of_last : IntegerSocket
        Index of Last
    """

    _name = "Group Info"
    _asset_name = "Group Info"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        group_id: IntegerSocket
        """Group ID"""

    class _Outputs(SocketAccessor):
        size: IntegerSocket
        """Size"""
        index_in_group: IntegerSocket
        """Index in Group"""
        index_of_first: IntegerSocket
        """Index of First"""
        index_of_last: IntegerSocket
        """Index of Last"""

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
        group_id = tree.inputs.integer("Group ID", 0, hide_value=True)
        size = tree.outputs.integer("Size")
        index_in_group = tree.outputs.integer("Index in Group")
        index_of_first = tree.outputs.integer("Index of First")
        index_of_last = tree.outputs.integer("Index of Last")

        index = g.Index()
        accumulate_field = g.AccumulateField.point.integer(g.Integer(), group_id)
        switch = g.Compare.integer.equal(
            accumulate_field.o.trailing, 0
        ).o.result.switch.integer(true=index)
        switch.point.total(group_id) >> index_of_first
        switch_1 = g.Compare.integer.equal(
            accumulate_field.o.leading, accumulate_field.o.total
        ).o.result.switch.integer(true=index)
        switch_1.point.total(group_id) >> index_of_last

        accumulate_field.o.total >> size
        accumulate_field.o.trailing >> index_in_group


ASSET = GroupInfo

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
