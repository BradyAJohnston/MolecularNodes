# Node-group asset "URes ID" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from ._shared.attribute_at_index import AttributeAtIndex


class UResID(AssetGeometryGroup):
    """
    URes ID

    Parameters
    ----------
    index : InputInteger
        Index

    Inputs
    ------
    i.index : IntegerSocket
        Index

    Outputs
    -------
    o.ures_id : IntegerSocket
        ures_id
    o.size : IntegerSocket
        The total of all of the values in the corresponding group
    """

    _name = "URes ID"
    _asset_name = "URes ID"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        ures_id: IntegerSocket
        size: IntegerSocket
        """The total of all of the values in the corresponding group"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        index: InputInteger = 0,
    ):
        super().__init__(**{"Index": index})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        index = tree.inputs.integer(
            "Index", 0, min_value=0, hide_value=True, default_input="INDEX"
        )
        ures_id = tree.outputs.integer("ures_id")
        size = tree.outputs.integer(
            "Size",
            description="The total of all of the values in the corresponding group",
        )

        group = AttributeAtIndex(index=index, name="ures_id")
        accumulate_field = g.AccumulateField.point.integer(group_index=group)

        group >> ures_id
        accumulate_field.o.total >> size


ASSET = UResID

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
