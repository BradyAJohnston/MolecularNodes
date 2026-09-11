# Node-group asset "Secondary Structure" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from ._shared.attribute_at_index import AttributeAtIndex


class SecondaryStructure(AssetGeometryGroup):
    """
    Secondary Structure

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
    o.sec_struct : IntegerSocket
        sec_struct
    """

    _name = "Secondary Structure"
    _asset_name = "Secondary Structure"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.secondary_structure"}

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        sec_struct: IntegerSocket

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
        sec_struct = tree.outputs.integer("sec_struct")

        AttributeAtIndex(index=index, name="sec_struct") >> sec_struct


ASSET = SecondaryStructure

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
