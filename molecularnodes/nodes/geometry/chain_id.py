# Node-group asset "Chain ID" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class ChainID(AssetGeometryGroup):
    """
    Integer representation of the Chain IDs that were present in the structure

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
    o.chain_id : IntegerSocket
        The `chain_id` attribute read from the points, corresponding to the different chains that were read from the structure
    """

    _name = "Chain ID"
    _asset_name = "Chain ID"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Integer representation of the Chain IDs that were present in the structure",
        "node_tool_idname": "geometry.chain_id",
    }

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        chain_id: IntegerSocket
        """The `chain_id` attribute read from the points, corresponding to the different chains that were read from the structure"""

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
        chain_id = tree.outputs.integer(
            "chain_id",
            description="The `chain_id` attribute read from the points, corresponding to the different chains that were read from the structure",
        )

        AttributeAtIndex(index=index, name="chain_id") >> chain_id


ASSET = ChainID

ASSET_METADATA = {
    "description": "Integer representation of the Chain IDs that were present in the structure",
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
