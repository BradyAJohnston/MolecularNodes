# Node-group asset "Entity ID" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class EntityID(AssetGeometryGroup):
    """
    The `entity_id` attribute read from the points, corresponding to the unique entities in the structure (that may appear several times as different chains)

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
    o.entity_id : IntegerSocket
        The `entity_id` attribute read from the points, corresponding to the unique entities in the structure (that may appear several times as different chains)
    """

    _name = "Entity ID"
    _asset_name = "Entity ID"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "The `entity_id` attribute read from the points, corresponding to the unique entities in the structure (that may appear several times as different chains)",
        "node_tool_idname": "geometry.entity_id",
    }

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        entity_id: IntegerSocket
        """The `entity_id` attribute read from the points, corresponding to the unique entities in the structure (that may appear several times as different chains)"""

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
        entity_id = tree.outputs.integer(
            "entity_id",
            description="The `entity_id` attribute read from the points, corresponding to the unique entities in the structure (that may appear several times as different chains)",
        )

        AttributeAtIndex(index=index, name="entity_id") >> entity_id


ASSET = EntityID

ASSET_METADATA = {
    "description": "The `entity_id` attribute read from the points, corresponding to the unique entities in the structure (that may appear several times as different chains)",
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
