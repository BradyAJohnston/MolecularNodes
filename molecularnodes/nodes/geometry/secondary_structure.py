# Node-group asset "Secondary Structure" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
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


class SecondaryStructure(AssetGeometryGroup):
    """
    Secondary Structure

    Outputs
    -------
    o.sec_struct : IntegerSocket
        sec_struct
    """

    _name = "Secondary Structure"
    _asset_name = "Secondary Structure"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.secondary_structure"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        sec_struct: IntegerSocket

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self) -> None:
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        sec_struct = tree.outputs.integer("sec_struct")

        named_attribute = g.NamedAttribute.integer("sec_struct")

        named_attribute >> sec_struct


ASSET = SecondaryStructure

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
