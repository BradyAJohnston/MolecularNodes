# Node-group asset "B Factor" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from .fallback_float import FallbackFloat


class BFactor(AssetGeometryGroup):
    """
    B Factor

    Outputs
    -------
    o.b_factor : FloatSocket
        b_factor
    """

    _name = "B Factor"
    _asset_name = "B Factor"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.b_factor"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        b_factor: FloatSocket

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self) -> None:
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        b_factor = tree.outputs.float("b_factor")

        FallbackFloat(name="b_factor") >> b_factor


ASSET = BFactor

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
