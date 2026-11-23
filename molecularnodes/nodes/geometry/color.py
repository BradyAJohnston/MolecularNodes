# Node-group asset "Color" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    PackageLibrary,
    SocketAccessor,
)


class Color(AssetGeometryGroup):
    """
    Color

    Outputs
    -------
    o.color : ColorSocket
        Read the `Color` attribute from the geometry
    """

    _name = "Color"
    _asset_name = "Color"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.color"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """Read the `Color` attribute from the geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self) -> None:
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        color = tree.outputs.color(
            "Color",
            (0.0, 0.0, 0.0, 1.0),
            description="Read the `Color` attribute from the geometry",
        )

        named_attribute = g.NamedAttribute.color("Color")

        named_attribute >> color


ASSET = Color

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
