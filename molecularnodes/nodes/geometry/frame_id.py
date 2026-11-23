# Node-group asset "Frame ID" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class FrameID(AssetGeometryGroup):
    """
    Read the `frame_id` attribute, created when multiple frames from a trajectory are merged into a single structure

    Outputs
    -------
    o.frame_id : IntegerSocket
        Read the `frame_id` attribute from the geometry
    """

    _name = "Frame ID"
    _asset_name = "Frame ID"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Read the `frame_id` attribute, created when multiple frames from a trajectory are merged into a single structure",
        "node_tool_idname": "geometry.residue_id",
    }

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        frame_id: IntegerSocket
        """Read the `frame_id` attribute from the geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self) -> None:
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        frame_id = tree.outputs.integer(
            "frame_id", description="Read the `frame_id` attribute from the geometry"
        )

        named_attribute = g.NamedAttribute.integer("frame_id")

        named_attribute >> frame_id


ASSET = FrameID

ASSET_METADATA = {
    "description": "Read the `frame_id` attribute, created when multiple frames from a trajectory are merged into a single structure",
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
