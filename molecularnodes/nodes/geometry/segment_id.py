# Node-group asset "Segment ID" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class SegmentID(AssetGeometryGroup):
    """
    Integer representation of the segment IDs present in the structure (MDAnalysis `segid`), assigned in order of the segments read from the topology

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
    o.segid : IntegerSocket
        The `segid` attribute read from the points, an integer index into the segments read from the topology
    """

    _name = "Segment ID"
    _asset_name = "Segment ID"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Integer representation of the segment IDs present in the structure (MDAnalysis `segid`), assigned in order of the segments read from the topology",
        "node_tool_idname": "geometry.segment_id",
    }

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        segid: IntegerSocket
        """The `segid` attribute read from the points, an integer index into the segments read from the topology"""

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
        segid = tree.outputs.integer(
            "segid",
            description="The `segid` attribute read from the points, an integer index into the segments read from the topology",
        )

        AttributeAtIndex(index=index, name="segid") >> segid


ASSET = SegmentID

ASSET_METADATA = {
    "description": "Integer representation of the segment IDs present in the structure (MDAnalysis `segid`), assigned in order of the segments read from the topology",
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
